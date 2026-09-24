from __future__ import absolute_import, division, print_function
from __future__ import annotations

from typing import Any
from typing import TYPE_CHECKING

from compas.datastructures.mesh.mesh import Mesh
from compas.geometry import discrete_coons_patch
from compas.itertools import pairwise
from compas.tolerance import TOL

from compas_singular.datastructures.mesh import meshes_join_and_weld
from compas_singular.datastructures.mesh_quad_coarse import CoarseQuadMesh
from compas_singular.datastructures.mesh_quad_coarse.coarse_network import (
    DEFAULT_PRECISION,
    check_faces,
    check_network,
    faces_from_network,
    weld_network,
)
from compas_singular.datastructures.mesh_quad_pseudo import PseudoQuadMesh
from compas_singular.datastructures.mesh_quad_coarse.patterns import (
    PATTERNS,
    create_pattern,
    patch_divisions,
    pattern_morph,
    pattern_morph_triangle,
    reconcile_strip_densities,
)

if TYPE_CHECKING:
    from compas_singular.datastructures import QuadMesh

__all__ = [	'CoarsePseudoQuadMesh']


class CoarsePseudoQuadMesh(PseudoQuadMesh, CoarseQuadMesh):

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super(CoarsePseudoQuadMesh, self).__init__(*args, **kwargs)

    # --------------------------------------------------------------------------
    # constructors
    # --------------------------------------------------------------------------

    @classmethod
    def from_coarse_polylines(cls, polylines: list[list[list[float]]], poles: list[list[float]] | None = None, holes: list[list[float]] | None = None,
                              precision: int | None = None, tol: float | None = None, collect_strips: bool = True) -> "CoarsePseudoQuadMesh":
        """**A coarse quad layout from a drawn network of edge-curves.**

        One polyline is one coarse EDGE: its two ENDS are corners of the layout, and
        the points between them are that edge's shape. Nothing else is needed -- no
        domain boundary, no separate face list, no field. The layout's own boundary
        is whatever ends up with one adjacent patch.

        This is the inverse of ``coarse_edges_to_curves``, so a layout baked out one
        curve per edge comes back in through here, curvature included::

            coarse = CoarsePseudoQuadMesh.from_coarse_polylines(curves)
            coarse.set_strips_density_target(t=0.5)
            dense = coarse.densification()

        **It validates rather than repairs.** A dangling curve, a corner sitting on
        another curve's interior, two curves crossing with no corner between them,
        or a patch that is not a triangle or a quadrilateral all raise, naming the
        curve or the corner. Each is a drawing error with a one-second fix in the
        CAD session, and repairing it here would silently hand back a layout nobody
        drew. That is the whole difference from ``Mesh.from_polylines``, which
        additionally discards any patch with every corner on a boundary polyline and
        so returns NOTHING for a square with one cut across it.

        **A triangular patch is a pseudo-quad, not a defect.** It is registered in
        ``attributes['face_pole']`` with a collapsed corner, which is what makes it
        densifiable and what ``collect_strips`` needs; ``poles`` only chooses WHICH
        of its corners collapses.

        Parameters
        ----------
        polylines : list[list[[x, y, z]]]
            One polyline per coarse edge. A straight edge is two points.
        poles : list[[x, y, z]], optional
            Preferred pole positions. A triangular patch is a pseudo-quad with or
            without one; a position given here decides which corner collapses.
        holes : list[[x, y, z]], optional
            One point inside each hole of the domain. A hole is a bounded cycle of
            the network and is otherwise indistinguishable from a patch, so without
            this it is filled.
        precision : int, optional
            Decimals the endpoints are welded at. Default is COMPAS's own, 3.
        tol : float, optional
            Distance below which a point is ON a curve, for the validity checks.
            Defaults to the weld resolution.
        collect_strips : bool, optional
            Collect strip data before returning, as
            :meth:`CoarseQuadMesh.from_quad_mesh` does.

        Returns
        -------
        CoarsePseudoQuadMesh
            Carrying the drawn shape of every edge -- see
            :meth:`CoarseQuadMesh.edges_to_curves`.

        Raises
        ------
        ValueError
            On a network that is not a clean planar arrangement of triangular and
            quadrilateral patches.
        """
        vertices, edges = weld_network(polylines, precision=precision)
        if tol is None:
            tol = 10.0 ** -(DEFAULT_PRECISION if precision is None else precision)
        check_network(vertices, edges, tol)

        faces = faces_from_network(vertices, edges, holes=holes)
        if not faces:
            raise ValueError(
                'no patch was recovered from {} curve(s). Every cycle in the '
                'network was either the outside or a hole.'.format(len(polylines)))
        check_faces(faces, vertices)

        # ``vertices`` is a LIST, so the mesh keys are its indices, 0..n-1, and the
        # edge list below can be read straight off the network. Every 3-sided face
        # is registered as a pseudo-quad here; ``poles`` only moves the collapse.
        mesh = cls.from_vertices_and_faces_with_poles(vertices, faces, poles or [])

        orphans = [i for i, (u, v, _points) in enumerate(edges)
                   if mesh.halfedge[u].get(v) is None
                   and mesh.halfedge[v].get(u) is None]
        if orphans:
            raise ValueError(
                'polyline(s) {} bound no patch on either side -- they lie between '
                'the outside and a hole, or between two holes.'.format(orphans))

        if not mesh.is_manifold():
            raise ValueError(
                'the recovered layout is not manifold. More than two patches meet '
                'along one edge, or a corner is pinched between patches that do '
                'not share an edge.')

        # The drawn shape of each edge, so ``densification`` reproduces the curves
        # rather than chording them. Without this the layout keeps the topology of
        # the drawing and loses its geometry.
        mesh.set_edges_to_curves({(u, v): points for u, v, points in edges})

        if collect_strips:
            mesh.collect_strips()
        return mesh

    # --------------------------------------------------------------------------
    # densification
    # --------------------------------------------------------------------------

    def quad_mesh(self, boundary_curvature: bool = True, skeleton_curvature: bool = True,
                 overwrite_edges_to_curves: dict[tuple[int, int], list[list[float]]] | None = None, field: Any = None, pattern_overwrite: "str | dict[int, str] | None" = None) -> "QuadMesh":
        """Generate a dense quad mesh from this layout, in a chosen pattern.

        Parameters
        ----------
        boundary_curvature : bool, optional
            Use the shape :meth:`edges_to_curves` has stored for edges on the
            layout's own boundary, instead of chording them. Defaults to True.
            Ignored -- treated as True -- when ``overwrite_edges_to_curves`` is
            given. See :meth:`densify`.
        skeleton_curvature : bool, optional
            Same, for the edges that are NOT on the boundary -- the interior
            edges a skeleton or field decomposition traced as separatrices.
            Defaults to True. Ignored -- treated as True -- when
            ``overwrite_edges_to_curves`` is given.
        overwrite_edges_to_curves : dict, optional
            Coarse edges ``(u, v)`` pointing to a curve to densify them along,
            each curve a list of XYZ points, overriding whatever
            :meth:`edges_to_curves` has stored -- for every edge, regardless of
            ``boundary_curvature`` / ``skeleton_curvature``.
        field : optional
            A ``CrossField``. Only valid with the ``ortho`` pattern -- see
            :meth:`densification`.
        pattern : str, optional
            One of ``compas_singular...patterns.PATTERNS``: ``'ortho'`` (the
            plain grid), ``'diagonal'`` (the grid with both patch diagonals cut
            through it) or ``'fan'`` (four polar fans meeting at the centre).

        Returns
        -------
        QuadMesh
            The dense mesh, also stored on this one -- ``get_quad_mesh()``.

        Notes
        -----
        ``'diagonal'`` and ``'fan'`` constrain the layout: a patch using either
        needs the same even density on both strips crossing it -- for the
        diagonal, that is what keeps it one unbroken line. Densities are
        raised to that in place before any patch is built (see
        ``reconcile_strip_densities``). On a pseudo-quad face ``'fan'`` is
        three polar fans instead of four, one per corner of the triangle -- the
        pole included -- meeting at its centre (see ``_pattern_fan_triangle``).
        It divides its sides like every other pattern, so it welds against any
        neighbour.
        """
        return self.densify(pattern_overwrite=pattern_overwrite,
                            boundary_curvature=boundary_curvature,
                            skeleton_curvature=skeleton_curvature,
                            overwrite_edges_to_curves=overwrite_edges_to_curves, field=field)

    def densify(self, pattern_overwrite: "str | dict[int, str] | None" = None, boundary_curvature: bool = True, skeleton_curvature: bool = True,
               overwrite_edges_to_curves: dict[tuple[int, int], list[list[float]]] | None = None, field: Any = None) -> "QuadMesh":
        """Build one dense patch per coarse face and weld them together.

        The loop is the whole of it: take a face's four sides, ask the pattern
        for a template at the divisions those sides imply, morph the template
        onto them. Everything a pattern differs in lives in the template;
        everything the patterns share -- curves, poles, welding, pole
        bookkeeping -- lives here and is written once.

        Parameters
        ----------
        pattern : str, optional
            A key of ``PATTERNS``.
        boundary_curvature : bool, optional
            Use the stored :meth:`edges_to_curves` shape for edges on the
            layout's own boundary. Defaults to True. Ignored -- treated as True
            -- when ``overwrite_edges_to_curves`` is given.
        skeleton_curvature : bool, optional
            Same, for the interior edges. Defaults to True. Ignored -- treated
            as True -- when ``overwrite_edges_to_curves`` is given.
        overwrite_edges_to_curves : dict, optional
            Overrides whatever :meth:`edges_to_curves` has stored, for every
            edge.
        field : optional

        Returns
        -------
        QuadMesh
        """
        if overwrite_edges_to_curves is not None:
            edges_to_curves = overwrite_edges_to_curves
        else:
            edges_to_curves = self._filtered_edges_to_curves(boundary_curvature, skeleton_curvature)

        if pattern_overwrite is None:
            pattern_dict = self.dense_patterns()
        else:
            if isinstance(pattern_overwrite, dict):
                pattern_dict = pattern_overwrite
            elif isinstance(pattern_overwrite, str):
                pattern_dict = {fkey: pattern_overwrite for fkey in self.faces()}

        if field is not None:
            # The field owns this: it carries its own background and builds its
            # own point locator, so a layout from ANY source -- a skeleton
            # decomposition, a hand-built mesh, one read back out of a document
            # -- can be densified with patch interiors that follow it, instead of
            # the bilinear blend of its own four sides that the Coons morph
            # gives and that never consults a field.
            #
            # Nothing about ``field`` is type-checked, so any object offering
            # ``densify(coarse, edges_to_curves=...)`` works.

            patterns = list(pattern_dict.values())
            if patterns.count('ortho') != len(patterns):
                print(
                    'a field integrates the patch interiors itself, so it cannot '
                    'be combined with the patterns other than ortho. Turn the field off to use the patterns.'
                    )
            dense, _stats = field.densify(self, edges_to_curves=edges_to_curves)
            self.set_quad_mesh(dense)
            return self.get_quad_mesh()
            
        # a pattern may constrain the layout -- fan to tile at all, diagonal to
        # stay continuous -- and since a density is shared along a whole strip,
        # it has to do so before any patch is built
        # Modifies the densities in place; what it RETURNS is {skey: (old, new)},
        # not a pattern map, so it must not be assigned back to pattern_dict.
        reconcile_strip_densities(self, pattern_dict)
        
        edge_strip = {}
        for strip, edges in self.strips(data=True):
            for u, v in edges:
                edge_strip[u, v] = strip
                edge_strip[v, u] = strip

        pole_map = {TOL.geometric_key(self.vertex_coordinates(pole)) for pole in self.poles()}

        meshes = []
        for fkey in self.faces():
            sides = self._patch_sides(fkey, edge_strip, edges_to_curves)
            kind = pattern_dict[fkey]

            if kind not in PATTERNS:
                raise ValueError('unknown pattern {!r} -- pick one of {}'.format(pattern_overwrite, sorted(PATTERNS)))

            nu, nw = patch_divisions(sides)
            if kind == 'fan' and any(side is None for side in sides):
                # the fan's centres are the patch corners, and a collapsed
                # quad has three: its own template, and a morph that treats the
                # three alike
                lam, template_faces = create_pattern('fan_triangle', nu, nw)
                vertices, faces = pattern_morph_triangle(lam, template_faces, sides)
            else:
                uw, template_faces = create_pattern(kind, nu, nw)
                vertices, faces = pattern_morph(uw, template_faces, sides)
            meshes.append(PseudoQuadMesh.from_vertices_and_faces_with_face_poles(vertices, faces))

        # a face with two coincident corners is a pole; remember which, by the
        # geometric keys of its corners, because welding is about to renumber
        # everything. The weld keys vertices the same way, so the corners of a
        # welded face are its corners before the weld with the repeat dropped.
        # Not by its centre: rounded to the same digits, a centre can land on
        # either side of one -- measured, 0.3125 became 0.313 before the weld
        # and 0.312 after, and the pole face was silently dropped.
        face_pole_map = {}
        for mesh in meshes:
            for fkey in mesh.faces():
                corners = [TOL.geometric_key(mesh.vertex_coordinates(v)) for v in mesh.face_vertices(fkey)]
                for gu, gv in pairwise(corners + corners[:1]):
                    if gu in pole_map and gu == gv:
                        face_pole_map[frozenset(corners)] = gu
                        break

        self.set_quad_mesh(meshes_join_and_weld(meshes))
        dense = self.get_quad_mesh()

        face_pole = {}
        for fkey in dense.faces():
            corners = {TOL.geometric_key(dense.vertex_coordinates(vkey)): vkey for vkey in dense.face_vertices(fkey)}
            gkey = face_pole_map.get(frozenset(corners))
            if gkey is not None:
                face_pole[fkey] = corners[gkey]

        dense.attributes['face_pole'] = face_pole
        return dense

    def _patch_sides(self, fkey: int, edge_strip: dict[tuple[int, int], int], edges_to_curves: dict[tuple[int, int], list[list[float]]] | None = None) -> list[list[list[float]] | None]:
        """The four side polylines of one coarse face.

        Each side is densified to its own strip's density, along the curve
        ``edges_to_curves`` gave it or along the straight chord otherwise. The
        four come back in the order ``discrete_coons_patch`` documents --
        ``[ab, bc, dc, ad]``, with the last two reversed -- so that ``u`` runs
        a->b and ``w`` runs a->d.

        A pseudo-quad has only three sides; the missing one is put back as
        ``None`` at the pole's position, which is how the morph learns which
        corner collapsed.

        Parameters
        ----------
        fkey : hashable
            A coarse face key.
        edge_strip : dict
            ``(u, v) -> strip key``, both ways round.
        edges_to_curves : dict, optional
            As in :meth:`quad_mesh`.

        Returns
        -------
        list
            ``[ab, bc, dc, ad]``, one of which may be None.
        """
        polylines = []
        for u, v in self.face_halfedges(fkey):
            d = self.get_strip_density(edge_strip[u, v])
            polylines.append(self._create_patch_edge(u, v, d, edges_to_curves))

        if self.is_face_pseudo_quad(fkey):
            pole = self.attributes['face_pole'][fkey]
            polylines.insert(self.face_vertices(fkey).index(pole), None)

        ab, bc, cd, da = polylines
        return [ab, bc, cd[::-1] if cd else None, da[::-1] if da else None]

    # NOTE old densification has to be phased out
    def densification(self, boundary_curvature: bool = True, skeleton_curvature: bool = True,
                      overwrite_edges_to_curves: dict[tuple[int, int], list[list[float]]] | None = None, field: Any = None) -> "QuadMesh":
        """Generate a denser quad mesh from the coarse quad mesh and its strip densities.

        Parameters
        ----------
        boundary_curvature : bool, optional
            Use the stored :meth:`edges_to_curves` shape for edges on the
            layout's own boundary. Defaults to True. Ignored -- treated as True
            -- when ``overwrite_edges_to_curves`` is given.
        skeleton_curvature : bool, optional
            Same, for the interior edges. Defaults to True. Ignored -- treated
            as True -- when ``overwrite_edges_to_curves`` is given.
        overwrite_edges_to_curves : dict, optional
            A dictionary with edges (u, v) pointing to a curve for
            densification, overriding whatever :meth:`edges_to_curves` has
            stored -- for every edge. The curves are lists of XYZ points.
        field : optional
            A ``CrossField`` (from ``FieldDecomposition.get_field()`` or
            ``CrossField.from_boundary(...)``). With it, each patch INTERIOR is
            integrated from the field instead of blended from its own four
            sides -- which is the only way a guide curve reaches a patch it
            forced no topology in. The layout does not have to have come from
            that field: a skeleton layout and a field solved on the same walls
            work together. Boundaries stay fixed either way, so the result still
            welds into a mesh whose strips can be collected and edited.

        Returns
        -------
        QuadMesh
            The dense mesh, also stored on this one -- ``get_quad_mesh()``.
        """
        if overwrite_edges_to_curves is not None:
            edges_to_curves = overwrite_edges_to_curves
        else:
            edges_to_curves = self._filtered_edges_to_curves(boundary_curvature, skeleton_curvature)

        if field is not None:
            # The field owns this: it carries its own background and builds its
            # own point locator, so a layout from ANY source -- a skeleton
            # decomposition, a hand-built mesh, one read back out of a document
            # -- can be densified with patch interiors that follow it, instead of
            # the bilinear blend of its own four sides that ``discrete_coons_patch``
            # gives and that never consults a field.
            #
            # Imported here and not at module scope: ``framefield.densify``
            # imports ``PseudoQuadMesh`` and ``meshes_join_and_weld`` from this
            # package, so a top-level import is a circular one -- this module is
            # reached while ``compas_singular.datastructures`` is still being
            # initialised. Nothing about ``field`` is type-checked, so any object
            # offering ``densify(coarse, edges_to_curves=...)`` works.
            dense, _stats = field.densify(self, edges_to_curves=edges_to_curves)
            self.set_quad_mesh(dense)
            return self.get_quad_mesh()

        edge_strip = {}
        for strip, edges in self.strips(data=True):
            for u, v in edges:
                edge_strip[u, v] = strip
                edge_strip[v, u] = strip

        pole_map = [TOL.geometric_key(self.vertex_coordinates(pole)) for pole in self.poles()]

        meshes = []
        for fkey in self.faces():
            polylines = []
            for u, v in self.face_halfedges(fkey):
                d = self.get_strip_density(edge_strip[u, v])
                polylines.append(self._create_patch_edge(u, v, d, edges_to_curves))

            if self.is_face_pseudo_quad(fkey):
                pole = self.attributes['face_pole'][fkey]
                idx = self.face_vertices(fkey).index(pole)
                polylines.insert(idx, None)

            ab, bc, cd, da = polylines

            if cd:
                dc = cd[::-1]
            else:
                dc = None

            if da:
                ad = da[::-1]
            else:
                ad = None

            vertices, faces = discrete_coons_patch(ab, bc, dc, ad)
            faces = [[u for u, v in pairwise(face + face[:1]) if u != v] for face in faces]
            mesh = PseudoQuadMesh.from_vertices_and_faces_with_face_poles(vertices, faces)
            meshes.append(mesh)

        face_pole_map = {}
        for mesh in meshes:
            for fkey in mesh.faces():
                for u, v in pairwise(mesh.face_vertices(fkey) + mesh.face_vertices(fkey)[: 1]):
                    if TOL.geometric_key(mesh.vertex_coordinates(u)) in pole_map and TOL.geometric_key(mesh.vertex_coordinates(u)) == TOL.geometric_key(mesh.vertex_coordinates(v)):
                        face_pole_map[TOL.geometric_key(mesh.face_center(fkey))] = TOL.geometric_key(mesh.vertex_coordinates(u))
                        break

        self.set_quad_mesh(meshes_join_and_weld(meshes))

        face_pole = {}
        for fkey in self.get_quad_mesh().faces():
            if TOL.geometric_key(self.get_quad_mesh().face_center(fkey)) in face_pole_map:
                for vkey in self.get_quad_mesh().face_vertices(fkey):
                    if TOL.geometric_key(self.get_quad_mesh().vertex_coordinates(vkey)) == face_pole_map[TOL.geometric_key(self.get_quad_mesh().face_center(fkey))]:
                        face_pole[fkey] = vkey
                        break

        self.get_quad_mesh().attributes['face_pole'] = face_pole
        return self.get_quad_mesh()

    # NOTE old tests

    # def quad_mesh(self, edges_to_curves=None, field=None, pattern="ortho"):
    #     # A layout that KNOWS the shape of its edges -- one from
    #     # ``from_coarse_polylines``, say -- does not make the caller hand them back.
    #     # A layout from any other constructor stores none, so this is ``{}``, which
    #     # is falsy and falls through to the straight-chord branch exactly as before.
    #     if edges_to_curves is None:
    #         edges_to_curves = self.edges_to_curves()

    #     if pattern == "ortho":
    #         return self.densification(edges_to_curves=edges_to_curves, field=field)

    #     elif pattern == "diagonal":
    #         edge_strip = {}
    #         for strip, edges in self.strips(data=True):
    #             for u, v in edges:
    #                 edge_strip[u, v] = strip
    #                 edge_strip[v, u] = strip
    #         meshes = []
    #         for fkey in self.faces():
    #             vkeys = self.face_vertices(fkey)
    #             face_corners = [self.vertex_coordinates(vkey) for vkey in vkeys]
    #             polylines = []
    #             for u, v in self.face_halfedges(fkey):
    #                 d = self.get_strip_density(edge_strip[u, v])
    
    #                 if edges_to_curves:
    #                     # d + 1 points (not d) to match the straight-chord branch below --
    #                     # linspace(0, 1, d) samples one point short and raises for d == 1.
    #                     polyline = []
    #                     if (u, v) in edges_to_curves:
    #                         curve = Polyline(edges_to_curves[u, v])
    #                         polyline = [curve.point_at(t) for t in linspace(0, 1, d + 1)]
    #                     else:
    #                         curve = Polyline(edges_to_curves[v, u])
    #                         polyline = [curve.point_at(t) for t in linspace(0, 1, d + 1)]
    #                         polyline[:] = polyline[::-1]
    #                 else:
    #                     polyline = []
    #                     curve = Polyline([self.vertex_coordinates(u), self.vertex_coordinates(v)])
    #                     for i in range(0, d + 1):
    #                         point = curve.point_at(float(i) / float(d))
    #                         polyline.append(point)
    
    #                 polylines.append(polyline)
    
    #             if self.is_face_pseudo_quad(fkey):
    #                 pole = self.attributes['face_pole'][fkey]
    #                 idx = self.face_vertices(fkey).index(pole)
    #                 polylines.insert(idx, None)
    
    #             ab, bc, cd, da = polylines
    
    #             if cd:
    #                 dc = cd[::-1]
    #             else:
    #                 dc = None
    
    #             if da:
    #                 ad = da[::-1]
    #             else:
    #                 ad = None


    #             nx, ny = 4, 4
    #             mesh = Mesh.from_meshgrid(dx=1.0, nx=nx, ny=ny)  # n x n quads on unit square

    #             # split the quads that sit on either diagonal of the grid into triangles,
    #             # tracing both diagonals as straight lines scaled by the grid's own
    #             # aspect ratio (nx : ny) so they stay continuous even if nx != ny
    #             n_cells, m_cells = nx, ny
    #             steps = max(n_cells, m_cells)

    #             main_diagonal = set()   # bottom-left -> top-right
    #             anti_diagonal = set()   # bottom-right -> top-left
    #             for k in range(steps + 1):
    #                 t = k / steps if steps else 0
    #                 i = round(t * (n_cells - 1))
    #                 main_diagonal.add((i, round(t * (m_cells - 1))))
    #                 anti_diagonal.add((i, round((1 - t) * (m_cells - 1))))

    #             for idx, fkey in enumerate(list(mesh.faces())):
    #                 i, j = divmod(idx, m_cells)
    #                 v0, v1, v2, v3 = mesh.face_vertices(fkey)
    #                 if (i, j) in main_diagonal:
    #                     mesh.split_face(fkey, v0, v2)
    #                 elif (i, j) in anti_diagonal:
    #                     mesh.split_face(fkey, v1, v3)

    #             A, B, C, D = face_corners[0], face_corners[1], face_corners[2], face_corners[3]

    #             # find grid bounds once
    #             xs = [mesh.vertex_attribute(v,'x') for v in mesh.vertices()]
    #             ys = [mesh.vertex_attribute(v,'y') for v in mesh.vertices()]
    #             x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)

    #             for v in mesh.vertices():
    #                 u = (mesh.vertex_attribute(v,'x')-x0)/(x1-x0) if x1!=x0 else 0
    #                 w = (mesh.vertex_attribute(v,'y')-y0)/(y1-y0) if y1!=y0 else 0
    #                 x = (1-u)*(1-w)*A[0]+u*(1-w)*B[0]+u*w*C[0]+(1-u)*w*D[0]
    #                 y = (1-u)*(1-w)*A[1]+u*(1-w)*B[1]+u*w*C[1]+(1-u)*w*D[1]
    #                 mesh.vertex_attributes(v, "xyz", [x,y,0])
    #             meshes.append(mesh)

    #             """
    #             vertices, faces = discrete_coons_patch(ab, bc, dc, ad)

    #             n = len(ab) if ab else len(dc)
    #             m = len(bc) if bc else len(ad)
    #             n_cells, m_cells = n - 1, m - 1
    #             steps = max(n_cells, m_cells)

    #             # trace both diagonals as straight lines across the (n_cells x m_cells)
    #             # grid, scaled by the patch's own aspect ratio so they stay continuous
    #             # even when the patch is not square
    #             main_diagonal = set()   # a -> c
    #             anti_diagonal = set()   # b -> d
    #             for k in range(steps + 1):
    #                 t = k / steps if steps else 0
    #                 i = round(t * (n_cells - 1))
    #                 main_diagonal.add((i, round(t * (m_cells - 1))))
    #                 anti_diagonal.add((i, round((1 - t) * (m_cells - 1))))

    #             new_faces = []
    #             for idx, face in enumerate(faces):
    #                 i, j = divmod(idx, m_cells)
    #                 v0, v1, v2, v3 = face
    #                 if (i, j) in main_diagonal:
    #                     new_faces.append([v0, v1, v2])
    #                     new_faces.append([v0, v2, v3])
    #                 elif (i, j) in anti_diagonal:
    #                     new_faces.append([v1, v2, v3])
    #                     new_faces.append([v1, v3, v0])
    #                 else:
    #                     new_faces.append(face)
    #             faces = new_faces

    #             faces = [[u for u, v in pairwise(face + face[:1]) if u != v] for face in faces]
    #             mesh = PseudoQuadMesh.from_vertices_and_faces_with_face_poles(vertices, faces)
    #             meshes.append(mesh)"""

    #         self.set_quad_mesh(meshes_join_and_weld(meshes))

    #         #NOTE implement pole check part

    #         return self.get_quad_mesh()

    #     elif pattern == "fan":

    #         edge_strip = {}
    #         for strip, edges in self.strips(data=True):
    #             for u, v in edges:
    #                 edge_strip[u, v] = strip
    #                 edge_strip[v, u] = strip

    #         meshes = []
    #         for fkey in self.faces():
    #             face_corners = [self.vertex_coordinates(vkey) for vkey in self.face_vertices(fkey)]
    #             A, B, C, D = face_corners[0], face_corners[1], face_corners[2], face_corners[3]

    #             # the fan needs the same number of divisions on all four sides of the
    #             # patch (a ring hits the two sides meeting at its pole at the same
    #             # parameter), so take one density for the whole patch
    #             densities = [self.get_strip_density(edge_strip[u, v]) for u, v in self.face_halfedges(fkey)]
    #             densities = [8 for u, v in self.face_halfedges(fkey)]
    #             n_rings = max(1, max(densities) // 2)   # rings per quadrant = divisions per HALF side
    #             n_rays = 2 * n_rings                    # rays per quadrant; even, so the diagonal is a ray

    #             centre = (0.5, 0.5)
    #             # (pole, first mid-edge, second mid-edge) per quadrant, ordered so the
    #             # angular sweep runs counter-clockwise around the pole in every quadrant
    #             quadrants = [
    #                 ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5)),
    #                 ((1.0, 0.0), (1.0, 0.5), (0.5, 0.0)),
    #                 ((1.0, 1.0), (0.5, 1.0), (1.0, 0.5)),
    #                 ((0.0, 1.0), (0.0, 0.5), (0.5, 1.0)),
    #             ]

    #             vertices = []
    #             vertex_index = {}

    #             def add_vertex(u, w):
    #                 # one vertex per (u, w) parameter pair, so the four quadrants share
    #                 # their seams, their poles and the centre instead of duplicating them
    #                 key = (round(u, 9), round(w, 9))
    #                 if key not in vertex_index:
    #                     x = (1 - u) * (1 - w) * A[0] + u * (1 - w) * B[0] + u * w * C[0] + (1 - u) * w * D[0]
    #                     y = (1 - u) * (1 - w) * A[1] + u * (1 - w) * B[1] + u * w * C[1] + (1 - u) * w * D[1]
    #                     z = (1 - u) * (1 - w) * A[2] + u * (1 - w) * B[2] + u * w * C[2] + (1 - u) * w * D[2]
    #                     vertex_index[key] = len(vertices)
    #                     vertices.append([x, y, z])
    #                 return vertex_index[key]

    #             faces = []
    #             for pole, m1, m2 in quadrants:
    #                 grid = {}
    #                 for i in range(n_rays + 1):
    #                     # a: angular parameter, walking the L-shaped far boundary
    #                     # m1 -> centre -> m2 at constant speed (half the rays per leg)
    #                     a = float(i) / n_rays
    #                     if a <= 0.5:
    #                         t = a / 0.5
    #                         lu = m1[0] + t * (centre[0] - m1[0])
    #                         lw = m1[1] + t * (centre[1] - m1[1])
    #                     else:
    #                         t = (a - 0.5) / 0.5
    #                         lu = centre[0] + t * (m2[0] - centre[0])
    #                         lw = centre[1] + t * (m2[1] - centre[1])
    #                     for j in range(n_rings + 1):
    #                         # b: radial parameter, sliding from the pole out to the L
    #                         b = float(j) / n_rings
    #                         grid[i, j] = add_vertex(pole[0] + b * (lu - pole[0]),
    #                                                 pole[1] + b * (lw - pole[1]))
    #                 for i in range(n_rays):
    #                     for j in range(n_rings):
    #                         face = [grid[i, j], grid[i + 1, j], grid[i + 1, j + 1], grid[i, j + 1]]
    #                         # the whole j = 0 ring sits on the pole, so the innermost
    #                         # quads have a repeated vertex -- collapse them to triangles
    #                         face = [a_ for a_, b_ in pairwise(face + face[:1]) if a_ != b_]
    #                         faces.append(face)

    #             meshes.append(PseudoQuadMesh.from_vertices_and_faces_with_face_poles(vertices, faces))

    #         self.set_quad_mesh(meshes_join_and_weld(meshes))

    #         return self.get_quad_mesh()


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':
    pass
