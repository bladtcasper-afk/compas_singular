from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from math import floor
# from math import ceil
from math import pi
from operator import itemgetter

from compas.geometry import Polyline
# from compas.geometry import length_vector
# from compas.geometry import length_vector_xy
from compas.geometry import subtract_vectors
from compas.geometry import angle_points
from compas.geometry import angle_vectors
from compas.geometry import angle_vectors_signed
# from compas.geometry import cross_vectors
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas.datastructures.graph.operations.join import graph_polylines
from compas.datastructures.mesh.operations.insert import mesh_insert_vertex_on_edge
from compas.datastructures.mesh.operations.substitute import mesh_substitute_vertex_in_faces
from compas.datastructures.mesh.operations.weld import mesh_unweld_edges
from compas.itertools import pairwise
from compas.itertools import window
from compas.tolerance import TOL

from compas_singular.algorithms import boundary_triangulation

from ..datastructures import CoarsePseudoQuadMesh
from ..datastructures import Network
from ..datastructures import Skeleton
from ..datastructures import mesh_weld
from ..datastructures import split_quad_in_pseudo_quads
from ..datastructures import trimesh_face_circle
from ..geometry import bounding_box_diagonal
from ..geometry import discretise_boundary
from ..geometry import discretise_line
from ..utilities import list_split

from .propagation import quadrangulate_mesh


__all__ = ['SkeletonDecomposition']


def _wall_spacing(outer, inners, target_length, alpha):
    """The segment length the walls were discretised at, for the features.

    ``discretise_boundary`` resolves ``alpha`` against the bounding-box diagonal
    of the loops; the features must be measured against that same D, not against
    their own extent, or a short feature comes back far denser than the wall
    beside it.
    """
    if target_length is not None:
        return float(target_length)
    if alpha is None:
        return None
    diagonal = bounding_box_diagonal(outer, *inners)
    return alpha * diagonal if diagonal > 0.0 else None


class SkeletonDecomposition(Skeleton):
    """Coarse quad mesh from a topological skeleton and its singularities.

    Implements Oval's thesis chapter 4. The correspondence, method by method:

    ==========================================  ==================================
    thesis                                      here
    ==========================================  ==================================
    S4.2.1  end / singular / regular faces      :meth:`corner_faces`, ``singular_faces``
    S4.2.2  PRUNING                             :meth:`branches_singularity_to_singularity`
    S4.2.2  GRAFTING                            :meth:`branches_singularity_to_boundary`
    S4.2.2  CLOSING                             :meth:`branches_boundary`
    S4.2.3.1  missed concavities                :meth:`branches_splitting_boundary_kinks`
    S4.2.3.2  unwanted triangles                :meth:`solve_triangular_faces`
    S4.2.3.3  flipped patches                   :meth:`branches_splitting_flipped_faces`
    S4.2.3.4  collapsed boundaries              :meth:`branches_splitting_collapsed_boundaries`
    S4.3.1  point features become poles         :meth:`split_quads_with_poles`, :meth:`store_pole_data`
    S4.3.2  seam propagation (Fig 4.20d)        :meth:`quadrangulate_polygonal_faces`
    ==========================================  ==================================

    Known deviations from the thesis are listed in ``HOW_IT_WORKS.md``; the
    benchmarks that measure them are ``examples/11_thesis_curve_features.py``.

    """

    def __init__(self, *args, **kwargs):
        super(SkeletonDecomposition, self).__init__(*args, **kwargs)
        self.mesh = None
        self.polylines = None
        self.relative_kink_angle_limit = pi / 8.
        self.flip_angle_limit = pi / 2.
        #: How far each copy of a duplicated vertex moves toward its own
        #: neighbour centroid, as a fraction of that distance. See
        #: :meth:`solve_triangular_faces`.
        self.collapsed_edge_opening = 0.5
        self.repair_notes = []

        self.origin = None #skeleton, mesh or boundary
        self.inputs = {}
        #: The ``poles`` :attr:`mesh` was built with, as a rounded tuple. The
        #: cache key of :meth:`coarse_mesh` -- see there.
        self._mesh_poles = None

    @classmethod
    def from_skeleton(cls, skeleton):
        """Construct a SkeletonDecomposition object from a Skeleton.

        Returns
        -------
        Skeleton
            A skeleton object.

        """
        return cls.from_vertices_and_faces(*skeleton.to_vertices_and_faces())

    @classmethod
    def from_mesh(cls, mesh):
        """Construct a SkeletonDecomposition object from a Mesh.

        Returns
        -------
        Skeleton
            A skeleton object.

        """
        skeleton = cls.from_vertices_and_faces(*mesh.to_vertices_and_faces())
        skeleton.feature_edges = frozenset(mesh.attributes.get('feature_edges') or ())
        skeleton.feature_points = list(mesh.attributes.get('feature_points') or ())
        return skeleton

    @classmethod
    def from_boundary(cls, outer_boundary, inner_boundaries=None, polyline_features=[], point_features=[],
                      target_length=None, alpha=0.04, d_min=5):
        """Triangulate a domain given by its walls, ready to decompose.

        Parameters
        ----------
        outer_boundary : list[[x, y, z]]
            The outer wall, as an open loop -- the last point is not the first.
        inner_boundaries : list[list[[x, y, z]]], optional
            The holes, same convention.
        polyline_features : list[list[[x, y, z]]], optional
            Feature curves the decomposition must follow. See
            ``examples/10_curve_features.py`` for what these do and where they
            still break.
        point_features : list[[x, y, z]], optional
            Points the decomposition must pass through. They become the POLES of
            the layout, and :meth:`coarse_mesh` takes them from here.
        target_length : float, optional
            **Background spacing, not the quad size.** The medial axis is read
            off a Delaunay triangulation of the boundary POINTS, so the walls are
            discretised to this before they are triangulated: a four-point square
            gives a two-triangle Delaunay mesh, which has no interior structure
            at all, and the decomposition that comes out of it is a caricature.
            ``None`` falls back to the thesis rule, ``alpha`` times the
            bounding-box diagonal -- the same default the field route uses, so
            the two front ends discretise the same walls the same way. Pass
            ``alpha=None`` as well to triangulate the loops exactly as handed in.
        alpha : float, optional
            Fraction of the bounding-box diagonal to use as the target length
            when none is given. Thesis eq. 4.1; see
            :func:`~compas_singular.geometry.discretise_boundary`.
        d_min : int, optional
            Fewest points per boundary loop, whatever the target length says.

        Returns
        -------
        SkeletonDecomposition
        """
        # Discretised for the TRIANGULATION only. What is stored below is what
        # the caller handed in, because that is the most faithful description of
        # the walls there is: ``edges_to_curves`` samples them again, finer, and
        # would otherwise be re-sampling this function's own approximation.
        inner_boundaries = list(inner_boundaries or [])
        outer_dense, inners_dense = discretise_boundary(
            outer_boundary, inner_boundaries, target_length=target_length,
            alpha=alpha, d_min=d_min)
        # The features get the SAME sampling as the walls. Thesis eq. 4.1 is
        # stated for every curve, and a feature is cut into the Delaunay along
        # its OWN segments: one longer than the surrounding sampling is not a
        # Delaunay edge, so the cut does not happen there at all. A straight
        # guide drawn in Rhino arrives as TWO points -- ``curve_points`` takes a
        # polyline at its own vertices -- and was triangulated as one long
        # segment the layout ignored completely.
        spacing = _wall_spacing(outer_dense, inners_dense, target_length, alpha)
        features_dense = [discretise_line(feature, spacing)
                          for feature in (polyline_features or [])]
        trimesh = boundary_triangulation(
            outer_boundary=outer_dense, inner_boundaries=inners_dense,
            polyline_features=features_dense, point_features=point_features)

        decomposition = cls.from_vertices_and_faces(*trimesh.to_vertices_and_faces())
        decomposition.feature_edges = frozenset(trimesh.attributes.get('feature_edges') or ())
        decomposition.feature_points = list(trimesh.attributes.get('feature_points') or ())
        decomposition.origin = "boundary"
        # ``inputs``, matching ``__init__``. It used to be written as ``input``,
        # so ``coarse_mesh`` read the empty dict and every point feature was
        # silently dropped instead of becoming a pole.
        # The RESOLVED length is stored, not the argument: ``edges_to_curves``
        # derives its wall sampling from it, and a stored ``None`` would leave
        # the walls unsampled for every caller who took the thesis default.
        if target_length is None and alpha is not None:
            diagonal = bounding_box_diagonal(outer_dense, *inners_dense)
            if diagonal > 0.0:
                target_length = alpha * diagonal
        decomposition.inputs = {"outer_boundary": outer_boundary, "inner_boundaries": inner_boundaries,
                                "polyline_features": polyline_features, "point_features": point_features,
                                "target_length": target_length}

        return decomposition

    # --------------------------------------------------------------------------
    # key elements
    # --------------------------------------------------------------------------

    def corner_faces(self):
        """Get the indices of the corner faces in the Delaunay mesh, i.e. the ones with one neighbour.

        Thesis S4.2.1 calls these END faces: "end faces have one adjacent face".
        Regular faces have two and singular faces three.

        Returns
        -------
        list
            List of face keys.

        """
        return [fkey for fkey in self.faces() if len(self.real_neighbors(fkey)) == 1]

    def corner_vertices(self):
        """Get the indices of the corner vertices of the topological skeleton, i.e. the two-valent boundary vertices in the Delaunay mesh.

        The "two-valent boundary vertices" the CLOSING operation of thesis
        S4.2.2 splits the boundary at.

        Returns
        -------
        list
            List of vertex keys.

        """
        return [vkey for bdry in self.vertices_on_boundaries() for vkey in bdry if len(self.vertex_neighbors(vkey)) == 2]

    def split_vertices(self):
        """Get the indices of the boundary split vertices, i.e. the vertices of the singular faces.

        The other half of what CLOSING splits the boundary at (thesis S4.2.2):
        "the vertices of the singular faces".

        Returns
        -------
        list
            List of vertex keys.

        """
        return [vkey for fkey in self.singular_faces() for vkey in self.face_vertices(fkey)]

    # --------------------------------------------------------------------------
    # branches
    # --------------------------------------------------------------------------

    def branches_singularity_to_singularity(self):
        """Get the branch polylines of the topological skeleton between singularities only, not corners.

        **PRUNING**, the first operation of thesis S4.2.2: "pruning removes the
        branches connected to end faces". Here that is a filter -- a branch whose
        end is an end-face circumcentre is dropped.

        Returns
        -------
        list
            List of polylines as list of point XYZ-coordinates.

        """
        map_corners = [TOL.geometric_key(trimesh_face_circle(self, corner)[0]) for corner in self.corner_faces()]
        return [
            branch for branch in self.branches()
            if TOL.geometric_key(branch[0]) not in map_corners and TOL.geometric_key(branch[-1]) not in map_corners]

    def branches_singularity_to_boundary(self):
        """Get new branch polylines between singularities and boundaries, at the location fo the split vertices. Not part of the topological skeleton.

        **GRAFTING**, the second operation of thesis S4.2.2: "grafting adds
        branches between the singular face circumcentres and their three
        vertices".

        Returns
        -------
        list
            List of polylines as list of point XYZ-coordinates.

        """
        grafts = [[trimesh_face_circle(self, fkey)[0], self.vertex_coordinates(vkey)]
                  for fkey in self.singular_faces() for vkey in self.face_vertices(fkey)]
        return self.merge_graft_targets(grafts)

    def merge_graft_targets(self, grafts):
        """Grafts landing on one feature at ADJACENT samples share a node.

        Two singular faces either side of a curve feature each graft to their own
        nearest sample of it, and the sampling puts those one step apart rather
        than at the same place. The patch spanning them comes out a triangle.

        Measured on Fig 4.18 -- the singular faces at (2.75, 2.75) and
        (6.90, 6.90) graft to (4.77, 4.63) and (4.63, 4.77), 0.2 apart, when both
        project to the segment's midpoint. Merging them takes that layout from 6
        interior singularities to 4 with no change in patch count.

        Adjacency along the chain is the test rather than a distance: the
        discretisation decides what counts as the same place, so there is nothing
        to tune. Targets two samples apart are left alone -- Fig 4.19's
        (3.36, 6.04) and (3.08, 6.32) -- and widening the rule to reach them
        makes that case worse, not better.

        Parameters
        ----------
        grafts : list[[[x, y, z], [x, y, z]]]
            Branches as ``[circumcentre, target]``.

        Returns
        -------
        list[[[x, y, z], [x, y, z]]]

        """
        if not self.feature_points:
            return grafts

        sample_of = {}
        for c, chain in enumerate(self.feature_points):
            for i, sample in enumerate(chain):
                sample_of.setdefault(TOL.geometric_key(sample), (c, i))

        on_feature = {}
        for _, target in grafts:
            key = TOL.geometric_key(target)
            if key in sample_of:
                on_feature[key] = sample_of[key]
        if len(on_feature) < 2:
            return grafts

        # walk each chain in order, grouping runs of samples that touch
        move = {}
        for c in range(len(self.feature_points)):
            hits = sorted(((i, key) for key, (chain, i) in on_feature.items() if chain == c))
            run = []
            for i, key in hits + [(None, None)]:
                if run and (i is None or i - run[-1][0] > 1):
                    if len(run) > 1:
                        shared = self.feature_points[c][run[0][0]]
                        for _, member in run:
                            move[member] = shared
                    run = []
                if i is not None:
                    run.append((i, key))

        if not move:
            return grafts
        return [[centre, move.get(TOL.geometric_key(target), target)] for centre, target in grafts]

    def branches_boundary(self):
        """Get new branch polylines from the Delaunay mesh boundaries split at the corner and plit vertices. Not part of the topological skeleton.

        **CLOSING**, the third operation of thesis S4.2.2: "closing adds boundary
        branches, split at the vertices of the singular faces and the two-valent
        boundary vertices".

        Returns
        -------
        list
            List of polylines as list of point XYZ-coordinates.

        """
        boundaries = [bdry + bdry[0:] for bdry in self.boundaries()]
        splits = self.corner_vertices() + self.split_vertices()
        split_boundaries = [split_boundary for boundary in boundaries for split_boundary in list_split(boundary, [boundary.index(split) for split in splits if split in boundary])]
        return [[self.vertex_coordinates(vkey) for vkey in boundary] for boundary in split_boundaries]

    # --------------------------------------------------------------------------
    # decomposition
    # --------------------------------------------------------------------------

    def decomposition_polylines(self):
        """Get all the branch polylines to form a decomposition of the Delaunay mesh.

        The three operations of thesis S4.2.2 -- pruning, grafting, closing --
        followed by the corrections of S4.2.3.
        These branches include the ones between singularities, between singularities and boundaries and along boundaries.
        Additional branches for some fixes: the ones to further split boundaries that only have two splits.

        Returns
        -------
        list
            List of polylines as list of point XYZ-coordinates.

        """
        branches = self.branches_singularity_to_singularity() + self.branches_singularity_to_boundary() + self.branches_boundary()
        # Thesis S4.2.3 orders the corrections, and says why:
        # concavities -> unwanted triangles -> flipped patches -> collapsed.
        # "Adding branches for other corrections can solve collapsed
        # boundaries. Therefore, correcting collapsed boundaries occurs
        # last." The unwanted-triangle correction is not here at all --
        # solve_triangular_faces runs on the built mesh, after this.
        branches += self.branches_splitting_boundary_kinks()
        branches += self.branches_splitting_flipped_faces()
        branches += self.branches_splitting_collapsed_boundaries()
        self.polylines = graph_polylines(Network.from_lines([(u, v) for polyline in branches for u, v in pairwise(polyline)]),
                                         splits=[self.vertex_coordinates(vkey) for vkey in self.corner_vertices()])
        return self.polylines

    def decomposition_polyline(self, geom_key_1, geom_key_2):
        """Retrieve the decomposition polyline with extremities corresponding to two geoemtric keys.

        Parameters
        ----------
        geom_key_1 : float
            Geometric key of one extremity.
        geom_key_2 : float
            Geometric key of the other extremity.

        Returns
        -------
        list, None
            A polyline as a list of point XYZ-coordinates if a polyline corresponds to the geometric keys, None otherwise.
        """
        polylines = {(TOL.geometric_key(polyline[0]), TOL.geometric_key(polyline[-1])): polyline for polyline in self.polylines}
        return polylines.get((geom_key_1, geom_key_2), polylines.get((geom_key_2, geom_key_1), None))

    def decomposition_mesh(self, poles):
        """Return a quad mesh based on the decomposition polylines.
        Some fixes are added to convert the mesh formed by the decomposition polylines into a (coarse) quad mesh.

        Returns
        -------
        mesh
            A coarse quad mesh based on the topological skeleton from a Delaunay mesh.

        """
        polylines = self.decomposition_polylines()
        boundary_keys = set([TOL.geometric_key(self.vertex_coordinates(vkey)) for bdry in self.vertices_on_boundaries() for vkey in bdry])
        boundary_polylines = [polyline for polyline in polylines if TOL.geometric_key(polyline[0]) in boundary_keys and TOL.geometric_key(polyline[1]) in boundary_keys]
        other_polylines = [polyline for polyline in polylines if TOL.geometric_key(polyline[0]) not in boundary_keys or TOL.geometric_key(polyline[1]) not in boundary_keys]
        self.repair_notes = []
        self.mesh = CoarsePseudoQuadMesh.from_polylines(boundary_polylines, other_polylines)
        self.solve_triangular_faces()
        self.quadrangulate_polygonal_faces()
        self.repair_polygonal_faces(poles)
        self.split_quads_with_poles(poles)
        self.store_pole_data(poles)
        return self.mesh

    def coarse_mesh(self, poles=None, force=False):
        """**The coarse quad layout. The same object every time you ask.**

        The entry point for the workflow::

            decomposition = SkeletonDecomposition.from_boundary(...)
            coarse = decomposition.coarse_mesh()
            coarse.set_strips_density_target(t=0.5)
            dense = coarse.densification()

        :meth:`decomposition_mesh` is the implementation and still works; this
        adds the two things that make it usable as a step of that sequence.

        **The poles come from the domain.** ``point_features`` handed to
        :meth:`from_boundary` are what ``split_quads_with_poles`` and
        ``store_pole_data`` turn into pseudo-quads, so there is no reason to make
        the caller pass them a second time.

        **It is idempotent.** ``decomposition_mesh`` rebuilds from scratch on
        every call, which makes it the one method here you cannot call twice: set
        a strip density on the layout, ask for the layout again, and you get a
        different object with none of your work on it -- no error, nothing to
        see. This returns :attr:`mesh` unchanged whenever it was built from the
        same poles.

        Parameters
        ----------
        poles : list[[x, y, z]], optional
            Preferred pole positions. ``None`` takes them from the
            ``point_features`` of :meth:`from_boundary`, or none at all for a
            decomposition built from a mesh or a skeleton. An explicit list wins.
            Different poles miss the cache.
        force : bool, optional
            Rebuild even on a cache hit -- for a caller that wants a clean layout
            back after mutating the one it was given.

        Returns
        -------
        CoarsePseudoQuadMesh
        """
        if poles is None:
            poles = self.inputs.get("point_features") or []
        poles_key = tuple(tuple(round(float(c), 6) for c in point) for point in poles)
        if not force and self.mesh is not None and self._mesh_poles == poles_key:
            return self.mesh
        mesh = self.decomposition_mesh(poles)
        self._mesh_poles = poles_key
        return mesh

    def edges_to_curves(self, coarse=None, wall_sampling=None, snap=True):
        """**The shape of every coarse edge**, for ``densification``.

        A coarse edge is a straight chord and has to be -- the layout is a
        topological quad graph, and strips, densities and poles are all defined
        on it -- so the SHAPE of each edge lives here instead. Without this a
        curved domain loses the area between chord and wall outright: a round
        hole comes out as a polygon.

        Each edge gets the best of three, and the tally says which: the piece of
        WALL between its two corners, the decomposition BRANCH whose ends match
        it, or a straight CHORD. See
        :func:`~compas_singular.datastructures.coarse_edges_to_curves`.

        The walls come from the domain :meth:`from_boundary` was given, resampled
        finer than the background: those points ARE the boundary from here on, so
        an arc handed over as two points is a chord.

        Parameters
        ----------
        coarse : CoarseQuadMesh, optional
            The layout to describe. Defaults to :attr:`mesh`.
        wall_sampling : float, optional
            How finely to resample the walls. Defaults to a quarter of the
            background ``target_length``, the ratio ``CMD_coarse_mesh`` uses.
        snap : bool, optional
            Put the layout's boundary corners ON the walls first. **This mutates
            the layout, and the order matters**: a corner is placed by the
            background triangulation, so on a curved wall it sits on a CHORD of
            it, and since an arc is pinned to its two ends those corners would be
            the only points of the dense boundary not on the curve -- measured at
            15.1% of the radius on a real plate. The arcs are anchored on corner
            positions, so snapping afterwards is too late.

        Returns
        -------
        (dict, dict)
            ``{(u, v): polyline}`` to hand to ``densification``, and the tally.
            A non-zero ``chord`` count on a boundary edge is a layout that
            quietly lost its curvature.

        Raises
        ------
        ValueError
            If this decomposition was not built by :meth:`from_boundary`, so
            there are no walls to derive the curves from.
        """
        from ..datastructures import coarse_edges_to_curves
        from ..datastructures import snap_corners_to_walls

        if coarse is None:
            coarse = self.mesh
        if coarse is None:
            raise ValueError('no coarse layout -- call coarse_mesh() first')
        if not self.inputs.get('outer_boundary'):
            raise ValueError(
                'edges_to_curves needs the domain walls, which only '
                'from_boundary records. Call '
                'compas_singular.datastructures.coarse_edges_to_curves directly '
                'with the loops you have.')

        if wall_sampling is None:
            target = self.inputs.get('target_length')
            wall_sampling = target * 0.25 if target else None
        # ``alpha``/``d_min`` off: this is output smoothing at an explicitly
        # given sampling, not the thesis discretisation of the input.
        outer_wall, inner_walls = discretise_boundary(
            self.inputs['outer_boundary'],
            self.inputs.get('inner_boundaries') or [],
            target_length=wall_sampling, alpha=None, d_min=None)
        walls = [outer_wall] + inner_walls

        if snap:
            snap_corners_to_walls(coarse, loops=walls)
        return coarse_edges_to_curves(coarse, loops=walls,
                                      polylines=self.polylines or [])

    # --------------------------------------------------------------------------
    # corrections
    # --------------------------------------------------------------------------

    def branches_splitting_collapsed_boundaries(self):
        """Add new branches to fix the problem of boundaries with less than three splits that would be collapsed in the decomposition mesh.

        Thesis S4.2.3.4, COLLAPSED BOUNDARIES: "If less than three branches
        represent a boundary, the boundary collapses in the coarse quad mesh."
        It runs LAST, because "adding branches for other corrections can solve
        collapsed boundaries".

        Returns
        -------
        new_branches : list
            List of polylines as list of point XYZ-coordinates.

        """
        new_branches = []

        all_splits = set(list(self.corner_vertices()) + list(self.split_vertices()))

        for polyedge in [bdry + bdry[:1] for bdry in self.boundaries()]:

            splits = set([vkey for vkey in polyedge if vkey in all_splits])
            new_splits = []

            if len(splits) == 0:
                new_splits += [vkey for vkey in list(itemgetter(0, int(floor(len(polyedge) / 3)), int(floor(len(polyedge) * 2 / 3)))(polyedge))]

            elif len(splits) == 1:
                i = polyedge.index(splits[0])
                new_splits += list(itemgetter(i - int(floor(len(polyedge) * 2 / 3)), i - int(floor(len(polyedge) / 3)))(polyedge))

            elif len(splits) == 2:
                one, two = list_split(polyedge, [polyedge.index(vkey) for vkey in splits])
                half = one if len(one) > len(two) else two
                new_splits.append(half[int(floor(len(half) / 2))])

            for vkey in new_splits:
                fkey = list(self.vertex_faces(vkey))[0]
                for edge in self.face_halfedges(fkey):
                    if vkey in edge and not self.is_edge_on_boundary(*edge):
                        new_branches += [[trimesh_face_circle(self, fkey)[0], self.vertex_coordinates(vkey_2)] for vkey_2 in edge]
                        all_splits.update(edge)
                        break

        return new_branches

    def branches_splitting_flipped_faces(self):
        """Add new branches to fix the problem of polyline patches that would form flipped faces in the decomposition mesh.

        Thesis S4.2.3.3, FLIPPED PATCHES: a patch whose face normal opposes its
        patch normal. The thesis puts this AFTER the unwanted-triangle
        correction, because it "requires a quad decomposition"; here
        :meth:`solve_triangular_faces` runs later still, on the built mesh, so
        this correction sees a decomposition that is not yet all quads.

        Returns
        -------
        new_branches : list
            List of polylines as list of point XYZ-coordinates.

        """
        new_branches = []
        centre_to_fkey = {TOL.geometric_key(trimesh_face_circle(self, fkey)[0]): fkey for fkey in self.faces()}

        # compute total rotation of polyline
        for polyline in self.branches_singularity_to_singularity():
            angles = [angle_vectors_signed(subtract_vectors(v, u), subtract_vectors(w, v), [0., 0., 1.]) for u, v, w in window(polyline, n=3)]
            # subdivide once per angle limit in rotation
            if abs(sum(angles)) > self.flip_angle_limit:
                # the step between subdivision points in polylines (+ 2 for the extremities, which will be discarded)
                alone = len(self.singular_faces()) == 0
                n = floor(abs(sum(angles)) / self.flip_angle_limit) + 1
                step = int(floor(len(polyline) / n))
                # add new branches from corresponding face in Delaunay mesh
                seams = polyline[:: step]
                if polyline[-1] != seams[-1]:
                    if len(seams) == n + 1:
                        del seams[-1]
                    seams.append(polyline[-1])
                if alone:
                    seams = seams[0:-1]
                else:
                    seams = seams[1:-1]
                for point in seams:
                    fkey = centre_to_fkey[TOL.geometric_key(point)]
                    for edge in self.face_halfedges(fkey):
                        if not self.is_edge_on_boundary(*edge):
                            new_branches += [[trimesh_face_circle(self, fkey)[0], self.vertex_coordinates(vkey)] for vkey in edge]
                            break

        return new_branches

    def boundary_interior_angle(self, vkey):
        """The interior angle of the domain at a boundary vertex, in radians.

        The sum of the incident face angles at ``vkey``. Convex is below pi,
        concave above: a square corner reads 90 deg, the reentrant corner of an
        L-plate 270 deg, and a vertex the domain wraps completely around -- the
        end of a curve-feature cut -- 360 deg.

        Measured this way rather than from the turn of the boundary walk, whose
        sign flips between the outer loop and a hole, and is undefined on the
        loop around a cut, which encloses no area.

        """
        total = 0.
        for fkey in self.vertex_faces(vkey):
            face_vertices = self.face_vertices(fkey)
            i = face_vertices.index(vkey)
            total += angle_points(self.vertex_coordinates(vkey),
                                  self.vertex_coordinates(face_vertices[i - 1]),
                                  self.vertex_coordinates(face_vertices[(i + 1) % len(face_vertices)]))
        return total

    def branches_splitting_boundary_kinks(self):
        """Add new branches to fix the problem of boundary kinks not marked by the skeleton

        Thesis S4.2.3.1, MISSED CONCAVITIES: "The skeleton marks convex but not
        concave kinks, and the added branches may not spot the concavities.
        Therefore, branches are added to include them." Only a CONCAVITY is
        corrected -- see the gate on :meth:`boundary_interior_angle` below.
        Due to a low density that did not spot the change of curvature at the kink.
        Does not modify the singularites on the contrarty to increasing the density.

        Returns
        -------
        new_branches : list
            List of polylines as list of point XYZ-coordinates.

        """
        new_branches = []

        singular_faces = set(self.singular_faces())
        for boundary in self.boundaries():
            angles = {(u, v, w): angle_vectors(subtract_vectors(self.vertex_coordinates(v), self.vertex_coordinates(u)), subtract_vectors(
                self.vertex_coordinates(w), self.vertex_coordinates(v))) for u, v, w in window(boundary + boundary[: 2], n=3)}
            for u, v, w, x, y in list(window(boundary + boundary[: 4], n=5)):

                # check if not a corner
                if self.vertex_degree(w) == 2:
                    continue

                # Thesis S4.2.3.1: "The skeleton marks convex but not concave
                # kinks". Only a CONCAVITY needs a branch adding; a convex kink
                # is already carried by the skeleton, and correcting one puts a
                # branch point a fraction of the discretisation off the corner.
                #
                # The original test is unsigned, so it cannot tell the two apart.
                # It shows up as soon as a curve feature lands on a wall: the
                # corner stops being two-valent, slips past the guard above, and
                # is corrected as though it were a concavity. Measured on
                # Fig 4.17, a diagonal corner to corner -- both of the square's
                # own corners, interior angle 90 deg, were being corrected.
                if self.boundary_interior_angle(w) <= pi:
                    continue

                angle = angles[(v, w, x)]
                adjacent_angles = (angles[(u, v, w)] + angles[(w, x, y)]) / 2

                if angle - adjacent_angles > self.relative_kink_angle_limit:
                    # check if not already marked via an adjacent singular face
                    if all([fkey not in singular_faces for fkey in self.vertex_faces(w)]):
                        fkeys = list(self.vertex_faces(w, ordered=True))
                        fkey = fkeys[int(floor(len(fkeys) / 2))]
                        for edge in self.face_halfedges(fkey):
                            if w in edge and not self.is_edge_on_boundary(*edge):
                                new_branches += [[trimesh_face_circle(self, fkey)[0], self.vertex_coordinates(vkey)] for vkey in edge]
                                break

        return new_branches

    def solve_triangular_faces(self):
        """Modify the decomposition mesh from polylines to make it a quad mesh by converting the degenerated quad faces that appear as triangular faces.

        Thesis S4.2.3.2, UNWANTED TRIANGLES: "If two adjacent singular faces have
        one or several coinciding vertices or a coinciding circumcentre, the
        missing branch yields triangular faces in the coarse mesh."

        **Deviation.** The thesis inserts the missing BRANCH, before the mesh is
        built: "If two of the three patch corners are on the boundary, the branch
        is inserted at the other corner, and reciprocally if two of the three
        patch corners are off the boundary." This does vertex surgery on the
        built mesh instead -- duplicate a corner (case 1), merge two (case 2).
        The two consequences: it cannot act on a triangle whose corners are all
        interior, and it never fires at a curve feature, because
        ``from_polylines`` has already welded the topological cut shut so no
        corner there answers ``is_vertex_on_boundary``. That is why a free curve
        extremity stays a pole instead of collapsing to the two-valent
        singularity of Fig 4.22b.
        """
        mesh = self.mesh

        for fkey in list(mesh.faces()):
            if len(mesh.face_vertices(fkey)) == 3:

                boundary_vertices = [vkey for vkey in mesh.face_vertices(fkey) if mesh.is_vertex_on_boundary(vkey)]
                case = sum(mesh.is_vertex_on_boundary(vkey) for vkey in mesh.face_vertices(fkey))

                if case == 1:
                    # convert triangular face to quad by duplicating the boundary vertex
                    # due to singular face vertices at the same location
                    u = boundary_vertices[0]
                    v = mesh.add_vertex(attr_dict={attr: xyz for attr, xyz in zip(['x', 'y', 'z'], mesh.vertex_coordinates(u))})

                    # modify adjacent faces
                    vertex_faces = mesh.vertex_faces(u, ordered=True)
                    mesh_substitute_vertex_in_faces(mesh, u, v, vertex_faces[: vertex_faces.index(fkey)])

                    # modify triangular face
                    mesh_insert_vertex_on_edge(mesh, (u, mesh.face_vertex_ancestor(fkey, u)), v)

                elif case == 2:
                    # remove triangular face and merge the two boundary vertices
                    # due to singularities at the same location
                    #
                    # Thesis S4.2.3.2 does the opposite here -- "If two of the
                    # three patch corners are on the boundary, the branch is
                    # inserted at the OTHER corner" -- i.e. split the interior
                    # corner, the mirror of case 1. That was implemented and
                    # measured: across a square, an L-plate, a hole, a point
                    # feature and all three curve-feature benchmarks, **case 2
                    # never fires once** (case 0 fires 12 times, case 1 four
                    # times), so the two rules cannot be told apart. Left as
                    # Robin wrote it rather than shipping an unexercised branch.
                    #
                    # The gap that matters is case 0 -- all three corners
                    # interior -- which neither this code nor S4.2.3.2 covers,
                    # and which is what a free curve extremity produces.
                    polyline = Polyline(self.decomposition_polyline(*map(lambda x: TOL.geometric_key(mesh.vertex_coordinates(x)), boundary_vertices)))
                    point = polyline.point_at(t=.5, snap=True)
                    new_vkey = mesh.add_vertex(attr_dict={'x': point.x, 'y': point.y, 'z': point.z})

                    # modify triangular face
                    mesh.delete_face(fkey)

                    # modify adjacent faces
                    for old_vkey in boundary_vertices:
                        mesh_substitute_vertex_in_faces(mesh, old_vkey, new_vkey, mesh.vertex_faces(old_vkey))
                        mesh.delete_vertex(old_vkey)

        to_move = {}
        # Give some length to the new edge. Case 1 duplicates a vertex, which
        # leaves a zero-length edge, and densification divides by it -- without
        # this the free-tip case raises ZeroDivisionError.
        #
        # Each copy moves toward its OWN neighbour centroid, so the pair opens up
        # and the step scales with the local patch size. The fraction was 0.1,
        # which is too small to produce an edge: measured on a 10-unit plate with
        # a free guide it left a 0.14 edge, and on a real Rhino plate a 0.019 one
        # -- aspect ratio 376, a 178.6 degree corner, and a singularity that
        # reads as a tiny edge rather than a point. At 0.5 the same case opens to
        # 0.70 with the layout unchanged.
        for edge in mesh.edges():
            threshold = 1e-6
            if mesh.edge_length(*edge) < threshold:
                for vkey in edge:
                    xyz = centroid_points([mesh.vertex_coordinates(nbr) for nbr in mesh.vertex_neighbors(vkey)])
                    xyz0 = mesh.vertex_coordinates(vkey)
                    to_move[vkey] = [self.collapsed_edge_opening * (a - a0)
                                     for a, a0 in zip(xyz, xyz0)]

        for vkey, xyz in to_move.items():
            attr = mesh.vertex[vkey]
            attr['x'] += xyz[0]
            attr['y'] += xyz[1]
            attr['z'] += xyz[2]

    def quadrangulate_polygonal_faces(self):
        """Turn the polygonal faces left by a curve feature into quad faces.

        Thesis S4.3.2, Fig 4.20d: "Pentagonal or higher-valency faces ... become
        quad faces by propagating the seams of the discrepancies on the curve
        features."

        A curve feature is a topological cut (see
        :func:`~compas_singular.algorithms.boundary_triangulation`), and the two
        sides of a cut do not receive the same skeleton branches. Where one side
        gets a branch landing that the other does not, the two sides disagree,
        and welding them back together in ``from_polylines`` turns that
        disagreement into a face with five or more vertices.

        Such a face is repaired by propagating the seam: the odd vertex is
        carried across the patch to the opposite side and the patch is relaid as
        a discrete Coons patch of quads. This is the step of Fig. 4.20d in Oval's
        thesis, and :func:`~compas_singular.algorithms.propagation.quadrangulate_mesh`
        is its implementation.

        Notes
        -----
        Faces of three vertices are NOT touched here. A triangle is a legitimate
        pseudo-quad with a pole, and is the business of
        :meth:`solve_triangular_faces` and :meth:`store_pole_data`.

        This method was disabled between 2021 and 2026. See ``HOW_IT_WORKS.md``;
        the short version is that it kept one un-migrated call to COMPAS's
        ``vertices_on_boundary()`` after that returned loops rather than
        vertices, and that it welded inside the ``exploded()`` loop and so threw
        away every component but the last.

        """
        supermesh = self.mesh

        # Seam propagation only has something to do where a face is bigger than
        # a quad. Returning early matters: mesh_weld below REBUILDS the mesh and
        # renumbers its keys, and the editing layers address faces by key.
        if not any(len(supermesh.face_vertices(fkey)) > 4 for fkey in supermesh.faces()):
            return

        delaunay_vertex_map = tuple(TOL.geometric_key(self.vertex_coordinates(vkey)) for vkey in self.vertices())
        # newly added vertices in mesh that were not in the Delaunay are missing...

        edges_to_unweld = [edge for edge in supermesh.edges() if sum([TOL.geometric_key(supermesh.vertex_coordinates(i)) in delaunay_vertex_map for i in edge]) == 2]
        mesh_unweld_edges(supermesh, edges_to_unweld)

        # A "discrepancy" is a position whose coincident copies, one per side of
        # the cut, do not agree on their boundary valency: a branch landed on one
        # side only. Collect them over every component BEFORE welding -- welding
        # inside this loop is what used to reduce self.mesh to the last piece.
        source_map = []
        for mesh in supermesh.exploded():
            candidate_map = {TOL.geometric_key(mesh.vertex_coordinates(vkey)): [] for vkey in mesh.vertices()}
            for boundary in mesh.vertices_on_boundaries():
                for vkey in boundary:
                    candidate_map[TOL.geometric_key(mesh.vertex_coordinates(vkey))].append(mesh.vertex_degree(vkey))
            source_map += [geom_key for geom_key, valencies in candidate_map.items() if len(list(set(valencies))) > 1]
        source_map = tuple(source_map)

        self.mesh = mesh_weld(supermesh)
        mesh = self.mesh

        sources = [vkey for vkey in mesh.vertices() if TOL.geometric_key(mesh.vertex_coordinates(vkey)) in source_map]

        # A seam vertex is a source for the polygonal face on one side of the cut
        # and a genuine corner of the quad on the other. One global source list
        # is still safe, because quadrangulate_mesh only rewrites faces whose
        # length is not 4 and so never reaches the quad.
        quadrangulate_mesh(mesh, sources)

    def quadrangulate_polygonal_faces_wip(self):
        pass
        # mesh = self.mesh

        # delaunay_vertex_map = tuple(TOL.geometric_key(self.vertex_coordinates(vkey)) for vkey in self.vertices())

        # for fkey in mesh.faces():
        # 	face_vertices = mesh.face_vertices(fkey)
        # 	if len(face_vertices) > 4:

    def repair_polygonal_faces(self, poles=()):
        """Last resort for a face seam propagation could not turn into quads.

        :meth:`quadrangulate_polygonal_faces` repairs a polygonal face by
        carrying its seam across to the opposite side, which needs the face to
        have exactly four real corners. A face with more than one seam on the
        same side -- a closed curve feature produces them -- has fewer, and comes
        through untouched. ``collect_strips`` cannot walk such a face, so rather
        than leave the whole densification to fail on it, hand it to the general
        repair the frame-field front end already carries: it splits an n-gon
        along its own diagonals, adding no vertex and so leaving every neighbour
        conforming.

        Inert unless a face larger than a quad survived, which is also what keeps
        the import out of the common path.

        Parameters
        ----------
        poles : list, optional
            Point features, passed on as preferred pole positions.

        """
        mesh = self.mesh
        if not any(len(mesh.face_vertices(fkey)) > 4 for fkey in mesh.faces()):
            return

        try:
            from ..framefield.repair import solve_non_quad_faces
        except Exception as exc:      # the frame-field extras are optional
            self.repair_notes.append(
                'polygonal faces left and the fallback repair is unavailable ({})'.format(exc))
            return

        self.mesh, note = solve_non_quad_faces(mesh, cls=type(mesh), poles=list(poles))
        if note:
            self.repair_notes.append('fallback repair: {}'.format(note))

    def split_quads_with_poles(self, poles):
        new_lines = []

        mesh = self.mesh
        pole_map = tuple([TOL.geometric_key(pole) for pole in poles])

        faces = list(mesh.faces())
        for fkey in faces:
            fv = mesh.face_vertices(fkey)
            if len(fv) == 4:
                for vkey in fv:
                    if TOL.geometric_key(mesh.vertex_coordinates(vkey)) in pole_map:
                        idx = fv.index(vkey)
                        xkey = fv[idx + 2 - len(fv)]
                        new_lines.append([mesh.vertex_coordinates(vkey), mesh.vertex_coordinates(xkey)])
                        split_quad_in_pseudo_quads(mesh, fkey, vkey)
                        break

        self.polylines += new_lines
        return new_lines

    def store_pole_data(self, poles):
        """Record, for every triangular face, which of its corners is the pole.

        A triangular face of the coarse mesh is not a defect: it is the
        pseudo-quad ``(p, a, b, p)`` with one side collapsed at the pole ``p``,
        and ``face_opposite_edge``, ``collect_strips`` and ``densification`` all
        handle it. What IS a defect is a triangle with no pole recorded --
        ``face_opposite_edge`` then raises ``KeyError`` on the first strip that
        walks into it.

        A point feature is the obvious source of a pole, but not the only one. A
        curve feature whose extremity lies off the boundary is a singularity too
        (Oval, thesis p. 100: "a pole if adjacent to several singular faces"), and
        so is a medial-axis degeneracy. Those triangles used to print
        ``pole missing`` and take the whole densification down with them, so any
        triangle without a point feature at a corner gets one chosen for it.

        Parameters
        ----------
        poles : list
            Point features, as XYZ coordinates. A triangle with one of these at a
            corner collapses there, exactly as before.

        """
        mesh = self.mesh
        pole_map = tuple([TOL.geometric_key(pole) for pole in poles])

        face_poles = {}
        for fkey in mesh.faces():
            if len(mesh.face_vertices(fkey)) == 3:
                for vkey in mesh.face_vertices(fkey):
                    if TOL.geometric_key(mesh.vertex_coordinates(vkey)) in pole_map:
                        face_poles[fkey] = vkey
                        break
                if fkey not in face_poles:
                    face_poles[fkey] = self._choose_pole(fkey)
                    self.repair_notes.append(
                        'face {} is a triangle with no point feature at a corner; '
                        'pole set to vertex {}'.format(fkey, face_poles[fkey]))

        mesh.attributes['face_pole'] = face_poles

    def _choose_pole(self, fkey):
        """Which corner of a triangular face to collapse the pseudo-quad at.

        A pseudo-quad ``(p, a, b)`` is the quad ``(p, a, b, p)``: ``(p, a)`` and
        ``(b, p)`` are one strip, and ``(a, b)`` is the side opposite the
        collapsed corner. Densification fans the face from ``p`` onto ``(a, b)``
        with the two rails at one density, so the fan is even when the two rails
        are the same length -- pick the corner whose two edges are closest to
        equal.

        Same rule, and same reasoning, as ``_choose_pole`` in
        ``compas_singular.framefield.repair``.

        """
        mesh = self.mesh
        fv = mesh.face_vertices(fkey)
        best = None
        for i in range(3):
            p, a, b = fv[i], fv[(i + 1) % 3], fv[(i + 2) % 3]
            ra = distance_point_point(mesh.vertex_coordinates(p), mesh.vertex_coordinates(a))
            rb = distance_point_point(mesh.vertex_coordinates(b), mesh.vertex_coordinates(p))
            if ra + rb <= 0.0:
                continue
            imbalance = abs(ra - rb) / (ra + rb)
            if best is None or imbalance < best[0]:
                best = (imbalance, p)
        return best[1] if best else fv[0]

# ==============================================================================
# Main
# ==============================================================================


if __name__ == '__main__':
    pass
