"""Step 1 -- a triangulation a field can actually live on.

``compas_singular``'s own ``boundary_triangulation`` puts vertices ONLY on the
boundaries (see ``algorithms/triangulation.py``); the interior is spanned by a
few large, often sliver Delaunay triangles whose circumcentres approximate the
medial axis. That is exactly right for the skeleton front end and useless for a
field: there is no interior resolution to carry one, and slivers wreck the
Laplacian.

So the two front ends want opposite triangulations from the same boundary, and
this module builds the other one: boundary densified to ``target_length``, plus
an interior grid at the same spacing.

Planar domain, so every tangent space is world XY -- no parallel transport, no
per-edge connection. ``face_basis`` returns the identity and exists only so the
surface case does not need an API change.
"""
from math import ceil

from compas_singular.datastructures import Mesh   # has .boundaries(); compas core does not
from compas_singular.geometry.polyline import bounding_box_diagonal
from compas_singular.geometry.polyline import discretise_boundary
from compas.geometry import delaunay_triangulation
from compas.geometry import is_point_in_polygon_xy
from compas.geometry import distance_point_point
from compas.geometry import normalize_vector
from compas.geometry import subtract_vectors
from compas.geometry import cross_vectors
from compas.geometry import length_vector
from compas.itertools import pairwise


# ``discretise_boundary`` is re-exported, not reimplemented: the skeleton front
# end discretises its walls with the SAME function, so the two routes cannot
# drift apart. See ``geometry/polyline.py``.
__all__ = ['BackgroundMesh', 'discretise_boundary', 'interior_grid']


def _as_open_loop(points):
    """Drop a repeated closing point and any consecutive duplicates."""
    pts = [[float(p[0]), float(p[1]), 0.0] for p in points]
    out = [pts[0]]
    for p in pts[1:]:
        if distance_point_point(p, out[-1]) > 1e-9:
            out.append(p)
    if len(out) > 1 and distance_point_point(out[0], out[-1]) < 1e-9:
        out.pop()
    return out


def interior_grid(outer, inners=(), target_length=None, margin=0.45, symmetry=None):
    """The INTERIOR points of the background mesh, at ``target_length`` spacing.

    The half of this module the skeleton front end has no use for. Its
    triangulation puts vertices only on the boundaries, which is right for
    reading a medial axis off the circumcentres and useless for a field: there
    is no interior resolution to carry one. These are the points that give it
    one. The boundary half comes from :func:`discretise_boundary`, shared with
    the skeleton route.

    Parameters
    ----------
    outer : list[[x, y, z]]
        The domain outline, open and already discretised. Its bounding box is
        the extent the grid is laid over.
    inners : list[list[[x, y, z]]], optional
        The holes, same convention.
    target_length : float
        Grid spacing. The same number the boundary was discretised with, so the
        interior and the wall meet at one element size.
    margin : float, optional
        Points closer to a boundary than ``margin * target_length`` are dropped,
        so the grid does not crowd the wall and produce slivers.
    symmetry : :class:`framefield.symmetry.Symmetry`, optional
        When given and non-trivial, the grid below is replaced by one that is
        EXACTLY invariant under the group. **This is where a symmetric domain's
        symmetry is won or lost.** The grid below is asymmetric in two
        independent ways -- ``_jitter`` hashes the grid indices, and the grid is
        anchored at ``min(xs)`` rather than centred -- and neither is
        recoverable downstream: measured on a symmetric square, the interior
        vertices are 0% invariant under a quarter turn with a worst deviation of
        0.65 at a 0.5 spacing, which is more than a whole triangle. See
        ``symmetry.py``.

    Returns
    -------
    list[[x, y, z]]
    """
    inners = list(inners or [])
    limit = margin * target_length

    if symmetry is not None and symmetry.enabled('background'):
        # Imported here rather than at module scope: symmetry.py imports
        # _jitter and _distance_to_loop from this module, so a top-level
        # import either way round is a cycle.
        from .symmetry import interior_points

        return list(interior_points(symmetry, target_length, outer, inners,
                                    margin=margin))

    xs = [p[0] for p in outer]
    ys = [p[1] for p in outer]
    points = []
    nx = int(ceil((max(xs) - min(xs)) / target_length))
    ny = int(ceil((max(ys) - min(ys)) / target_length))
    for i in range(nx + 1):
        for j in range(ny + 1):
            p = [min(xs) + (i + 0.5) * target_length + _jitter(i, j, 0.2 * target_length),
                 min(ys) + (j + 0.5) * target_length + _jitter(j, i, 0.2 * target_length),
                 0.0]
            if not is_point_in_polygon_xy(p, outer):
                continue
            if _distance_to_loop(p, outer) < limit:
                continue
            if any(is_point_in_polygon_xy(p, loop) for loop in inners):
                continue
            if any(_distance_to_loop(p, loop) < limit for loop in inners):
                continue
            points.append(p)
    return points


def _distance_to_loop(p, loop):
    """Shortest distance from a point to a closed polyline, in XY."""
    best = float('inf')
    for a, b in pairwise(loop + loop[:1]):
        ab = subtract_vectors(b, a)
        length2 = ab[0] * ab[0] + ab[1] * ab[1]
        if length2 == 0.0:
            continue
        t = ((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1]) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + ab[0] * t, a[1] + ab[1] * t, 0.0]
        best = min(best, distance_point_point(p, q))
    return best


def _jitter(i, j, amount):
    """Deterministic sub-cell offset.

    A perfectly regular grid makes every square's four corners cocircular, so
    Qhull picks a diagonal arbitrarily and the triangulation is not reproducible
    across runs or platforms. A fixed pseudo-random nudge removes the degeneracy
    without introducing randomness.
    """
    h = (i * 73856093) ^ (j * 19349663)
    return (((h >> 8) & 0xFFFF) / 65535.0 - 0.5) * amount


class BackgroundMesh(object):
    """A well-shaped triangulation of a planar domain, with boundary tangents.

    Attributes
    ----------
    mesh : :class:`compas.datastructures.Mesh`
        The triangulation. All faces have positive area in XY.
    outer : list[[x, y, z]]
        The densified outer boundary, as an open loop (no repeated last point).
    inners : list[list[[x, y, z]]]
        The densified inner boundaries (holes), same convention.
    target_length : float
        The spacing the triangulation was built at.
    """

    def __init__(self, mesh, outer, inners, target_length):
        self.mesh = mesh
        self.outer = outer
        self.inners = inners
        self.target_length = target_length
        self._boundary_tangents = None

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------

    @classmethod
    def from_boundary(cls, outer_boundary, inner_boundaries=None, target_length=None,
                      margin=0.45, symmetry=None, alpha=0.04, d_min=5):
        """Triangulate the region inside ``outer_boundary`` and outside the inners.

        Parameters
        ----------
        outer_boundary : list[[x, y, z]]
            The domain outline, closed or open (a repeated last point is dropped).
        inner_boundaries : list[list[[x, y, z]]], optional
            Holes.
        target_length : float, optional
            Edge-length target, for the boundary and the interior grid alike.
            Defaults to ``alpha`` times the bounding-box diagonal -- thesis
            eq. 4.1, see :func:`discretise_boundary` -- which at the default
            alpha gives a few thousand triangles, enough for a smooth field
            without making the solve slow.
        margin : float, optional
            Interior grid points closer to a boundary than ``margin *
            target_length`` are dropped, so the grid does not crowd the wall and
            produce slivers.
        symmetry : :class:`framefield.symmetry.Symmetry`, optional
            Passed to :func:`interior_grid`, which is where a symmetric domain's
            symmetry is won or lost -- see there.
        alpha : float, optional
            Fraction of the bounding-box diagonal to use as ``target_length``
            when none is given.
        d_min : int, optional
            Fewest points per boundary loop, whatever the target length says.

        Returns
        -------
        BackgroundMesh
        """
        # The target length is resolved HERE rather than left to
        # ``discretise_boundary``, because the interior grid has to be laid out
        # at the same spacing as the wall and the number is stored on the
        # instance. ``_as_open_loop`` first so the scale is read off the real
        # points; ``discretise_boundary`` cleans the loops again, idempotently.
        if target_length is None:
            loops = [_as_open_loop(outer_boundary)]
            loops += [_as_open_loop(loop) for loop in (inner_boundaries or [])]
            target_length = alpha * bounding_box_diagonal(*loops)

        outer, inners = discretise_boundary(outer_boundary, inner_boundaries,
                                            target_length=target_length,
                                            d_min=d_min)

        points = list(outer)
        for loop in inners:
            points.extend(loop)
        points.extend(interior_grid(outer, inners, target_length,
                                    margin=margin, symmetry=symmetry))

        faces = delaunay_triangulation(points)
        mesh = Mesh.from_vertices_and_faces(points, faces)

        # drop zero-area faces, then faces whose centroid is outside the domain.
        # centroid rather than circumcentre: a sliver's circumcentre can land far
        # outside a face that is perfectly inside the domain.
        for fkey in list(mesh.faces()):
            a, b, c = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
            if length_vector(cross_vectors(subtract_vectors(b, a), subtract_vectors(c, a))) < 1e-12:
                mesh.delete_face(fkey)
                continue
            centre = mesh.face_centroid(fkey)
            if not is_point_in_polygon_xy(centre, outer):
                mesh.delete_face(fkey)
            elif any(is_point_in_polygon_xy(centre, loop) for loop in inners):
                mesh.delete_face(fkey)

        for vkey in list(mesh.vertices()):
            if not mesh.vertex_faces(vkey):
                mesh.delete_vertex(vkey)

        return cls(mesh, outer, inners, target_length)

    # ------------------------------------------------------------------
    # tangent spaces
    # ------------------------------------------------------------------

    def face_basis(self, fkey):
        """Orthonormal tangent basis of a face.

        Planar domain, so this is the world XY basis for every face. Kept in the
        API because the surface case needs it to be per-face, and a caller
        written against it now will not have to change.
        """
        return ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])

    def boundary_tangents(self):
        """The adjacent boundary EDGE directions at every boundary vertex.

        Two per vertex (one each side), not one averaged tangent.

        Taking the chord between the two neighbours instead -- the obvious
        implementation -- is wrong at a corner, and wrong in the worst possible
        way. At a square's 90-degree corner the chord bisects, giving 45 degrees:
        exactly halfway between the two arms of the cross, the one direction
        equally far from both. The field then has to unwind that spurious 45
        degrees somewhere, and a square comes out with two interior
        singularities instead of the none it should have.

        The two edge directions are the honest input; ``constraints.from_boundary``
        combines them in the 4th-power representation, where a right-angle corner's
        two edges are the SAME cross and cancel no information at all.

        Returns
        -------
        dict[int, list[[x, y, z]]]
        """
        if self._boundary_tangents is None:
            tangents = {}
            for loop in self.mesh.boundaries():
                n = len(loop)
                for i, vkey in enumerate(loop):
                    here = self.mesh.vertex_coordinates(vkey)
                    prev = self.mesh.vertex_coordinates(loop[(i - 1) % n])
                    nxt = self.mesh.vertex_coordinates(loop[(i + 1) % n])
                    tangents[vkey] = [
                        normalize_vector(subtract_vectors(here, prev)),
                        normalize_vector(subtract_vectors(nxt, here)),
                    ]
            self._boundary_tangents = tangents
        return self._boundary_tangents

    def boundary_vertices(self):
        """Set of vertex keys on any boundary."""
        return set(self.boundary_tangents())

    # ------------------------------------------------------------------
    # validation
    # ------------------------------------------------------------------

    def validate(self):
        """Check the invariants Step 2 relies on.

        Returns
        -------
        dict
            ``ok`` plus the individual counts, so a caller can print them.

        Notes
        -----
        Euler characteristic is checked against ``1 - len(holes)`` for a planar
        domain, which catches a triangulation that has silently torn or kept a
        face bridging a hole.
        """
        mesh = self.mesh
        report = {
            'vertices': mesh.number_of_vertices(),
            'faces': mesh.number_of_faces(),
            'edges': mesh.number_of_edges(),
        }

        negative = []
        for fkey in mesh.faces():
            a, b, c = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
            area2 = ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
            if area2 <= 0.0:
                negative.append(fkey)
        report['negative_area_faces'] = len(negative)

        report['triangles_only'] = all(len(mesh.face_vertices(f)) == 3 for f in mesh.faces())

        loops = mesh.boundaries()
        report['boundary_loops'] = len(loops)
        report['expected_loops'] = 1 + len(self.inners)

        chi = report['vertices'] - report['edges'] + report['faces']
        report['euler'] = chi
        report['expected_euler'] = 1 - len(self.inners)

        report['ok'] = (
            report['negative_area_faces'] == 0
            and report['triangles_only']
            and report['boundary_loops'] == report['expected_loops']
            and report['euler'] == report['expected_euler']
        )
        return report
