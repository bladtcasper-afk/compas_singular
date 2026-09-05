from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from compas.geometry import is_point_in_polygon_xy
from compas.geometry import length_vector
from compas.geometry import subtract_vectors
from compas.geometry import cross_vectors
from compas.geometry import delaunay_triangulation as delaunay_from_points
from compas.geometry import distance_point_point
from compas.geometry import intersection_segment_segment_xy
from compas.datastructures.graph.operations.join import graph_polylines
from compas.datastructures.mesh.operations.weld import mesh_unweld_edges
from compas.itertools import pairwise
from compas.tolerance import TOL

from ..datastructures import Mesh
from ..datastructures import Network
from ..datastructures import trimesh_face_circle


__all__ = [
    'boundary_triangulation',
    'weld_polyline_features',
    'arrange_polyline_features',
    'as_points',
    'as_curves',
]


def as_points(curve, close=None):
    """Coerce one curve to a list of ``[x, y, z]``.

    Accepts a :class:`compas.geometry.Polyline`, anything else exposing
    ``.points``, a list of :class:`compas.geometry.Point`, or the list of
    ``[x, y, z]`` this module has always taken. A list of ``[x, y, z]`` is
    returned unchanged, so existing callers pay nothing and their Delaunay
    tie-breaks cannot move.

    Parameters
    ----------
    curve : Polyline | list
        The curve.
    close : bool, optional
        ``False`` drops a final point coincident with the first, which is the
        convention for a boundary loop. ``True`` adds one, which is the
        convention for a closed curve feature -- the cut has to come back round
        to where it started. ``None`` leaves the ends alone.

    Returns
    -------
    list[list[float]]

    """
    if hasattr(curve, 'points'):
        curve = list(curve.points)

    points = curve
    if any(not isinstance(point, list) for point in curve):
        points = [[float(point[0]), float(point[1]),
                   float(point[2]) if len(point) > 2 else 0.0] for point in curve]

    if close is not None and len(points) > 1:
        closed = TOL.geometric_key(points[0]) == TOL.geometric_key(points[-1])
        if close and not closed:
            points = list(points) + [list(points[0])]
        elif closed and not close:
            points = list(points[:-1])
    return points


def as_curves(curves, close=None):
    """Coerce one curve, or a list of curves, to a list of lists of ``[x, y, z]``.

    A single curve is wrapped in a list, so ``polyline_features=my_polyline``
    means the same as ``polyline_features=[my_polyline]``.

    Parameters
    ----------
    curves : Polyline | list
        One curve or a list of curves, in any of the forms :func:`as_points`
        accepts.
    close : bool, optional
        Passed through to :func:`as_points` for every curve.

    Returns
    -------
    list[list[list[float]]]

    """
    if curves is None:
        return []
    if hasattr(curves, 'points'):                    # a single Polyline
        return [as_points(curves, close)]
    curves = list(curves)
    if not curves:
        return []
    if hasattr(curves[0], 'points'):                 # a list of Polylines
        return [as_points(curve, close) for curve in curves]
    try:                                             # a single curve of points
        float(curves[0][0])
        float(curves[0][1])
        return [as_points(curves, close)]
    except (TypeError, IndexError, ValueError):
        return [as_points(curve, close) for curve in curves]


def weld_polyline_features(polyline_features):
    """Weld feature polylines into chains that share their junctions.

    This is what :meth:`compas_singular.rhino.RhinoSurface.discrete_mapping`
    already does to the curves selected in Rhino, and what a headless caller
    otherwise has to remember to do. Polylines that SHARE A POINT are joined
    into one network and split again at every node that is not two-valent, so
    the junction becomes one vertex of the Delaunay rather than two.

    Without it the topological cut leaks through the junction and the
    decomposition returns one large polygonal face instead of a layout: an X
    fed as two full diagonals comes back with an 11-gon.

    Note what this does NOT do, here or in Rhino: polylines that cross without
    sharing a point are left alone, because no node exists to weld. Supply a
    crossing as the chains that meet at it -- an X as four arms from the centre,
    which is also what four drawn lines give you in Rhino.

    Parameters
    ----------
    polyline_features : list
        List of planar polylines as lists of vertex coordinates.

    Returns
    -------
    list
        The welded polylines. A single polyline, or an empty list, is returned
        unchanged.

    """
    # One chain has no junction to weld, and round-tripping it through the
    # network would renumber its points and move the Delaunay's tie-breaks.
    if len(polyline_features) < 2:
        return polyline_features
    return graph_polylines(Network.from_lines(
        [(u, v) for polyline in polyline_features for u, v in pairwise(polyline)]))


def arrange_polyline_features(polyline_features, tol=1e-9):
    """Put a vertex where two curve features cross.

    A feature is embedded by cutting the Delaunay along it, and a cut can only
    run along edges the triangulation has. Two curves crossing in their interiors
    share no point, so the crossing is no vertex and **no triangulation can hold
    both crossing segments** -- one is simply absent, and the cut leaks straight
    through the crossing.

    That is not a shortcoming of the Delaunay. Chew's constrained Delaunay
    triangulation is defined for "a set of vertices together with a set of
    NONCROSSING edges", so a crossing breaks its precondition too. The fix is a
    planar arrangement, not a stronger triangulator.

    Measured on Fig 4.19, whose two features cross at (4.7, 4.7): exactly one
    segment is missing from the Delaunay, and the two skeleton branches that
    should be pruned instead run through the crossing, passing within 0.026 of
    it. Splitting there takes the layout from 20 patches with 2 triangles and 10
    singularities to 20 all-quad patches and 4.

    Parameters
    ----------
    polyline_features : list
        Planar polylines as lists of vertex coordinates.
    tol : float, optional
        Distance below which a crossing counts as an existing vertex.

    Returns
    -------
    list
        The polylines, with a vertex inserted at every crossing.

    """
    if len(polyline_features) < 2:
        return polyline_features

    chains = [[list(point) for point in polyline] for polyline in polyline_features]
    cuts = [[] for _ in chains]
    for i, one in enumerate(chains):
        for j, two in enumerate(chains):
            if j <= i:
                continue
            for a, (p1, p2) in enumerate(pairwise(one)):
                for b, (p3, p4) in enumerate(pairwise(two)):
                    point = intersection_segment_segment_xy((p1, p2), (p3, p4))
                    if point is None:
                        continue
                    point = [point[0], point[1], 0.0]
                    # a crossing AT a point one of them already has is a shared
                    # node, which weld_polyline_features handles on its own
                    for index, seg, ends in ((i, a, (p1, p2)), (j, b, (p3, p4))):
                        if all(distance_point_point(point, end) > tol for end in ends):
                            cuts[index].append((seg, point))

    if not any(cuts):
        return polyline_features

    out = []
    for chain, chain_cuts in zip(chains, cuts):
        if not chain_cuts:
            out.append(chain)
            continue
        by_segment = {}
        for seg, point in chain_cuts:
            by_segment.setdefault(seg, []).append(point)
        points = [chain[0]]
        for index, (a, b) in enumerate(pairwise(chain)):
            for point in sorted(by_segment.get(index, []),
                                key=lambda q: distance_point_point(a, q)):
                if distance_point_point(points[-1], point) > tol:
                    points.append(point)
            points.append(b)
        out.append(points)
    return out


def boundary_triangulation(outer_boundary, inner_boundaries, polyline_features=[], point_features=[], delaunay=None):
    """Generate Delaunay triangulation between a planar outer boundary and planar inner boundaries. All vertices lie the boundaries.

    Thesis S4.2.1: the surface is triangulated "using points on the boundaries as
    vertices", and faces outside the boundaries are deleted. The boundaries "must
    be discretised into points densely enough to capture the relevant curvature
    changes but not too densely to avoid unnecessary heavy computation" --
    d_i = max(ceil(l_i / (alpha * D)), d_min), with alpha 0.01 to 0.05 of the
    bounding-box diagonal D and d_min 5 to 10. That sampling is the caller's job;
    this function takes the points as given.

    Thesis S4.3.2: curve features are embedded by a TOPOLOGICAL CUT -- "Topological
    cuts are made in the Delaunay mesh along the curve features to consider them as
    boundaries" -- which is the unwelding at the end.

    Not implemented: the constrained Delaunay of Fig 4.20 (Chew 1989). Measured
    harmless at these sampling densities; see ``HOW_IT_WORKS.md``.

    Parameters
    ----------
    outer_boundary : Polyline | list
        Planar outer boundary, as a list of vertex coordinates or a
        :class:`compas.geometry.Polyline`. A final point coincident with the
        first is dropped.
    inner_boundaries : list
        Planar inner boundaries, as a list of the above.
    polyline_features : Polyline | list
        Planar curve features, as a list of the above -- or one curve on its
        own. A closed feature must repeat its first point at the end, and one
        given as a closed ``Polyline`` already does.
    point_features : list
        Planar point features, as a flat list of ``[x, y, z]`` or of
        :class:`compas.geometry.Point`. NOT a list of lists of points.
    delaunay : callable or proxy
        Delaunay triangulation function.

    Returns
    -------
    delaunay_mesh : Mesh
        The Delaunay mesh.

    """
    if not delaunay:
        delaunay = delaunay_from_points

    # accept Polylines and Points as readily as lists of [x, y, z]
    outer_boundary = as_points(outer_boundary, close=False)
    inner_boundaries = as_curves(inner_boundaries, close=False)
    polyline_features = as_curves(polyline_features)
    point_features = [as_points([point])[0] for point in (point_features or [])]

    # a planar arrangement of the features -- a vertex at every crossing --
    # then the chains split at every junction, as the Rhino front end does
    polyline_features = arrange_polyline_features(polyline_features)
    polyline_features = weld_polyline_features(polyline_features)

    # generate planar Delaunay triangulation
    vertices = [pt for boundary in [outer_boundary] + inner_boundaries + polyline_features for pt in boundary] + point_features
    faces = delaunay(vertices)

    delaunay_mesh = Mesh.from_vertices_and_faces(vertices, faces)

    # delete false faces with aligned vertices
    for fkey in list(delaunay_mesh.faces()):
        a, b, c = [delaunay_mesh.vertex_coordinates(vkey) for vkey in delaunay_mesh.face_vertices(fkey)]
        ab = subtract_vectors(b, a)
        ac = subtract_vectors(c, a)
        if length_vector(cross_vectors(ab, ac)) == 0:
            delaunay_mesh.delete_face(fkey)

    # delete faces outisde the borders
    for fkey in list(delaunay_mesh.faces()):
        centre = trimesh_face_circle(delaunay_mesh, fkey)[0]
        if not is_point_in_polygon_xy(centre, outer_boundary) or any([is_point_in_polygon_xy(centre, inner_boundary) for inner_boundary in inner_boundaries]):
            delaunay_mesh.delete_face(fkey)

    # topological cut along the feature polylines through unwelding
    vertex_map = {TOL.geometric_key(delaunay_mesh.vertex_coordinates(vkey)): vkey for vkey in delaunay_mesh.vertices()}
    edges = [edge for polyline in polyline_features for edge in pairwise([vertex_map[TOL.geometric_key(point)] for point in polyline])]
    mesh_unweld_edges(delaunay_mesh, edges)

    # Record the cut. The unwelding above does not separate the first and last
    # segment of a chain -- their end vertices are never split -- so the faces
    # either side of those segments stay adjacent and the skeleton would run
    # straight across the feature there (thesis Fig 4.20a). ``Skeleton``
    # discounts those adjacencies; see ``Skeleton.real_neighbors``.
    # The features as point chains. The grafting needs to know which samples
    # belong to the same curve, and in what order -- see
    # ``SkeletonDecomposition.branches_singularity_to_boundary``.
    delaunay_mesh.attributes['feature_points'] = [
        [list(point) for point in polyline] for polyline in polyline_features]
    delaunay_mesh.attributes['feature_edges'] = frozenset(
        (TOL.geometric_key(a), TOL.geometric_key(b))
        for polyline in polyline_features for a, b in pairwise(polyline))

    return delaunay_mesh


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':
    pass
