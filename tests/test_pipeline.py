"""Headless regression tests for the core compas_singular pipelines.

These exercise the decomposition / densification / grammar logic against the
current COMPAS 2.x API (no visualisation, so they run anywhere).
"""
import os
import json

import pytest

from math import pi

from compas.tolerance import TOL

from compas_singular import blocks
from compas_singular.datastructures import CoarseQuadMesh
from compas_singular.datastructures import CoarsePseudoQuadMesh
from compas_singular.datastructures import QuadMesh
from compas_singular.datastructures.lizard import Lizard
from compas_singular.algorithms import boundary_triangulation
from compas_singular.algorithms import SkeletonDecomposition

HERE = os.path.dirname(__file__)
DATA = os.path.abspath(os.path.join(HERE, '..', 'examples', 'data'))


def test_coarse_quad_densification():
    mesh = CoarseQuadMesh.from_json(os.path.join(DATA, 'coarse_quad_mesh_british_museum.json'))
    assert mesh.number_of_faces() == 12
    mesh.collect_strips()
    mesh.set_strips_density(3)
    mesh.densification()
    assert mesh.get_quad_mesh().number_of_faces() > mesh.number_of_faces()


def test_pseudo_quad_densification_preserves_poles():
    mesh = CoarsePseudoQuadMesh.from_json(os.path.join(DATA, 'coarse_quad_mesh_british_museum_poles.json'))
    # the pole map must survive JSON (de)serialisation with integer face keys
    assert all(isinstance(fkey, int) for fkey in mesh.attributes['face_pole'])
    assert sorted(mesh.poles()) == [0, 3, 16, 19]
    mesh.collect_strips()
    mesh.set_strips_density_target(t=.5)
    mesh.densification()
    assert mesh.get_quad_mesh().number_of_faces() > 0


def test_skeleton_decomposition_planar():
    with open(os.path.join(DATA, '01_decomposition.json')) as f:
        outer_boundary, inner_boundaries, polyline_features, point_features = json.load(f)
    trimesh = boundary_triangulation(outer_boundary, inner_boundaries, polyline_features, point_features)
    decomposition = SkeletonDecomposition.from_mesh(trimesh)
    coarsemesh = decomposition.decomposition_mesh(point_features)
    coarsemesh.collect_strips()
    coarsemesh.set_strips_density_target(0.5)
    coarsemesh.densification()
    densemesh = coarsemesh.get_quad_mesh()
    assert densemesh.number_of_faces() > coarsemesh.number_of_faces()


def _ngon(n, radius, phase=0.0):
    from math import cos, sin
    return [[radius * cos(phase + 2 * pi * i / n),
             radius * sin(phase + 2 * pi * i / n), 0.0] for i in range(n)]


@pytest.mark.parametrize('radii,spacing', [((4.0, 2.0), 0.5), ((5.0, 2.0), 0.25)])
def test_polygonal_walls_do_not_divide_by_zero(radii, spacing):
    """A polygonal domain triangulates without a ZeroDivisionError.

    A straight wall discretised at ``spacing`` becomes a RUN of collinear
    points, and Qhull answers those with flat simplices. A flat triangle has no
    circumcircle, and ``trimesh_face_circle`` divides by its doubled squared
    area -- so ``boundary_triangulation`` must delete every one of them before
    the next loop asks for a circumcentre.

    The rounding is the whole difficulty, which is why these two shapes are
    pinned rather than any polygon: a 12-gon around a hexagon puts three flat
    faces in the Delaunay, and an exact ``== 0`` area test caught two of them.
    The third's area came out 0.0 from one pair of its edge vectors and 8e-17
    from another, so it survived the filter and raised in ``trimesh_face_circle``.
    """
    outer = _ngon(12, radii[0], phase=pi / 12)
    inner = _ngon(6, radii[1], phase=pi / 12)
    decomposition = SkeletonDecomposition.from_boundary(
        outer, inner_boundaries=[inner], target_length=spacing)
    coarse = decomposition.coarse_mesh()
    assert coarse.number_of_faces() > 0
    assert coarse.is_manifold()


def _l_plate_dense(density=3, resolution=6):
    """An L-plate, whose reentrant corner gives the layout one valence-5 joint.

    Deliberately not a square with a point feature: that comes out a POLE, and
    a pole carries index 4, which no single face can carry.
    """
    from compas.geometry import Line, Point

    corners = [Point(0, 0, 0), Point(10, 0, 0), Point(10, 4, 0),
               Point(4, 4, 0), Point(4, 9, 0), Point(0, 9, 0)]
    corners = corners + [corners[0]]
    outer = []
    for i in range(len(corners) - 1):
        outer.extend(Line(corners[i], corners[i + 1]).to_polyline(n=resolution).points)
    outer = [[p.x, p.y, p.z] for p in outer]

    trimesh = boundary_triangulation(outer, [], [], [])
    coarse = SkeletonDecomposition.from_mesh(trimesh).decomposition_mesh([])
    coarse.collect_strips()
    coarse.set_strips_density(density)
    coarse.densification()
    return coarse.get_quad_mesh()


def test_block_from_point_conserves_index():
    """Truncating a joint TRADES it for a face -- it must not invent one."""
    mesh = _l_plate_dense()
    joints = blocks.joints(mesh)
    assert len(joints) == 1, 'the L-plate should give exactly one joint'
    valence = mesh.vertex_degree(joints[0])
    x, y, z = mesh.vertex_coordinates(joints[0])

    before = blocks.index_sum(mesh)
    degrees_before = blocks.face_degrees(mesh)
    assert set(degrees_before) == {4}, 'the primal should be all quads'

    # click NEAR the joint, not exactly on it -- the point is the handle
    out, report = blocks.block_points(mesh, [[x + 0.12, y - 0.09, z]], ratio=0.45)

    assert report['built'] == 1, report['blocks'][0]['reason']
    assert blocks.index_sum(out) == before
    assert out.is_manifold()

    degrees_after = blocks.face_degrees(out)
    assert degrees_after.get(valence) == 1, 'one block, of the joint\'s degree'
    assert set(degrees_after) == {4, valence}
    assert report['blocks'][0]['degree'] == valence
    # the block landed on the point, within one element
    assert report['blocks'][0]['offset'] < 1.0


def test_block_points_refuses_a_pole_with_a_reason():
    """A pole cannot be truncated, and the refusal has to say so."""
    mesh = CoarsePseudoQuadMesh.from_json(
        os.path.join(DATA, 'coarse_quad_mesh_british_museum_poles.json'))
    mesh.collect_strips()
    mesh.set_strips_density(2)
    mesh.densification()
    dense = mesh.get_quad_mesh()

    pole = dense.vertex_coordinates(
        [v for v in dense.vertices()
         if any(len(dense.face_vertices(f)) != 4
                for f in dense.vertex_faces(v, ordered=True))
         and not dense.is_vertex_on_boundary(v)][0])

    out, report = blocks.block_points(dense, [pole])
    assert report['built'] == 0
    assert 'pole' in report['blocks'][0]['reason'].lower()
    assert out is dense, 'a refusal must not modify the mesh'


def _square_centre_pole_dense(density=2, size=10.0, resolution=20):
    """A square with one centre point feature -- a FULL pole, clear of the wall."""
    from compas.geometry import Line, Point

    h = size / 2.0
    corners = [Point(h, h, 0), Point(-h, h, 0), Point(-h, -h, 0), Point(h, -h, 0)]
    corners = corners + [corners[0]]
    outer = []
    for i in range(len(corners) - 1):
        outer.extend(Line(corners[i], corners[i + 1]).to_polyline(n=resolution).points)
    outer = [[p.x, p.y, p.z] for p in outer]

    point = [0.0, 0.0, 0.0]
    trimesh = boundary_triangulation(outer, [], [], [point])
    coarse = SkeletonDecomposition.from_mesh(trimesh).decomposition_mesh([point])
    coarse.collect_strips()
    coarse.set_strips_density(density)
    coarse.densification()
    return coarse.get_quad_mesh(), point


def test_pole_blocks_collapses_the_fan_and_conserves_index():
    """The other route: the fan becomes one face, and the ring pays for it."""
    dense, point = _square_centre_pole_dense(density=2)

    pole = [v for v in dense.vertices()
            if not dense.is_vertex_on_boundary(v)
            and any(len(dense.face_vertices(f)) != 4
                    for f in dense.vertex_faces(v, ordered=True))
            and dense.vertex_coordinates(v)[:2] == [0.0, 0.0]][0]
    fan = len(dense.vertex_faces(pole))
    # a FULL pole: every face in the fan is a triangle, and none of the fan
    # touches the wall, which is what makes the index conserve exactly
    assert all(len(dense.face_vertices(f)) == 3 for f in dense.vertex_faces(pole))

    before = blocks.index_sum(dense)
    faces_before = dense.number_of_faces()
    joints_before = len(blocks.joints(dense))
    out, report = blocks.pole_blocks(dense, [point])

    assert report['built'] == 1, report['blocks'][0]['reason']
    assert blocks.index_sum(out) == before
    assert out.is_manifold()
    # the whole fan became exactly one face, and nothing else moved
    assert out.number_of_faces() == faces_before - fan + 1
    assert report['blocks'][0]['degree'] == fan
    # the block sits exactly on the picked point, because the point WAS the pole
    assert report['blocks'][0]['offset'] < 1e-6
    assert report['strips'] == 0, 'the collapse must not propagate'
    # and the price: the ring vertices each drop to valence 3
    assert report['ring_joints'] == len(blocks.joints(out)) - joints_before


def test_pole_blocks_reports_index_loss_at_the_wall():
    """A fan touching the boundary still collapses -- and says what it cost.

    ``index_sum`` counts interior vertices only, so a ring vertex ON the wall
    loses a valence without its term being counted and the total drops. That is
    the bookkeeping, not a broken collapse: the result must still be manifold,
    and the drop must be visible in the report rather than silent.
    """
    mesh = CoarsePseudoQuadMesh.from_json(
        os.path.join(DATA, 'coarse_quad_mesh_british_museum_poles.json'))
    mesh.collect_strips()
    mesh.set_strips_density(2)
    mesh.densification()
    dense = mesh.get_quad_mesh()

    poles = [v for v in dense.vertices()
             if not dense.is_vertex_on_boundary(v)
             and any(len(dense.face_vertices(f)) != 4
                     for f in dense.vertex_faces(v, ordered=True))]
    assert poles, 'the fixture should have an interior pole'
    pole = poles[0]
    # this fixture's poles are PARTIAL -- quads in the fan as well as triangles
    assert any(len(dense.face_vertices(f)) == 4 for f in dense.vertex_faces(pole))
    fan = len(dense.vertex_faces(pole))

    out, report = blocks.pole_blocks(dense, [dense.vertex_coordinates(pole)])

    assert report['built'] == 1, report['blocks'][0]['reason']
    assert out.is_manifold()
    assert out.number_of_faces() == dense.number_of_faces() - fan + 1
    # the block is closed on the fan's boundary CYCLE, which is longer than the
    # neighbour ring when there are quads in the fan
    assert report['blocks'][0]['degree'] > dense.vertex_degree(pole)
    assert report['index_after'] <= report['index_before']
    assert any('INDEX CHANGED' in line for line in blocks.report_lines(report))


@pytest.mark.parametrize('string', ['ata', 'atta', 'attta', 'attpptta'])
def test_lizard_grammar_produces_manifold_mesh(string):
    vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0],
                [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0],
                [0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0]]
    faces = [[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]]
    mesh = QuadMesh.from_vertices_and_faces(vertices, faces)
    mesh.collect_strips()
    lizard = Lizard(mesh)
    lizard.initiate()
    lizard.from_string_to_rules(string)
    assert mesh.is_manifold()
    assert all(len(mesh.halfedge[vkey]) > 0 for vkey in mesh.vertices())


# ==============================================================================
# Curve features (Oval, thesis 4.3.2, Figs 4.17-4.22)
#
# Restored in 2026 after five years switched off; see HOW_IT_WORKS.md. The first
# test is the tripwire that keeps the restoration from moving anything else:
# quadrangulate_polygonal_faces welds, and welding renumbers keys.
# ==============================================================================

def _square(side=10.0, spacing=0.2):
    from compas.geometry import distance_point_point

    def sample(a, b):
        n = max(1, int(round(distance_point_point(a + [0.0], b + [0.0]) / spacing)))
        return [[a[0] + i / n * (b[0] - a[0]), a[1] + i / n * (b[1] - a[1]), 0.0] for i in range(n)]

    c = [[0, 0], [side, 0], [side, side], [0, side]]
    return [p for i in range(4) for p in sample(c[i], c[(i + 1) % 4])]


def _segment(a, b, spacing=0.2):
    from compas.geometry import distance_point_point
    n = max(1, int(round(distance_point_point(a, b) / spacing)))
    return [[a[0] + i / n * (b[0] - a[0]), a[1] + i / n * (b[1] - a[1]), 0.0] for i in range(n + 1)]


def _coarse(features=(), points=()):
    trimesh = boundary_triangulation(_square(), [], list(features), list(points))
    return SkeletonDecomposition.from_mesh(trimesh).decomposition_mesh(list(points))


def _face_sizes(mesh):
    sizes = {}
    for fkey in mesh.faces():
        n = len(mesh.face_vertices(fkey))
        sizes[n] = sizes.get(n, 0) + 1
    return sizes


@pytest.mark.parametrize('features, points, faces, vertices, sizes', [
    ((), (), 4, 9, {4: 4}),
    ((), ([5.0, 5.0, 0.0],), 12, 17, {3: 4, 4: 8}),
    ((), ([3.0, 3.0, 0.0], [7.0, 7.0, 0.0]), 14, 18, {3: 6, 4: 8}),
    # re-measured 2026-09-15: was 10 faces / 21 vertices, a layout that dropped
    # both free tips, covered 95.9% of the square and left two corner vertices
    # belonging to no face. Each tip is now a two-valent vertex (Fig 4.21c).
    ((_segment([1.4, 1.4, 0.0], [8.6, 8.6, 0.0]),), (), 12, 19, {4: 12}),
])
def test_curve_feature_repair_leaves_working_inputs_alone(features, points, faces, vertices, sizes):
    """Inputs that already worked must come out untouched, key for key.

    quadrangulate_polygonal_faces calls mesh_weld, which REBUILDS the mesh and
    renumbers every key -- and the editing and agent layers address faces by key.
    Its early return is what keeps that off the common path; this is the guard on
    the early return. The numbers are the pre-restoration output, except where
    noted.
    """
    coarse = _coarse(features, points)
    assert coarse.number_of_faces() == faces
    assert coarse.number_of_vertices() == vertices
    assert _face_sizes(coarse) == sizes


def test_curve_feature_extremity_on_boundary_densifies():
    """Fig 4.17: a curve feature with both extremities on the boundary."""
    coarse = _coarse([_segment([0.0, 5.0, 0.0], [10.0, 5.0, 0.0])])
    assert max(_face_sizes(coarse)) <= 4, 'seam propagation left a polygonal face'
    coarse.collect_strips()
    coarse.set_strips_density_target(0.5)
    coarse.densification()
    assert coarse.get_quad_mesh().number_of_faces() > coarse.number_of_faces()


def test_curve_feature_free_extremities_densifies():
    """Fig 4.18: a curve feature with both extremities off the boundary."""
    coarse = _coarse([_segment([2.0, 5.0, 0.0], [8.0, 5.0, 0.0])])
    assert max(_face_sizes(coarse)) <= 4
    coarse.collect_strips()
    coarse.set_strips_density_target(0.5)
    coarse.densification()
    assert coarse.get_quad_mesh().number_of_faces() > 0


def test_multiple_curve_features_densify():
    """Fig 4.19: several curve features, meeting at a shared node."""
    arms = [_segment([5.0, 5.0, 0.0], [x, y, 0.0])
            for x, y in ((0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0))]
    coarse = _coarse(arms)
    assert max(_face_sizes(coarse)) <= 4
    coarse.collect_strips()
    coarse.set_strips_density_target(0.5)
    coarse.densification()
    assert coarse.get_quad_mesh().number_of_faces() > 0


def test_closed_curve_feature_densifies():
    """A closed curve feature needs the fallback repair -- seam propagation
    cannot reach a face carrying more than one seam on the same side."""
    import math
    n = 80
    loop = [[5.0 + 2.5 * math.cos(2 * math.pi * i / n), 5.0 + 2.5 * math.sin(2 * math.pi * i / n), 0.0]
            for i in range(n)]
    coarse = _coarse([loop + loop[:1]])
    assert max(_face_sizes(coarse)) <= 4
    coarse.collect_strips()
    coarse.set_strips_density_target(0.5)
    coarse.densification()
    assert coarse.get_quad_mesh().number_of_faces() > 0


def test_every_triangular_coarse_face_carries_a_pole():
    """A triangle is a pseudo-quad, but only if its pole is recorded.

    An unrecorded one used to print 'pole missing' and take collect_strips down
    with a KeyError, which is the whole crash the curve-feature path died on.
    """
    coarse = _coarse([_segment([2.0, 5.0, 0.0], [8.0, 5.0, 0.0])])
    face_pole = coarse.attributes['face_pole']
    for fkey in coarse.faces():
        if len(coarse.face_vertices(fkey)) == 3:
            assert fkey in face_pole
            assert face_pole[fkey] in coarse.face_vertices(fkey)


def test_crossing_features_are_welded_at_their_junction():
    """weld_polyline_features is what discrete_mapping already does in Rhino."""
    from compas_singular.algorithms import weld_polyline_features

    arms = [_segment([5.0, 5.0, 0.0], [x, y, 0.0]) for x, y in ((0.0, 5.0), (10.0, 5.0))]
    welded = weld_polyline_features(arms)
    assert len(welded) == 1, 'two arms meeting end to end should join into one chain'
    # a single feature, or none, is passed through untouched
    assert weld_polyline_features([arms[0]]) == [arms[0]]
    assert weld_polyline_features([]) == []


def test_guide_feature_mesh_embeds_the_guide():
    """A guide as a curve feature lands ON the mesh edges by construction."""
    from compas.geometry import closest_point_on_segment, distance_point_point
    from compas_singular.guide_lines import guide_feature_mesh

    guide = [[2.0, 5.0, 0.0], [8.0, 5.0, 0.0]]
    target = 0.5
    dense, info = guide_feature_mesh(_square(), [guide], target_length=target)
    assert dense.number_of_faces() > 0

    edges = [(dense.vertex_coordinates(u), dense.vertex_coordinates(v)) for u, v in dense.edges()]
    offsets = []
    for i in range(41):
        point = [2.0 + 6.0 * i / 40.0, 5.0, 0.0]
        offsets.append(min(distance_point_point(point, closest_point_on_segment(point, e)) for e in edges))

    # Away from the extremities the guide IS a run of mesh edges, exactly -- that
    # is what the topological cut buys, and it is not a snapping tolerance.
    interior = offsets[4:-4]
    assert max(interior) < 1e-6, 'guide leaves the mesh edges by {}'.format(max(interior))

    # At a free extremity it can drift, because the medial axis does not always
    # put a branch point on the tip. Fig 4.22's unwanted-triangle operation is
    # what closes this and is not implemented; see HOW_IT_WORKS.md.
    assert max(offsets) < 0.1 * target, 'extremity drifts by {}'.format(max(offsets))


def test_boundary_triangulation_accepts_polylines_and_points():
    """A Polyline, a list of Point, or a list of [x, y, z] all mean the same."""
    from compas.geometry import Point, Polyline
    from compas_singular.algorithms import as_curves, as_points

    diagonal = _segment([1.4, 1.4, 0.0], [8.6, 8.6, 0.0])
    reference = _face_sizes(_coarse([diagonal]))

    polyline = Polyline([Point(*p) for p in diagonal])
    for form in (diagonal,                       # a bare list of points
                 [diagonal],                     # the historical form
                 polyline,                       # a bare Polyline
                 [polyline],                     # one Polyline in a list
                 [[Point(*p) for p in diagonal]]):
        assert _face_sizes(_coarse(form)) == reference

    # a Polyline outer boundary, and a Point feature
    trimesh = boundary_triangulation(Polyline([Point(*p) for p in _square()]), [], [], [Point(5, 5, 0)])
    coarse = SkeletonDecomposition.from_mesh(trimesh).decomposition_mesh([[5.0, 5.0, 0.0]])
    assert _face_sizes(coarse) == {3: 4, 4: 8}

    # a list of [x, y, z] is handed back as-is, so nothing downstream can shift
    assert as_points(diagonal) is diagonal

    # closed conventions differ by role: a boundary drops the repeat, a feature keeps it
    ring = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
    assert len(as_curves([ring], close=False)[0]) == 3
    assert len(as_curves([ring[:-1]], close=True)[0]) == 4


def _plate(corners, spacing=0.2):
    """A closed polygon, sampled, as boundary_triangulation wants it."""
    pts = []
    for i in range(len(corners)):
        a, b = corners[i], corners[(i + 1) % len(corners)]
        n = max(1, int(round(((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2) ** 0.5 / spacing)))
        pts += [[a[0] + k / n * (b[0] - a[0]), a[1] + k / n * (b[1] - a[1]), 0.0] for k in range(n)]
    return pts


def test_boundary_interior_angle_separates_convex_from_concave():
    """Thesis 4.2.3.1: only CONCAVE kinks get a correction branch.

    The turn of the boundary walk cannot be used for this -- its sign flips
    between the outer loop and a hole, and it is undefined on the loop around a
    curve-feature cut, which encloses no area. The interior angle can.
    """
    from math import degrees

    square = SkeletonDecomposition.from_mesh(
        boundary_triangulation(_plate([[0, 0], [10, 0], [10, 10], [0, 10]]), [], [], []))
    ell = SkeletonDecomposition.from_mesh(
        boundary_triangulation(_plate([[0, 0], [10, 0], [10, 4], [4, 4], [4, 9], [0, 9]]), [], [], []))

    def angle_at(decomposition, xy):
        for vkey in decomposition.vertices():
            x, y, _ = decomposition.vertex_coordinates(vkey)
            if abs(x - xy[0]) < 1e-6 and abs(y - xy[1]) < 1e-6:
                return degrees(decomposition.boundary_interior_angle(vkey))
        raise AssertionError('no vertex at %s' % (xy,))

    assert angle_at(square, (0, 0)) == pytest.approx(90.0, abs=1e-6)     # convex
    assert angle_at(ell, (10, 4)) == pytest.approx(90.0, abs=1e-6)       # convex
    assert angle_at(ell, (4, 4)) == pytest.approx(270.0, abs=1e-6)       # concave


def test_convex_corner_carrying_a_curve_feature_is_not_corrected():
    """Fig 4.17: a curve extremity landing on a wall must not add a singularity.

    The cut runs through the corner the feature lands on, so the corner becomes
    two copies, one per side, each an ordinary two-valent convex corner of its
    half. Before the cut reached the wall the corner was one THREE-valent vertex:
    it slipped past the corner test, the end of the uncut slit next to it read as
    a concavity, and the correction put a spurious pole beside the corner.
    """
    outer = _plate([[0, 0], [10, 0], [10, 10], [0, 10]])
    diagonal = _segment([0.0, 0.0, 0.0], [10.0, 10.0, 0.0])

    decomposition = SkeletonDecomposition.from_mesh(
        boundary_triangulation(outer, [], [diagonal], []))

    # the diagonal cuts the square into two halves
    assert len(decomposition.vertices_on_boundaries()) == 2

    for corner in ([0.0, 0.0], [10.0, 10.0]):
        copies = [v for v in decomposition.vertices() if decomposition.vertex_neighbors(v)
                  and abs(decomposition.vertex_coordinates(v)[0] - corner[0]) < 1e-9
                  and abs(decomposition.vertex_coordinates(v)[1] - corner[1]) < 1e-9]
        assert len(copies) == 2, 'the cut did not reach the corner %s' % corner
        for vkey in copies:
            assert vkey in decomposition.corner_vertices()
            assert decomposition.boundary_interior_angle(vkey) <= pi

    # no slit end is left, so there is no concavity to correct
    assert decomposition.branches_splitting_boundary_kinks() == []


def test_no_skeleton_branch_crosses_a_curve_feature():
    """Thesis 4.20a vs 4.20b: the cut exists so no branch crosses a feature.

    The faces either side of a cut edge are no longer adjacent, so the skeleton
    cannot run across it. ``Skeleton.real_neighbors`` discounts any adjacency
    across a feature edge that survives anyway.
    """
    outer = _plate([[0, 0], [10, 0], [10, 10], [0, 10]])
    diagonal = _segment([0.0, 0.0, 0.0], [10.0, 10.0, 0.0])
    trimesh = boundary_triangulation(outer, [], [diagonal], [])
    decomposition = SkeletonDecomposition.from_mesh(trimesh)

    assert decomposition.feature_edges, 'the cut was not recorded'

    # every corner of the square is an END face, including the two the feature
    # lands on -- a triangle has three, and the cut makes two triangles
    assert len(decomposition.corner_faces()) == 6

    # so every branch runs to an end face and PRUNING removes all of them
    assert len(decomposition.branches()) == 6
    assert decomposition.branches_singularity_to_singularity() == []

    # which is the thesis figure: three quads per half, one singularity each
    coarse = decomposition.decomposition_mesh([])
    assert coarse.number_of_faces() == 6
    assert all(len(coarse.face_vertices(f)) == 4 for f in coarse.faces())
    interior = [v for v in coarse.vertices()
                if not coarse.is_vertex_on_boundary(v) and len(coarse.vertex_neighbors(v)) != 4]
    assert len(interior) == 2


def test_real_neighbors_is_inert_without_curve_features():
    """No features -> no discounted adjacency -> every layout bit-identical."""
    outer = _plate([[0, 0], [10, 0], [10, 10], [0, 10]])
    decomposition = SkeletonDecomposition.from_mesh(boundary_triangulation(outer, [], [], []))
    assert decomposition.feature_edges == frozenset()
    for fkey in decomposition.faces():
        assert decomposition.real_neighbors(fkey) == decomposition.face_neighbors(fkey)


def test_grafts_at_adjacent_samples_of_a_feature_share_a_node():
    """Thesis 4.2.2 grafting: two singular faces either side of a feature each
    graft to their own nearest sample, and the sampling puts those one step
    apart. The patch between them comes out a triangle.

    Fig 4.18: the singular faces at (2.75, 2.75) and (6.90, 6.90) grafted to
    (4.77, 4.63) and (4.63, 4.77) -- 0.2 apart -- when both project to the
    segment's midpoint.
    """
    outer = _plate([[0, 0], [10, 0], [10, 10], [0, 10]])
    segment = _segment([2.8, 6.6, 0.0], [6.6, 2.8, 0.0])
    decomposition = SkeletonDecomposition.from_mesh(
        boundary_triangulation(outer, [], [segment], []))

    assert decomposition.feature_points, 'the feature chains were not recorded'

    targets = {TOL.geometric_key(t) for _, t in decomposition.branches_singularity_to_boundary()}
    on_feature = {k for k in targets
                  if any(TOL.geometric_key(s) == k for s in decomposition.feature_points[0])}
    # four singular faces, each grafting once to the segment; the two at the
    # midpoint now share a node, the two at the tips stay distinct
    assert len(on_feature) == 3

    coarse = decomposition.decomposition_mesh([])
    assert all(len(coarse.face_vertices(f)) == 4 for f in coarse.faces())
    interior = [v for v in coarse.vertices()
                if not coarse.is_vertex_on_boundary(v) and len(coarse.vertex_neighbors(v)) != 4]
    # the four singular points, and the two free tips as two-valent vertices
    assert len(interior) == 6


def test_graft_merge_only_joins_ADJACENT_samples():
    """Two samples apart is a different place, and merging there makes it worse.

    Guards the rule against being widened into a distance threshold: on 4.19 the
    pair (3.36, 6.04) / (3.08, 6.32) is two steps apart and must survive.
    """
    outer = _plate([[0, 0], [10, 0], [10, 10], [0, 10]])
    chain = _segment([0.0, 0.0, 0.0], [4.0, 0.0, 0.0], spacing=1.0)
    decomposition = SkeletonDecomposition.from_mesh(boundary_triangulation(outer, [], [], []))
    decomposition.feature_points = [chain]

    centre = [0.0, 5.0, 0.0]
    two_apart = [[centre, list(chain[0])], [centre, list(chain[2])]]
    assert decomposition.merge_graft_targets(two_apart) == two_apart

    adjacent = [[centre, list(chain[0])], [centre, list(chain[1])]]
    merged = decomposition.merge_graft_targets(adjacent)
    assert merged[0][1] == merged[1][1] == chain[0]


def test_graft_merge_keeps_a_chain_extremity():
    """Merged grafts land on the chain's END, whichever way the chain runs.

    Chain order would pick the sample nearer the start. Where the chain runs
    towards a free tip that is one sample inward, and the tip is left with no
    branch -- both free arms of Fig 4.19 lost theirs.
    """
    outer = _plate([[0, 0], [10, 0], [10, 10], [0, 10]])
    chain = _segment([0.0, 0.0, 0.0], [4.0, 0.0, 0.0], spacing=1.0)
    decomposition = SkeletonDecomposition.from_mesh(boundary_triangulation(outer, [], [], []))
    decomposition.feature_points = [chain]

    centre = [0.0, 5.0, 0.0]
    at_the_end = [[centre, list(chain[-2])], [centre, list(chain[-1])]]
    merged = decomposition.merge_graft_targets(at_the_end)
    assert merged[0][1] == merged[1][1] == chain[-1]


def test_a_feature_landing_on_walls_gives_an_all_quad_layout():
    """``examples/000_testing.py``: a rectangle cut corner to corner by a line.

    Each half is a triangle -- one singularity, three quads -- and each half's
    graft on the line is a discrepancy on the other side. Seam propagation
    (Fig 4.20d) carries each across its half to the wall. It never ran: the
    sources were collected per component, so none was found, and a single global
    list made each pentagon see two sources and three corners.
    """
    rectangle = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [5.0, 10.0, 0.0], [0.0, 10.0, 0.0]]
    decomposition = SkeletonDecomposition.from_boundary(rectangle, polyline_features=[[[0.0, 0.0, 0.0], [5.0, 10.0, 0.0]]])
    coarse = decomposition.coarse_mesh()

    assert _face_sizes(coarse) == {4: 10}
    assert decomposition.repair_notes == []
    for corner in rectangle:
        assert any(coarse.vertex_coordinates(v) == corner for v in coarse.vertices()), 'corner %s lost' % corner


def test_free_feature_extremities_stay_layout_vertices():
    """Fig 4.21: a free extremity is a singularity of the layout.

    Both sides of the slit join into one segment, so a grafted tip has only two
    branches and was merged straight through -- the tip disappeared and the
    coarse edges crossed the feature.
    """
    coarse = _coarse([_segment([2.8, 6.6, 0.0], [6.6, 2.8, 0.0])])
    for tip in ([2.8, 6.6], [6.6, 2.8]):
        vkeys = [v for v in coarse.vertices()
                 if abs(coarse.vertex_coordinates(v)[0] - tip[0]) < 1e-6
                 and abs(coarse.vertex_coordinates(v)[1] - tip[1]) < 1e-6]
        assert len(vkeys) == 1, 'tip %s is not a layout vertex' % tip
        assert len(coarse.vertex_neighbors(vkeys[0])) == 2


def test_crossing_seam_strips_leave_no_polygon():
    """A line cutting off a corner: two seam strips cross in the corner patch.

    A vertex made by the first strip is a real corner once its face is split;
    counted as a source for the second strip, the face showed two sources and
    three corners and was left a pentagon for the fallback repair.
    """
    square = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]]
    decomposition = SkeletonDecomposition.from_boundary(square, polyline_features=[[[0.0, 5.41, 0.0], [5.03, 0.0, 0.0]]])
    coarse = decomposition.coarse_mesh()

    assert _face_sizes(coarse) == {4: 17}
    assert decomposition.repair_notes == []


def test_seam_propagation_that_never_ends_is_abandoned():
    """Around two crossing features a seam strip can close on itself.

    Unguarded, this layout grew past 1000 faces. It is now restored and left to
    the fallback repair.
    """
    rectangle = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [5.0, 10.0, 0.0], [0.0, 10.0, 0.0]]
    crossing = [[[4.052, 5.821, 0.0], [1.391, 9.395, 0.0]], [[3.654, 9.331, 0.0], [2.072, 6.471, 0.0]]]
    decomposition = SkeletonDecomposition.from_boundary(rectangle, polyline_features=crossing, target_length=0.2236)
    coarse = decomposition.coarse_mesh()

    assert any('did not terminate' in note for note in decomposition.repair_notes)
    assert coarse.number_of_faces() < 200
    assert max(_face_sizes(coarse)) <= 4


def test_seam_propagation_skips_a_face_through_a_vertex_twice():
    """A bent line near a corner produced a face passing through one vertex twice,
    and ``quadrangulate_face`` raised ``ValueError: not enough values to unpack``.
    """
    rectangle = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [5.0, 10.0, 0.0], [0.0, 10.0, 0.0]]
    bent = [[[5.0, 2.589, 0.0], [3.505, 1.12, 0.0], [3.824, 0.0, 0.0]]]
    decomposition = SkeletonDecomposition.from_boundary(rectangle, polyline_features=bent)
    coarse = decomposition.coarse_mesh()
    assert max(_face_sizes(coarse)) <= 4


def test_crossing_features_get_a_vertex_at_the_crossing():
    """Two features crossing in their interiors share no point, so the crossing
    is no vertex -- and no triangulation can hold both crossing segments. The cut
    leaks straight through, and the branches that should be pruned run through it
    instead.

    Chew's CDT cannot help: it is defined for a set of NONCROSSING edges. The fix
    is a planar arrangement.
    """
    from compas_singular.algorithms import arrange_polyline_features

    one = _segment([1.0, 1.0, 0.0], [9.0, 9.0, 0.0])
    two = _segment([9.0, 1.0, 0.0], [1.0, 9.0, 0.0])
    assert not {TOL.geometric_key(p) for p in one} & {TOL.geometric_key(p) for p in two}

    arranged = arrange_polyline_features([one, two])
    crossing = TOL.geometric_key([5.0, 5.0, 0.0])
    assert all(any(TOL.geometric_key(p) == crossing for p in chain) for chain in arranged)

    # a single feature, and features that do not cross, are returned untouched
    assert arrange_polyline_features([one]) == [one]
    apart = [_segment([0.0, 1.0, 0.0], [4.0, 1.0, 0.0]), _segment([0.0, 3.0, 0.0], [4.0, 3.0, 0.0])]
    assert arrange_polyline_features(apart) == apart


def test_two_crossing_features_decompose_to_all_quads():
    """Fig 4.19. Before the arrangement its two branches ran within 0.026 of the
    crossing instead of being pruned, and the layout carried two triangles.
    """
    outer = _plate([[0, 0], [10, 0], [10, 10], [0, 10]])
    diagonal = _segment([0.0, 0.0, 0.0], [10.0, 10.0, 0.0])
    crossing = _segment([2.8, 6.6, 0.0], [6.6, 2.8, 0.0])

    decomposition = SkeletonDecomposition.from_mesh(
        boundary_triangulation(outer, [], [diagonal, crossing], []))
    coarse = decomposition.decomposition_mesh([])

    assert all(len(coarse.face_vertices(f)) == 4 for f in coarse.faces())
    interior = [v for v in coarse.vertices()
                if not coarse.is_vertex_on_boundary(v) and len(coarse.vertex_neighbors(v)) != 4]
    assert len(interior) == 4


def test_a_two_point_guide_is_discretised_like_the_walls():
    """A straight polyline drawn in Rhino arrives as TWO points.

    ``curve_points`` takes a polyline at its own vertices, and ``from_boundary``
    resampled the walls but handed the features through untouched. A feature is
    cut into the Delaunay along its OWN segments, so one long segment is not a
    Delaunay edge and the cut does not happen at all -- the layout came back as
    though there were no feature.
    """
    from compas_singular.geometry import discretise_line

    line = discretise_line([[0.0, 0.0, 0.0], [10.0, 10.0, 0.0]], 0.2)
    assert len(line) == 72
    # the EXTREMITIES must survive: whether an end lands on a wall decides
    # whether the layout gets a node there
    assert line[0] == [0.0, 0.0, 0.0]
    assert line[-1] == [10.0, 10.0, 0.0]
    assert max(((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2) ** 0.5
               for a, b in zip(line, line[1:])) <= 0.2 + 1e-9

    # every input point survives, so a kink stays a kink
    kinked = discretise_line([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0]], 1.0)
    assert [2.0, 0.0, 0.0] in kinked

    # None is the explicit opt-out
    assert discretise_line([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], None) == [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]

    # and end to end: the two-point guide now gives the Fig 4.17 layout
    square = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]]
    decomposition = SkeletonDecomposition.from_boundary(
        square, polyline_features=[[[0.0, 0.0, 0.0], [10.0, 10.0, 0.0]]], target_length=0.2)
    coarse = decomposition.decomposition_mesh([])
    assert len(decomposition.corner_faces()) == 6
    assert coarse.number_of_faces() == 6
    assert all(len(coarse.face_vertices(f)) == 4 for f in coarse.faces())


def test_a_collapsed_edge_is_opened_enough_to_be_an_edge():
    """``solve_triangular_faces`` case 1 duplicates a vertex, leaving a
    zero-length edge that densification divides by. The tail of the method opens
    the pair up again.

    The fraction was 0.1, too small to produce an edge: on a real Rhino plate it
    left a 0.019 edge -- aspect ratio 376, a 178.6 degree corner, a singularity
    reading as a tiny edge rather than a point.

    Built by hand: the free guide that used to reach case 1 no longer does, now
    its tips are layout vertices. Two triangles, each with ONE boundary corner.
    On this unit layout the old fraction opens the pair to 0.09, the current one
    to 0.43.
    """
    vertices = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [0.7, 1, 0], [1.3, 1, 0],
                [2, 1, 0], [0, 2, 0], [1, 2, 0], [2, 2, 0]]
    faces = [[0, 1, 4, 3], [1, 2, 6, 5], [1, 5, 4], [3, 4, 8, 7], [4, 5, 8], [5, 6, 9, 8]]
    decomposition = SkeletonDecomposition()
    decomposition.mesh = CoarsePseudoQuadMesh.from_vertices_and_faces(vertices, faces)

    decomposition.solve_triangular_faces()
    coarse = decomposition.mesh

    assert _face_sizes(coarse) == {4: 6}, 'case 1 did not turn the triangles into quads'
    shortest = min(coarse.edge_length(*e) for e in coarse.edges())
    assert shortest > 0.25, 'collapsed edge left at %.4f' % shortest


def test_case_two_survives_a_missing_decomposition_polyline():
    """``solve_triangular_faces`` case 2 merges two coincident boundary
    singularities. It consulted ``decomposition_polyline`` for WHERE to put the
    merged vertex and passed the result straight to ``Polyline`` -- so when no
    branch joined the pair it raised ``TypeError: 'NoneType' object is not
    iterable``, seen on a real plate with no curve features at all.

    The merge is the point of the case; the position is secondary, and the two
    vertices are coincident anyway.
    """
    square = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]]
    decomposition = SkeletonDecomposition.from_boundary(square, target_length=0.4)

    # force the lookup to fail the way the real plate did
    decomposition.decomposition_polyline = lambda *a, **k: None
    coarse = decomposition.decomposition_mesh([])

    assert coarse.number_of_faces() > 0
    for fkey in coarse.faces():
        for vkey in coarse.face_vertices(fkey):
            x, y, z = coarse.vertex_coordinates(vkey)
            assert x == x and y == y and z == z, 'NaN in the merged vertex'


def _arc(p, q, bulge, n=24):
    """A chord from p to q bowed sideways by `bulge`."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    length = (dx * dx + dy * dy) ** 0.5
    nx, ny = -dy / length, dx / length
    return [[p[0] + dx * (i / n) + nx * bulge * 4 * (i / n) * (1 - i / n),
             p[1] + dy * (i / n) + ny * bulge * 4 * (i / n) * (1 - i / n), 0.0]
            for i in range(n + 1)]


def test_symmetric_curve_features_do_not_crash_the_collapsed_boundary_fix():
    """``branches_splitting_collapsed_boundaries`` built ``splits`` as a SET and
    then indexed it: ``splits[0]`` raises ``TypeError: 'set' object is not
    subscriptable`` whenever a boundary loop has exactly ONE split.

    Reached by a square with two symmetric curve features -- 10 of 48
    configurations crashed. A set would also have made the two-split branch's
    output depend on iteration order.
    """
    square = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]]
    failures = []
    for spacing in (0.3, 0.4):
        for bulge in (0.5, 1.5, 2.5):
            features = [_arc([0, 5, 0], [10, 5, 0], bulge),
                        _arc([0, 5, 0], [10, 5, 0], -bulge)]
            try:
                decomposition = SkeletonDecomposition.from_boundary(
                    square, polyline_features=features, target_length=spacing)
                coarse = decomposition.decomposition_mesh([])
                assert coarse.number_of_faces() > 0
            except Exception as exc:                     # noqa: BLE001
                failures.append('spacing %s bulge %s: %s' % (spacing, bulge, exc))
    assert not failures, failures
