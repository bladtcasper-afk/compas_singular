"""Headless tests for the COMPAS-only constrained smoothing.

No Rhino: every projection goes through :mod:`compas.geometry`.
"""
import random
from math import pi

import pytest

from compas.datastructures import Mesh
from compas.geometry import Bezier
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import Polyline
from compas.geometry import distance_point_point

from compas_singular.datastructures import automated_boundary_constraints
from compas_singular.datastructures import boundary_constrained_smoothing
from compas_singular.datastructures import closest_point_on_constraint
from compas_singular.datastructures import constrained_smoothing
from compas_singular.datastructures import mesh_boundary_corners
from compas_singular.datastructures import mesh_boundary_loops
from compas_singular.datastructures import mesh_boundary_polylines
from compas_singular.datastructures import smoothing_region


SIDE = 10.0


def _distance_to_square_border(xyz, side=SIDE):
    x, y, _ = xyz
    return min(abs(x), abs(y), abs(x - side), abs(y - side))


def _perturbed_grid(nx=8, amplitude=0.35, seed=0):
    """A square grid whose interior vertices have been shaken."""
    mesh = Mesh.from_meshgrid(dx=SIDE, nx=nx)
    rng = random.Random(seed)
    for vertex in mesh.vertices():
        if mesh.is_vertex_on_boundary(vertex):
            continue
        x, y, z = mesh.vertex_coordinates(vertex)
        mesh.vertex_attributes(vertex, 'xyz', [x + rng.uniform(-amplitude, amplitude),
                                               y + rng.uniform(-amplitude, amplitude), z])
    return mesh


def _grid_with_hole(nx=8):
    """A square grid with the four central faces removed, so it has two boundaries."""
    mesh = Mesh.from_meshgrid(dx=SIDE, nx=nx)
    centre = [SIDE / 2, SIDE / 2, 0]
    faces = sorted(mesh.faces(), key=lambda face: distance_point_point(mesh.face_centroid(face), centre))
    for face in faces[:4]:
        mesh.delete_face(face)
    mesh.remove_unused_vertices()
    return mesh


# ==============================================================================
# Closest point
# ==============================================================================

def test_closest_point_on_constraint_dispatch():
    polyline = Polyline([[0, 0, 0], [10, 0, 0], [10, 10, 0]])

    assert closest_point_on_constraint(None, [5, 3, 0]) is None
    assert closest_point_on_constraint(Point(1, 2, 3), [5, 3, 0]) == [1.0, 2.0, 3.0]
    assert closest_point_on_constraint([1, 2, 3], [5, 3, 0]) == [1.0, 2.0, 3.0]
    assert closest_point_on_constraint(polyline, [5, 3, 0]) == [5.0, 0.0, 0.0]
    assert closest_point_on_constraint(list(polyline), [5, 3, 0]) == [5.0, 0.0, 0.0]
    assert closest_point_on_constraint(Line([0, 0, 0], [10, 0, 0]), [5, 3, 0]) == [5.0, 0.0, 0.0]

    # analytical closest point, from the geometry itself
    on_circle = closest_point_on_constraint(Circle(5.0, Frame.worldXY()), [10, 0, 0])
    assert distance_point_point(on_circle, [5, 0, 0]) < 1e-9

    # no analytical closest point: falls back on the polyline discretisation
    bezier = Bezier([[0, 0, 0], [5, 5, 0], [10, 0, 0]])
    on_bezier = closest_point_on_constraint(bezier, [5, 10, 0])
    assert distance_point_point(on_bezier, [5, 2.5, 0]) < 1e-2


def test_closest_point_on_constraint_rejects_unsupported():
    with pytest.raises(TypeError):
        closest_point_on_constraint('not geometry', [0, 0, 0])


# ==============================================================================
# Boundary
# ==============================================================================

def test_mesh_boundary_loops_are_not_repeating():
    mesh = Mesh.from_meshgrid(dx=SIDE, nx=8)
    loops = mesh_boundary_loops(mesh)
    assert len(loops) == 1
    assert loops[0][0] != loops[0][-1]
    assert len(loops[0]) == len(set(loops[0])) == 4 * 8


def test_mesh_boundary_polylines_are_closed():
    polylines = mesh_boundary_polylines(_grid_with_hole())
    assert len(polylines) == 2
    assert all(polyline[0] == polyline[-1] for polyline in polylines)


def test_mesh_boundary_corners_finds_the_square_corners():
    mesh = Mesh.from_meshgrid(dx=SIDE, nx=8)
    corners = mesh_boundary_corners(mesh, corner_angle=pi / 6)
    assert len(corners) == 4
    for vertex in corners:
        x, y, _ = mesh.vertex_coordinates(vertex)
        assert x in (0.0, SIDE) and y in (0.0, SIDE)


# ==============================================================================
# Smoothing
# ==============================================================================

def test_boundary_constrained_smoothing_keeps_the_outline_and_the_corners():
    mesh = _perturbed_grid()
    corners_before = {vertex: mesh.vertex_coordinates(vertex) for vertex in mesh_boundary_corners(mesh)}

    constraints = boundary_constrained_smoothing(mesh, kmax=50, damping=0.5, algorithm='area')

    # every boundary vertex is constrained, and the corners are pinned to a point
    assert set(constraints) == set(mesh.vertices_on_boundary())
    assert sum(1 for constraint in constraints.values() if isinstance(constraint, Point)) == 4

    # the outline is kept exactly ...
    for vertex in mesh.vertices_on_boundary():
        assert _distance_to_square_border(mesh.vertex_coordinates(vertex)) < 1e-9

    # ... and so are the corners
    for vertex, xyz in corners_before.items():
        assert distance_point_point(mesh.vertex_coordinates(vertex), xyz) < 1e-9


def test_boundary_constrained_smoothing_lets_the_boundary_slide():
    start = _perturbed_grid()

    fixed = start.copy()
    fixed.smooth_area(fixed=list(fixed.vertices_on_boundary()), kmax=50, damping=0.5)

    constrained = start.copy()
    boundary_constrained_smoothing(constrained, kmax=50, damping=0.5, algorithm='area')

    def slide(mesh):
        return max(distance_point_point(mesh.vertex_coordinates(vertex), start.vertex_coordinates(vertex))
                   for vertex in mesh.vertices_on_boundary() if vertex not in mesh_boundary_corners(mesh))

    assert slide(fixed) < 1e-9
    assert slide(constrained) > 1e-3


def test_boundary_constrained_smoothing_can_be_told_not_to_pin_corners():
    mesh = _perturbed_grid()
    constraints = automated_boundary_constraints(mesh, fix_corners=False)
    assert not any(isinstance(constraint, Point) for constraint in constraints.values())


def test_automated_boundary_constraints_matches_each_boundary_to_its_curve():
    mesh = _grid_with_hole()
    outer, inner = mesh_boundary_polylines(mesh)

    constraints = automated_boundary_constraints(mesh, curves=[outer, inner])

    for vertex, constraint in constraints.items():
        if isinstance(constraint, Point):
            continue
        xyz = mesh.vertex_coordinates(vertex)
        on_outer = _distance_to_square_border(xyz) < 1e-9
        assert constraint is (outer if on_outer else inner)


def test_boundary_constrained_smoothing_holds_every_boundary_of_a_mesh_with_a_hole():
    """The hole must be held too.

    ``mesh.vertices_on_boundary()`` returns the LONGEST boundary only, so the usual
    ``smooth(fixed=mesh.vertices_on_boundary())`` idiom leaves a hole free to blow up.
    The automated constraints go through every boundary loop instead.
    """
    mesh = _grid_with_hole()
    before = {vertex: mesh.vertex_coordinates(vertex)
              for loop in mesh_boundary_loops(mesh) for vertex in loop}
    hole = mesh_boundary_polylines(mesh)[1]  # the hole as it is BEFORE smoothing

    boundary_constrained_smoothing(mesh, kmax=50, damping=0.5, algorithm='area')

    for vertex, xyz in before.items():
        moved = mesh.vertex_coordinates(vertex)
        if _distance_to_square_border(xyz) < 1e-9:
            assert _distance_to_square_border(moved) < 1e-9
        else:
            assert distance_point_point(moved, closest_point_on_constraint(hole, moved)) < 1e-9

    # and the free version really does lose the hole, which is what this guards against
    loose = _grid_with_hole()
    loose.smooth_area(fixed=list(loose.vertices_on_boundary()), kmax=50, damping=0.5)
    drift = max(distance_point_point(loose.vertex_coordinates(vertex), closest_point_on_constraint(hole, loose.vertex_coordinates(vertex)))
                for vertex in mesh_boundary_loops(loose)[1])
    assert drift > 1e-3


def test_constrained_smoothing_pulls_a_free_mesh_onto_its_constraints():
    mesh = _perturbed_grid()
    circle = Circle(4.0, Frame([SIDE / 2, SIDE / 2, 0], [1, 0, 0], [0, 1, 0]))

    constraints = {vertex: circle for vertex in mesh.vertices_on_boundary()}
    constrained_smoothing(mesh, kmax=30, damping=0.5, constraints=constraints, algorithm='centroid')

    for vertex in mesh.vertices_on_boundary():
        xyz = mesh.vertex_coordinates(vertex)
        assert abs(distance_point_point(xyz, [SIDE / 2, SIDE / 2, 0]) - 4.0) < 1e-6


def test_constrained_smoothing_honours_fixed_vertices():
    mesh = _perturbed_grid()
    interior = [vertex for vertex in mesh.vertices() if not mesh.is_vertex_on_boundary(vertex)]
    frozen = interior[0]
    before = mesh.vertex_coordinates(frozen)

    constrained_smoothing(mesh, kmax=20, damping=0.5, algorithm='centroid', fixed=[frozen])

    assert distance_point_point(mesh.vertex_coordinates(frozen), before) < 1e-9


# ==============================================================================
# Region
# ==============================================================================

def _core_around(mesh, centre, radius):
    """The vertices within `radius` of a point."""
    return [vertex for vertex in mesh.vertices()
            if distance_point_point(mesh.vertex_coordinates(vertex), centre) < radius]


def _rings(mesh, core, count):
    """The `count` rings of vertices around a core, as a list of sets."""
    out = []
    seen = set(core)
    frontier = set(core)
    for _ in range(count):
        frontier = {n for vertex in frontier for n in mesh.vertex_neighbors(vertex)} - seen
        seen |= frontier
        out.append(frontier)
    return out


def _relax_region_as_it_was(mesh, weight, kmax, damping, constraints):
    """The original `relax_region` from the Rhino command, kept verbatim.

    This is what `smoothing_region` was lifted out of. It is here so the restructuring
    can be checked against it rather than asserted to be equivalent.
    """
    for _ in range(kmax):
        target = {}
        for vertex, w in weight.items():
            if w <= 0.0:
                continue
            areas, centroids = [], []
            for face in mesh.vertex_faces(vertex, ordered=True):
                if face is None:
                    continue
                areas.append(mesh.face_area(face))
                centroids.append(mesh.face_centroid(face))
            total = sum(areas)
            if not areas or total <= 0.0:
                continue
            target[vertex] = [
                sum(a * c[i] for a, c in zip(areas, centroids)) / total for i in range(3)]

        for vertex, goal in target.items():
            xyz = mesh.vertex_coordinates(vertex)
            step = damping * weight[vertex]
            mesh.vertex_attributes(vertex, 'xyz', [
                xyz[i] + step * (goal[i] - xyz[i]) for i in range(3)])

        for vertex, constraint in constraints.items():
            moved = closest_point_on_constraint(constraint, mesh.vertex_coordinates(vertex))
            if moved is not None:
                mesh.vertex_attributes(vertex, 'xyz', moved)


def test_smoothing_region_matches_the_original_relax_region():
    """The restructuring must not have changed a single coordinate."""
    start = _perturbed_grid(nx=12, seed=2)
    core = _core_around(start, [SIDE / 2, SIDE / 2, 0], 3.0)

    new = start.copy()
    weights = smoothing_region(new, core, kmax=20, damping=0.5, blend=3)

    # the old call site: taper outside, boundary constraints assembled and filtered there
    old = start.copy()
    boundary = automated_boundary_constraints(old)
    constraints = {v: c for v, c in boundary.items() if weights.get(v, 0.0) > 0.0}
    _relax_region_as_it_was(old, weights, 20, 0.5, constraints)

    for vertex in start.vertices():
        assert distance_point_point(new.vertex_coordinates(vertex),
                                    old.vertex_coordinates(vertex)) < 1e-12


def test_smoothing_region_only_moves_the_region():
    start = _perturbed_grid(nx=12, seed=3)
    core = _core_around(start, [SIDE / 2, SIDE / 2, 0], 2.0)

    mesh = start.copy()
    weights = smoothing_region(mesh, core, kmax=30, damping=0.5, blend=3)

    moved = [vertex for vertex in mesh.vertices()
             if distance_point_point(mesh.vertex_coordinates(vertex),
                                     start.vertex_coordinates(vertex)) > 1e-9]
    assert moved                            # something happened
    assert set(moved) <= set(weights)       # and only inside the region

    # the core really relaxed, and the ring beyond the blend is untouched
    assert max(distance_point_point(mesh.vertex_coordinates(vertex),
                                    start.vertex_coordinates(vertex)) for vertex in core) > 1e-3
    for vertex in _rings(start, core, 4)[3]:
        assert vertex not in weights


def test_smoothing_region_tapers_the_damping_over_the_blend():
    mesh = _perturbed_grid(nx=12, seed=4)
    core = _core_around(mesh, [SIDE / 2, SIDE / 2, 0], 2.0)

    weights = smoothing_region(mesh, core, kmax=1, damping=0.5, blend=3)

    assert all(weights[vertex] == 1.0 for vertex in core)
    previous = 1.0
    for ring in _rings(mesh, core, 3):
        weight = {weights[vertex] for vertex in ring}
        assert len(weight) == 1                 # one weight per ring
        assert 0.0 < weight.pop() < previous    # strictly decreasing, never zero
        previous = min(weights[vertex] for vertex in ring)

    # no blend at all: the core and nothing else
    assert set(smoothing_region(mesh, core, kmax=1, damping=0.5, blend=0)) == set(core)


def test_smoothing_region_keeps_a_region_on_the_boundary():
    """A zone touching the outline slides along it instead of being dragged inward."""
    start = _perturbed_grid(nx=12, seed=5)
    core = _core_around(start, [SIDE / 2, 0, 0], 2.0)
    assert any(start.is_vertex_on_boundary(vertex) for vertex in core)

    mesh = start.copy()
    weights = smoothing_region(mesh, core, kmax=30, damping=0.5, blend=3)

    on_boundary = [vertex for vertex in weights if mesh.is_vertex_on_boundary(vertex)]
    assert on_boundary
    for vertex in on_boundary:
        assert _distance_to_square_border(mesh.vertex_coordinates(vertex)) < 1e-9

    # and they slid along it rather than staying pinned
    assert max(distance_point_point(mesh.vertex_coordinates(vertex),
                                    start.vertex_coordinates(vertex)) for vertex in on_boundary) > 1e-3


def test_smoothing_region_accepts_ready_made_weights():
    start = _perturbed_grid(nx=12, seed=6)
    core = _core_around(start, [SIDE / 2, SIDE / 2, 0], 2.0)

    tapered = start.copy()
    weights = smoothing_region(tapered, core, kmax=20, damping=0.5, blend=2)

    given = start.copy()
    assert smoothing_region(given, weights, kmax=20, damping=0.5, blend=99) == weights

    for vertex in start.vertices():
        assert distance_point_point(tapered.vertex_coordinates(vertex),
                                    given.vertex_coordinates(vertex)) < 1e-12


def test_smoothing_region_fixes_vertices_without_a_positive_weight():
    start = _perturbed_grid(nx=12, seed=7)
    core = _core_around(start, [SIDE / 2, SIDE / 2, 0], 2.0)
    frozen = core[0]

    mesh = start.copy()
    weights = {vertex: (0.0 if vertex == frozen else 1.0) for vertex in core}
    smoothing_region(mesh, weights, kmax=20, damping=0.5)

    assert distance_point_point(mesh.vertex_coordinates(frozen), start.vertex_coordinates(frozen)) < 1e-9


def test_windowed_projection_matches_the_full_search():
    """The projector only searches around the previous result: same answer, less work."""
    import sys

    # not `import compas_singular.datastructures.mesh.smoothing`: the package exports a
    # `mesh` module of its own, which shadows the `mesh` subpackage on attribute lookup
    smoothing = sys.modules['compas_singular.datastructures.mesh.smoothing']

    def run(window):
        default = smoothing._PolylineProjector.__init__.__defaults__
        smoothing._PolylineProjector.__init__.__defaults__ = (window,)
        try:
            mesh = _perturbed_grid(nx=12, seed=1)
            boundary_constrained_smoothing(mesh, kmax=50, damping=0.5, algorithm='area')
            return [mesh.vertex_coordinates(vertex) for vertex in sorted(mesh.vertices())]
        finally:
            smoothing._PolylineProjector.__init__.__defaults__ = default

    windowed = run(10)
    exhaustive = run(10 ** 6)
    assert all(distance_point_point(a, b) < 1e-12 for a, b in zip(windowed, exhaustive))
