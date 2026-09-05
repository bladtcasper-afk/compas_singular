"""Constrained smoothing of a quad mesh, without Rhino.

The old constrained smoothing (``compas_singular.rhino.constrained_smoothing``) needed
Rhino to project a vertex onto its constraint. This one projects with
:mod:`compas.geometry` instead, so it runs headless, and it can constrain the boundary
vertices automatically:

* ``constrained_smoothing``            -- smooth, snapping constrained vertices back onto
                                          their point / polyline / curve after every iteration.
* ``automated_boundary_constraints``   -- find the boundary vertices and hand each one the
                                          piece of boundary it has to stay on, pinning kinks.
* ``boundary_constrained_smoothing``   -- the two above in one call.

Freezing the boundary (``fixed=boundary_vertices``) keeps the outline but leaves the rows
next to it distorted, because those vertices cannot follow. Constraining the boundary lets
the vertices *slide along* it: the outline is kept exactly, and the pattern relaxes through
to the edge.
"""
from math import degrees
from math import pi

from compas.colors import Color
from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import Polyline
from compas.geometry import Translation
from compas.geometry import angle_points
from compas.geometry import distance_point_point

from compas_singular.algorithms import SkeletonDecomposition
from compas_singular.algorithms import boundary_triangulation
from compas_singular.datastructures import automated_boundary_constraints
from compas_singular.datastructures import boundary_constrained_smoothing
from compas_singular.datastructures import closest_point_on_constraint
from compas_singular.datastructures import constrained_smoothing
from compas_singular.datastructures import mesh_boundary_loops


# -----------------------------------------------------------------------------
# Define boundary curves
# -----------------------------------------------------------------------------

outer = [Point(5, 5, 0), Point(5, -5, 0), Point(-5, -5, 0), Point(-5, 5, 0), Point(0, 7, 0)]
inner = [Point(1, 1, 0), Point(1, -1, 0), Point(-1, -1, 0), Point(-1, 1, 0)]

outer_curve = Polyline(outer + [outer[0]])
inner_curve = Polyline(inner + [inner[0]])


def densify(curve, resolution):
    # Important to discretise the sides and never the entire polyline to make sure that
    # hard corners are kept as a point.
    if isinstance(curve, Polyline):
        line_pts = curve.split_at_corners(pi / 10)
    else:
        line_pts = list(curve) + [curve[0]]

    densified_pts = []
    for start, end in zip(line_pts[:-1], line_pts[1:]):
        densified_pts.extend(Line(start, end).to_polyline(n=resolution).points)
    return densified_pts


outer_pts = densify(outer, 5)
inner_pts = densify(inner, 5)


def to_coordinates(pts):
    return [[pt[0], pt[1], pt[2]] for pt in pts]


# -----------------------------------------------------------------------------
# Pattern: triangulation -> skeleton decomposition -> densification
# -----------------------------------------------------------------------------

point_features = []

trimesh = boundary_triangulation(
    outer_boundary=to_coordinates(outer_pts),
    inner_boundaries=[to_coordinates(inner_pts)],
    point_features=point_features)

decomposition = SkeletonDecomposition.from_mesh(trimesh)
coarse_mesh = decomposition.decomposition_mesh(point_features)

coarse_mesh.collect_strips()
coarse_mesh.set_strips_density(d=5)
coarse_mesh.densification()
dense_mesh = coarse_mesh.get_quad_mesh()


# -----------------------------------------------------------------------------
# Smoothing
# -----------------------------------------------------------------------------

def all_boundary_vertices(mesh):
    # NOT mesh.vertices_on_boundary(): in COMPAS 2 that returns the LONGEST boundary only,
    # so on a mesh with a hole the hole is left out -- and left free to blow up when it is
    # used as the `fixed` set of a smoothing algorithm.
    return [vertex for loop in mesh_boundary_loops(mesh) for vertex in loop]


boundary_vertices = all_boundary_vertices(dense_mesh)
input_curves = [outer_curve, inner_curve]


def smoothed_copy(mesh, smoother=None, kmax=50, damping=0.5):
    copy = mesh.copy()
    if smoother:
        smoother(copy, kmax=kmax, damping=damping)
    return copy


def fixed_boundary(mesh, kmax, damping):
    """The old way: the boundary vertices cannot move at all."""
    mesh.smooth_area(fixed=all_boundary_vertices(mesh), kmax=kmax, damping=damping)


def constrained_to_mesh_boundary(mesh, kmax, damping):
    """Auto: every boundary vertex slides along the boundary of the mesh itself."""
    boundary_constrained_smoothing(mesh, kmax=kmax, damping=damping, algorithm='area')


def constrained_to_input_curves(mesh, kmax, damping):
    """Auto, but snapping onto the input curves instead of their mesh discretisation."""
    boundary_constrained_smoothing(
        mesh, curves=input_curves, kmax=kmax, damping=damping, algorithm='area')


def constrained_by_hand(mesh, kmax, damping):
    """The same machinery, driven by a constraint dictionary written by hand."""
    constraints = {vertex: outer_curve for vertex in all_boundary_vertices(mesh)}
    for vertex in all_boundary_vertices(mesh):
        xyz = mesh.vertex_coordinates(vertex)
        if distance_point_point(xyz, closest_point_on_constraint(inner_curve, xyz)) < \
           distance_point_point(xyz, closest_point_on_constraint(outer_curve, xyz)):
            constraints[vertex] = inner_curve
    for corner in outer + inner:
        closest = min(mesh.vertices(), key=lambda v: distance_point_point(corner, mesh.vertex_coordinates(v)))
        constraints[closest] = Point(*corner)
    constrained_smoothing(mesh, kmax=kmax, damping=damping, constraints=constraints, algorithm='area')


smoothing_variants = [
    ("raw (unsmoothed)", smoothed_copy(dense_mesh)),
    ("area, boundary fixed", smoothed_copy(dense_mesh, fixed_boundary)),
    ("area, boundary constrained (mesh boundary)", smoothed_copy(dense_mesh, constrained_to_mesh_boundary)),
    ("area, boundary constrained (input curves)", smoothed_copy(dense_mesh, constrained_to_input_curves)),
    ("area, constraints by hand", smoothed_copy(dense_mesh, constrained_by_hand)),
]


# -----------------------------------------------------------------------------
# Report: are the corners kept, is the boundary kept, is the pattern better?
# -----------------------------------------------------------------------------

def boundary_drift(mesh):
    """Largest distance from a boundary vertex to the input curves."""
    return max(
        min(distance_point_point(mesh.vertex_coordinates(vertex), closest_point_on_constraint(curve, mesh.vertex_coordinates(vertex)))
            for curve in input_curves)
        for vertex in all_boundary_vertices(mesh))


def boundary_slide(mesh):
    """Largest distance a boundary vertex travelled from where it started."""
    return max(
        distance_point_point(mesh.vertex_coordinates(vertex), dense_mesh.vertex_coordinates(vertex))
        for vertex in all_boundary_vertices(mesh))


def angle_deviation(mesh):
    """Mean deviation of the face corner angles from 90 degrees, in degrees."""
    deviations = []
    for face in mesh.faces():
        vertices = mesh.face_vertices(face)
        count = len(vertices)
        for index in range(count):
            a = mesh.vertex_coordinates(vertices[index - 1])
            b = mesh.vertex_coordinates(vertices[index])
            c = mesh.vertex_coordinates(vertices[(index + 1) % count])
            deviations.append(abs(degrees(angle_points(b, a, c)) - 90.0))
    return sum(deviations) / len(deviations)


constraints = automated_boundary_constraints(dense_mesh, curves=input_curves)
corners = [vertex for vertex, constraint in constraints.items() if isinstance(constraint, Point)]
print("{} boundary vertices, {} of them pinned as corners".format(len(boundary_vertices), len(corners)))
print()
print("max drift  : how far a boundary vertex ended up OFF the input curves (must stay 0)")
print("max slide  : how far a boundary vertex travelled ALONG them")
print("mean |90-a|: mean deviation of the face angles from a right angle")
print()
print("{:<45} {:>10} {:>10} {:>14}".format("variant", "max drift", "max slide", "mean |90-a|"))
for name, mesh in smoothing_variants:
    print("{:<45} {:>10.4f} {:>10.4f} {:>13.2f}d".format(
        name, boundary_drift(mesh), boundary_slide(mesh), angle_deviation(mesh)))


# -----------------------------------------------------------------------------
# Visualisation
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    from compas_viewer import Viewer

    gap = dense_mesh.aabb().xsize + 2.0
    for index, (name, mesh) in enumerate(smoothing_variants):
        mesh.transform(Translation.from_vector([index * gap, 0, 0]))

    colors = [
        Color.grey(),
        Color.from_hex("#4C72B0"),
        Color.from_hex("#DD8452"),
        Color.from_hex("#55A868"),
        Color.from_hex("#C44E52"),
    ]

    viewer = Viewer()

    # The variants sit in a row, and the pattern is flat: look straight down on it,
    # from far enough back to hold the whole row (the default camera frames one object).
    width = (len(smoothing_variants) - 1) * gap + dense_mesh.aabb().xsize
    middle = 0.5 * (len(smoothing_variants) - 1) * gap
    viewer.renderer.camera.target = [middle, 1.0, 0.0]
    # a small offset in y, so that "up" is not ambiguous in a straight top view
    viewer.renderer.camera.position = [middle, -0.04 * width, 0.85 * width]

    step1 = viewer.scene.add_group(name="input boundary")
    step1.add(outer_curve)
    step1.add(inner_curve)

    step2 = viewer.scene.add_group(name="constrained smoothing")
    for (name, mesh), color in zip(smoothing_variants, colors):
        step2.add(mesh, name=name, facecolor=color)

    step3 = viewer.scene.add_group(name="pinned corners")
    for vertex in corners:
        step3.add(Point(*dense_mesh.vertex_coordinates(vertex)))

    viewer.show()
