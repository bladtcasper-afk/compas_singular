"""Recreate the "square with a central pole" quad pattern.

A square domain filled with quads, with a single high-valence singularity at the
centre: radial fans converging on it and four coarse patch seams running out to
the four corners. That pattern is what a *pole* (a point feature) at the centre
of a square gives you -- four pseudo-quad patches meeting at the pole, densified
with Coons patches.

Two routes are built side by side:

ROUTE A -- the idiomatic workflow (Delaunay -> skeleton decomposition -> coarse
    pseudo-quad mesh), the same pipeline as ``02_decomposition_discrete_planar.py``
    but with the geometry built in code. The patch layout is *derived* from the
    medial axis, so whether the pole survives as a coarse vertex depends on how
    the Delaunay comes out on a perfectly symmetric square.

ROUTE B -- the coarse pseudo-quad mesh written out by hand: five vertices, four
    triangular faces, pole declared. Ten lines, and it reproduces the target
    topology exactly by construction.

Measured on a 10-unit square, RESOLUTION 20, TARGET 0.5:

    route A   601 V,  576 F -- pole (valence 36) + four valence-5 vertices at
                               (+-2.071, +-2.071), one on each diagonal seam
    route B  1201 V, 1200 F -- pole (valence 80) and nothing else

Route A is the closer match to the reference: the decomposition wraps a ring of
patches around the pole, so the fan is shallow and blends into the surrounding
grid, at the cost of four extra valence-5 singularities on the diagonals. Route
B is topologically minimal -- a single singularity -- but with only four patches
the pole absorbs the *whole* boundary discretisation (valence = 4 x the side
density), giving a very steep fan. Pick on which of those two costs matters.

Run with the ``singular312`` env. ``--no-view`` for the numbers only.
"""

import sys
from math import pi

from compas.geometry import Line, Point, Polyline, Translation

from compas.tolerance import TOL
from compas.datastructures.mesh.smoothing import mesh_smooth_area
from compas_singular.algorithms import SkeletonDecomposition, boundary_triangulation
from compas_singular.datastructures import CoarsePseudoQuadMesh


# -----------------------------------------------------------------------------
# parameters -- nudge these to match a reference image
# -----------------------------------------------------------------------------

SIZE = 10.0                     # side length of the square
POLE = [0.0, 0.0, 0.0]          # where the singularity sits
RESOLUTION = 20                 # boundary discretisation per side, keep it even
TARGET = 0.5                    # target quad edge length for densification
SMOOTH_ITERS = 0                # >0 relaxes the dense mesh, boundary + pole pinned

VIEW = '--no-view' not in sys.argv


# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------

def densify(curve, resolution):
    """Discretise a closed outline, side by side.

    Important to discretise the sides and never the entire polyline, so that the
    hard corners are kept as points. (Same helper as the other examples in this
    folder -- duplicated rather than shared, per this folder's convention.)
    """
    if isinstance(curve, list) and isinstance(curve[0], Point):
        line_pts = curve
    elif isinstance(curve, Polyline):
        line_pts = curve.split_at_corners(pi / 10)
    else:
        raise TypeError(f"Wrong input type, expected Polyline or list[Point], got {type(curve)}")

    lines = [Line(line_pts[i], line_pts[i + 1]) for i in range(len(line_pts) - 1)]

    densified_pts = []
    for line in lines:
        polyline = line.to_polyline(n=resolution)
        densified_pts.extend(polyline.points)
    return densified_pts


def to_coordinates(pts):
    return [[pt.x, pt.y, pt.z] for pt in pts]


def has_pole(mesh, pole):
    """True if ``pole`` survived into the mesh as a declared pole vertex."""
    face_pole = mesh.attributes.get('face_pole') or {}
    gkey = TOL.geometric_key(pole)
    return bool(face_pole) and any(
        TOL.geometric_key(mesh.vertex_coordinates(v)) == gkey for v in face_pole.values())


def irregular_vertices(mesh):
    """Interior vertices that are not 4-valent -- the singularities."""
    boundary = {v for ring in mesh.vertices_on_boundaries() for v in ring}
    return [v for v in mesh.vertices()
            if v not in boundary and len(mesh.vertex_neighbors(v)) != 4]


def densify_coarse(coarse, target):
    coarse.collect_strips()
    coarse.set_strips_density_target(t=target)
    coarse.densification()
    return coarse.get_quad_mesh()


def relax(dense, kmax):
    """Pin the outline and the pole, relax everything else."""
    if not kmax:
        return dense
    boundary = {v for ring in dense.vertices_on_boundaries() for v in ring}
    poles = set((dense.attributes.get('face_pole') or {}).values())
    mesh_smooth_area(dense, fixed=list(boundary | poles), kmax=kmax, damping=0.5)
    return dense


def report(name, mesh, pole):
    irregular = irregular_vertices(mesh)
    at_pole = [v for v in irregular
               if TOL.geometric_key(mesh.vertex_coordinates(v)) == TOL.geometric_key(pole)]
    print(f"  {name:<28} V {mesh.number_of_vertices():>5}  F {mesh.number_of_faces():>5}"
          f"  irregular {len(irregular)}  at pole {len(at_pole)}")
    for v in irregular:
        x, y, _ = mesh.vertex_coordinates(v)
        print(f"      valence {len(mesh.vertex_neighbors(v)):>2}  at ({x:7.3f}, {y:7.3f})")


# -----------------------------------------------------------------------------
# geometry: the square outline
# -----------------------------------------------------------------------------

h = SIZE / 2.0
corners = [Point(h, h, 0), Point(-h, h, 0), Point(-h, -h, 0), Point(h, -h, 0)]
outline = corners + [corners[0]]

point_features = [POLE]


# -----------------------------------------------------------------------------
# ROUTE A -- skeleton decomposition from the outline + one point feature
# -----------------------------------------------------------------------------

print("\nROUTE A -- Delaunay -> SkeletonDecomposition")

trimesh = None
coarse_a = None
dense_a = None

# The Delaunay of a perfectly symmetric square is cocircular-prone, so the pole
# does not always survive as a coarse vertex. Different boundary discretisations
# give different triangulations -- try a few before giving up.
for resolution in (RESOLUTION, RESOLUTION + 4, RESOLUTION + 10):
    outer = to_coordinates(densify(outline, resolution))

    trimesh = boundary_triangulation(
        outer_boundary=outer,
        inner_boundaries=[],
        polyline_features=[],
        point_features=point_features)

    decomposition = SkeletonDecomposition.from_mesh(trimesh)
    candidate = decomposition.decomposition_mesh(point_features)

    if has_pole(candidate, POLE):
        coarse_a = candidate
        print(f"  pole survived at boundary resolution {resolution}")
        break
    print(f"  resolution {resolution}: pole did not survive the decomposition")

if coarse_a is None:
    print("  -> route A could not place the pole; route B below is the answer.")
else:
    dense_a = relax(densify_coarse(coarse_a, TARGET), SMOOTH_ITERS)


# -----------------------------------------------------------------------------
# ROUTE B -- the coarse pseudo-quad mesh written out by hand
# -----------------------------------------------------------------------------
#
# Five vertices, four triangular faces, all wound the same way (inconsistent
# winding does not raise -- it silently fails to weld and the dense mesh comes
# out as four disconnected patches).
#
# ``from_vertices_and_faces_with_poles`` walks every 3-vertex face and records
# the vertex whose geometric key matches a pole in ``attributes['face_pole']``,
# which is what ``densification`` reads to insert the degenerate ``None`` side
# into the Coons patch.

print("\nROUTE B -- explicit coarse pseudo-quad mesh")

vertices = to_coordinates(corners) + [POLE]
p = len(vertices) - 1
faces = [[p, 0, 1], [p, 1, 2], [p, 2, 3], [p, 3, 0]]

coarse_b = CoarsePseudoQuadMesh.from_vertices_and_faces_with_poles(
    vertices, faces, poles=[POLE])

assert has_pole(coarse_b, POLE), "pole was not registered on the coarse mesh"

# The four seams to the corners form one closed strip around the pole, so they
# densify uniformly (equal lengths on a square). Each square side is its own
# strip, terminating against the degenerate pole edge -- so the boundary grid
# and the radial fan are sized independently and both land near TARGET.
dense_b = relax(densify_coarse(coarse_b, TARGET), SMOOTH_ITERS)

print(f"  coarse: {coarse_b.number_of_vertices()} vertices, "
      f"{coarse_b.number_of_faces()} faces, {len(coarse_b.poles())} pole(s)")


# -----------------------------------------------------------------------------
# report
# -----------------------------------------------------------------------------

print("\nDENSE MESHES")
if dense_a is not None:
    report("route A (decomposition)", dense_a, POLE)
else:
    print("  route A (decomposition)      -- not built, pole did not survive")
report("route B (explicit coarse)", dense_b, POLE)
print()


# -----------------------------------------------------------------------------
# view
# -----------------------------------------------------------------------------

if VIEW:
    from compas_viewer import Viewer

    viewer = Viewer()

    if dense_a is not None:
        viewer.scene.add(dense_a, show_points=True, show_lines=True, show_faces=True)

    shown_b = dense_b.copy()
    if dense_a is not None:
        shown_b.transform(Translation.from_vector([1.5 * SIZE, 0, 0]))
    viewer.scene.add(shown_b, show_points=True, show_lines=True, show_faces=True)

    viewer.show()
