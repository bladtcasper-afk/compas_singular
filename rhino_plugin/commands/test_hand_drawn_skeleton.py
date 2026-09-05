#! python3

# r: compas

"""**Step: quad mesh, with the coarse layout rebuilt from the POLYLINE network.**

    reads   TopologyProblem::Skeleton::Polylines   the network
            TopologyProblem::Skeleton::Poles       collapsed corners, if any
            TopologyProblem::InputBoundaries::*    outer + inner loops
    writes  TopologyProblem::QuadMesh

The layout is built with ``CoarsePseudoQuadMesh.from_polylines``, so the
network on the layer IS the layout -- ``Skeleton::Mesh`` is not read at all.

WHY THE BOUNDARIES HAVE TO GO IN
--------------------------------
``from_polylines(boundary_polylines, other_polylines)`` makes a mesh vertex at
every polyline ENDPOINT, and uses the boundary list to decide which recovered
faces are real. Handing it the interior curves alone does not merely lose the
wall -- it returns one meaningless face. So the walls are rebuilt here as a
true partition of each loop, cut at:

* every polyline endpoint that lands on that loop, and
* every CORNER of the loop, meaning a tangent jump above 45 degrees.

The corner marks are not optional. Without them a wall arc runs straight past a
square's corner, and the patch on it comes back 3-sided: measured on a square
with a 4-arm cross, 4 triangles instead of 4 quads. The 45 degree threshold is
the one ``framefield.repair.SHARP_TURN`` uses, and for the same reason -- a
turn that sharp cannot be a sample of a curve, so a circle discretised into 36
segments still yields 4 arcs, not 36.

Rebuilding rather than trusting the arcs already on the layer also guarantees
the wall is tiled exactly ONCE. Two overlapping arcs are a boundary drawn
twice, and nothing downstream recovers from that.

THE ONE CASE THIS CANNOT MESH
-----------------------------
``from_polylines`` keeps a face only when at least one of its corners is NOT on
a boundary polyline. A network whose cuts all run wall-to-wall therefore
produces NOTHING -- measured, a square with a single cut across it and a square
plate with a hole and four radial cuts both come back with 0 faces. That is a
property of the constructor, not of the input, so it is reported here rather
than worked around. Draw at least one interior node, or switch to
``framefield.arrangement.faces_from_arrangement``, which was written for
exactly this case.
"""
from math import atan2, pi

from CMD_start import get_settings, import_compas_singular
import_compas_singular()

import rhinoscriptsyntax as rs

import compas_rhino as cr
from compas_rhino.conversions import curve_to_compas_polyline, point_to_compas

from compas.geometry import distance_point_point
from compas.itertools import pairwise

from compas_singular.datastructures import CoarsePseudoQuadMesh
from compas_singular.rhino.helpers.helpers import bake_mesh, read_boundaries

#: Two points this close are the same node. Tied to the background spacing, the
#: same scale every other stage snaps at.
settings = get_settings()
TOL = settings["triangulation_spacing"] * 0.5
#: A tangent jump above this is a corner of the domain, whatever its neighbours do.
SHARP_TURN = pi / 4.0
#: Anything shorter than this rounds to a single node in a geometric key, so an
#: arc below it is a self-loop rather than an edge.
GKEY_RESOLUTION = 1e-3


# ----------------------------------------------------------------------
# geometry
# ----------------------------------------------------------------------

def to_xyz(curve, close=False):
    """Point list, flattened to z = 0. ``close`` drops a repeated closing point.

    A closed Rhino curve converts with its first point repeated at the end.
    Left in, that is a zero-length segment and a duplicate mark, so a loop is
    always read with ``close=True`` and an open polyline never is.
    """
    points = [[float(p[0]), float(p[1]), 0.0] for p in curve]
    if close and len(points) > 1 and distance_point_point(points[0], points[-1]) < 1e-9:
        points = points[:-1]
    return points


def cumulative(loop):
    """Arc-length position of each loop vertex, plus the loop's total length."""
    cum = [0.0]
    for a, b in pairwise(list(loop) + list(loop[:1])):
        cum.append(cum[-1] + distance_point_point(a, b))
    return cum


def loop_param(loop, cum, p):
    """``(distance, arc length)`` of the loop point closest to ``p``."""
    best_d, best_s = float("inf"), 0.0
    for i, (a, b) in enumerate(pairwise(list(loop) + list(loop[:1]))):
        abx, aby = b[0] - a[0], b[1] - a[1]
        length2 = abx * abx + aby * aby
        if length2 == 0.0:
            continue
        t = ((p[0] - a[0]) * abx + (p[1] - a[1]) * aby) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + abx * t, a[1] + aby * t, 0.0]
        d = distance_point_point(p, q)
        if d < best_d:
            best_d, best_s = d, cum[i] + t * (cum[i + 1] - cum[i])
    return best_d, best_s


def loop_corners(loop):
    """Loop vertices where the tangent jumps by more than :data:`SHARP_TURN`."""
    out = []
    n = len(loop)
    for i in range(n):
        a, b, c = loop[i - 1], loop[i], loop[(i + 1) % n]
        v1 = (b[0] - a[0], b[1] - a[1])
        v2 = (c[0] - b[0], c[1] - b[1])
        if (v1[0] or v1[1]) and (v2[0] or v2[1]):
            turn = abs(atan2(v1[0] * v2[1] - v1[1] * v2[0],
                             v1[0] * v2[0] + v1[1] * v2[1]))
            if turn > SHARP_TURN:
                out.append([float(b[0]), float(b[1]), 0.0])
    return out


def partition_loop(loop, marks, tol):
    """Cut a closed loop into arcs at ``marks``: a partition, exactly once round.

    Marks are ordered by ARC LENGTH, not by vertex index. Index order is not
    total -- two marks on the same segment sort arbitrarily -- and a walk that
    takes them in the wrong order runs the long way round, so the arcs it emits
    OVERLAP and the wall ends up drawn one and a half times.
    """
    cum = cumulative(loop)
    total = cum[-1]
    if total <= 0.0:
        return []

    placed = []
    for p in marks:
        d, s = loop_param(loop, cum, p)
        if d <= tol:
            placed.append((s % total, [float(p[0]), float(p[1]), 0.0]))

    # Fewer than two cuts leaves a closed arc, which bounds a degenerate patch.
    # Three arbitrary cuts is what SkeletonDecomposition does about the same
    # case in branches_splitting_collapsed_boundaries.
    if len(placed) < 2:
        placed = [(total * k / 3.0, None) for k in range(3)]

    placed.sort(key=lambda m: m[0])
    unique = [placed[0]]
    for s, p in placed[1:]:
        if s - unique[-1][0] > max(tol * 0.25, GKEY_RESOLUTION):
            unique.append((s, p))
    if len(unique) > 1 and (total - unique[-1][0]) + unique[0][0] <= max(tol * 0.25, GKEY_RESOLUTION):
        unique.pop()
    if len(unique) < 2:
        unique = [(total * k / 3.0, None) for k in range(3)]

    def point_at(s):
        s = s % total
        for i, (c0, c1) in enumerate(zip(cum, cum[1:])):
            if c0 - 1e-12 <= s <= c1 + 1e-12:
                a, b = loop[i], loop[(i + 1) % len(loop)]
                span = c1 - c0
                t = 0.0 if span == 0 else (s - c0) / span
                return [a[k] + (b[k] - a[k]) * t for k in range(3)]
        return [float(c) for c in loop[0]]

    arcs = []
    n = len(unique)
    for k in range(n):
        s0, p0 = unique[k]
        s1, p1 = unique[(k + 1) % n]
        a0 = list(p0) if p0 is not None else point_at(s0)
        a1 = list(p1) if p1 is not None else point_at(s1)
        s1w = s1 if s1 > s0 else s1 + total

        # the loop's own vertices in between, so a curved wall keeps its shape
        between = []
        for i, sv in enumerate(cum[:-1]):
            for s in (sv, sv + total):
                if s0 + 1e-9 < s < s1w - 1e-9:
                    between.append((s, [float(c) for c in loop[i]]))
        between.sort(key=lambda m: m[0])

        arc = [a0]
        for _, p in between:
            if distance_point_point(p, arc[-1]) > GKEY_RESOLUTION:
                arc.append(p)
        if distance_point_point(a1, arc[-1]) <= GKEY_RESOLUTION and len(arc) > 1:
            arc.pop()
        arc.append(a1)
        if len(arc) >= 2:
            arcs.append(arc)
    return arcs


def split_network(polylines, loops, tol):
    """``(wall arcs, interior polylines)`` -- the two arguments from_polylines wants.

    A polyline every one of whose points lies on a loop is a piece of wall and
    is dropped: the walls are rebuilt from the loops instead, so an arc already
    baked on the layer cannot end up in the network twice.
    """
    interior = []
    for points in polylines:
        on_wall = False
        for loop in loops:
            cum = cumulative(loop)
            if all(loop_param(loop, cum, p)[0] <= tol for p in points):
                on_wall = True
                break
        if not on_wall:
            interior.append(points)

    walls = []
    for loop in loops:
        cum = cumulative(loop)
        # Corners first, then every landing -- including the ends of the wall
        # arcs already on the layer, which is how the splits the decomposition
        # chose are preserved rather than re-derived.
        marks = list(loop_corners(loop))
        for points in polylines:
            for p in (points[0], points[-1]):
                if loop_param(loop, cum, p)[0] <= tol:
                    marks.append(p)
        walls.extend(partition_loop(loop, marks, tol))
    return walls, interior


# ----------------------------------------------------------------------
# read
# ----------------------------------------------------------------------

outer, inners, guides, point_features = read_boundaries()
loops = [to_xyz(outer, close=True)] + [to_xyz(loop, close=True) for loop in inners]

guids = rs.ObjectsByLayer("TopologyProblem::Skeleton::Polylines")
if not guids:
    raise RuntimeError("No polylines on 'Skeleton::Polylines' -- run the coarse "
                       "mesh command first, or draw the network there.")
polylines = [to_xyz(curve_to_compas_polyline(cr.objects.find_object(guid).Geometry))
             for guid in guids]

poles = [list(point_to_compas(cr.objects.find_object(guid).Geometry))
         for guid in rs.ObjectsByLayer("TopologyProblem::Skeleton::Poles") or []]


# ----------------------------------------------------------------------
# layout
# ----------------------------------------------------------------------

walls, interior = split_network(polylines, loops, TOL)
print("network: {} polyline(s) in -> {} wall arc(s) + {} interior".format(
    len(polylines), len(walls), len(interior)))

coarse_mesh = CoarsePseudoQuadMesh.from_polylines(walls, interior)

if coarse_mesh.number_of_faces() == 0:
    raise RuntimeError(
        "from_polylines recovered no face. It keeps a face only when at least "
        "one corner is OFF the boundary, so a network whose cuts all run "
        "wall-to-wall yields nothing -- a square with one cut across it, or a "
        "plate with a hole and radial cuts, both give 0 faces. Draw at least "
        "one interior node, or use framefield.arrangement.faces_from_arrangement.")

# Every 3-sided face is registered as a pseudo-quad here, which is what makes it
# densifiable; a pole point on the layer only decides WHICH corner collapses.
# Note the two repairs SkeletonDecomposition applies after its own from_polylines
# -- solve_triangular_faces and split_quads_with_poles -- are NOT run, so a
# 5-sided patch survives to break collect_strips. The histogram is the warning.
vertices, faces = coarse_mesh.to_vertices_and_faces()
coarse_mesh = CoarsePseudoQuadMesh.from_vertices_and_faces_with_poles(vertices, faces, poles)

sides = {}
for fkey in coarse_mesh.faces():
    n = len(coarse_mesh.face_vertices(fkey))
    sides[n] = sides.get(n, 0) + 1
print("layout: {} patch(es), {} corner(s), sides {}".format(
    coarse_mesh.number_of_faces(), coarse_mesh.number_of_vertices(), sides))
odd = sorted(n for n in sides if n not in (3, 4))
if odd:
    print("  patches with {} side(s) -- densification needs quads; split them in "
          "the network and re-run.".format(odd))


# ----------------------------------------------------------------------
# density -> quad mesh
# ----------------------------------------------------------------------

coarse_mesh.collect_strips()
coarse_mesh.set_strips_density_target(t=settings["target_length"])
coarse_mesh.densification()
dense_mesh = coarse_mesh.get_quad_mesh()

print("quad mesh: {} face(s) at target length {}".format(
    dense_mesh.number_of_faces(), settings["target_length"]))

rs.AddLayer("QuadMesh", parent="TopologyProblem")
bake_mesh(dense_mesh, "QuadMesh")