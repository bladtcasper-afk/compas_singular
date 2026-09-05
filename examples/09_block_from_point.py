"""**A picked point that becomes a BLOCK instead of a joint.**

A point feature becomes a POLE: a vertex where the mesh lines gather. A block
system cannot build that -- blocks meet along faces, never in a point. This
example is the other half of the pair: a point that becomes the CENTRE of a
non-quad FACE with quads all round it.

There are two ways to get one, they charge in different currency, and the whole
point of this file is to put their prices side by side.

    COLLAPSE   ``pole_blocks``   drop the pole and the fan of triangles round
                                 it, close the hole with ONE face. Lands
                                 EXACTLY on the picked point, costs zero seams
                                 -- but every ring vertex loses its edge to the
                                 pole and drops to valence 3, so the block
                                 arrives ringed by n new joints. This is what
                                 ``mesh_dual_conway`` does to a pole, done at
                                 one vertex instead of to the whole plate.

    TRUNCATE   ``block_points``  cut the corner off a genuine valence-n JOINT.
                                 The ring stays regular and the singularity is
                                 properly TRADED for the face -- but it needs a
                                 joint rather than a pole to start from, and
                                 each of the n repair splits runs to a wall as
                                 a seam.

WHY A POINT CANNOT SIMPLY BE TOLD "BE A PENTAGON"
--------------------------------------------------

``index_sum`` is conserved: a valence-5 joint and a pentagonal face both carry
-1, which is why one can be traded for the other and why neither can be made to
disappear. A full POLE carries +4, and no single face carries +4 -- that would
be a face of degree 0. So a pole can never become one clean block; its +4 has
to be spread somewhere, and in the collapse route it is spread over the ring.

Measured here, on a 10-unit square with one centre point (route 1):

    strip density   block degree   joints after
              1          4 (a quad, not a block)          0
              2          8                               12
              3         12                               16
              5         20                               24

The block's degree follows the fan, i.e. 4 x the strip density. That is the one
handle this route gives you, and at density 1 the "block" is just a quad.

Route 3 shows the truncation route doing the clean trade on an L-shaped plate,
where the reentrant corner gives the layout a real valence-5 joint to spend.

Run with the ``singular312`` env. ``--no-view`` for the numbers only.
"""

import sys

from compas.geometry import Line, Point, Translation, Vector

from compas_singular import blocks as B
from compas.tolerance import TOL
from compas_singular.algorithms import SkeletonDecomposition, boundary_triangulation


# -----------------------------------------------------------------------------
# parameters
# -----------------------------------------------------------------------------

SIZE = 10.0                     # side of the square in route 1
POINT = [0.0, 0.0, 0.0]         # the picked point -- the block's centre
RESOLUTION = 20                 # boundary discretisation per side
DENSITY = 3                     # strips per coarse patch; sets the block degree
RATIO = 0.45                    # block size for the truncation route
RELAX = 20                      # smoothing passes once the block is an element

VIEW = '--no-view' not in sys.argv


# -----------------------------------------------------------------------------
# helpers -- the same discretisation helper the other examples in this folder
# use, duplicated rather than shared, per this folder's convention
# -----------------------------------------------------------------------------

def densify(pts, resolution):
    """Discretise a closed outline side by side, so hard corners are kept."""
    lines = [Line(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
    out = []
    for line in lines:
        out.extend(line.to_polyline(n=resolution).points)
    return [[p.x, p.y, p.z] for p in out]


def decompose(outline, resolution, poles=()):
    """Outline -> Delaunay -> skeleton -> coarse pseudo-quad mesh."""
    outer = densify(outline, resolution)
    trimesh = boundary_triangulation(outer, [], [], list(poles))
    decomposition = SkeletonDecomposition.from_mesh(trimesh)
    return decomposition.decomposition_mesh(list(poles))


def to_dense(coarse, density):
    coarse.collect_strips()
    coarse.set_strips_density(density)
    coarse.densification()
    return coarse.get_quad_mesh()


def describe(name, mesh):
    lo, hi, asp = B.quality(mesh)
    print('  {:<26} F {:>5}  V {:>5}  faces {:<22} joints {:<3} index {:>3}  '
          'angles {:.1f}-{:.1f} deg'.format(
              name, mesh.number_of_faces(), mesh.number_of_vertices(),
              str(B.face_degrees(mesh)), len(B.joints(mesh)),
              B.index_sum(mesh), lo, hi))


def outline_deviation(before, after):
    """Max distance from a relaxed boundary vertex to the outline it started on.

    The honest check on sliding: NOT the area, which also moves when relaxation
    folds a face, but whether the boundary is still where it was.
    """
    from compas.geometry import closest_point_on_polyline, distance_point_point

    outlines = []
    for ring in before.vertices_on_boundaries():
        pts = [before.vertex_coordinates(v) for v in ring]
        outlines.append(pts if pts[0] == pts[-1] else pts + pts[:1])
    worst = 0.0
    for ring in after.vertices_on_boundaries():
        for v in set(ring):
            p = after.vertex_coordinates(v)
            worst = max(worst, min(
                distance_point_point(p, closest_point_on_polyline(p, o))
                for o in outlines))
    return worst


def valence_at(mesh, point):
    gk = TOL.geometric_key(point)
    for v in mesh.vertices():
        if TOL.geometric_key(mesh.vertex_coordinates(v)) == gk:
            return mesh.vertex_degree(v)
    return None


# -----------------------------------------------------------------------------
# geometry
# -----------------------------------------------------------------------------

h = SIZE / 2.0
square = [Point(h, h, 0), Point(-h, h, 0), Point(-h, -h, 0), Point(h, -h, 0)]
square = square + [square[0]]

ell = [Point(0, 0, 0), Point(10, 0, 0), Point(10, 4, 0),
       Point(4, 4, 0), Point(4, 9, 0), Point(0, 9, 0)]
ell = ell + [ell[0]]

L_PICK = [2.42, 2.22, 0.0]      # near, not exactly on, the reentrant joint


# -----------------------------------------------------------------------------
# ROUTE 1 -- the pole, as it is today
# -----------------------------------------------------------------------------

print('\nROUTE 1 -- point_features -> POLE (what happens today)')

coarse_pole = decompose(square, RESOLUTION, poles=[POINT])
print('  coarse faces {} -- the decomposition puts {} TRIANGLES on an interior'
      ' point'.format(B.face_degrees(coarse_pole),
                      B.face_degrees(coarse_pole).get(3, 0)))
dense_pole = to_dense(coarse_pole, DENSITY)
describe('pole', dense_pole)
print('  the point is a vertex of valence {} -- lines gathering, not a face'
      .format(valence_at(dense_pole, POINT)))

spin, why = B.blockable(dense_pole, [
    v for v in dense_pole.vertices()
    if TOL.geometric_key(dense_pole.vertex_coordinates(v)) == TOL.geometric_key(POINT)][0])
print('  truncating it: REFUSED -- {}'.format(why))
print('  and it could not work anyway: a full pole carries index +4, and no '
      'single face carries +4.')


# -----------------------------------------------------------------------------
# ROUTE 2 -- collapse the pole's fan into one face
# -----------------------------------------------------------------------------

print('\nROUTE 2 -- pole_blocks: collapse the fan, exact placement, no seams')

dense_collapse, rep2 = B.pole_blocks(dense_pole, [POINT])
for line in B.report_lines(rep2):
    print('  ' + line)
describe('collapsed', dense_collapse)

print('\n  the degree follows the fan, i.e. the strip density:')
for density in (1, 2, 3, 5):
    trial = to_dense(decompose(square, RESOLUTION, poles=[POINT]), density)
    out, rep = B.pole_blocks(trial, [POINT])
    degree = rep['blocks'][0]['degree']
    print('    density {}  ->  block degree {:>2}   joints after {:>2}{}'.format(
        density, degree, len(B.joints(out)),
        '   (a quad -- not a block)' if degree == 4 else ''))

# The collapsed ring is a near-circle of points that were a fan a moment ago,
# so some of its corners come out almost straight. Relaxing is what turns the
# block from a hole that was closed into an element.
collapsed_relaxed = B.relax(dense_collapse.copy(), RELAX, slide=True)
describe('collapsed, relaxed', collapsed_relaxed)


# -----------------------------------------------------------------------------
# ROUTE 3 -- truncate a real joint, the clean trade
# -----------------------------------------------------------------------------

print('\nROUTE 3 -- block_points: truncate a valence-n JOINT on an L-plate')

coarse_l = decompose(ell, 6)
dense_l = to_dense(coarse_l, 3)
describe('L-plate primal', dense_l)
for j in B.joints(dense_l):
    x, y, _ = dense_l.vertex_coordinates(j)
    print('    joint valence {} at ({:6.3f},{:6.3f}) : {}'.format(
        dense_l.vertex_degree(j), x, y, B.blockable(dense_l, j)[1]))

blocked_l, rep3 = B.block_points(dense_l, [L_PICK], ratio=RATIO)
for line in B.report_lines(rep3):
    print('  ' + line)
describe('L-plate blocked', blocked_l)

seams = B.seam_edges(dense_l, blocked_l)
print('  {} seam edges -- what the clean trade cost'.format(len(seams)))
# Pinning every boundary vertex keeps the outline but freezes the spacing of
# the boundary row, so the seams cannot pull it into line. Sliding lets a
# boundary vertex move to its neighbours' centroid and then projects it back
# onto the outline, keeping only the tangential part of the move -- same
# outline, to the same tolerance, but the row is free to even out.
pinned_relaxed = B.relax(blocked_l.copy(), RELAX)
blocked_relaxed = B.relax(blocked_l.copy(), RELAX, slide=True)
describe('L-plate relaxed, pinned', pinned_relaxed)
describe('L-plate relaxed, sliding', blocked_relaxed)
print('  boundary left the outline by {:.9f} pinned, {:.9f} sliding -- sliding '
      'moves ALONG it, so it gives up nothing'.format(
          outline_deviation(blocked_l, pinned_relaxed),
          outline_deviation(blocked_l, blocked_relaxed)))


# -----------------------------------------------------------------------------
# the comparison this file exists for
# -----------------------------------------------------------------------------

print('\nTHE TRADE, SIDE BY SIDE')
print('  collapse   block exactly on the point, 0 seams, degree {} -- but {:+d} '
      'joints on the ring'.format(rep2['blocks'][0]['degree'],
                                  rep2['ring_joints']))
print('  truncate   block {:.3f} from the point, {} seams, degree {} -- ring '
      'stays regular, {} joints -> {}'.format(
          rep3['blocks'][0]['offset'], rep3['strips'],
          rep3['blocks'][0]['degree'], len(B.joints(dense_l)),
          len(B.joints(blocked_l))))
print('  index conserved on both: {} -> {} and {} -> {}\n'.format(
    rep2['index_before'], rep2['index_after'],
    rep3['index_before'], rep3['index_after']))


# -----------------------------------------------------------------------------
# view -- row 1 the square (pole / collapsed / relaxed),
#         row 2 the L-plate (primal / blocked / relaxed)
# -----------------------------------------------------------------------------

if VIEW:
    from compas_viewer import Viewer

    viewer = Viewer()
    pitch = 1.4 * SIZE

    # Look straight down on the two rows. The default oblique view puts the
    # panels behind one another, and the whole point of the layout is to read
    # them side by side.
    camera = viewer.renderer.camera
    camera.rotation = Vector(0.0, 0.0, 0.0)          # straight down -Z
    camera.target = Vector(pitch, -0.55 * SIZE, 0.0)
    camera.position = Vector(pitch, -0.55 * SIZE, 4.8 * SIZE)

    row1 = [dense_pole, dense_collapse, collapsed_relaxed]
    for i, mesh in enumerate(row1):
        shown = mesh.copy()
        shown.transform(Translation.from_vector([i * pitch, 0, 0]))
        viewer.scene.add(shown, show_points=True, show_lines=True,
                         show_faces=True)

    row2 = [dense_l, blocked_l, blocked_relaxed]
    for i, mesh in enumerate(row2):
        shown = mesh.copy()
        shown.transform(Translation.from_vector([i * pitch - h, -1.6 * SIZE, 0]))
        viewer.scene.add(shown, show_points=True, show_lines=True,
                         show_faces=True)

    viewer.show()
