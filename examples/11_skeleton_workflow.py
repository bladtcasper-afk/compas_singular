from __future__ import print_function

import math

from compas.datastructures import Mesh
from compas.geometry import Polyline

from compas_singular.algorithms import SkeletonDecomposition
from compas_singular.datastructures.mesh.smoothing import boundary_constrained_smoothing
from compas_singular.editing import CoarseLayoutEditor
from compas_singular.framefield.quality import mesh_quality
from compas_singular.rhino.dual_mesh import dual_mesh


# =============================================================================
# 1  SETTINGS                                                    (CMD_start)
# =============================================================================


SPACING = 0.5

WALL_SAMPLING = SPACING * 0.25

TARGET_LENGTH = 0.5


# =============================================================================
# 2  INPUT                                                (CMD_boundary_selection)
# =============================================================================

OUTER = [[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]]

HOLE_CENTRE, HOLE_RADIUS = [2.0, 0.0, 0.0], 1.6

POINT_FEATURES = [[-4.0,0.0,0.0]]

def sample_circle(centre, radius, spacing):
    """A circle as a point list. Sample it as finely as the mesh should follow it
    -- these points ARE the hole from here on, so an arc given as two is a chord.
    Straight walls need no such treatment: ``from_boundary`` subdivides them."""
    n = max(8, int(round(2.0 * math.pi * radius / spacing)))
    return [[centre[0] + radius * math.cos(2.0 * math.pi * i / n),
             centre[1] + radius * math.sin(2.0 * math.pi * i / n), 0.0]
            for i in range(n)]


HOLES = [sample_circle(HOLE_CENTRE, HOLE_RADIUS, WALL_SAMPLING)]


# =============================================================================
# 3  COARSE LAYOUT                              (CMD_coarse_mesh, "Skeleton")
# =============================================================================

# ``target_length`` is the BACKGROUND spacing, not the quad size: the medial axis
# is read off a Delaunay triangulation of the boundary POINTS, so the walls are
# resampled to it first. A four-point square handed in raw gives a two-triangle
# mesh and a caricature layout.
decomposition = SkeletonDecomposition.from_boundary(
    OUTER, inner_boundaries=HOLES, point_features=POINT_FEATURES,
    target_length=SPACING)

# The medial axis. Read BEFORE the layout: ``coarse_mesh`` runs the
# ``branches_splitting_*`` fix-ups, which insert vertices into the triangulation.
branches = decomposition.branches()
singular_points = decomposition.singular_points()

coarse = decomposition.coarse_mesh()

polylines = list(decomposition.polylines)


# =============================================================================
# 4  EDIT                                          (CMD_edit_coarse_mesh)
# =============================================================================

"""
The editor is not very intuitive to use in code yet.
"""

# =============================================================================
# 5  DENSITIES                                             (CMD_densities)
# =============================================================================

coarse.collect_strips()

coarse.set_strips_density_target(TARGET_LENGTH)

# =============================================================================
# 6  QUAD MESH                                             (CMD_quad_mesh)
# =============================================================================

# A coarse edge is a straight chord and has to be, so the SHAPE of each edge is
# handed in separately: the wall arc between its two corners, the branch polyline
# whose ends match it, or a chord. Without it the round hole comes out a polygon
# -- 54.6% of the radius off, against 0.08% with it.
edges_to_curves, tally = decomposition.edges_to_curves()

dense = coarse.densification(edges_to_curves=edges_to_curves)

raw = dense.copy()                                   # kept for the viewer only

# =============================================================================
# 7  SMOOTHING
# =============================================================================

walls = [decomposition.inputs['outer_boundary']] + decomposition.inputs['inner_boundaries']

constraints = boundary_constrained_smoothing(
    dense, curves=[Polyline(loop + loop[:1]) for loop in walls],
    kmax=50, damping=0.5, algorithm='area')

# =============================================================================
# 8  DUAL
# =============================================================================

plain = Mesh.from_vertices_and_faces(*dense.to_vertices_and_faces(keep_keys=False))
dual = dual_mesh(plain, redistribute=True)

# =============================================================================
# VIEW
# =============================================================================

from compas.colors import Color
from compas.geometry import Point
from compas_viewer import Viewer

PITCH = 15.0

def shifted(mesh, dx):
    clone = mesh.copy()
    for vkey in clone.vertices():
        x, y, z = clone.vertex_coordinates(vkey)
        clone.vertex_attributes(vkey, 'xyz', [x + dx, y, z])
    return clone

viewer = Viewer()

group = viewer.scene.add_group(name='3  skeleton')
# the decomposition IS the Delaunay mesh -- SkeletonDecomposition is a Mesh
group.add(shifted(decomposition, 0.0), show_points=False, opacity=0.25,
            linecolor=Color(0.6, 0.6, 0.6), name='Delaunay background')
for branch in branches:
    group.add(Polyline([[p[0], p[1], 0.02] for p in branch]),
                linewidth=3, linecolor=Color.red(), name='skeleton branch')
for p in singular_points:
    group.add(Point(p[0], p[1], 0.04), pointsize=12,
                pointcolor=Color.blue(), name='singular point')

group = viewer.scene.add_group(name='4/5  coarse layout')
group.add(shifted(coarse, PITCH), show_points=True, opacity=0.5,
            name='layout')
for polyline in polylines:
    group.add(Polyline([[p[0] + PITCH, p[1], 0.02] for p in polyline]),
                linecolor=Color.green(), name='decomposition polyline')

viewer.scene.add(shifted(raw, 2 * PITCH), show_points=False, name='6  quad mesh')
viewer.scene.add(shifted(dense, 3 * PITCH), show_points=False, name='7  smoothed')
viewer.scene.add(shifted(dual, 4 * PITCH), show_points=False, name='8  dual')

viewer.show()