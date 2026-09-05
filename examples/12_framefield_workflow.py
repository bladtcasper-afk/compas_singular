"""THE COMPLETE FRAME-FIELD WORKFLOW -- the plug-in's eight steps, in plain Python.

The twin of ``11_skeleton_workflow.py``. Same eight Rhino commands, same order,
same result type -- only section 3 differs. Where the skeleton route reads a
medial axis off a Delaunay triangulation, this one solves a CROSS FIELD on a
background mesh, traces its separatrices, and cuts the domain into patches along
them. Everything from section 4 down is shared between the two routes.

    Rhino command                    section here
    ------------------------------   -----------------------------------------
    1  CMD_start                      SETTINGS
    2  CMD_boundary_selection         INPUT
    3  CMD_coarse_mesh (FrameField)   COARSE LAYOUT
    4  CMD_edit_coarse_mesh           EDIT
    5  CMD_densities                  DENSITIES
    6  CMD_quad_mesh                  QUAD MESH
    7  CMD_smoothen                   SMOOTHING
    8  CMD_dual                       DUAL

**To use it: change the numbers in SETTINGS and the curves in INPUT.** Nothing
below those two sections is specific to this domain.

Why you would pick this route over the skeleton:

* it takes GUIDES -- cables, force lines, anything the mesh should run along or
  across. A medial axis is equidistant from the walls by definition and cannot
  be steered; a separatrix follows a direction field and can. On a plain square
  a cable forces no topology at all, and the field still reaches the mesh
  through the patch interiors (section 6);
* it keeps the SYMMETRY of the input. ``symmetry='auto'`` detects the group of
  outer boundary + holes + guides and holds it all the way to the layout;
* measured over the 18-domain baseline it wins 15 to 1.

And where it does not: a domain of uniform curvature loses quality at the patch
CORNERS, and a sampled arc that turns more than 45 degrees between two points IS
a corner to a cross field -- an oval hole sampled too coarsely sprays
singularities. ``SPACING`` is the knob for both.

Six things are load-bearing, and each is marked WHY where it happens:

* ``SPACING`` is the BACKGROUND, not the quad size (1);
* ``relax=True`` is what you want with guides (3);
* ``decomposition.route()`` and ``.warnings()`` have to be read (3/6) -- the
  front end always returns a mesh, including when it gave up on the field;
* the layout is densified by :meth:`densify`, not by ``densification()`` (6) --
  that is what puts the field inside the patches as well as on their edges;
* ``smooth_quad_mesh`` only keeps a result that improves min angle, max angle
  and aspect together (7);
* the dual wants a plain COMPAS mesh (8).

Run:  python 12_framefield_workflow.py
"""
from __future__ import print_function

import math

from compas.datastructures import Mesh
from compas.geometry import Line

from compas_singular.framefield.decomposition import FieldDecomposition
from compas_singular.rhino.dual_mesh import dual_mesh


# =============================================================================
# 1  SETTINGS                                                    (CMD_start)
# =============================================================================

# The BACKGROUND triangulation spacing -- the mesh the field is solved on, not
# the quad size. Finer means a better field and a slower solve.
SPACING = 0.5

TARGET_LENGTH = 0.5

# 'tangent' makes elements run ALONG the guides, 'perpendicular' across.
MODE = 'tangent'

# Diffusion + normalisation instead of a single Dirichlet solve. The default
# solver is measurably wrong about where singularities sit (disc radius 0.45
# against a published 0.85); with guides, use this.
RELAX = True


# =============================================================================
# 2  INPUT                                                (CMD_boundary_selection)
# =============================================================================

OUTER = [[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]]

HOLE_CENTRE, HOLE_RADIUS = [2.0, 0.0, 0.0], 1.6

# Guides: cables, force lines, anything the mesh should follow. A LIST of
# curves, like the holes -- wrap a lone guide in a list.
GUIDES = [[[-6.0, -2.0, 0.0], [-1.0, -1.0, 0.0], [3.0, 2.5, 0.0], [6.0, 3.0, 0.0]]]

POINT_FEATURES = []


def densify_loop(points, spacing):
    """Resample a closed loop so no segment is longer than ``spacing``."""
    loop = list(points)
    dense = []
    for a, b in zip(loop, loop[1:] + loop[:1]):
        line = Line(a, b)
        n = max(1, int(round(line.length / spacing)))
        dense.extend([list(point) for point in line.to_polyline(n=n).points[:-1]])
    return dense


def sample_circle(centre, radius, spacing):
    n = max(8, int(round(2.0 * math.pi * radius / spacing)))
    return [[centre[0] + radius * math.cos(2.0 * math.pi * i / n),
             centre[1] + radius * math.sin(2.0 * math.pi * i / n), 0.0]
            for i in range(n)]


# Sampled FINER than the background. A sampled arc that turns more than 45
# degrees between two points reads as a corner to a cross field, and every
# corner launches a separatrix. The background is rebuilt from these anyway.
outer = densify_loop(OUTER, SPACING)
holes = [sample_circle(HOLE_CENTRE, HOLE_RADIUS, SPACING * 0.5)]


# =============================================================================
# 3  COARSE LAYOUT                            (CMD_coarse_mesh, "FrameField")
# =============================================================================

# Solves the field and traces its separatrices. This is the expensive call;
# ``CMD_start.get_decomposition`` caches it across steps 3, 4 and 6.
decomposition = FieldDecomposition.from_boundary(
    outer,
    inner_boundaries=holes,
    guides=GUIDES,
    mode=MODE,
    target_length=SPACING,
    relax=RELAX,
    symmetry='auto')

# The separatrices and the walls, as polylines -- the field route's answer to
# the skeleton's branches. ``decomposition_mesh`` runs this itself; called first
# only so the polylines are available before the layout is built.
polylines = decomposition.decomposition_polylines()

coarse = decomposition.coarse_mesh(POINT_FEATURES)

# The field travels on its own: it carries its own background and builds its own
# point locator, so it can densify a layout that did not come from it -- see
# ``11_skeleton_workflow.py``. Equally, leave it out of section 6 and this
# layout densifies with plain Coons interiors.
field = decomposition.get_field()

# Read these. The front end always returns a layout, including when the
# separatrix network did not close and it fell back to covering the shape.
print('field   ', decomposition.field.report())
print('symmetry', None if decomposition.symmetry is None else decomposition.symmetry.names())
for warning in decomposition.warnings():
    print('warning ', warning)


# =============================================================================
# 4  EDIT                                          (CMD_edit_coarse_mesh)
# =============================================================================

"""
The editor is not very intuitive to use in code yet.

On this route it does have a working commit, unlike the skeleton one:

    editor = CoarseLayoutEditor(decomposition)     # move_vertex / insert_curve
    ok, notes = editor.commit()                    # re-warps the separatrices

and a layout edited elsewhere goes back in through ``quad_mesh(coarse=...)``.
"""

# =============================================================================
# 5  DENSITIES                                             (CMD_densities)
# =============================================================================

coarse.collect_strips()

coarse.set_strips_density_target(TARGET_LENGTH)

# Override individual strips on top of that. These survive into section 6 --
# both ``densify`` and a no-size ``quad_mesh()`` read them rather than reset
# them. Strip KEYS renumber whenever the layout is rebuilt, so anything that
# outlives one run should key its densities geometrically and resolve them just
# before the call; ``agent.core.rebuild`` does that on top of ``strip_address``.
skey = sorted(coarse.strips())[0]
coarse.set_strip_density(skey, coarse.get_strip_density(skey) * 2)

# =============================================================================
# 6  QUAD MESH                                             (CMD_quad_mesh)
# =============================================================================

# ``field=`` is what makes this different from a bare ``densification()``. That
# fills each patch with a bilinear blend of its own four sides and never consults
# the field; with the field, the patch INTERIOR is integrated from it too -- the
# only way a guide reaches a patch it forced no topology in. ``edges_to_curves``
# gives the coarse EDGES the curvature of the separatrices and walls they came
# from; the two are independent and both are worth having.
edges_to_curves = decomposition.edges_to_curves()

dense = coarse.densification(edges_to_curves=edges_to_curves, field=field)

raw = dense.copy()                                   # kept for the viewer only

# ``decomposition.quad_mesh()`` is sections 5 + 6 in one call, plus two
# fallbacks for a layout whose separatrix network did not close. Use it in
# anything that has to return a mesh no matter what; then read ``route()``. Its
# density rule: a size argument wins everywhere, no size argument keeps what is
# already set.
#
#     decomposition.quad_mesh()                       # keeps section 5
#     decomposition.quad_mesh(target_length=0.3)      # 0.3 on every strip
#     decomposition.quad_mesh(densities={skey: 8})    # those strips, then base
print('route   ', decomposition.route())
print('quality ', decomposition.quality(dense))

# =============================================================================
# 7  SMOOTHING                                              (CMD_smoothen)
# =============================================================================

# ``densify`` solves each patch with its boundary held fixed, so the seams are
# out of its reach at any stiffness; this is the pass that reaches them. It keeps
# the result only if min angle, max angle and aspect all improve, so a mesh with
# nothing to gain comes back untouched. In place.
#
# ``datastructures.mesh.smoothing`` is the route-independent alternative --
# ``boundary_constrained_smoothing`` there is what ``11_skeleton_workflow.py``
# uses, and it takes the walls as explicit curves.
report = decomposition.smooth_quad_mesh(dense)

print('relaxed ', report['accepted'], report['after'])

# =============================================================================
# 8  DUAL                                                       (CMD_dual)
# =============================================================================

plain = Mesh.from_vertices_and_faces(*dense.to_vertices_and_faces(keep_keys=False))
dual = dual_mesh(plain, redistribute=True)

# =============================================================================
# VIEW
# =============================================================================

from compas_viewer import Viewer

from compas_singular.framefield.viz import Grid
from compas_singular.framefield.viz import add_background
from compas_singular.framefield.viz import add_dense
from compas_singular.framefield.viz import add_field
from compas_singular.framefield.viz import add_guides
from compas_singular.framefield.viz import add_layout
from compas_singular.framefield.viz import add_singularities

viewer = Viewer()
grid = Grid(pitch=16.0, cols=5)

dx, dy = grid.cell(0)
group = viewer.scene.add_group(name='3  cross field')
add_background(group, decomposition, dx, dy)
add_field(group, decomposition, dx, dy)
add_singularities(group, decomposition, dx, dy)
add_guides(group, GUIDES, dx, dy)

dx, dy = grid.cell(1)
group = viewer.scene.add_group(name='3/5  coarse layout')
add_layout(group, decomposition, dx, dy)

dx, dy = grid.cell(2)
group = viewer.scene.add_group(name='6  quad mesh')
add_dense(group, raw, dx, dy)

dx, dy = grid.cell(3)
group = viewer.scene.add_group(name='7  smoothed')
add_dense(group, dense, dx, dy)

dx, dy = grid.cell(4)
group = viewer.scene.add_group(name='8  dual')
add_dense(group, dual, dx, dy, name='blocks')

grid.frame(viewer)
viewer.show()