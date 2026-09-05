#! python3

# r: compas

"""**Step 4 -- drag the layout's corners.** Optional, and the point of the tool.

    reads   TopologyProblem::Skeleton::Mesh, ::Poles, InputBoundaries::*
    writes  TopologyProblem::Skeleton::Mesh, ::Poles
    scratch TopologyProblem::Skeleton::Edit    draggable corners, cleared on exit

Type ``move`` to drag corners, ``commit`` to keep the result, ``reset`` to put
them all back, ``finish`` to leave without committing.

**The layout is edited as a mesh held in Python, not as the polylines.**
Polylines duplicate every shared corner -- an interior corner where four patches
meet is four coincident points -- so dragging one of them tears the layout, and
the 3-decimal weld turns the moved copy into a NEW vertex without complaining.
The mesh shares that corner, so moving it moves all four patches, which is what
the user meant.

**Committing goes through** ``FieldDecomposition.edit_coarse``, which is where
four things that have to happen already live: coincident corners are welded,
boundary corners are projected onto the actual wall, whatever is no longer
four-sided is repaired, and a layout that will not densify is REFUSED -- while
you are still in Rhino and can fix it. It also warps each traced separatrix onto
the moved corners, so a nudged corner costs the nudge and not the curvature.

Two consequences worth knowing:

* the snap is not cosmetic. ``_on_loop`` tests membership of a wall at 1e-6, so
  a corner dragged to what LOOKS like the boundary is not on it; the densifier
  then finds no boundary arc and cuts the corner off with a chord. That is the
  failure that took an ellipse to 65% coverage. Boundary corners are projected
  as you drag, so the preview shows where the corner will actually land;
* ``edit_coarse`` RENUMBERS vertex keys -- measured 2 of 8 unchanged on a disc.
  Nothing may be keyed by index across it, which is why step 5 keys densities on
  geometry.
"""
from ff_common import EDIT_LAYER
from ff_common import MESH_LAYER
from ff_common import bake_coarse
from ff_common import ensure_layer
from ff_common import get_decomposition
from ff_common import read_boundaries
from ff_common import read_coarse

import rhinoscriptsyntax as rs

from rhino_plugin.coarse_edit import CoarseLayoutObject


outer, inners, guides = read_boundaries()
coarse, poles = read_coarse()
print("layout: {} patch(es), {} corner(s), {} pole(s)".format(
    coarse.number_of_faces(), coarse.number_of_vertices(), len(poles)))

# Needed to commit: ``edit_coarse`` warps the separatrices onto the moved
# corners and snaps to the domain walls, and only the decomposition has either.
decomposition = get_decomposition()

ensure_layer(EDIT_LAYER, (0, 120, 200))

obj = CoarseLayoutObject(decomposition, coarse=coarse, loops=[outer] + inners)
# The draggable corners get their own layer: Skeleton::Mesh holds the mesh
# object and points sharing it would be deleted by the rebake.
obj.settings = dict(CoarseLayoutObject.settings, layer=EDIT_LAYER)

# The mesh is in the way while dragging -- it still shows the corners where they
# were. Hidden rather than deleted, so a cancelled edit leaves the document as
# it was found.
mesh_guids = rs.ObjectsByLayer(MESH_LAYER)
rs.HideObjects(mesh_guids)

# A cached decomposition still carries the notes from an edit committed earlier
# this session. Left in place, a run the user CANCELS would look committed and
# the stale layout would be baked over the current one.
decomposition.edit_notes = {}

try:
    obj.update()
finally:
    obj.clear()
    rs.ShowObjects(mesh_guids)

if not decomposition.edit_notes:
    print("nothing committed -- '{}' is unchanged.".format(MESH_LAYER))
else:
    bake_coarse(decomposition.mesh)
    notes = decomposition.edit_notes
    print("committed: {} patch(es) in, {} out, {} corner(s) snapped to a wall"
          .format(notes.get("faces_in"), notes.get("faces_out"),
                  notes.get("snapped")))
    if notes.get("off_wall"):
        print("  {} corner(s) sit on no wall -- expected only where a patch was "
              "deleted".format(notes["off_wall"]))
    # Densities stored in step 5 are keyed on edge midpoints, and this edit moved
    # some of those. Step 5 reports how many overrides survived.
    print("layout written back to '{}'".format(MESH_LAYER))

print("next: ff_05_densities")
