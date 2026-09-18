#! python3
# r: compas

# Development bootstrap -- delete once compas_singular is installed into Rhino's
# Python. MUST run before any compas_singular import: Rhino resets sys.path between
# runs but keeps sys.modules, so put the source on the path and drop a stale copy
# (see CMD_start for why every module, framefield included, has to go).
import sys
SINGULAR_SRC = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular\src"
if SINGULAR_SRC not in sys.path:
    sys.path.insert(0, SINGULAR_SRC)
if not getattr(sys, "compas_singular_keep_modules", False):  # set by headless tests
    for _mod in list(sys.modules):
        if _mod == "compas_singular" or _mod.startswith("compas_singular."):
            del sys.modules[_mod]

import rhinoscriptsyntax as rs
from compas_singular.rhino.helpers import clear_layer




items = ("InputBoundaries", "No", "Yes"), ("FrameField", "No", "Yes"), ("Skeleton", "No", "Yes"), ("QuadMesh", "No", "Yes"), ("TrashBin", "No", "Yes")
defaults = (False, True, True, True, True)

to_clear = rs.GetBoolean("Clear?", items, defaults)

objects_deleted_count = 0
cleared = []
for i in range(len(items)):
    if to_clear[i]:
        objects_deleted_count += clear_layer(items[i][0], True)
        cleared.append(items[i][0])


print(f"Deleted {objects_deleted_count} object(s) from: {', '.join(cleared) or 'nothing'}.")
print("note: clearing an earlier layer (e.g. InputBoundaries or Skeleton) leaves any "
      "later layer (Skeleton, QuadMesh) describing a domain or layout that no longer exists.")
print("next: CMD_boundary_selection to start a new domain, or CMD_coarse_mesh / "
      "CMD_read_coarse_mesh if the domain is untouched.")