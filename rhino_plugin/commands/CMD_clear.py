#! python3
# r: compas

#Temporary import of compas_singular development library
from CMD_start import import_compas_singular
import_compas_singular()

import rhinoscriptsyntax as rs
from compas_singular.rhino.helpers.helpers import clear_layer




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