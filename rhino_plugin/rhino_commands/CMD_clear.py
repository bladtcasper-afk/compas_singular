#! python3
# r: compas

#Temporary import of compas_singular development library
from CMD_start import import_compas_singular
import_compas_singular()

import rhinoscriptsyntax as rs
from compas_singular.rhino.helpers.helpers import clean_layer




items = ("InputBoundaries", "No", "Yes"), ("FrameField", "No", "Yes"), ("Skeleton", "No", "Yes"), ("QuadMesh", "No", "Yes"), ("TrashBin", "No", "Yes")
defaults = (False, True, True, True, True)

to_clear = rs.GetBoolean("Clear?", items, defaults)

objects_deleted_count = 0
for i in range(len(items)):
    if to_clear[i]:
        objects_deleted_count += clean_layer(items[i][0], True)


print(f"Deleted {objects_deleted_count} object(s).")