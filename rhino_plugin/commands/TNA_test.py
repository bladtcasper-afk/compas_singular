#! python3
# r: compas


import rhinoscriptsyntax as rs
import compas_rhino as cr
from compas_rhino.conversions import mesh_to_compas

from compas.scene import Scene

mesh_id = rs.GetObject(message="Pick a mesh for the formdiagram of the TNA.", filter=rs.filter.mesh)
mesh = mesh_to_compas(cr.objects.find_object(mesh_id).Geometry)

from CMD_start import get_settings, import_compas_singular
import_compas_singular()
from compas_singular.rhino.helpers.tna import tna_analysis

scene = Scene()

result, form, force = tna_analysis(mesh)


scene.add(result)

scene.draw()