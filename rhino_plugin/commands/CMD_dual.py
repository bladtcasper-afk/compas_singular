#! python3
# r: compas
# r: pydantic

from collections import Counter

import rhinoscriptsyntax as rs
import compas_rhino as cr

from compas.scene import Scene
from compas_singular.rhino.dual_mesh import dual_mesh
from compas_singular.rhino.helpers import read_mesh, bake_mesh, clear_layer
from compas.datastructures import mesh_conway_dual
from compas_singular.rhino.helpers import mesh_from_rhino

if not rs.IsLayer("QuadMesh"):
    raise RuntimeError("No dense quad mesh has been generated yet. The layer does not exist yet.")

def quad_mesh_filter(rhobj, geometry, component_index):
    layer = rs.ObjectLayer(rhobj)
    quad_layer = rs.LayerName("QuadMesh", fullpath=True)
    return layer == quad_layer or rs.IsLayerChildOf(quad_layer, layer)

mesh_id = rs.GetObject(message="Pick a dense quad mesh to smoothen from the layer QuadMesh and its sublayers", filter=rs.filter.mesh, preselect=False, select=False, custom_filter=quad_mesh_filter, subobjects=False)
if not mesh_id:
    print("Cancelled -- nothing was drawn or changed.")
mesh = mesh_from_rhino(cr.objects.find_object(mesh_id).Geometry)
dual = dual_mesh(mesh, redistribute=True)

#scene = Scene()
#scene.add(mesh.dual(include_boundary=True))
#scene.draw()

clear_layer("Dual", clean_sublayers=True)

rs.AddLayer("Dual", parent="QuadMesh")
bake_mesh(dual, "Dual")

print("dual: {} faces from {} primal vertices, baked to 'QuadMesh::Dual'".format(
    dual.number_of_faces(), mesh.number_of_vertices()))
print("next: this is the end of the pipeline for this mesh -- re-run CMD_quad_mesh "
      "or edit in CMD_edit_quad_mesh if further changes are needed, then dual again.")


"""
def spread(m):
    areas = [m.face_area(f) for f in m.faces()]
    return max(areas) / min(areas)


# a stray non-quad here is a singularity, which is expected -- but a boundary
# ring count above 2 means the input was UNWELDED and every crease read as an
# edge of the mesh.
degrees = dict(Counter(len(dual.face_vertices(f)) for f in dual.faces()))
rings = len(dual.vertices_on_boundaries())
print("dual: {} faces from {} primal vertices, degrees {}, {} boundary ring(s)".format(
    dual.number_of_faces(), mesh.number_of_vertices(), degrees, rings))
print("block spread max/min: {:.2f} raw -> {:.2f} redistributed  (primal {:.2f})".format(
    spread(raw), spread(dual), spread(mesh)))


conway_dual = mesh_conway_dual(mesh)
raw = dual_mesh(mesh, redistribute=False)
scene = Scene()
scene.clear()
scene.add(mesh.dual(include_boundary=True))
scene.add(conway_dual)
scene.add(raw)
scene.draw()
"""