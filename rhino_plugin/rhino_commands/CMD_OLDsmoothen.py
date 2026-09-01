#! python3
# r:compas

from CMD_start import get_settings, import_compas_singular
import_compas_singular()
from compas.scene import Scene

from compas.datastructures.mesh.smoothing import mesh_smooth_centroid, mesh_smooth_area, mesh_smooth_centerofmass
from compas.datastructures.mesh.conway import *
from compas_singular.rhino.helpers.helpers import read_mesh, bake_mesh, clean_layer
from compas_rhino.conversions import mesh_to_compas

import rhinoscriptsyntax as rs
import compas_rhino as cr

if not rs.IsLayer("QuadMesh"):
    raise RuntimeError("No dense quad mesh has been generatedy yet. The layer does not exist yet.")

def quad_mesh_filter(rhobj, geometry, component_index):
    layer = rs.ObjectLayer(rhobj)
    quad_layer = rs.LayerName("QuadMesh", fullpath=True)
    return layer == quad_layer or rs.IsLayerChildOf(quad_layer, layer)

mesh = rs.GetObject(message="Pick a dense quad mesh to smoothen from the layer QuadMesh and its sublayers", filter=rs.filter.mesh, preselect=False, select=False, custom_filter=quad_mesh_filter, subobjects=False)
if mesh is None:
    raise RuntimeError("Pick a dense quad mesh to smoothen.")
mesh = mesh_to_compas(cr.objects.find_object(mesh).Geometry)


rs.AddLayer("Smoothened", parent="QuadMesh")
clean_layer("Smoothened", clean_sublayers=True)

mode = rs.GetString(message="Use the guides and points features as constraints.", defaultString="No", strings=["Yes", "No"])
constrains = {}
if mode=="yes":
    pass

mode = rs.GetString(message="Select smoothing method", defaultString="Area", strings=["Area", "Centroid", "CenterOfMass"])
mode = (mode or "Exit").lower()
if mesh is None:
    mode = "Exit".lower()

boundary = [v for ring in mesh.vertices_on_boundaries() for v in ring]

if mode=="Area".lower():
    smesh_area = mesh.copy()
    mesh_smooth_area(smesh_area, fixed=boundary, kmax=100, damping=0.5) 
    rs.AddLayer("Area", parent="Smoothened")
    smesh_area = bake_mesh(smesh_area, "Area")
elif mode=="Centroid".lower():
    rs.AddLayer("Centroid", parent="Smoothened")
    smesh_centroid = mesh.copy()
    mesh_smooth_centroid(smesh_centroid, fixed=boundary, kmax=100, damping=0.5)
    smesh_centroid = bake_mesh(smesh_centroid, "Centroid")
elif mode=="CenterOfMass".lower():
    rs.AddLayer("CenterOfMass", parent="Smoothened")
    smesh_center = mesh.copy()
    mesh_smooth_centerofmass(smesh_center, fixed=boundary, kmax=100, damping=0.5)
    smesh_centermass = bake_mesh(smesh_center, "CenterOfMass")
elif mode=="All".lower():
    rs.AddLayer("Centroid", parent="Smoothened")
    smesh_centroid = mesh.copy()
    mesh_smooth_centroid(smesh_centroid, fixed=boundary, kmax=100, damping=0.5)
    rs.AddLayer("Area", parent="Smoothened")
    smesh_area = mesh.copy()
    mesh_smooth_area(smesh_area, fixed=boundary, kmax=100, damping=0.5) 
    rs.AddLayer("CenterOfMass", parent="Smoothened")
    smesh_center = mesh.copy()
    mesh_smooth_centerofmass(smesh_center, fixed=boundary, kmax=100, damping=0.5)

    smesh_centroid = bake_mesh(smesh_centroid, "Centroid")
    smesh_area = bake_mesh(smesh_area, "Area")
    smesh_center = bake_mesh(smesh_center, "CenterOfMass")
else:
    pass

#rs.AddLayer("ForceDensity", parent="Smoothened")


"""conway_meshes = []
conway_meshes.append(mesh_conway_ambo(mesh))
conway_meshes.append(mesh_conway_bevel(mesh))
conway_meshes.append(mesh_conway_dual(mesh))
conway_meshes.append(mesh_conway_expand(mesh))
conway_meshes.append(mesh_conway_gyro(mesh))
conway_meshes.append(mesh_conway_join(mesh))
conway_meshes.append(mesh_conway_kis(mesh))
conway_meshes.append(mesh_conway_meta(mesh))
conway_meshes.append(mesh_conway_needle(mesh))
conway_meshes.append(mesh_conway_ortho(mesh))
conway_meshes.append(mesh_conway_zip(mesh))
conway_meshes.append(mesh_conway_truncate(mesh))
conway_meshes.append(mesh_conway_snub(mesh))


scene = Scene()
scene.clear()
for mesh in conway_meshes:
    scene.add(mesh)


scene.draw()"""
