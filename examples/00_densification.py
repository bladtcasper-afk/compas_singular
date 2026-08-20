import os

from compas_singular.datastructures import CoarseQuadMesh
from compas_viewer import Viewer

HERE = os.path.dirname(__file__)
FILE = os.path.join(HERE, 'data/coarse_quad_mesh_british_museum.json')

# read input data
coarse_quad_mesh = CoarseQuadMesh.from_json(FILE)

# view the coarse quad mesh
viewer = Viewer()
viewer.scene.add(coarse_quad_mesh, show_points=True, show_lines=True, show_faces=True)
viewer.show()

# collect strip data
coarse_quad_mesh.collect_strips()

# densification with uniform density
coarse_quad_mesh.set_strips_density(3)
coarse_quad_mesh.densification()

# view the dense quad mesh
viewer = Viewer()
viewer.scene.add(coarse_quad_mesh.get_quad_mesh(), show_points=True, show_lines=True, show_faces=True)
viewer.show()

coarse_quad_mesh.set_strips_density(2)
coarse_quad_mesh.densification()

# view the dense quad mesh
viewer = Viewer()
viewer.scene.add(coarse_quad_mesh.get_quad_mesh(), show_points=True, show_lines=True, show_faces=True)
viewer.show()

# densification with target length
coarse_quad_mesh.set_strips_density_target(t=.5)
coarse_quad_mesh.densification()

# view the dense quad mesh
viewer = Viewer()
viewer.scene.add(coarse_quad_mesh.get_quad_mesh(), show_points=True, show_lines=True, show_faces=True)
viewer.show()

# change density of one strip
skey = list(coarse_quad_mesh.strips())[0]
coarse_quad_mesh.set_strip_density(skey, 10)
coarse_quad_mesh.densification()

# view the dense quad mesh
viewer = Viewer()
viewer.scene.add(coarse_quad_mesh.get_quad_mesh(), show_points=True, show_lines=True, show_faces=True)
viewer.show()
