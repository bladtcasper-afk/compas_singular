import os

from compas_singular.datastructures import CoarsePseudoQuadMesh
from compas_viewer import Viewer

HERE = os.path.dirname(__file__)
FILE = os.path.join(HERE, 'data/coarse_quad_mesh_british_museum_poles.json')

# read input data
coarse_pseudo_quad_mesh = CoarsePseudoQuadMesh.from_json(FILE)

# view the coarse pseudo quad mesh
viewer = Viewer()
viewer.scene.add(coarse_pseudo_quad_mesh, show_points=True, show_lines=True, show_faces=True)
viewer.show()

# collect strip data
coarse_pseudo_quad_mesh.collect_strips()

# densification with target length
coarse_pseudo_quad_mesh.set_strips_density_target(t=.5)
coarse_pseudo_quad_mesh.densification()

# view the dense quad mesh
viewer = Viewer()
viewer.scene.add(coarse_pseudo_quad_mesh.get_quad_mesh(), show_points=True, show_lines=True, show_faces=True)
viewer.show()
