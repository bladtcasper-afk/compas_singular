from compas_tna.diagrams import ForceDiagram, FormDiagram
from compas_tna.equilibrium import horizontal_nodal, vertical_from_zmax, horizontal_nodal_numpy, horizontal_numpy
from compas_tna.envelope import DomeEnvelope


def tna_analysis(mesh):
    form = FormDiagram.from_mesh(mesh)
    # form = FormDiagram.create_circular_radial()
    form.assign_support_type("all")
    force = ForceDiagram.from_formdiagram(form)
    try:
        eq_h = horizontal_nodal_numpy(form, force)
    except:
        print("failed")
    eq_v = vertical_from_zmax(form, zmax=5.0)
    return eq_v[0], form, force

from compas.datastructures.mesh.mesh import Mesh
Mesh.to