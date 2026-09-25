from compas_tna.diagrams import ForceDiagram, FormDiagram
from compas_tna.equilibrium import horizontal_nodal, vertical_from_zmax, horizontal_nodal_numpy, horizontal_numpy
from compas_tna.envelope import DomeEnvelope

form = FormDiagram.create_circular_spiral()
form.update_boundaries()
form.assign_support_type("all")
force = ForceDiagram.from_formdiagram(form)
eq_h = horizontal_nodal(form, force)
eq_v = vertical_from_zmax(form, zmax=5.0)

def tna_analysis(mesh):
    form = FormDiagram.create_circular_spiral()
    form.update_boundaries()
    form.assign_support_type("all")
    force = ForceDiagram.from_formdiagram(form)
    eq_h = horizontal_nodal(form, force)
    eq_v = vertical_from_zmax(form, zmax=5.0)
    return eq_v[0]