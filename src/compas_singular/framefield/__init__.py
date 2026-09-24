"""Frame-field front end for compas_singular.

Replaces ``SkeletonDecomposition``'s medial-axis front end with a
frame-field -> separatrix -> coarse-layout one, leaving everything downstream
(``collect_strips``, ``densification``, ``add_strip``)
untouched.

See ``frame-field-front-end-module-layout.md`` in the parent folder for the
design and the step-by-step plan this implements.
"""
from compas_singular.framefield.background import BackgroundMesh                     # noqa: F401
from compas_singular.framefield.constraints import Constraint                        # noqa: F401
from compas_singular.framefield.constraints import from_boundary                     # noqa: F401
from compas_singular.framefield.constraints import from_curves                       # noqa: F401
from compas_singular.framefield.densify import field_densification                   # noqa: F401
from compas_singular.framefield.field import CrossField                              # noqa: F401
from compas_singular.framefield.guides import guide_metrics                          # noqa: F401
from compas_singular.framefield.quality import curve_alignment                       # noqa: F401
from compas_singular.framefield.quality import hard_floor                            # noqa: F401
from compas_singular.framefield.quality import mesh_quality                          # noqa: F401
from compas_singular.framefield.relax import relax_mesh                              # noqa: F401
from compas_singular.framefield.symmetry import Symmetry                             # noqa: F401

__all__ = [
    'BackgroundMesh',
    'Constraint',
    'from_boundary',
    'from_curves',
    'CrossField',
    'guide_metrics',
    'relax_mesh',
    'field_densification',
    'mesh_quality',
    'curve_alignment',
    'hard_floor',
    'Symmetry',
]
