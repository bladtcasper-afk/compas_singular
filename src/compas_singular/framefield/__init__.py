"""Frame-field front end for compas_singular.

Replaces ``SkeletonDecomposition``'s medial-axis front end with a
frame-field -> separatrix -> coarse-layout one, leaving everything downstream
(``collect_strips``, ``densification``, ``add_strip``, ``guide_lines``)
untouched.

See ``frame-field-front-end-module-layout.md`` in the parent folder for the
design and the step-by-step plan this implements.
"""
from .background import BackgroundMesh                     # noqa: F401
from .cache import SolveCache                              # noqa: F401
from .cache import solve                                   # noqa: F401
from .constraints import Constraint                        # noqa: F401
from .constraints import from_boundary                     # noqa: F401
from .constraints import from_curves                       # noqa: F401
from .densify import field_densification                   # noqa: F401
from .edit import coarse_from_skeleton                     # noqa: F401
from .edit import face_polylines                           # noqa: F401
from .edit import warp_polyline                            # noqa: F401
from .field import CrossField                              # noqa: F401
from .guides import guide_metrics                          # noqa: F401
from .quality import curve_alignment                       # noqa: F401
from .quality import hard_floor                            # noqa: F401
from .quality import mesh_quality                          # noqa: F401
from .relax import relax_mesh                              # noqa: F401
from .symmetry import Symmetry                             # noqa: F401

__all__ = [
    'BackgroundMesh',
    'SolveCache',
    'solve',
    'Constraint',
    'from_boundary',
    'from_curves',
    'CrossField',
    'guide_metrics',
    'relax_mesh',
    'field_densification',
    'coarse_from_skeleton',
    'face_polylines',
    'warp_polyline',
    'mesh_quality',
    'curve_alignment',
    'hard_floor',
    'Symmetry',
]
