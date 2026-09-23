"""
********************************************************************************
compas_singular.algorithms
********************************************************************************

.. currentmodule:: compas_singular.algorithms


Coloring
========

Algorithm for two-colouring of quad meshes.

.. autosummary::
    :toctree: generated/
    :nosignatures:

    TwoColourableProjection


Decomposition
=============

Algorithm for decomposition of surfaces into coarse quad meshes. Optional point and curve features.

.. autosummary::
    :toctree: generated/
    :nosignatures:

    surface_discrete_mapping
    boundary_triangulation
    Skeleton
    SkeletonDecomposition
    DecompositionRemap


Interpolation
=============

Algorithm for combination and interpolation of quad meshes.

.. autosummary::
    :toctree: generated/
    :nosignatures:


"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import compas

if not compas.IPY:
    from compas_singular.algorithms.isomorphism import *  # noqa: F401 F403

from compas_singular.algorithms.layout import *  # noqa: F401 F403
from compas_singular.algorithms.propagation import *  # noqa: F401 F403
from compas_singular.algorithms.triangulation import *  # noqa: F401 F403
from compas_singular.algorithms.twocoloring import *  # noqa: F401 F403

from compas_singular.algorithms.skeleton_decomposition import *  # noqa: F401 F403


import types  # noqa: E402

# Only re-export names bound in this namespace, never the submodules themselves:
# ``from .foo import *`` also binds ``foo`` as an attribute of this package, and
# re-exporting that module object shadows any function or subpackage of the same
# name further up the chain.
__all__ = [name for name, obj in list(globals().items())
           if not name.startswith('_') and not isinstance(obj, types.ModuleType)]
