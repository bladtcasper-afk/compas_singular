"""
********************************************************************************
compas_singular.geometry
********************************************************************************

.. currentmodule:: compas_singular.geometry


Array
=====

Array functions.

.. autosummary::
    :toctree: generated/
    :nosignatures:

    circle_evaluate
    archimedean_spiral_evaluate
    line_array
    rectangular_array
    circular_array
    spiral_array


Polyline
========

Polyline class with additional methods.

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Polyline


Projection
==========

Projection functions.

.. autosummary::
    :toctree: generated/
    :nosignatures:

    closest_point_on_circle
    closest_point_on_line
    closest_point_on_segment
    closest_point_on_polyline
    closest_point_on_polylines

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas_singular.geometry.array import *  # noqa: F401 F403
from compas_singular.geometry.polyline import *  # noqa: F401 F403
from compas_singular.geometry.projection import *  # noqa: F401 F403

import types  # noqa: E402

# Only re-export names bound in this namespace, never the submodules themselves:
# ``from .foo import *`` also binds ``foo`` as an attribute of this package, and
# re-exporting that module object shadows any function or subpackage of the same
# name further up the chain.
__all__ = [name for name, obj in list(globals().items())
           if not name.startswith('_') and not isinstance(obj, types.ModuleType)]
