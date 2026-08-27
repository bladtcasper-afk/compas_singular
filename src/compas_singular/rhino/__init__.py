"""
********************************************************************************
compas_singular.rhino
********************************************************************************

.. currentmodule:: compas_singular.rhino


Curve
=====

Curve class for Rhino with additional methods

.. autosummary::
    :toctree: generated/
    :nosignatures:

    RhinoCurve


Surface
=======

Surface class for Rhino with additional methods

.. autosummary::
    :toctree: generated/
    :nosignatures:

    RhinoSurface


Artist
======

Artist to select mesh elements.

.. autosummary::
    :toctree: generated/
    :nosignatures:

    select_mesh_polyedge
    select_quad_mesh_polyedge
    select_quad_mesh_strip
    select_quad_mesh_strip


Draw
====

Drawing specific objects.

.. autosummary::
    :toctree: generated/
    :nosignatures:

    draw_graph


"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .artists import *  # noqa: F401 F403
from .geometry import *  # noqa: F401 F403
from .objects import *  # noqa: F401 F403
from .constraints import *  # noqa: F401 F403

import types  # noqa: E402

# Only re-export names bound in this namespace, never the submodules themselves:
# ``from .foo import *`` also binds ``foo`` as an attribute of this package, and
# re-exporting that module object shadows any function or subpackage of the same
# name further up the chain.
__all__ = [name for name, obj in list(globals().items())
           if not name.startswith('_') and not isinstance(obj, types.ModuleType)]
