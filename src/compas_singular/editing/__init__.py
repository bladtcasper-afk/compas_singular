"""
********************************************************************************
editing
********************************************************************************

.. currentmodule:: compas_singular.editing

**Hand-editing a layout, with no CAD in it.**

Every operation a user performs on a coarse layout in Rhino -- dragging a
corner, cutting patches with a drawn line, deleting a strip, committing the
result back to the decomposition -- is defined here, on the mesh, and can be
run from a plain script with no Rhino open. ``compas_singular.rhino`` and the
``CMD_`` commands supply the picks, the previews and the prompts; they supply
nothing else.

The split is not tidiness. This is the half that can be wrong in ways a user
cannot see -- a cut that leaves a five-sided patch, a curve map that no longer
matches its edges -- so it has to be testable without a CAD session.

Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    CoarseLayoutEditor
    GuideCurve

Functions
=========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    guide_chain
    attach_chain
    chain_quality
    collect_polyedges
    mean_edge_length

"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .coarse_layout import *  # noqa: F401 F403
from .guide_chain import *  # noqa: F401 F403

__all__ = [name for name in dir() if not name.startswith('_')]
