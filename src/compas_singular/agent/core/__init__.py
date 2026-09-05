"""
********************************************************************************
agent.core
********************************************************************************

.. currentmodule:: compas_singular.agent.core

**The session, with no model API in it.**

Everything an automated caller needs to hold the pipeline as one job: how to
NAME things that renumber (:mod:`address`), how to densify without destroying
the caller's work (:mod:`rebuild`), and the state machine that says which of the
three meshes is in hand and what a backward move would cost (:mod:`session`).

Nothing here imports an SDK, and nothing here imports Rhino. That is what makes
the whole layer testable offline: a scripted plan drives exactly the same
functions a model would, with no key and no CAD session.

Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    MeshEditSession
    AddressBook
    DensifyRefused
    StageError

Functions
=========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    point_address
    vertex_address
    edge_address
    face_address
    strip_address
    apply_densities
    read_densities
    densify_preserving

"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .address import *  # noqa: F401 F403
from .digest import *  # noqa: F401 F403
from .inputs import *  # noqa: F401 F403
from .rebuild import *  # noqa: F401 F403
from .render import *  # noqa: F401 F403
from .session import *  # noqa: F401 F403
from .tools import *  # noqa: F401 F403

__all__ = [name for name in dir() if not name.startswith('_')]
