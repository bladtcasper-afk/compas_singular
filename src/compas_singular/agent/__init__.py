"""
********************************************************************************
agent
********************************************************************************

.. currentmodule:: compas_singular.agent

**Driving the meshing pipeline from something that is not a person.**

The pipeline has six stages a user works through by hand -- solve a field, edit
the coarse layout, densify, set densities, smooth, edit the dense mesh -- and
each one already has a backend with no CAD in it. This package is the layer that
holds those backends as ONE session with a state machine over it, so a caller
that is not sitting at a mouse can drive them: an LLM tool loop, a batch script,
a regression sweep.

**The dependency runs one way and never back.** ``compas_singular.agent`` imports
the rest of the library; nothing in the library imports ``agent``. Delete this
package and everything else still works.

**This module deliberately does NOT star-import** :mod:`~compas_singular.agent.runner`,
unlike every other ``__init__`` in this library. ``runner`` is the half that
talks to a model API and imports its SDK; ``core`` is pure Python and compas.
Star-importing ``runner`` here would make that SDK a hard dependency of
``compas_singular`` and break every headless and Rhino import that never asks
for an agent. Reach the loop explicitly::

    from compas_singular.agent.runner.loop import run_session

Sub-packages
============

``core``
    The session, the addressing, and the rebuild. No model API, no CAD.
``runner``
    The tool schemas, the prompt and the loop. Imports an SDK.

"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .core import *  # noqa: F401 F403

__all__ = [name for name in dir() if not name.startswith('_')]
