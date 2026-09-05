"""
********************************************************************************
agent.runner
********************************************************************************

.. currentmodule:: compas_singular.agent.runner

**The half that talks to a model.**

:mod:`~compas_singular.agent.core` is the session, the tools and the state
machine, and knows nothing about any model. This package adds the loop that
turns an instruction into a sequence of tool calls, and the backend that carries
those calls to a provider.

**Importing this package does NOT import an SDK.** :class:`AnthropicBackend`
imports ``anthropic`` inside its constructor, so the whole layer stays
importable -- and testable, through :class:`ScriptedBackend` -- in an
interpreter that has never installed it. That matters here: neither the
``singular312`` env nor Rhino's own CPython has the SDK, and the geometry half
must not stop working because of it.

    from compas_singular.agent.core import MeshEditSession
    from compas_singular.agent.runner import run_session

    session = MeshEditSession(decomposition, outer=outline)
    result = run_session(session, 'even out the mesh around the hole',
                         confirm=True, on_event=lambda e: print(e['message']))
    print(result.summary())

Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    SessionResult
    Backend
    AnthropicBackend
    GeminiBackend
    ScriptedBackend
    KeyNotFound

Functions
=========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    run_session
    system_prompt
    read_key_file
    resolve_api_key
    repair_ssl_env

"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .backend import *  # noqa: F401 F403
from .backend_gemini import *  # noqa: F401 F403
from .keys import *  # noqa: F401 F403
from .loop import *  # noqa: F401 F403
from .prompt import system_prompt  # noqa: F401

__all__ = [name for name in dir() if not name.startswith('_')]
