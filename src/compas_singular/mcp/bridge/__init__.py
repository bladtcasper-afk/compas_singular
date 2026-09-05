"""
********************************************************************************
bridge
********************************************************************************

.. currentmodule:: compas_singular.mcp.bridge

**The wire to Rhino: a directory of JSON files, and nothing else.**

The server runs in its own process and Rhino runs in its own, and neither one
imports the other's runtime. What passes between them is plain JSON in a spool
directory: the server writes a request, Rhino's ``CMD_mcp_link`` picks it up on
the main thread and writes a response.

**Why a spool and not a socket.** A Rhino document may only be touched from the
main thread, so anything arriving on a background thread has to be handed across
anyway. A spool does that handoff with no thread, no port and no listener to
leak: the request simply sits on disk until Rhino is idle enough to take it, and
a request posted while Rhino is closed waits rather than failing.

**Both modules here are imported by BOTH sides**, so they are standard library
plus compas and nothing more. In particular they must stay importable inside
Rhino's interpreter after ``CMD_start.ensure_paths()``, which is why neither one
imports anything from the rest of :mod:`compas_singular.mcp`.

Modules
=======

``spool``
    The request/response file protocol.
``wire``
    Meshes and curves as plain, index-based JSON.

"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


__all__ = []
