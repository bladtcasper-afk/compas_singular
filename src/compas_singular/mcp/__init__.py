"""
********************************************************************************
mcp
********************************************************************************

.. currentmodule:: compas_singular.mcp

**The meshing tools as an MCP server: the model lives outside, not inside.**

:mod:`compas_singular.agent` calls a model from inside Python -- it owns the tool
vocabulary, the state machine, the system prompt AND the transport, so changing
any one of the four means editing the package. This is the other arrangement.
The tools are exposed over the Model Context Protocol, an MCP client supplies the
model and the conversation, and this package supplies only the mesh operations
and the knowledge about them.

**It imports nothing from** :mod:`compas_singular.agent`, **and never should.**
The whole point of having two is that they can be compared, and two things that
share a session model, an addressing scheme or a rebuild path are not two things.
Where this package solves a problem ``agent`` has already solved -- geometric
addressing is the obvious one -- it solves it again, independently, so that a bug
in one cannot hide in the other. Delete either package and the other still works.

**Standalone from Rhino, but connected to it.** The server runs as its own
process, with Rhino closed, on meshes loaded from disk. When Rhino IS open and
``CMD_mcp_link`` is attached, :mod:`~compas_singular.mcp.bridge` lets the server
pull geometry out of the live document and push results back -- through a
directory of JSON files, so neither side ever imports the other's runtime.

**Three capabilities, not one.** Tools are the verbs. Resources carry the
guidance -- what good quality looks like, when to reach for which smoother -- as
Markdown you can edit without a release, where ``agent`` compiles the equivalent
into ``runner/prompt.py``. Prompts carry the workflows. And the library keeps
finished sessions as worked examples, which ``agent`` writes but never reads
back.

Scope
=====

**Reads and improves** an existing mesh: pull one, smooth or relax it, push it
back. **Builds and edits a coarse layout, too**: from a domain's boundary, line
and point features, via a topological-skeleton decomposition -- not a frame
field, and not a read of hand-drawn patch curves -- through to a densified
dense mesh, with strip add/remove as the only topology edit either mesh ever
gets. The DENSE mesh's own topology is never changed by any tool here.

No field is solved anywhere in this package, which is why it does not need
:mod:`~compas_singular.framefield`'s solver. It does now reach for numpy/scipy
(the skeleton decomposition's triangulation) and, for ``relax_fdm`` only,
``compas_fd`` -- a missing one of those is a clean tool refusal, not an
import-time failure of the server itself.

Modules
=======

``protocol``
    JSON-RPC 2.0 over stdio. Standard library only.
``server``
    Binds the registry, the session and the library to the protocol's methods.
``registry``
    The ``@tool`` decorator and the JSON schemas it emits.
``session``
    Mesh, walls, guides, history, remarks, and a position-map undo.
``handle``
    Geometric addressing. No raw vertex key ever crosses the protocol.
``describe``
    Turns a metric dict into a reading a person or a model can act on.
``library``
    Guidance, prompts and the corpus of worked examples.
``bridge``
    The Rhino wire: a spool of JSON files, and the mesh encoding on it.

"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


__all__ = []
