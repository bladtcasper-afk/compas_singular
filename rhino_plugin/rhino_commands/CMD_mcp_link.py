#! python3
# r: compas
"""**Attach this document to the standalone MCP server. Toggle on, toggle off.**

    reads   Outer / Inner / Guides                the domain
            TopologyProblem::QuadMesh             the dense mesh
            Mesh + Poles                          the coarse layout, on request
    writes  TopologyProblem::QuadMesh             only on a push
            ...::QuadMesh::MCP::Before            the mesh that was there
    prints  one line per request served

**RUN IT ONCE TO ATTACH, AGAIN TO DETACH.** While it is attached you keep using
Rhino normally: model, run any other CMD_ command, undo your own work, change
layers and views. The link only acts when Rhino is idle, and it never acts while
you are inside a command.

**This is a COURIER, not an implementation.** It understands three verbs --
``ping``, ``pull``, ``push`` -- and nothing about meshes. Every decision about a
mesh is made by ``compas_singular.mcp`` in a separate process, which is tested
without Rhino in ``examples/mcp_tests``. Nothing here decides anything.

**Why a spool of files and not a socket.** A Rhino document may only be touched
from the main thread, so a background listener would have to hand its work
across anyway. A directory of JSON files does that with no thread and no port to
leak, and a request posted while Rhino is closed simply waits instead of failing.

**FOUR THINGS KEEP THIS OUT OF YOUR WAY**, and they are the whole design:

1. It never runs inside a command. ``RhinoApp.Idle`` fires while you are in the
   middle of a pick too, so the first thing the handler does is check
   ``Command.InCommand()`` and give up the tick if you are busy. Requests queue
   and drain when you finish -- which is also why the link never fights
   ``CMD_densities`` or any other command for the document.
2. A write is one named undo record. Your Ctrl+Z steps over your own modelling;
   the link's bake is a single step called "MCP push" that you can undo on
   purpose.
3. It never changes your selection, your current layer, or your view.
4. It writes only under ``TopologyProblem``. Never to a layer you made.

**A purge is expected, not an error.** ``CMD_start.import_compas_singular``
empties every ``compas_singular`` module out of ``sys.modules``, so running any
other command in the family leaves this one holding an orphaned copy. The handler
notices and rebinds. That is safe here only because the wire format is plain
JSON with no compas types in it -- nothing is pickled and nothing is
``isinstance``-checked across the boundary, so two copies of the package alive at
once cannot produce the ``not the same object as ...Mesh`` failure.
"""

# Temporary import of the compas_singular development library. MUST come before
# any compas_singular import: it fixes sys.path and purges a stale copy.
from CMD_start import import_compas_singular
import_compas_singular()

import sys
import traceback

import Rhino
import rhinoscriptsyntax as rs
import scriptcontext as sc

from CMD_start import ensure_paths


ROOT = "TopologyProblem"
QUADMESH_LAYER = ROOT + "::QuadMesh"
BEFORE_LAYER = QUADMESH_LAYER + "::MCP::Before"

#: Layers that may be pulled from, by the short name the tool passes.
PULL_LAYERS = {
    "QuadMesh": QUADMESH_LAYER,
    "Mesh": "Mesh",
}

#: Where the toggle keeps its state, so a second run finds the first one's
#: handler. ``sc.sticky`` survives between script runs; a module global does not.
STICKY_KEY = "mcp.link.handler"

_BOUND = {}


def _bind():
    """Import everything the handler needs and hold the references.

    Called once at attach and again whenever a purge is detected. Uses
    ``ensure_paths`` -- safe to call anywhere -- rather than
    ``import_compas_singular``, which purges and is only safe in a module body.
    """
    ensure_paths()
    import compas_singular
    from compas_singular.mcp.bridge import spool
    from compas_singular.mcp.bridge import wire
    from compas_singular.rhino.helpers.helpers import bake_mesh
    from compas_singular.rhino.helpers.helpers import clean_layer
    from compas_singular.rhino.helpers.helpers import read_boundary_loops
    from compas_singular.rhino.helpers.helpers import read_mesh
    from compas_singular.rhino.helpers.helpers import read_polylines
    _BOUND.update({
        "package": compas_singular, "spool": spool, "wire": wire,
        "bake_mesh": bake_mesh, "clean_layer": clean_layer,
        "read_boundary_loops": read_boundary_loops, "read_mesh": read_mesh,
        "read_polylines": read_polylines,
    })
    return _BOUND


def _fresh():
    """Rebind if another command has purged the package out from under us."""
    if _BOUND.get("package") is not sys.modules.get("compas_singular"):
        _bind()
    return _BOUND


def ensure_layer(path, color=None):
    """Create a ``::`` layer path, parents first, and return the full path.

    ``rs.AddLayer`` does not create intermediate parents, and ``bake_mesh``'s own
    layer creation hard-codes ``parent="TopologyProblem"``, which is wrong for a
    layer three levels down. ``name`` is the LEAF, never the full path: passing
    the path as the name AND a parent makes Rhino create a layer literally
    called "TopologyProblem::QuadMesh" nested inside "TopologyProblem".

    Same implementation as ``CMD_ai_edit.ensure_layer`` and
    ``CMD_edit_quad_mesh.ensure_layer``.
    """
    parts = path.split("::")
    for i in range(len(parts)):
        name = "::".join(parts[:i + 1])
        if not rs.IsLayer(name):
            rs.AddLayer(name=parts[i],
                        parent="::".join(parts[:i]) if i else None,
                        color=color if i == len(parts) - 1 else None)
    return path


def document_name():
    try:
        return rs.DocumentName() or "(unsaved)"
    except Exception:
        return "(unknown)"


def _busy():
    """Whether Rhino is inside a command, i.e. whether the user is mid-action.

    ``Idle`` fires during a pick as well as between commands, so without this
    the link would bake geometry underneath somebody's cursor.
    """
    try:
        return bool(Rhino.Commands.Command.InCommand())
    except Exception:
        # If the API is not where we expect it, refuse to act rather than act
        # blind: a missed request costs a wait, a wrongly-timed one costs work.
        return True


# ======================================================================
# the three verbs
# ======================================================================

def _verb_ping(args):
    return {"document": document_name()}


def _verb_pull(args):
    """Read the domain and a mesh out of the document. Changes nothing."""
    bound = _fresh()
    layer = args.get("layer") or "QuadMesh"
    spacing = float(args.get("spacing") or 0.125)
    resolved = PULL_LAYERS.get(layer, layer)

    # NOT ``read_boundaries``: that one DELETES any inner boundary, guide or
    # pole lying outside the outer boundary. Correct for a person driving a
    # selection command, unacceptable for a read. These two only read.
    outer, inners = [], []
    try:
        outer, inners = bound["read_boundary_loops"](spacing)
    except RuntimeError:
        pass                      # no Outer layer; the tool warns about it
    try:
        guides = bound["read_polylines"]("Guides")
    except Exception:
        guides = []

    mesh = None
    if layer == "Mesh":
        # The coarse layout carries poles, which a plain mesh read would drop.
        from compas_singular.rhino.helpers.helpers import read_coarse
        try:
            mesh, _poles = read_coarse()
        except RuntimeError:
            mesh = None
    else:
        try:
            mesh = bound["read_mesh"](resolved)
        except RuntimeError:
            mesh = None

    return {
        "document": document_name(),
        "layer": layer,
        "outer": [list(p) for p in outer] if outer else [],
        "inners": [[list(p) for p in loop] for loop in inners],
        "guides": [[list(p) for p in g] for g in guides],
        "mesh": bound["wire"].mesh_to_wire(mesh) if mesh is not None else None,
    }


def _verb_push(args):
    """Bake a mesh onto a layer, keeping what was there. The only write."""
    bound = _fresh()
    layer = args.get("layer") or "QuadMesh"
    target = PULL_LAYERS.get(layer, layer)
    payload = args.get("mesh")
    if not payload:
        raise ValueError("the push carried no mesh")

    mesh = bound["wire"].mesh_from_wire(payload)

    # Your selection is yours. Restore it whatever happens below.
    try:
        selected = list(rs.SelectedObjects() or [])
    except Exception:
        selected = []

    # ONE named undo record, so Ctrl+Z steps over this as a single deliberate
    # action instead of unpicking it vertex by vertex among your own edits.
    serial = sc.doc.BeginUndoRecord("MCP push")
    try:
        ensure_layer(target)
        ensure_layer(BEFORE_LAYER)
        # Move what is there aside rather than deleting it, so the two can be
        # compared. The previous backup goes, not the previous mesh.
        bound["clean_layer"](BEFORE_LAYER)
        moved = 0
        for guid in rs.ObjectsByLayer(target) or []:
            try:
                rs.ObjectLayer(guid, BEFORE_LAYER)
                moved += 1
            except Exception:
                pass
        guid = bound["bake_mesh"](mesh, target, clear_existing=False)
        # VERIFY rather than trust: bake_mesh raises only on an empty guid.
        on_layer = len(rs.ObjectsByLayer(target) or [])
    finally:
        sc.doc.EndUndoRecord(serial)
        try:
            rs.UnselectAllObjects()
            alive = [g for g in selected if rs.IsObject(g)]
            if alive:
                rs.SelectObjects(alive)
        except Exception:
            pass

    return {
        "document": document_name(),
        "layer": target,
        "faces": mesh.number_of_faces(),
        "vertices": mesh.number_of_vertices(),
        "before_layer": BEFORE_LAYER if moved else None,
        "moved_aside": moved,
        "objects_on_layer": on_layer,
        "undo_record": "MCP push",
    }


VERBS = {"ping": _verb_ping, "pull": _verb_pull, "push": _verb_push}


def handle(request_id, verb, args):
    """Serve one request. Never raises: a fault becomes the answer."""
    bound = _fresh()
    function = VERBS.get(verb)
    if function is None:
        bound["spool"].answer(
            request_id, False,
            error="this link understands ping, pull and push, not "
                  "{!r}".format(verb))
        return
    try:
        result = function(args or {})
    except Exception as exc:
        traceback.print_exc()
        bound["spool"].answer(
            request_id, False,
            error="{} failed in Rhino: {}: {}".format(
                verb, type(exc).__name__, exc))
        print("MCP link: {} FAILED -- {}".format(verb, exc))
        return
    bound["spool"].answer(request_id, True, result=result)
    print("MCP link: served {}".format(verb))


# ======================================================================
# the pump
# ======================================================================

def on_idle(sender, e):
    """One request per tick, on the main thread, never inside a command.

    One per tick rather than draining the queue, so a burst of requests can
    never hold the UI for longer than a single operation.
    """
    try:
        bound = _fresh()
        if _busy():
            return                      # you are mid-command; try again later
        job = bound["spool"].take()
        if job is None:
            bound["spool"].heartbeat(document_name())
            return
        handle(*job)
    except Exception:
        # An exception escaping an event handler can detach it or spam the
        # command line. Swallow it here, having printed it once.
        traceback.print_exc()


def attached():
    return sc.sticky.get(STICKY_KEY) is not None


def attach():
    _bind()
    Rhino.RhinoApp.Idle += on_idle
    sc.sticky[STICKY_KEY] = on_idle
    _BOUND["spool"].heartbeat(document_name(), force=True)
    folder = _BOUND["spool"].default_directory()
    print("MCP link ATTACHED to {}".format(document_name()))
    print("  spool: {}".format(folder))
    print("  Keep using Rhino normally. The link only acts when Rhino is idle,")
    print("  and never while you are inside a command.")
    print("  Run CMD_mcp_link again to detach.")


def detach():
    handler = sc.sticky.pop(STICKY_KEY, None)
    if handler is not None:
        try:
            Rhino.RhinoApp.Idle -= handler
        except Exception:
            traceback.print_exc()
    try:
        _fresh()["spool"].clear_link()
    except Exception:
        pass
    print("MCP link DETACHED. The server will report that nothing is listening.")


def main():
    if attached():
        detach()
    else:
        attach()


if __name__ == "__main__":
    main()
