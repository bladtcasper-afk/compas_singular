"""CMD_mcp_link, with Rhino stubbed out. Verification step D.

Run under Rhino 8's own CPython -- the same interpreter the link runs on inside
Rhino, so a name that resolves here resolves there. This catches the class of
bug that only fails at use: a helper imported but never called, or an API that
moved.

The checks that matter most are the GATES, because they are what keeps the link
out of the way of somebody using Rhino: it must not act inside a command, it
must wrap a write in one undo record, and it must give the selection back.
"""
import sys
import types

from _harness import check, dense_mesh, report, spool_dir

spool_dir()

# Import compas BEFORE the Rhino stubs go in. ``compas.plugins`` imports .NET's
# ``System`` when it believes it is inside Rhino, and it decides that by looking
# for a ``Rhino`` module -- so stubbing first makes a real compas import fail on
# a module CPython does not have. Importing here caches it while compas.RHINO is
# still False.
import compas.geometry            # noqa: F401
import compas.plugins             # noqa: F401

# ----------------------------------------------------------------------
# a Rhino that records what was asked of it
# ----------------------------------------------------------------------
CALLS = []


class StubDoc(object):
    def __init__(self):
        self.undo_records = []
        self.open_record = None

    def BeginUndoRecord(self, name):
        self.undo_records.append(name)
        self.open_record = name
        return len(self.undo_records)

    def EndUndoRecord(self, serial):
        self.open_record = None
        return True


class StubRs(types.ModuleType):
    """Just enough rhinoscriptsyntax to drive the link."""

    def __init__(self):
        types.ModuleType.__init__(self, "rhinoscriptsyntax")
        self.layers = {"Outer", "Mesh", "TopologyProblem"}
        self.objects = {}
        self.selected = ["guid-a", "guid-b"]
        self.unselected = False

    def DocumentName(self):
        return "plate.3dm"

    def IsLayer(self, name):
        return name in self.layers

    def AddLayer(self, name=None, parent=None, color=None):
        full = "{}::{}".format(parent, name) if parent else name
        self.layers.add(full)
        CALLS.append(("AddLayer", full))
        return full

    def ObjectsByLayer(self, layer):
        return list(self.objects.get(layer, []))

    def ObjectLayer(self, guid, layer):
        CALLS.append(("ObjectLayer", guid, layer))
        self.objects.setdefault(layer, []).append(guid)
        return layer

    def SelectedObjects(self):
        return list(self.selected)

    def UnselectAllObjects(self):
        self.unselected = True
        CALLS.append(("UnselectAllObjects",))

    def SelectObjects(self, guids):
        CALLS.append(("SelectObjects", tuple(guids)))
        return len(guids)

    def IsObject(self, guid):
        return True

    def CurrentLayer(self, *a, **k):
        raise AssertionError("the link must never change the current layer")

    def ZoomExtents(self, *a, **k):
        raise AssertionError("the link must never change the view")


class StubCommand(object):
    active = 0

    @staticmethod
    def InCommand():
        return StubCommand.active


rhino = types.ModuleType("Rhino")
rhino.Commands = types.SimpleNamespace(Command=StubCommand)


class _Idle(object):
    def __init__(self):
        self.handlers = []

    def __iadd__(self, handler):
        self.handlers.append(handler)
        return self

    def __isub__(self, handler):
        self.handlers.remove(handler)
        return self


rhino.RhinoApp = types.SimpleNamespace(Idle=_Idle())

stub_rs = StubRs()
stub_sc = types.ModuleType("scriptcontext")
stub_sc.doc = StubDoc()
stub_sc.sticky = {}

sys.modules["Rhino"] = rhino
sys.modules["rhinoscriptsyntax"] = stub_rs
sys.modules["scriptcontext"] = stub_sc

# CMD_start runs on import and would touch the real document; stub the two
# functions CMD_mcp_link actually uses out of it.
# ``compas_rhino.conversions`` genuinely needs RhinoCommon types, and
# ``helpers.py`` imports it at module scope. The link binds the whole of helpers
# at attach time on purpose -- no late imports inside the idle handler -- so the
# test has to satisfy that import. The four names below are all it takes; the
# bake and read functions themselves are replaced where they are exercised.
conversions = types.ModuleType("compas_rhino.conversions")
for _name in ("point_to_rhino", "polyline_to_rhino", "vertices_and_faces_to_rhino",
              "mesh_to_compas", "point_to_compas"):
    setattr(conversions, _name, lambda *a, **k: None)
sys.modules["compas_rhino.conversions"] = conversions

cmd_start = types.ModuleType("CMD_start")
cmd_start.import_compas_singular = lambda: None
cmd_start.ensure_paths = lambda: None
sys.modules["CMD_start"] = cmd_start

import CMD_mcp_link as link

print("\n=== 1. the link imports inside a (stubbed) Rhino ===")
check("CMD_mcp_link imports", link is not None)
check("it binds its dependencies", bool(link._bind()))
check("it understands exactly three verbs",
      sorted(link.VERBS) == ["ping", "pull", "push"], sorted(link.VERBS))

print("\n=== 2. it never reads through the DESTRUCTIVE boundary reader ===")
source = open(link.__file__).read()
# The name appears in the comment saying why it is NOT used, so look for a
# CALL and an IMPORT rather than the bare name.
# The destructive reader's name appears in the comment saying why it is NOT
# used, so look for the two ways it could actually be reached: a direct call,
# and the bound-dict key the handler resolves its helpers through.
check("read_boundaries is never called", "read_boundaries(" not in source)
check("and never imported", "import read_boundaries" not in source)
check("and is not in the bound helpers", '"read_boundaries"' not in source)
check("read_boundary_loops is what is bound instead",
      '"read_boundary_loops"' in source)
check("and the reason is written down where the next reader will see it",
      "DELETES" in source)

print("\n=== 3. ping answers, and the spool round-trips through handle() ===")
spool = link._BOUND["spool"]
rid = spool.post("ping")
link.handle(*spool.take())
out = spool.wait(rid, timeout=2)
check("ping is answered", out.get("ok") is True, out)
check("naming the document", out["result"]["document"] == "plate.3dm")

print("\n=== 4. an unknown verb is refused, not crashed on ===")
rid = spool.post("frobnicate")
link.handle(*spool.take())
out = spool.wait(rid, timeout=2)
check("refused", out.get("ok") is False)
check("and says what it does understand", "ping, pull and push" in out["reason"])

print("\n=== 5. a verb that raises becomes a refusal carrying the cause ===")
rid = spool.post("push", {"layer": "QuadMesh"})          # no mesh in the payload
link.handle(*spool.take())
out = spool.wait(rid, timeout=2)
check("refused rather than raised", out.get("ok") is False)
check("naming the verb and the error", "push failed in Rhino" in out["reason"],
      out["reason"][:60])

print("\n=== 6. THE GATE: nothing happens while the user is in a command ===")
rid = spool.post("ping")
StubCommand.active = 1
link.on_idle(None, None)
check("the request is still queued", spool.pending() == [rid], spool.pending())
check("and no answer was written",
      spool.wait(rid, timeout=0.2).get("timed_out") is True)
StubCommand.active = 0
link.on_idle(None, None)
check("it drains once the command ends", spool.wait(rid, timeout=2).get("ok") is True)

print("\n=== 7. an idle tick with nothing to do only beats ===")
before = len(CALLS)
link.on_idle(None, None)
check("no document calls on an empty queue", len(CALLS) == before, CALLS[before:])
check("but the heartbeat is written", spool.link_state() is not None)
check("naming the document", spool.link_state()["document"] == "plate.3dm")

print("\n=== 8. a push is ONE undo record, and gives the selection back ===")
mesh, _walls = dense_mesh(density=2)
wire = link._BOUND["wire"]
stub_rs.objects["TopologyProblem::QuadMesh"] = ["old-mesh-guid"]
baked = {}
link._BOUND["bake_mesh"] = lambda m, layer, **kw: baked.setdefault("layer", layer)
link._BOUND["clear_layer"] = lambda layer, **kw: 0
CALLS[:] = []
records_before = len(stub_sc.doc.undo_records)

rid = spool.post("push", {"layer": "QuadMesh", "mesh": wire.mesh_to_wire(mesh)})
link.handle(*spool.take())
out = spool.wait(rid, timeout=2)
check("the push succeeds", out.get("ok") is True, out.get("reason"))
check("exactly one undo record was opened",
      len(stub_sc.doc.undo_records) == records_before + 1,
      stub_sc.doc.undo_records)
check("named so a person recognises it in Rhino",
      stub_sc.doc.undo_records[-1] == "MCP push")
check("and it was closed", stub_sc.doc.open_record is None)
check("the previous mesh was moved aside, not deleted",
      ("ObjectLayer", "old-mesh-guid", link.BEFORE_LAYER) in CALLS, CALLS)
check("the tool is told where it went",
      out["result"]["before_layer"] == link.BEFORE_LAYER)
check("the selection was restored",
      ("SelectObjects", ("guid-a", "guid-b")) in CALLS, CALLS)
check("the face count came through",
      out["result"]["faces"] == mesh.number_of_faces())

print("\n=== 9. it survives another command purging compas_singular ===")
import compas_singular
orphan = types.ModuleType("compas_singular")
sys.modules["compas_singular"] = orphan          # what import_compas_singular does
check("the link notices the package changed",
      link._BOUND["package"] is not sys.modules["compas_singular"])
sys.modules["compas_singular"] = compas_singular
rebound = link._fresh()
check("and rebinds rather than refusing",
      rebound["package"] is sys.modules["compas_singular"])
rid = spool.post("ping")
link.handle(*spool.take())
check("serving still works afterwards", spool.wait(rid, timeout=2).get("ok") is True)

print("\n=== 10. attach and detach toggle cleanly ===")
check("not attached to begin with", link.attached() is False)
link.attach()
check("attached", link.attached() is True)
check("the idle handler is registered", link.on_idle in rhino.RhinoApp.Idle.handlers)
check("the heartbeat is live", spool.link_state() is not None)
link.detach()
check("detached", link.attached() is False)
check("the handler is gone", link.on_idle not in rhino.RhinoApp.Idle.handlers)
check("and the server is told nothing is listening", spool.link_state() is None)

sys.exit(report("test_rhino_side"))
