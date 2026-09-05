"""Pull, improve, push -- the whole loop, with a fake link on the other end.

The closest thing to verification step E that can run without Rhino. A thread
plays ``CMD_mcp_link``: it takes requests off the same spool the real link uses
and answers them with the same payload shapes. What that exercises is the half
``test_rhino_side`` cannot -- ``tools_rhino``'s SUCCESS path, which is otherwise
only ever seen refusing.

The mesh work in between is real. Only the document is fake.
"""
import json
import sys
import threading
import time

from _harness import check, dense_mesh, library_dir, report, spool_dir

spool_dir()
library_dir()

from compas_singular.mcp.bridge import spool, wire
from compas_singular.mcp.server import build

MESH, WALLS = dense_mesh(density=4)
PUSHED = {}
STOP = threading.Event()


def fake_link():
    """What CMD_mcp_link does, minus Rhino."""
    while not STOP.is_set():
        job = spool.take()
        if job is None:
            spool.heartbeat("fake.3dm")
            time.sleep(0.01)
            continue
        request_id, verb, args = job
        if verb == "ping":
            spool.answer(request_id, True, result={"document": "fake.3dm"})
        elif verb == "pull":
            spool.answer(request_id, True, result={
                "document": "fake.3dm", "layer": args.get("layer"),
                "outer": WALLS[0], "inners": WALLS[1:], "guides": [],
                "mesh": wire.mesh_to_wire(MESH)})
        elif verb == "push":
            PUSHED["mesh"] = args["mesh"]
            spool.answer(request_id, True, result={
                "document": "fake.3dm", "layer": args.get("layer"),
                "before_layer": "TopologyProblem::QuadMesh::MCP::Before",
                "undo_record": "MCP push"})
        else:
            spool.answer(request_id, False, error="unknown verb")


worker = threading.Thread(target=fake_link, daemon=True)
worker.start()
time.sleep(0.1)

dispatcher = build()
session = dispatcher.handler.session


def call(_tool, **arguments):
    reply = dispatcher.handle({"jsonrpc": "2.0", "id": 1, "method": "tools/call",
                               "params": {"name": _tool, "arguments": arguments}})
    return json.loads(reply["result"]["content"][0]["text"])


print("\n=== 1. status finds the link ===")
out = call("rhino_status")
check("attached", out.get("attached") is True, out)
check("naming the document", out.get("document") == "fake.3dm")
check("with a fresh heartbeat", out.get("stale") is False)

print("\n=== 2. pull brings back the mesh AND the walls ===")
out = call("rhino_pull", layer="QuadMesh")
check("ok", out.get("ok") is True, out.get("reason"))
check("the faces arrive", out["faces"] == MESH.number_of_faces(), out["faces"])
check("the walls arrive", out["walls"] == len(WALLS), out["walls"])
check("so there is no missing-wall warning", "warning" not in out)
check("it reads the mesh straight away", len(out["reading"]) > 40)
check("the session holds them", len(session.walls) == len(WALLS))
check("coordinates survive the wire exactly",
      session.mesh.vertex_coordinates(list(session.mesh.vertices())[0])
      == MESH.vertex_coordinates(list(MESH.vertices())[0]))

print("\n=== 3. improve it for real ===")
start = session.triple()
out = call("relax")
check("relax ran", out.get("ok") is True, out.get("reason"))
end = session.triple()
if out["accepted"] is not None:
    check("and improved all three", out["all_improved"] is True, out)
    check("min angle rose", end[0] > start[0], (start, end))
else:
    check("or left the mesh alone", start == end)

print("\n=== 4. push sends the improved mesh back ===")
check("it refuses to push a mesh nothing has looked at",
      call("rhino_push", layer="QuadMesh").get("ok") is False)
looked = call("inspect", image=True)
check("the final look reports the inputs it was judged against",
      looked["inputs"]["walls"] == len(WALLS), looked.get("inputs"))
out = call("rhino_push", layer="QuadMesh")
check("ok", out.get("ok") is True, out.get("reason"))
check("the previous mesh was kept", "Before" in (out.get("kept_previous_on") or ""))
check("as one named undo step", out.get("undo_record") == "MCP push")
check("the link received a mesh", "mesh" in PUSHED)
check("with the same face count",
      len(PUSHED["mesh"]["faces"]) == MESH.number_of_faces())
check("carrying the IMPROVED positions, not the pulled ones",
      PUSHED["mesh"]["vertices"]
      == wire.mesh_to_wire(session.mesh)["vertices"])
if out.get("ok") and end != start:
    check("which really differ from what was pulled",
          PUSHED["mesh"]["vertices"] != wire.mesh_to_wire(MESH)["vertices"])

print("\n=== 5. the journal reads back as a record of all of it ===")
text = dispatcher.handler.read_resource("session://current")
check("it names the source layer", "QuadMesh" in text, text[:80])
check("and lists the steps", "`relax`" in text and "`rhino_push`" in text)

print("\n=== 6. a link that goes away is reported, not hung on ===")
STOP.set()
worker.join(timeout=2)
spool.clear_link()
out = call("rhino_pull", timeout=0.3)
check("pull refuses once the link is gone", out.get("ok") is False)
check("promptly", out.get("timed_out") is True)
check("naming the fix", "CMD_mcp_link" in out.get("reason", ""), out.get("reason"))

sys.exit(report("test_roundtrip"))
