"""Every tool, driven through the server, with no Rhino and no model.

Verification step C. This is the test that says the mesh work is right; the
protocol test says the wire is right, and they are deliberately separate.
"""
import json
import os
import sys
import tempfile

from _harness import check, dense_mesh, library_dir, report, spool_dir

spool_dir()
library_dir()

from compas_singular.mcp.server import build

dispatcher = build()
handler = dispatcher.handler
session = handler.session


def call(_tool, **arguments):
    """Call a tool the way a client does. The parameter is underscored because
    ``save_example`` itself takes a ``name``, and a plain ``name`` here would
    collide with it."""
    reply = dispatcher.handle({"jsonrpc": "2.0", "id": 1, "method": "tools/call",
                               "params": {"name": _tool, "arguments": arguments}})
    return json.loads(reply["result"]["content"][0]["text"])


print("\n=== 1. every tool refuses cleanly with nothing loaded ===")
for name in ("inspect", "relax", "smooth_boundary_constrained", "save_mesh"):
    out = call(name, **({"path": "x.json"} if name == "save_mesh" else {}))
    check("{} refuses".format(name), out.get("ok") is False, out.get("reason", "")[:44])
check("undo refuses with nothing to undo", call("undo").get("ok") is False)
check("smooth_region refuses",
      call("smooth_region", region={"kind": "all"}).get("ok") is False)

print("\n=== 2. load a mesh from disk: the server works with Rhino closed ===")
mesh, walls = dense_mesh(density=4)
tmp = tempfile.mkdtemp()
path = os.path.join(tmp, "dense.json")
mesh.to_json(path)

out = call("load_mesh", path=path)
check("loads without walls", out.get("ok") is True, out.get("reason"))
check("and warns the boundary cannot be corrected", "warning" in out)
out = call("load_mesh", path=path, walls=walls)
check("loads with walls", out.get("ok") is True and out["walls"] == len(walls),
      out.get("walls"))
check("no warning once walls are given", "warning" not in out)
check("faces survive the load", out["faces"] == mesh.number_of_faces(), out["faces"])
check("a missing file is a refusal", call("load_mesh", path="nope.json").get("ok") is False)

print("\n=== 3. inspect reports numbers AND a reading ===")
out = call("inspect")
check("ok", out.get("ok") is True)
for key in ("min_angle", "max_angle", "aspect_max", "share_below", "worst_face"):
    check("quality carries {}".format(key), key in out["quality"])
check("a reading in prose", len(out["reading"]) > 40, out["reading"][:44])
check("a verdict", out["verdict"] in ("good", "usable", "poor", "unusable"), out["verdict"])
check("the worst face is a HANDLE, not a key",
      out["worst_face_at"] is None or out["worst_face_at"].startswith("v:"),
      out["worst_face_at"])
check("no raw face key escapes", "worst_face" not in out["reading"])
check("singularities are listed as handles",
      all(h.startswith("v:") for h in out["singularities"]["at"]))

print("\n=== 4. relax is gated: it never makes the mesh worse ===")
before = session.triple()
out = call("relax")
after = session.triple()
check("ok", out.get("ok") is True, out.get("reason"))
if out["accepted"] is None:
    check("an unaccepted relax leaves the mesh untouched", before == after, (before, after))
    check("and reports no change rather than a failure", out["all_improved"] is None)
    check("and says a null accept is a pass", "pass" in out["note"].lower(),
          out["note"][:50])
else:
    check("an accepted relax improved all three", out["all_improved"] is True, out)
    check("min angle rose", after[0] >= before[0], (before, after))
    check("max angle fell", after[1] <= before[1], (before, after))
    check("aspect fell", after[2] <= before[2], (before, after))
check("a bad seams value is refused", call("relax", seams="sideways").get("ok") is False)

print("\n=== 5. an ungated pass snapshots itself, and undo restores exactly ===")
positions = dict(session.positions())
out = call("smooth_boundary_constrained", kmax=30)
check("ok", out.get("ok") is True, out.get("reason"))
check("it reports whether all three improved", "all_improved" in out)
check("an undo is offered", out["undo_available"] is True)
moved = sum(1 for k, v in session.positions().items() if v != positions[k])
check("it actually moved vertices", moved > 0, moved)
undone = call("undo")
check("undo succeeds", undone.get("ok") is True, undone.get("reason"))
check("every vertex is back exactly where it was", session.positions() == positions)
check("a bad algorithm is refused",
      call("smooth_boundary_constrained", algorithm="wishful").get("ok") is False)

print("\n=== 6. smooth_region: selectors, and the blend guard ===")
out = call("smooth_region", region={"kind": "worst_faces", "count": 3}, kmax=20)
check("worst_faces works", out.get("ok") is True, out.get("reason"))
check("it says what it selected", "worst face" in out.get("selected", ""),
      out.get("selected"))
call("undo")
check("blend=0 is refused, with the crease as the reason",
      "crease" in call("smooth_region", region={"kind": "all"}, blend=0).get("reason", ""))
check("an unknown region kind is refused",
      call("smooth_region", region={"kind": "banana"}).get("ok") is False)
check("a malformed region is refused",
      call("smooth_region", region={"kind": "near_point"}).get("ok") is False)

print("\n=== 7. remarks and history ===")
check("an empty remark is refused", call("remark", text="   ").get("ok") is False)
check("a remark on a nonexistent place is refused",
      call("remark", text="x", about="v:999.000,999.000,0.000").get("ok") is False)
out = call("remark", text="relax was enough; the rest is the layout.")
check("a remark is kept", out.get("ok") is True and out["total"] == 1, out)
hist = call("history")
check("history lists the steps", len(hist["steps"]) >= 2, len(hist["steps"]))
check("history carries the remark", len(hist["remarks"]) == 1)
check("each step records before and after",
      all("before" in s and "after" in s for s in hist["steps"]))

print("\n=== 8. save and recall a worked example ===")
check("an unknown verdict is refused",
      call("save_example", name="x", verdict="excellent").get("ok") is False)
out = call("save_example", name="British Museum dense", verdict="good",
           lesson="relax alone was enough.", instruction="improve this mesh")
check("saved", out.get("ok") is True, out.get("reason"))
check("it reports where", os.path.isfile(out.get("path", "")), out.get("path"))
found = call("recall_examples")
check("recall finds it", len(found["examples"]) == 1, found)
check("recall ranks by fingerprint", "score" in found["examples"][0])
check("the example is readable as a resource",
      "British" in (handler.read_resource(found["examples"][0]["uri"]) or ""))
second = call("save_example", name="British Museum dense", verdict="bad")
check("saving twice keeps both rather than overwriting",
      second.get("path") != out.get("path"), second.get("path"))
check("a bad example is recalled too, on purpose",
      len(call("recall_examples", limit=5)["examples"]) == 2)

print("\n=== 9. save_mesh round-trips in double precision ===")
check("save_mesh refuses until the mesh has been looked at",
      call("save_mesh", path=os.path.join(tmp, "out.json")).get("ok") is False)
call("inspect", image=True)
out = call("save_mesh", path=os.path.join(tmp, "out.json"))
check("and works once it has", out.get("ok") is True, out.get("reason"))
reloaded = call("load_mesh", path=os.path.join(tmp, "out.json"))
check("reloads with the same faces", reloaded["faces"] == out["faces"])
call("inspect", image=True)
check("writing into a missing directory is refused",
      call("save_mesh", path=os.path.join(tmp, "nope", "x.json")).get("ok") is False)

print("\n=== 10. the rhino tools refuse gracefully with no Rhino ===")
out = call("rhino_status")
check("status answers without Rhino", out.get("ok") is True, out)
check("and says nothing is attached", out["attached"] is False)
check("and names the command that fixes it", "CMD_mcp_link" in out["reading"])
out = call("rhino_pull", timeout=0.2)
check("pull refuses rather than hanging", out.get("ok") is False)
check("and says why", "not listening" in out.get("reason", ""), out.get("reason", "")[:40])
out = call("rhino_push", timeout=0.2)
check("push refuses rather than hanging", out.get("ok") is False)

print("\n=== 11. bad calls are refusals, not crashes ===")
check("an unknown tool", call("no_such_tool").get("ok") is False)
check("it lists what there is", "available" in call("no_such_tool"))
out = call("inspect", nonsense=1)
check("an unexpected argument is refused", out.get("ok") is False,
      out.get("reason", "")[:44])
check("the message names the TOOL, not the implementation",
      "inspect" in out.get("reason", "") and "_t_" not in out.get("reason", ""),
      out.get("reason", "")[:60])

sys.exit(report("test_tools_offline"))
