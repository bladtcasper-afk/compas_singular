"""The spool: round trip, atomic claim, ordering, timeouts, crash litter.

Verification step A. No Rhino, no model, no network.
"""
import os
import sys

from _harness import check, report, spool_dir

folder = spool_dir()
from compas_singular.mcp.bridge import spool

print("\n=== 1. a request goes out and an answer comes back ===")
check("no link before one attaches", spool.link_state() is None)
rid = spool.post("pull", {"layer": "QuadMesh"})
check("the request is queued", spool.pending() == [rid], spool.pending())
job = spool.take()
check("take returns (id, verb, args)", job == (rid, "pull", {"layer": "QuadMesh"}), job)
check("a claimed request leaves the queue", spool.pending() == [])
check("a claimed request is not handed out twice", spool.take() is None)
spool.answer(rid, True, result={"faces": 12})
got = spool.wait(rid, timeout=2)
check("wait returns the result", got.get("ok") and got["result"] == {"faces": 12}, got)

print("\n=== 2. ids sort into arrival order, and take is oldest-first ===")
ids = [spool.post("ping") for _ in range(5)]
check("ids sort lexically into arrival order", spool.pending() == ids)
check("take returns the oldest", spool.take()[0] == ids[0])
for _ in range(4):
    spool.take()
check("queue drains", spool.pending() == [])

print("\n=== 3. a timeout is a refusal that names the cause ===")
missing = spool.wait("nosuchid", timeout=0.2)
check("times out rather than hanging", missing.get("timed_out") is True)
check("with no link, it says so", "not listening" in missing["reason"], missing["reason"][:50])
spool.heartbeat("plate.3dm", force=True)
state = spool.link_state()
check("the heartbeat records the document", state["document"] == "plate.3dm")
check("a fresh heartbeat is not stale", state["stale"] is False)
busy = spool.wait("nosuchid", timeout=0.2)
check("with a link, the reason is different", "inside a command" in busy["reason"],
      busy["reason"][:50])

print("\n=== 4. a refusal from Rhino survives the wire ===")
rid = spool.post("pull")
spool.take()
spool.answer(rid, False, error="No mesh on that layer.")
bad = spool.wait(rid, timeout=2)
check("ok is false", bad.get("ok") is False)
check("the reason is Rhino's own", bad["reason"] == "No mesh on that layer.", bad)

print("\n=== 5. a response is consumed, not replayed ===")
again = spool.wait(rid, timeout=0.2)
check("a second wait does not see the old answer", again.get("timed_out") is True)

print("\n=== 6. crash litter is pruned; the heartbeat is not ===")
spool.post("ping")
rid = spool.post("ping")
spool.take()                                    # leave a .req.taken behind
removed = spool.prune(max_age=-1)
check("stale files are removed", removed >= 2, removed)
check("the heartbeat survives pruning", spool.link_state() is not None)
spool.clear_link()
check("clear_link detaches", spool.link_state() is None)

print("\n=== 7. writes are atomic: no reader sees a partial file ===")
import json
rid = spool.post("push", {"mesh": {"vertices": [[0, 0, 0]] * 500}})
path = os.path.join(folder, rid + ".req.json")
check("the request file parses completely", isinstance(json.load(open(path)), dict))
names = [n for n in os.listdir(folder) if n.endswith(".tmp")]
check("no temporary files are left behind", names == [], names)

sys.exit(report("test_spool"))
