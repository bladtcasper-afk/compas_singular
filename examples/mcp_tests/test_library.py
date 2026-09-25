"""Guidance, prompts, thresholds and the corpus. Verification step C2.

The library is what makes this approach different from ``agent``, where the
equivalent is Python string constants. So the things worth checking are that it
is really editable (the override works), that it degrades to nothing rather than
failing, and that a name off the wire cannot escape the folder.
"""
import json
import os
import shutil
import sys
import tempfile

from _harness import check, report

from compas_singular.mcp import library

print("\n=== 1. the shipped library loads ===")
shipped = library.library_root()
check("it is inside the package", os.path.isdir(shipped), shipped)
names = [r["uri"] for r in library.list_resources()]
for wanted in ("guidance://quality", "guidance://smoothing", "guidance://remarks"):
    check("ships {}".format(wanted), wanted in names, names)
prompts = [p["name"] for p in library.list_prompts()]
for wanted in ("improve_pulled_mesh", "diagnose_only", "compare_with_agent"):
    check("ships prompt {}".format(wanted), wanted in prompts, prompts)
check("every resource has a description",
      all(r["description"] for r in library.list_resources()))
body = library.read_resource("guidance://smoothing")
check("guidance reads back", body and len(body) > 500, len(body or ""))
check("and carries the measurement that matters",
      "accepted: null" in body and "pin_singularities" in body)

print("\n=== 2. thresholds are DATA, and merge per metric ===")
limits = library.thresholds()
check("min_angle bands load", limits["min_angle"]["good"] == 45.0, limits["min_angle"])
check("all four metrics are banded",
      set(limits) >= {"min_angle", "max_angle", "aspect_max", "share_below"})

folder = tempfile.mkdtemp(prefix="mcp-lib-")
os.environ["COMPAS_SINGULAR_MCP_LIBRARY"] = folder
check("the override redirects the root", library.library_root() == folder)
check("an empty library lists nothing rather than failing",
      library.list_resources() == [] and library.list_prompts() == [])
check("and thresholds fall back to the built-in defaults",
      library.thresholds()["min_angle"]["good"] == 45.0)

with open(os.path.join(folder, "thresholds.json"), "w") as stream:
    json.dump({"min_angle": {"good": 60.0}}, stream)
merged = library.thresholds()
check("a partial override is applied", merged["min_angle"]["good"] == 60.0)
check("without dropping the rest of that metric's bands",
      merged["min_angle"]["usable"] == 25.0, merged["min_angle"])
check("or the other metrics", merged["aspect_max"]["good"] == 2.0)

with open(os.path.join(folder, "thresholds.json"), "w") as stream:
    stream.write("{ this is not json")
check("a malformed thresholds file falls back rather than crashing",
      library.thresholds()["min_angle"]["good"] == 45.0)

print("\n=== 3. a name off the wire cannot escape the folder ===")
for bad in ("guidance://../../../secrets", "guidance://.hidden",
            "guidance://sub/dir"):
    check("refuses {}".format(bad), library.read_resource(bad) is None)
for bad in ("../../etc/passwd", ".hidden", "sub/dir"):
    check("refuses prompt {!r}".format(bad), library.get_prompt(bad) is None)
check("an unknown scheme is None", library.read_resource("evil://x") is None)

print("\n=== 4. examples save, render and recall ===")
record = {"name": "disc", "title": "A disc", "verdict": "good",
          "instruction": "improve it", "lesson": "relax was enough",
          "before": {"min_angle": 12.0, "max_angle": 170.0, "aspect_max": 9.0},
          "after": {"min_angle": 44.0, "max_angle": 130.0, "aspect_max": 2.0},
          "fingerprint": {"loops": 2, "faces_band": "small",
                          "dominant": "min_angle", "severity": "unusable"},
          "steps": [{"action": "relax", "arguments": {"seams": "free"}}],
          "remarks": ["the hole drove everything"]}
path = library.save_example(record)
check("saved", os.path.isfile(path), path)
check("it is listed", "example://disc" in [r["uri"] for r in library.list_resources()])
text = library.read_resource("example://disc")
check("it renders as markdown", text.startswith("# A disc"), text[:24])
check("with the verdict", "**Verdict: good**" in text)
check("a before/after table", "| before |" in text and "| after |" in text)
check("the steps", "`relax(seams=free)`" in text, text)
check("the remarks", "the hole drove everything" in text)
check("and the lesson", "relax was enough" in text)

print("\n=== 5. recall ranks by the shape of the PROBLEM ===")
near = dict(record, name="near", fingerprint=dict(record["fingerprint"]))
far = dict(record, name="far",
           fingerprint={"loops": 9, "faces_band": "large",
                        "dominant": "aspect_max", "severity": "good"})
library.save_example(near)
library.save_example(far)
ranked = library.recall(record["fingerprint"], limit=5)
check("the identically-shaped problem ranks first",
      ranked[0]["fingerprint"] == record["fingerprint"], ranked[0])
check("above the differently-shaped one",
      ranked[0]["score"] > ranked[-1]["score"], [r["score"] for r in ranked])
check("a verdict filter is honoured",
      library.recall(record["fingerprint"], verdicts=("bad",)) == [])
check("saving the same name twice keeps both",
      os.path.basename(library.save_example(record)) == "disc-2.json")

print("\n=== 6. a malformed example is skipped, not raised on ===")
with open(os.path.join(folder, "sessions", "broken.json"), "w") as stream:
    stream.write("{ nope")
check("load_examples skips it",
      "broken" not in [r["name"] for r in library.load_examples()])
check("and the rest still load", len(library.load_examples()) >= 3)

shutil.rmtree(folder, ignore_errors=True)
sys.exit(report("test_library"))
