"""Drive all six stages through the tool surface -- scripted, no model, no key.

Every call goes through ``tools.call(session, name, **kwargs)``, which is
exactly what the runner will do. If this passes, the runner has nothing left to
get wrong except the prompt.

Run:
    C:/Users/Casper/anaconda3/envs/singular312/python.exe test_tools_offline.py
"""
import os
import struct
import sys
from math import cos, pi, sin

REPO = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular"
for path in (os.path.join(REPO, "src"),):
    if path not in sys.path:
        sys.path.insert(0, path)

from compas_singular.framefield.decomposition import FieldDecomposition   # noqa: E402
from compas_singular.agent.core import MeshEditSession                   # noqa: E402
from compas_singular.agent.core import render_png, render_svg            # noqa: E402
from compas_singular.agent.core import tools                             # noqa: E402
from compas_singular.agent.core.digest import digest, render_text        # noqa: E402


FAILURES = []
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_out")


def check(name, condition, detail=""):
    print("  {} {}{}".format("ok  " if condition else "FAIL", name,
                             (" -- " + str(detail)) if detail else ""))
    if not condition:
        FAILURES.append(name)


def square(size=10.0):
    return [[0.0, 0.0, 0.0], [size, 0.0, 0.0], [size, size, 0.0],
            [0.0, size, 0.0], [0.0, 0.0, 0.0]]


def circle(cx, cy, r, n=32):
    pts = [[cx + r * cos(2 * pi * i / n), cy + r * sin(2 * pi * i / n), 0.0]
           for i in range(n)]
    return pts + pts[:1]


def run(session, name, **kwargs):
    """Call a tool and print a one-line result."""
    result = tools.call(session, name, **kwargs)
    flag = "ok  " if result.get("ok") else "--  "
    print("    {} {:<26} {}".format(flag, name, result.get("reason", "")[:64]))
    return result


# ======================================================================

print("\nregistry")

check("tools registered", len(tools.tool_names()) >= 25, len(tools.tool_names()))
check("every schema is an object with no extras",
      all(s["input_schema"]["type"] == "object"
          and s["input_schema"]["additionalProperties"] is False
          for s in tools.schemas()))
check("every required arg is declared",
      all(r in s["input_schema"]["properties"]
          for s in tools.schemas() for r in s["input_schema"]["required"]))
check("descriptions are substantial",
      all(len(s["description"]) > 60 for s in tools.schemas()),
      [s["name"] for s in tools.schemas() if len(s["description"]) <= 60])
check("the four destructive tools are marked confirm",
      set(n for n in tools.tool_names() if tools.TOOLS[n].confirm)
      == {"solve_field", "set_guides", "commit", "revert_to_coarse"},
      sorted(n for n in tools.tool_names() if tools.TOOLS[n].confirm))

check("unknown tool is a refusal, not a raise",
      tools.call(None, "no_such_tool")["ok"] is False)


# ======================================================================

print("\nstage 0 -- the field")

outer = square()
inners = [circle(5.0, 5.0, 2.0)]
d = FieldDecomposition.from_boundary(outer, inners, target_length=1.0)
s = MeshEditSession(d, outer=outer, inners=inners, target_length=2.0)

r = run(s, "check_inputs")
check("check_inputs ran", r["ok"], r.get("reason"))
check("  and passed a well-sampled input", r.get("faults") == 0, r.get("faults"))
check("  and counted the square's corners", r.get("corners", 0) >= 4,
      r.get("corners"))

r = run(s, "inspect_field")
check("inspect_field ran", r["ok"])
check("  reports poincare_hopf", "poincare_hopf" in r["field"])
check("  warns that ok is not a verdict",
      "poincare_hopf_means" in r["field"])
check("  counts singularities", "singularities" in r["field"],
      r["field"].get("singularities"))

# A badly sampled hole must be caught at the input, not blamed on the solver.
bad_inner = [[5 + 4.2 * cos(2 * pi * i / 48) * 0.4,
              5 + 0.5 * sin(2 * pi * i / 48) * 0.4, 0.0] for i in range(48)]
s_bad = MeshEditSession(d, outer=outer, inners=[bad_inner], target_length=2.0)
r = run(s_bad, "check_inputs")
check("an aliased hole is refused at the input", r["ok"] is False)
check("  with the resample advice, not a solver knob",
      "RESAMPLE" in r.get("advice", "") and "not relax" in r.get("advice", ""),
      r.get("advice", "")[:80])


# ======================================================================

print("\nstage gating")

r = run(s, "smooth")
check("a dense tool at the field stage is refused", r["ok"] is False)
check("  and names the stage", "only defined at stage" in r["reason"], r["reason"])


# ======================================================================

print("\nstage 1 -- the coarse layout")

r = run(s, "enter_coarse")
check("entered coarse", r["ok"], r.get("reason"))
check("  and reports the patch count", r["faces"] > 0, r["faces"])

r = run(s, "inspect")
check("inspect gives a coarse section", "coarse" in r, sorted(r))
check("  with labels while they are few enough",
      "labels" in r["coarse"] or "labels_note" in r["coarse"])

book = s.book()
v_labels = sorted(book.vertices)
interior = [lab for lab in v_labels
            if not s.mesh.is_vertex_on_boundary(book.vertex(lab))]
if interior:
    lab = interior[0]
    x, y, _z = s.mesh.vertex_coordinates(book.vertex(lab))
    r = run(s, "move_corner", label=lab, x=x + 0.25, y=y + 0.25)
    check("moved an interior corner", r["ok"], r.get("reason"))
else:
    print("    (no interior corner on this layout -- move_corner not exercised)")

r = run(s, "move_corner", label="v9999", x=0.0, y=0.0)
check("an unknown label is a refusal", r["ok"] is False)
check("  naming the label", "v9999" in r["reason"], r["reason"])

e_labels = sorted(s.book().edges)
r = run(s, "plan_strip_deletion", label=e_labels[0])
check("plan_strip_deletion ran", r["ok"], r.get("reason"))
plan = r["plan"]
check("  the plan performed a trial", "ok" in plan and "faces" in plan,
      sorted(plan))
check("  and reports boundaries_lost as a count",
      isinstance(plan["boundaries_lost"], int), plan["boundaries_lost"])

faces_before = s.mesh.number_of_faces()
r = run(s, "remove_strip", label=e_labels[0])
if plan["ok"]:
    check("remove_strip agreed with its plan", r["ok"], r.get("reason"))
    check("  and shrank the layout", s.mesh.number_of_faces() < faces_before,
          "{} -> {}".format(faces_before, s.mesh.number_of_faces()))
else:
    check("remove_strip refused where the plan said it would", r["ok"] is False)
    check("  leaving the layout alone",
          s.mesh.number_of_faces() == faces_before)

r = run(s, "reset_coarse")
check("reset_coarse restored the layout",
      r["ok"] and s.mesh.number_of_faces() == faces_before,
      "{} vs {}".format(s.mesh.number_of_faces(), faces_before))

r = run(s, "commit")
check("commit ran", r["ok"], r.get("reason"))
check("  and reissued labels", "relabelled" in r, sorted(r))
check("  reporting what moved and what was lost",
      "moved" in r["relabelled"] and "lost" in r["relabelled"],
      r["relabelled"])


# ======================================================================

print("\nstages 2 and 3 -- densify and densities")

r = run(s, "list_strips")
check("list_strips ran", r["ok"], r.get("reason"))
check("  and gave geometric addresses",
      r["count"] > 0 and len(r["strips"][0]["address"]) == 3,
      r["strips"][0] if r["strips"] else None)
strips = r["strips"]

r = run(s, "rebuild_dense", target_length=1.5)
check("rebuild_dense ran", r["ok"], r.get("reason"))
check("  at stage dense", s.stage == "dense", s.stage)
check("  on the field route", r.get("route") == "field", r.get("route"))
base_faces = r["faces"]

address = strips[0]["address"]
before = strips[0]["density"]
r = run(s, "set_strip_density", address=address, density=before + 4)
check("set_strip_density recorded", r["ok"], r.get("reason"))
check("  as pending", "next rebuild_dense" in r.get("note", ""), r.get("note"))

r = run(s, "rebuild_dense", target_length=1.5)
check("rebuild applied the per-strip density", r["ok"], r.get("reason"))
check("  and nothing was unmatched", r.get("unmatched") == 0, r.get("unmatched"))
check("  and the mesh changed size", r["faces"] != base_faces,
      "{} -> {}".format(base_faces, r["faces"]))

r = run(s, "set_strip_density", address=address, density=0)
check("a density of 0 is refused up front", r["ok"] is False, r.get("reason"))

r = run(s, "set_field_aware", enabled=False)
check("field_aware can be turned off", r["ok"] and r["field_aware"] is False)
check("  and says when it takes effect", "rebuild_dense" in r["note"])
run(s, "set_field_aware", enabled=True)

r = run(s, "set_target_length", target_length=-1)
check("a negative target_length is refused", r["ok"] is False)


# ======================================================================

print("\nstage 4 -- smoothing")

r = run(s, "smooth")
check("smooth ran", r["ok"], r.get("reason"))
check("  and reported before/after", "before" in r and "after" in r,
      "{} -> {}".format(r.get("before"), r.get("after")))
if r.get("accepted") is None:
    check("  a rejected pass says the mesh is unchanged",
          "unchanged" in r.get("note", ""), r.get("note"))
else:
    b, a = r["before"], r["after"]
    check("  an accepted pass improved all three",
          a[0] >= b[0] and a[1] <= b[1] and a[2] <= b[2],
          "{} -> {}".format(b, a))


# ======================================================================

print("\nstage 5 -- smoothing to a guide")

guide = [[0.4, 5.0, 0.0], [9.6, 5.0, 0.0]]
r = run(s, "select_guide_chain", guide=guide, max_angle=90)
if r["ok"]:
    check("a chain was selected", r["chain"] > 0, r["chain"])
    check("  measured before attaching", "quality" in r, sorted(r))
    check("  and the outline gate is reported",
          "boundary_interior" in r["quality"], sorted(r["quality"]))
    r2 = run(s, "attach_guide_chain", hold="fixed", kmax=20)
    check("attached", r2["ok"], r2.get("reason"))
    check("  boundary vertices slid rather than moved",
          r2.get("slid", 0) >= 0 and r2.get("moved", 0) <= r["chain"],
          "moved {} of {}".format(r2.get("moved"), r["chain"]))
else:
    print("    (no chain follows that guide here: {})".format(r["reason"][:60]))
    check("a failed selection is a refusal, not a raise", r["ok"] is False)
    r2 = run(s, "attach_guide_chain")
    check("attach without a selection is refused", r2["ok"] is False)


# ======================================================================

print("\nstage 6 -- the dense mesh")

e_labels = sorted(s.book().edges)
r = run(s, "plan_dense_strip_deletion", label=e_labels[0])
check("plan_dense_strip_deletion ran", r["ok"], r.get("reason"))
dense_plan = r["plan"]

faces_before = s.mesh.number_of_faces()
r = run(s, "remove_line", label=e_labels[0])
check("remove_line agreed with its plan", r["ok"] == dense_plan["ok"],
      "tool {} vs plan {}".format(r["ok"], dense_plan["ok"]))
if r["ok"]:
    check("  and shrank the mesh", s.mesh.number_of_faces() < faces_before,
          "{} -> {}".format(faces_before, s.mesh.number_of_faces()))

e_labels = sorted(s.book().edges)
faces_before = s.mesh.number_of_faces()
added = None
for lab in e_labels[:12]:
    rr = tools.call(s, "add_line", label=lab)
    if rr["ok"]:
        added = rr
        break
if added:
    check("add_line grew the mesh", s.mesh.number_of_faces() > faces_before,
          "{} -> {}".format(faces_before, s.mesh.number_of_faces()))
    pair = added["addition"].get("new_vertices") or []
    check("  and relaxed the new strip open",
          added["addition"].get("relaxed") is True, added["addition"].get("relaxed"))
else:
    check("every add_line refusal gave a reason",
          bool(s.dense_editor.last_reason), s.dense_editor.last_reason[:60])

v_labels = sorted(s.book().vertices)
lab = v_labels[0]
x, y, _z = s.mesh.vertex_coordinates(s.book().vertex(lab))
r = run(s, "move_vertex", label=lab, x=x + 0.05, y=y + 0.05)
check("move_vertex ran on the dense mesh", r["ok"], r.get("reason"))


# ======================================================================

print("\nsession control")

r = run(s, "discarded_by", transition="revert_to_coarse")
check("discarded_by prices the revert", r["ok"] and r["final"] is True, r)

r = run(s, "snapshot", label="before revert")
check("snapshot pushed", r["ok"] and r["depth"] >= 1, r.get("depth"))

r = run(s, "revert_to_coarse")
check("reverted", r["ok"] and s.stage == "coarse", s.stage)

r = run(s, "undo")
check("undo came back to dense", r["ok"] and s.stage == "dense",
      "{} / {}".format(r.get("restored"), s.stage))


# ======================================================================

print("\ndigest and pictures")

data = digest(s)
text = render_text(data)
check("digest renders as text", len(text) > 200, "{} chars".format(len(text)))
check("  and mentions the stage", "stage" in text)
check("  and carries the dense quality", "min_angle" in text)
check("  with hard_floor beside the route",
      "hard_floor" in text and "route" in text)

if not os.path.isdir(OUT):
    os.makedirs(OUT)

png = render_png(s, path=os.path.join(OUT, "dense.png"), width=700, height=700)
check("render_png produced a PNG", png[:8] == b"\x89PNG\r\n\x1a\n", png[:8])
w, h = struct.unpack(">II", png[16:24])
check("  of the requested size", (w, h) == (700, 700), (w, h))
check("  and a plausible size on disk", 1000 < len(png) < 4_000_000, len(png))

labels = sorted(s.book().vertices)[:6]
png2 = render_png(s, path=os.path.join(OUT, "labelled.png"), labels=labels)
check("labels change the picture", png2 != png)

svg = render_svg(s, path=os.path.join(OUT, "dense.svg"))
check("render_svg produced SVG", svg.startswith("<svg") and svg.endswith("</svg>"))
check("  containing the mesh", svg.count("<line") > 10, svg.count("<line"))

s2 = MeshEditSession(FieldDecomposition.from_boundary(square(), target_length=1.0),
                     outer=square())
png3 = render_png(s2, path=os.path.join(OUT, "field.png"), width=500, height=500)
check("the field stage renders too", png3[:8] == b"\x89PNG\r\n\x1a\n")

print("    wrote {}".format(OUT))


# ======================================================================

print("\n{} check(s) failed{}".format(
    len(FAILURES), (": " + ", ".join(FAILURES)) if FAILURES else ""))
sys.exit(1 if FAILURES else 0)
