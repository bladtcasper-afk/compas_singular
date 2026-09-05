"""Headless exercise of the agent core -- addressing, rebuild, state machine.

No model API, no key, no Rhino. A scripted plan drives exactly the functions a
model would.

Run:
    C:/Users/Casper/anaconda3/envs/singular312/python.exe test_session.py
"""
import os
import sys
from math import cos, pi, sin

REPO = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular"
for path in (os.path.join(REPO, "src"),):
    if path not in sys.path:
        sys.path.insert(0, path)

from compas_singular.framefield.decomposition import FieldDecomposition   # noqa: E402
from compas_singular.agent.core import AddressBook                       # noqa: E402
from compas_singular.agent.core import DensifyRefused                    # noqa: E402
from compas_singular.agent.core import MeshEditSession                   # noqa: E402
from compas_singular.agent.core import StageError                        # noqa: E402
from compas_singular.agent.core import edge_address                      # noqa: E402
from compas_singular.agent.core import point_address                     # noqa: E402
from compas_singular.agent.core.rebuild import densify_preserving        # noqa: E402


FAILURES = []


def check(name, condition, detail=""):
    print("  {} {}{}".format("ok  " if condition else "FAIL", name,
                             (" -- " + str(detail)) if detail else ""))
    if not condition:
        FAILURES.append(name)


def square(size=10.0):
    return [[0.0, 0.0, 0.0], [size, 0.0, 0.0], [size, size, 0.0],
            [0.0, size, 0.0], [0.0, 0.0, 0.0]]


def circle(cx, cy, r, n=24):
    pts = [[cx + r * cos(2 * pi * i / n), cy + r * sin(2 * pi * i / n), 0.0]
           for i in range(n)]
    return pts + pts[:1]


def build(inners=None, target_length=1.0):
    outer = square()
    d = FieldDecomposition.from_boundary(outer, inners, target_length=target_length)
    return d, outer


# ======================================================================
# addressing
# ======================================================================

print("\naddressing")

check("point_address rounds", point_address([1.00049, 2.0, 0.0]) == (1.0, 2.0, 0.0),
      point_address([1.00049, 2.0, 0.0]))
check("point_address normalises -0.0",
      point_address([-0.0001, 0.0, 0.0]) == (0.0, 0.0, 0.0),
      point_address([-0.0001, 0.0, 0.0]))
check("point_address pads 2D", point_address([1.0, 2.0]) == (1.0, 2.0, 0.0))

d, outer = build()
d.decomposition_mesh()
coarse = d.mesh

book = AddressBook(coarse)
check("book covers the mesh",
      book.counts()["vertices"] == coarse.number_of_vertices()
      and book.counts()["faces"] == coarse.number_of_faces(), book.counts())

again = AddressBook(coarse)
check("same geometry gives the same labels", again.vertices == book.vertices)

label = sorted(book.vertices)[0]
check("label resolves to a live vertex", book.vertex(label) in coarse.vertex,
      "{} -> {}".format(label, book.vertex(label)))
check("round trip label -> key -> label",
      book.label_of_vertex(book.vertex(label)) == label)

e_label = sorted(book.edges)[0]
u, v = book.edge(e_label)
check("edge address is order-independent",
      edge_address(coarse, (u, v)) == edge_address(coarse, (v, u)))

# carry onto an identical mesh: everything kept, nothing moved.
carried, report = book.carry(coarse.copy())
check("carry onto an identical mesh keeps everything", not report["lost"],
      "{} lost".format(len(report["lost"])))
check("  and nothing moved", not report["moved"], report["moved"][:3])

# carry onto a mesh with a corner moved: exactly that one is lost.
moved = coarse.copy()
vkey = list(moved.vertices())[0]
p = moved.vertex_coordinates(vkey)
moved.vertex_attributes(vkey, "xyz", [p[0] + 2.0, p[1] + 2.0, p[2]])
_, report = book.carry(moved)
check("carry reports a moved corner as lost", len(report["lost"]) >= 1,
      "{} lost".format(len(report["lost"])))
_, near = book.carry(moved, tol=5.0)
check("  and tol can recover it", len(near["lost"]) < len(report["lost"]),
      "{} lost with tol".format(len(near["lost"])))


# ======================================================================
# rebuild -- densities preserved, no fallback
# ======================================================================

print("\nrebuild")

d, outer = build()
d.decomposition_mesh()
dense, report = densify_preserving(d, target_length=2.0)
check("densified", dense.number_of_faces() > 0, dense.number_of_faces())
check("route is field", report["route"] == "field", report["route"])
check("report carries metrics", "min_angle" in report["metrics"],
      sorted(report["metrics"])[:4])
check("report carries densities by address", len(report["densities"]) > 0,
      "{} strip(s)".format(len(report["densities"])))

# The whole point: change ONE strip and rebuild -- the others keep their value.
densities = dict(report["densities"])
target = sorted(densities)[0]
densities[target] = densities[target] + 3
dense2, report2 = densify_preserving(d, densities=densities, target_length=2.0)
check("per-strip density applied",
      report2["densities"][target] == densities[target],
      "{} -> {}".format(report["densities"][target], report2["densities"][target]))
check("  and the mesh grew", dense2.number_of_faces() > dense.number_of_faces(),
      "{} -> {}".format(dense.number_of_faces(), dense2.number_of_faces()))
check("  with nothing unmatched", not report2["density_report"]["unmatched"],
      report2["density_report"]["unmatched"][:2])

others_held = all(report2["densities"][a] == report["densities"][a]
                  for a in report["densities"] if a != target)
check("every other strip kept its density", others_held)

stale = dict(densities)
stale[(999.0, 999.0, 0.0)] = 4
_, report3 = densify_preserving(d, densities=stale, target_length=2.0)
check("an address naming no strip is reported, not raised",
      len(report3["density_report"]["unmatched"]) == 1,
      report3["density_report"]["unmatched"])

layout_before = d.mesh
edited_before = d._edited
try:
    densify_preserving(d, density=0, target_length=None)
    refused = False
    reason = ""
except DensifyRefused as exc:
    refused = True
    reason = exc.reason
except Exception as exc:                                    # noqa: BLE001
    refused = False
    reason = "{}: {}".format(type(exc).__name__, exc)
check("a density of 0 is refused rather than falling back", refused, reason)
check("  and the layout object is untouched", d.mesh is layout_before)
check("  and _edited is unchanged",
      d._edited == edited_before, "{} -> {}".format(edited_before, d._edited))


# ----------------------------------------------------------------------
# the claim rebuild.py is FOR: quad_mesh would have destroyed the layout.
# Same decomposition, same bad density, both routes.
# ----------------------------------------------------------------------

print("\nrebuild vs quad_mesh on the same bad density")

d_a, _ = build()
d_a.decomposition_mesh()
d_a._edited = True                 # stand in for a hand edit
layout_a = d_a.mesh
try:
    d_a.quad_mesh(density=0)
    qm_raised = False
except Exception as exc:                                    # noqa: BLE001
    qm_raised = True
    qm_reason = "{}: {}".format(type(exc).__name__, exc)

print("    quad_mesh: raised={}, route={!r}, _edited={}, mesh is layout={}".format(
    qm_raised, d_a.route() if not qm_raised else "-",
    d_a._edited, d_a.mesh is layout_a))

d_b, _ = build()
d_b.decomposition_mesh()
d_b._edited = True
layout_b = d_b.mesh
try:
    densify_preserving(d_b, density=0)
    dp_refused = False
except DensifyRefused:
    dp_refused = True

check("densify_preserving refused", dp_refused)
check("  and kept the layout object", d_b.mesh is layout_b)
check("  and kept the _edited flag", d_b._edited is True)
check("quad_mesh did NOT refuse (it fell back)", not qm_raised,
      qm_reason if qm_raised else "")
if not qm_raised:
    check("  and that is what costs the layout: route or _edited changed",
          d_a.route() != "field" or d_a._edited is False or d_a.mesh is not layout_a,
          "route={!r}, _edited={}, same layout={}".format(
              d_a.route(), d_a._edited, d_a.mesh is layout_a))


# ======================================================================
# the state machine
# ======================================================================

print("\nstate machine")

d, outer = build()
s = MeshEditSession(d, outer=outer, target_length=2.0)
check("starts at field", s.stage == "field", s.stage)
check("no mesh at field", s.mesh is None)
check("no editor at field", s.editor is None)
check("inputs complete", s.inputs()["complete"])

s.enter_coarse()
check("entered coarse", s.stage == "coarse")
check("coarse editor exists", s.coarse_editor is not None)
check("mesh is the layout", s.mesh is s.coarse_editor.mesh)
check("book built for the layout",
      s.book().counts()["faces"] == s.mesh.number_of_faces())

editor = s.coarse_editor
s.enter_coarse()
check("enter_coarse is idempotent", s.coarse_editor is editor)

dense, rep = s.rebuild_dense()
check("entered dense", s.stage == "dense")
check("dense editor exists", s.dense_editor is not None)
check("mesh is the dense mesh", s.mesh is s.dense_editor.mesh)
check("densities captured", len(s.densities) > 0, len(s.densities))
check("history records the rebuild",
      any(e["action"] == "rebuild_dense" for e in s.history))

try:
    s.enter_coarse()
    raised = False
except StageError as exc:
    raised = True
    msg = str(exc)
check("enter_coarse from dense is refused", raised, "" if raised else "no StageError")
check("  and points at revert_to_coarse", "revert_to_coarse" in msg)

cost = s.discarded_by("revert_to_coarse")
check("revert is priced before it happens", cost["allowed"] and cost["final"],
      cost)
check("  and names the dense face count",
      cost["dense_faces"] == dense.number_of_faces(), cost["dense_faces"])

s.revert_to_coarse()
check("reverted", s.stage == "coarse")
check("dense editor gone", s.dense_editor is None)
check("decomposition.dense cleared", s.decomposition.dense is None)

cost = s.discarded_by("adopt_decomposition")
check("adopt is priced too", cost["allowed"] and cost["final"])
check("  and names the coarse face count", cost["coarse_faces"] > 0,
      cost["coarse_faces"])

d2, _ = build()
s.adopt_decomposition(d2)
check("adopting resets to field", s.stage == "field")
check("  layout dropped", s.coarse_editor is None)
check("  densities dropped", s.densities == {})
check("  and the new decomposition is in place", s.decomposition is d2)

check("unknown transition is refused",
      s.discarded_by("nonsense")["allowed"] is False)


# ======================================================================
# undo
# ======================================================================

print("\nundo")

d, outer = build()
s = MeshEditSession(d, outer=outer, target_length=2.0)
s.enter_coarse()
faces_before = s.mesh.number_of_faces()

check("nothing to undo yet", s.can_undo() is False)
ok, why = s.undo()
check("undo on an empty stack is a refusal, not a crash", ok is False, why)

s.snapshot("before dense")
check("can undo now", s.can_undo())

s.rebuild_dense()
check("at dense", s.stage == "dense")

ok, label = s.undo()
check("undo restored", ok, label)
check("  the stage", s.stage == "coarse", s.stage)
check("  the layout", s.mesh.number_of_faces() == faces_before,
      "{} vs {}".format(s.mesh.number_of_faces(), faces_before))
check("  and the label", label == "before dense", label)
check("editor came back wired to its decomposition",
      s.coarse_editor.decomposition is s.decomposition)


# ======================================================================
# a hole, and a coarse edit carried across a commit
# ======================================================================

print("\ncommit renumbers -- carry reports it")

d, outer = build(inners=[circle(5.0, 5.0, 2.0)], target_length=1.0)
s = MeshEditSession(d, outer=outer, inners=[circle(5.0, 5.0, 2.0)],
                    target_length=2.0)
s.enter_coarse()
before = s.book()
# All three kinds: carry() reports them in ONE dict, and the prefixes keep the
# labels unique across kinds, so the accounting below has to span all three.
n_labels = (len(before.vertices) + len(before.edges) + len(before.faces))

vkey = next(v for v in s.mesh.vertices() if not s.mesh.is_vertex_on_boundary(v))
p = s.mesh.vertex_coordinates(vkey)
moved_ok = s.coarse_editor.move_vertex(vkey, [p[0] + 0.3, p[1] + 0.3, p[2]])
check("moved an interior corner", moved_ok, s.coarse_editor.last_reason)

ok, notes = s.coarse_editor.commit()
check("committed", ok, notes.get("error", ""))

carry = s.refresh_book()
check("carry ran after the commit", carry != {})
kept = len(carry["kept"])
print("    {} of {} labels kept, {} moved, {} lost".format(
    kept, n_labels, len(carry["moved"]), len(carry["lost"])))
check("carry accounts for every old label exactly once",
      kept + len(carry["lost"]) == n_labels,
      "kept {} + lost {} != {}".format(kept, len(carry["lost"]), n_labels))
check("every kept label names a label that exists now",
      all(new in s.book().vertices or new in s.book().edges
          or new in s.book().faces for new in carry["kept"].values()))
check("moved is a subset of kept",
      set(carry["moved"]) <= set(carry["kept"]))
check("last_carry is readable afterwards", s.last_carry == carry)

dense, rep = s.rebuild_dense()
check("edited layout still densifies", dense.number_of_faces() > 0,
      "{} faces, route {}".format(dense.number_of_faces(), rep["route"]))
check("state() reports the stage", s.state()["stage"] == "dense", s.state())


# ======================================================================

# ======================================================================
# adopt_dense -- editing a mesh that came from somewhere else
# ======================================================================

print("\nadopt_dense")

from compas_singular.datastructures import QuadMesh                      # noqa: E402

d, outer = build()
s = MeshEditSession(d, outer=outer, target_length=2.0)
s.enter_coarse()
made, _rep = s.rebuild_dense(target_length=2.0)

# Round-trip it the way a document would: vertices and faces, nothing else.
verts, faces = made.to_vertices_and_faces()
foreign = QuadMesh.from_vertices_and_faces(verts, faces)

d2, outer2 = build()
s2 = MeshEditSession(d2, outer=outer2, target_length=2.0)
s2.adopt_dense(foreign)
check("adopting goes straight to dense", s2.stage == "dense", s2.stage)
check("  with an editor", s2.dense_editor is not None)
check("  and the mesh is the one handed in", s2.mesh is foreign)
check("  with no layout built", s2.coarse_editor is None)
check("  recorded in the history",
      any(e["action"] == "adopt_dense" for e in s2.history))
check("  and the decomposition knows about it", s2.decomposition.dense is foreign)

cost = s2.discarded_by("revert_to_coarse")
check("reverting is REFUSED -- no layout reproduces an adopted mesh",
      cost["allowed"] is False, cost)
check("  and says why", "no layout to go back to" in cost["reason"],
      cost["reason"][:70])

try:
    s2.revert_to_coarse()
    raised = False
except StageError:
    raised = True
check("calling it anyway raises rather than corrupting the stage", raised)
check("  leaving the session at dense", s2.stage == "dense", s2.stage)

# The tool must turn that into a refusal, not let it escape into the loop.
from compas_singular.agent.core import tools as _tools                   # noqa: E402
r = _tools.call(s2, "revert_to_coarse")
check("the tool reports it as a refusal, not a crash", r["ok"] is False)
check("  carrying the reason", "no layout" in r["reason"], r["reason"][:60])

# Editing still works on the adopted mesh -- that is the point of it.
book = s2.book()
edge = sorted(book.edges)[0]
plan = _tools.call(s2, "plan_dense_strip_deletion", label=edge)
check("dense tools work on an adopted mesh", plan["ok"], plan.get("reason"))

# A dense-start session is genuinely the end of the line: revert has no layout
# to reach, and enter_coarse refuses at the dense stage. So the refusal must NOT
# advise enter_coarse as a way out -- that would send a caller round a loop.
check("the refusal does not advise something that also refuses",
      "enter_coarse" not in cost["reason"], cost["reason"][-60:])
try:
    s2.enter_coarse()
    raised = False
except StageError:
    raised = True
check("  and enter_coarse from dense does indeed refuse", raised)


print("\n{} check(s) failed{}".format(
    len(FAILURES), (": " + ", ".join(FAILURES)) if FAILURES else ""))
sys.exit(1 if FAILURES else 0)
