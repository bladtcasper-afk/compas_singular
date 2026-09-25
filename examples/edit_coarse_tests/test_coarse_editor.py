"""Headless exercise of CoarseLayoutEditor -- no Rhino, no prompts.

The whole point of the split: every operation the Rhino command offers has to
be reachable, and checkable, from a plain script.
"""
import os
import sys
from math import cos, pi, sin

REPO = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular"
for path in (os.path.join(REPO, "src"), os.path.join(REPO, "examples", "New approach")):
    if path not in sys.path:
        sys.path.insert(0, path)

from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402
from compas_singular.editing import CoarseLayoutEditor                   # noqa: E402


FAILURES = []


def check(name, condition, detail=""):
    print("  {} {}{}".format("ok  " if condition else "FAIL", name,
                             (" -- " + str(detail)) if detail else ""))
    if not condition:
        FAILURES.append(name)


def all_quad(mesh):
    """Every face 4-sided, except a pseudo-quad which is 3 and has a pole."""
    face_pole = mesh.attributes.get("face_pole") or {}
    bad = []
    for fkey in mesh.faces():
        n = len(mesh.face_vertices(fkey))
        if n == 4 and fkey not in face_pole:
            continue
        if n == 3 and fkey in face_pole:
            continue
        bad.append((fkey, n))
    return bad


def square(size=10.0):
    return [[0.0, 0.0, 0.0], [size, 0.0, 0.0], [size, size, 0.0],
            [0.0, size, 0.0], [0.0, 0.0, 0.0]]


def circle(cx, cy, r, n=24):
    pts = [[cx + r * cos(2 * pi * i / n), cy + r * sin(2 * pi * i / n), 0.0]
           for i in range(n)]
    return pts + pts[:1]


def build(inners=None):
    d = FieldDecomposition.from_boundary(square(), inners, target_length=1.0)
    d.decomposition_mesh()
    return d


def mesh_digest(mesh):
    """Enough of a mesh to prove a refusal changed nothing."""
    return (mesh.number_of_vertices(), mesh.number_of_faces(),
            sorted(tuple(round(c, 9) for c in mesh.vertex_coordinates(v))
                   for v in mesh.vertices()))


def half_square():
    """A square cut in half -- the smallest layout with an INTERIOR edge."""
    e = CoarseLayoutEditor(build())
    assert e.insert_curve([[0.0, 5.0, 0.0], [10.0, 5.0, 0.0]], extend=True)
    return e


# ----------------------------------------------------------------------
print("\n=== 1. move a boundary corner ===")
ed = CoarseLayoutEditor(build())
print("  layout: {} patches, {} corners".format(
    ed.mesh.number_of_faces(), ed.mesh.number_of_vertices()))

boundary = [v for v in ed.mesh.vertices() if ed.is_vertex_on_boundary(v)]
check("layout has boundary corners", len(boundary) > 0, len(boundary))
vkey = boundary[0]
before = ed.mesh.vertex_coordinates(vkey)
ed.move_vertex(vkey, [before[0] + 2.0, before[1] + 2.0, 0.0])
after = ed.mesh.vertex_coordinates(vkey)
dist = min(w.project(after)[0] for w in ed._walls)
check("moved corner lands ON a wall", dist < 1e-9, "d={:.3e}".format(dist))
check("it actually moved", sum((a - b) ** 2 for a, b in zip(after, before)) > 1e-9)
ed.reset()
check("reset restores the corner", ed.mesh.vertex_coordinates(vkey) == before)

# ----------------------------------------------------------------------
print("\n=== 2. a cut from wall to wall ===")
ed = CoarseLayoutEditor(build())
n_in = ed.mesh.number_of_faces()
ok = ed.insert_curve([[0.0, 5.0, 0.0], [10.0, 5.0, 0.0]], extend=True)
check("wall-to-wall cut accepted", ok, ed.last_reason)
check("patches went up", ed.mesh.number_of_faces() > n_in,
      "{} -> {}".format(n_in, ed.mesh.number_of_faces()))
check("still all-quad", not all_quad(ed.mesh), all_quad(ed.mesh))

# ----------------------------------------------------------------------
print("\n=== 3/4. a cut that stops on an INTERIOR edge ===")
# Aimed at a real interior edge. A one-patch square has none, and a test using
# one measures the wrong refusal: not on the layout, rather than stops short.
ed = half_square()
interior = [(u, v) for u, v in ed.mesh.edges()
            if not ed.mesh.is_edge_on_boundary((u, v))]
check("half-cut layout has an interior edge", len(interior) > 0, len(interior))
u, v = interior[0]
pu, pv = ed.mesh.vertex_coordinates(u), ed.mesh.vertex_coordinates(v)
mid = [(pu[0] + pv[0]) / 2.0, (pu[1] + pv[1]) / 2.0, 0.0]
start = [mid[0], 0.0, 0.0]          # bottom wall, straight up to that edge

digest = mesh_digest(ed.mesh)
ok = ed.insert_curve([start, mid], extend=None)
check("short cut refused with extend=None", not ok)
check("refusal names the five-sided patch", "five sides" in ed.last_reason,
      ed.last_reason[:70])
check("mesh untouched by the refusal", mesh_digest(ed.mesh) == digest)

asked = []
ed = half_square()
n_in = ed.mesh.number_of_faces()
ok = ed.insert_curve([start, mid], extend=lambda n: asked.append(n) or True)
check("extended cut accepted", ok, ed.last_reason)
check("the extend callable WAS asked", asked == [1], asked)
check("patches went up", ed.mesh.number_of_faces() > n_in,
      "{} -> {}".format(n_in, ed.mesh.number_of_faces()))
check("still all-quad", not all_quad(ed.mesh), all_quad(ed.mesh))

# ----------------------------------------------------------------------
print("\n=== 5. a cut that starts nowhere -> refused ===")
ed = CoarseLayoutEditor(build())
digest = mesh_digest(ed.mesh)
ok = ed.insert_curve([[-5.0, 5.0, 0.0], [-5.0, 6.0, 0.0]], extend=True)
check("off-layout cut refused", not ok, ed.last_reason[:60])
check("mesh untouched", mesh_digest(ed.mesh) == digest)

# ----------------------------------------------------------------------
print("\n=== 6. commit a cut, then densify ===")
ed = CoarseLayoutEditor(build())
ed.insert_curve([[0.0, 5.0, 0.0], [10.0, 5.0, 0.0]], extend=True)
faces_before = ed.mesh.number_of_faces()
ok, notes = ed.commit()
check("commit accepted", ok, notes.get("error"))
check("faces_in matches what went in", notes.get("faces_in") == faces_before,
      "{} vs {}".format(notes.get("faces_in"), faces_before))
check("faces_out reported", notes.get("faces_out") is not None, notes.get("faces_out"))
check("lost_curves reported", "lost_curves" in notes, notes.get("lost_curves"))

ed.mesh.collect_strips()      # edit_coarse returns a mesh with none collected
dense = ed.decomposition.quad_mesh(coarse=ed.mesh, target_length=1.0)
check("committed layout densifies",
      dense is not None and dense.number_of_faces() > 0,
      dense.number_of_faces() if dense else None)

# ----------------------------------------------------------------------
print("\n=== 7. delete a strip ===")
ed = CoarseLayoutEditor(build([circle(5.0, 5.0, 1.5)]))
ed.mesh.collect_strips()
print("  layout: {} patches, {} strips, {} poles".format(
    ed.mesh.number_of_faces(), len(list(ed.mesh.strips())), len(ed.mesh.poles())))

plans = [(edge, ed.plan_strip_deletion(edge)) for edge in ed.mesh.edges()]
check("every edge resolves to a strip",
      all(p[1]["skey"] is not None for p in plans), len(plans))

refused = [p for p in plans if not p[1]["ok"]]
allowed = [p for p in plans if p[1]["ok"]]
print("  {} deletable, {} refused".format(len(allowed), len(refused)))
check("the plan separates deletable from not",
      len(allowed) > 0 and len(refused) > 0,
      "{} / {}".format(len(allowed), len(refused)))
check("every refusal carries a reason", all(p[1]["reason"] for p in refused),
      refused[0][1]["reason"][:60] if refused else "")

edge, info = refused[0]
digest = mesh_digest(ed.mesh)
ok = ed.delete_strip(edge)
check("delete refuses what the plan refused", not ok, ed.last_reason[:70])
check("mesh untouched by the refusal", mesh_digest(ed.mesh) == digest)
check("plan and delete give the SAME reason", ed.last_reason == info["reason"])

edge, info = allowed[0]
n_in = ed.mesh.number_of_faces()
ok = ed.delete_strip(edge)
check("allowed strip deleted", ok, ed.last_reason)
check("patches went down", ed.mesh.number_of_faces() < n_in,
      "{} -> {}".format(n_in, ed.mesh.number_of_faces()))
check("still all-quad after delete", not all_quad(ed.mesh), all_quad(ed.mesh))
check("result is manifold", ed.mesh.is_manifold())
check("deletion reported", ed.last_deletion.get("skey") is not None, ed.last_deletion)

ok2, notes2 = ed.commit()
check("commit after delete", ok2, notes2.get("error"))

# ----------------------------------------------------------------------
print("\n=== 8. reset undoes a delete (topology, not just coordinates) ===")
ed = CoarseLayoutEditor(build([circle(5.0, 5.0, 1.5)]))
ed.mesh.collect_strips()
n_in = ed.mesh.number_of_faces()
for edge in ed.mesh.edges():
    if ed.plan_strip_deletion(edge)["ok"] and ed.delete_strip(edge):
        break
check("something was deleted", ed.mesh.number_of_faces() != n_in,
      "{} -> {}".format(n_in, ed.mesh.number_of_faces()))
ed.reset()
check("reset restores the face count", ed.mesh.number_of_faces() == n_in,
      ed.mesh.number_of_faces())

# ----------------------------------------------------------------------
print("\n=== 9. the module imports no Rhino ===")
rhino_modules = [name for name in sys.modules
                 if name.split(".")[0] in ("rhinoscriptsyntax", "Rhino",
                                           "scriptcontext", "System")]
check("no Rhino module was imported", not rhino_modules, rhino_modules)

print("\n" + "=" * 60)
if FAILURES:
    print("FAILED ({}): {}".format(len(FAILURES), ", ".join(FAILURES)))
    sys.exit(1)
print("ALL CHECKS PASSED")
