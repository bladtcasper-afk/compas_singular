"""Headless exercise of DenseMeshEditor -- no Rhino, no prompts.

The whole point of the split: every operation CMD_edit_quad_mesh offers has to
be reachable, and checkable, from a plain script.

Run:
    C:/Users/Casper/anaconda3/envs/singular312/python.exe test_dense_editor.py
"""
import os
import sys

REPO = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular"
for path in (os.path.join(REPO, "src"),):
    if path not in sys.path:
        sys.path.insert(0, path)

from compas_singular.datastructures import QuadMesh                  # noqa: E402
from compas_singular.editing.dense_mesh import DenseMeshEditor       # noqa: E402


FAILURES = []


def check(name, condition, detail=""):
    print("  {} {}{}".format("ok  " if condition else "FAIL", name,
                             (" -- " + str(detail)) if detail else ""))
    if not condition:
        FAILURES.append(name)


def all_quad(mesh):
    """Every face 4-sided. The rule both line operations are checked against."""
    return [(fkey, len(mesh.face_vertices(fkey))) for fkey in mesh.faces()
            if len(mesh.face_vertices(fkey)) != 4]


def grid(nx=4, ny=3, size=1.0):
    """An nx by ny quad grid. Vertex (i, j) at index j * (nx + 1) + i."""
    vertices = [[i * size, j * size, 0.0]
                for j in range(ny + 1) for i in range(nx + 1)]

    def index(i, j):
        return j * (nx + 1) + i

    faces = [[index(i, j), index(i + 1, j), index(i + 1, j + 1), index(i, j + 1)]
             for j in range(ny) for i in range(nx)]
    return QuadMesh.from_vertices_and_faces(vertices, faces)


def square_walls(nx=4, ny=3, size=1.0):
    """The outline of the grid, as one closed loop."""
    w, h = nx * size, ny * size
    return [[[0.0, 0.0, 0.0], [w, 0.0, 0.0], [w, h, 0.0], [0.0, h, 0.0],
             [0.0, 0.0, 0.0]]]


def interior_edge(mesh):
    """An edge with two faces -- something a strip actually runs through."""
    for u, v in mesh.edges():
        if not mesh.is_edge_on_boundary(u, v):
            return (u, v)
    return None


def boundary_edge(mesh):
    for u, v in mesh.edges():
        if mesh.is_edge_on_boundary(u, v):
            return (u, v)
    return None


# ======================================================================
# construction
# ======================================================================

print("\nconstruction")

mesh = grid()
editor = DenseMeshEditor(mesh, walls=square_walls())
check("takes a QuadMesh", editor.mesh is mesh)
check("no walls is legal", DenseMeshEditor(grid()).walls == [])
check("walls parsed", len(editor.walls) == 1, "{} loop(s)".format(len(editor.walls)))
check("starts clean", editor.last_reason == "")
check("grid is all quads", not all_quad(editor.mesh))
check("strips available", editor.strips_available())

try:
    from compas.datastructures import Mesh
    DenseMeshEditor(Mesh.from_vertices_and_faces(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], [[0, 1, 2, 3]]))
    check("plain Mesh refused", False, "no TypeError raised")
except TypeError as e:
    check("plain Mesh refused", "QuadMesh" in str(e), str(e)[:60])


# ======================================================================
# move_vertex
# ======================================================================

print("\nmove_vertex")

editor = DenseMeshEditor(grid(), walls=square_walls())

interior = [v for v in editor.mesh.vertices()
            if not editor.mesh.is_vertex_on_boundary(v)]
vkey = interior[0]
before = editor.mesh.vertex_coordinates(vkey)
check("moved an interior vertex",
      editor.move_vertex(vkey, [before[0] + 0.2, before[1] + 0.1, 0.0]))
after = editor.mesh.vertex_coordinates(vkey)
check("interior vertex went exactly where it was told",
      abs(after[0] - before[0] - 0.2) < 1e-9 and abs(after[1] - before[1] - 0.1) < 1e-9,
      after)

bnd = [v for v in editor.mesh.vertices() if editor.mesh.is_vertex_on_boundary(v)]
bvkey = next(v for v in bnd
             if abs(editor.mesh.vertex_coordinates(v)[1]) < 1e-9
             and 0.5 < editor.mesh.vertex_coordinates(v)[0] < 3.5)
start = editor.mesh.vertex_coordinates(bvkey)
editor.move_vertex(bvkey, [start[0], start[1] - 0.5, 0.0])
held = editor.mesh.vertex_coordinates(bvkey)
check("boundary vertex held on the wall", abs(held[1]) < 1e-6,
      "y = {:.4f}, was pushed to -0.5".format(held[1]))

editor.move_vertex(bvkey, [start[0], start[1] - 0.5, 0.0], project=False)
loose = editor.mesh.vertex_coordinates(bvkey)
check("project=False lets it leave the wall", abs(loose[1] + 0.5) < 1e-9, loose)

check("unknown vertex refused", editor.move_vertex(9999, [0, 0, 0]) is False)
check("  and says why", "not in the mesh" in editor.last_reason,
      editor.last_reason)

check("move_vertex works with a 2-item point",
      editor.move_vertex(vkey, [1.0, 1.0]))

no_walls = DenseMeshEditor(grid())
bv = next(v for v in no_walls.mesh.vertices()
          if no_walls.mesh.is_vertex_on_boundary(v))
p = no_walls.mesh.vertex_coordinates(bv)
no_walls.move_vertex(bv, [p[0], p[1] - 0.5, 0.0])
check("no walls: boundary vertex left where it was put",
      abs(no_walls.mesh.vertex_coordinates(bv)[1] - (p[1] - 0.5)) < 1e-9)


# ======================================================================
# add_line
# ======================================================================

print("\nadd_line")

editor = DenseMeshEditor(grid(), walls=square_walls())
faces_before = editor.mesh.number_of_faces()
edge = interior_edge(editor.mesh)
pkey, polyedge = editor.polyedge_through(edge)
check("polyedge found through an interior edge", polyedge is not None,
      "{} vertices".format(len(polyedge) if polyedge else 0))

check("added a line", editor.add_line(edge), editor.last_reason)
check("faces grew", editor.mesh.number_of_faces() > faces_before,
      "{} -> {}".format(faces_before, editor.mesh.number_of_faces()))
check("still all quads", not all_quad(editor.mesh), all_quad(editor.mesh)[:3])
check("addition reported", editor.last_addition.get("strip") is not None,
      editor.last_addition)

# The zero-width trap: the pair add_strip created must have separated.
pair = editor.last_addition["new_vertices"]
if len(pair) >= 2:
    a = editor.mesh.vertex_coordinates(pair[0])
    b = editor.mesh.vertex_coordinates(pair[1])
    spread = max(abs(a[0] - b[0]), abs(a[1] - b[1]))
    check("relax opened the new strip", spread > 1e-6,
          "widest new pair separation {:.6f}".format(spread))

unrelaxed = DenseMeshEditor(grid(), walls=square_walls())
e2 = interior_edge(unrelaxed.mesh)
unrelaxed.add_line(e2, relax=False)
check("relax=False records that it did not relax",
      unrelaxed.last_addition.get("relaxed") is False)

check("unknown edge refused", editor.add_line((9998, 9999)) is False)
check("  and says why", "no polyedge" in editor.last_reason, editor.last_reason)

before_mesh = editor.mesh
editor.add_line((9998, 9999))
check("a refusal leaves the mesh untouched", editor.mesh is before_mesh)


# ======================================================================
# plan_strip_deletion / remove_line
# ======================================================================

print("\nplan_strip_deletion")

editor = DenseMeshEditor(grid(), walls=square_walls())
edge = interior_edge(editor.mesh)
plan = editor.plan_strip_deletion(edge)
check("plan mutates nothing",
      editor.mesh.number_of_faces() == 12, editor.mesh.number_of_faces())
check("plan names the strip", plan["strip"] is not None, plan["strip"])
check("plan counts the faces", plan["faces"] > 0, plan["faces"])
check("plan reports collateral as a list", isinstance(plan["collateral"], list))
check("plan reports boundaries_lost as a COUNT",
      isinstance(plan["boundaries_lost"], int), plan["boundaries_lost"])
check("  and which boundaries, separately",
      isinstance(plan["boundaries_lost_vertices"], list)
      and len(plan["boundaries_lost_vertices"]) == plan["boundaries_lost"],
      plan["boundaries_lost_vertices"])
check("plan predicts the result size",
      plan["ok"] is False or plan["faces_after"] is not None, plan)

bad_plan = editor.plan_strip_deletion((9998, 9999))
check("plan on an unknown edge is not ok", bad_plan["ok"] is False)
check("  and says why", "no strip" in bad_plan["reason"], bad_plan["reason"])

print("\nremove_line")

editor = DenseMeshEditor(grid(), walls=square_walls())
edge = interior_edge(editor.mesh)
plan = editor.plan_strip_deletion(edge)
faces_before = editor.mesh.number_of_faces()

if plan["ok"]:
    check("removed a line", editor.remove_line(edge), editor.last_reason)
    check("faces shrank", editor.mesh.number_of_faces() < faces_before,
          "{} -> {}".format(faces_before, editor.mesh.number_of_faces()))
    check("still all quads", not all_quad(editor.mesh), all_quad(editor.mesh)[:3])
    check("deletion reported", editor.last_deletion.get("strip") is not None,
          editor.last_deletion)
    check("plan predicted the face count",
          plan["faces_after"] == editor.mesh.number_of_faces(),
          "planned {}, got {}".format(plan["faces_after"],
                                      editor.mesh.number_of_faces()))
else:
    check("plan refused, so remove_line refuses too",
          editor.remove_line(edge) is False, plan["reason"])
    check("  with the same reason", editor.last_reason == plan["reason"])

check("unknown edge refused", editor.remove_line((9998, 9999)) is False)
check("  and says why", "no strip" in editor.last_reason, editor.last_reason)


# ======================================================================
# the sweep: every strip of a grid, planned and performed
# ======================================================================

print("\nsweep -- every strip of a 4x3 grid")

base = grid()
probe = DenseMeshEditor(base.copy())
edges = list(probe.mesh.edges())
planned_ok = 0
performed_ok = 0
disagreements = []

for e in edges:
    ed = DenseMeshEditor(base.copy())
    p = ed.plan_strip_deletion(e)
    if p["strip"] is None:
        continue
    got = ed.remove_line(e)
    if p["ok"]:
        planned_ok += 1
    if got:
        performed_ok += 1
    if bool(got) != bool(p["ok"]):
        disagreements.append((e, p["ok"], got, p["reason"], ed.last_reason))
    if got and all_quad(ed.mesh):
        disagreements.append((e, "non-quad result", all_quad(ed.mesh)))

check("plan and perform agree on every edge", not disagreements,
      disagreements[:2])
check("some strips are deletable", performed_ok > 0,
      "{} of {} edges".format(performed_ok, len(edges)))
print("    planned ok {}, performed ok {}, edges {}".format(
    planned_ok, performed_ok, len(edges)))


# ======================================================================
# poles disable the line operations
# ======================================================================

print("\npoles disable the line operations")

tri = QuadMesh.from_vertices_and_faces(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0],
     [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0]],
    [[0, 1, 4, 3], [1, 2, 4]])
poled = DenseMeshEditor(tri)
check("strips_available is False with a triangle", poled.strips_available() is False)
check("add_line refused", poled.add_line((0, 1)) is False)
check("  and names the pole rule", "not defined" in poled.last_reason,
      poled.last_reason[:70])
check("remove_line refused", poled.remove_line((0, 1)) is False)
check("plan is not ok either",
      poled.plan_strip_deletion((0, 1))["ok"] is False)
tri_v = poled.mesh.vertex_coordinates(0)
check("move_vertex still works on a poled mesh",
      poled.move_vertex(0, [tri_v[0] - 0.1, tri_v[1] - 0.1, 0.0]))


# ======================================================================
# reset / snapshot
# ======================================================================

print("\nreset / snapshot")

editor = DenseMeshEditor(grid(), walls=square_walls())
start_faces = editor.mesh.number_of_faces()
editor.add_line(interior_edge(editor.mesh))
grown = editor.mesh.number_of_faces()
editor.reset()
check("reset restores the topology", editor.mesh.number_of_faces() == start_faces,
      "{} -> {} -> {}".format(start_faces, grown, editor.mesh.number_of_faces()))

editor.add_line(interior_edge(editor.mesh))
editor.snapshot()
kept = editor.mesh.number_of_faces()
editor.add_line(interior_edge(editor.mesh))
editor.reset()
check("snapshot moves the restore point", editor.mesh.number_of_faces() == kept,
      "{} vs {}".format(editor.mesh.number_of_faces(), kept))

editor.reset()
check("reset clears last_reason", editor.last_reason == "")


# ======================================================================

print("\n{} check(s) failed{}".format(
    len(FAILURES), (": " + ", ".join(FAILURES)) if FAILURES else ""))
sys.exit(1 if FAILURES else 0)
