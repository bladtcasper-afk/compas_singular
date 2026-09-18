#! python3

# r: compas

"""**Step 7 -- edit the FINAL mesh by hand.**

    reads   TopologyProblem::QuadMesh (or a sublayer)  the mesh you pick
            TopologyProblem::InputBoundaries::*         the walls, optional
    writes  TopologyProblem::QuadMesh::Edited           on 'save' only
            TopologyProblem::QuadMesh::Unedited         this session's starting mesh
    scratch TopologyProblem::QuadMesh::Edit             pickable points and lines

**This is the end of the workflow, so the mesh no longer has to be all quads.**
Removing an edge leaves a hexagon, drawing a diagonal leaves two triangles, and
that is allowed: what the result is for is the user's call. The operations:

* **move_vertex** -- drag a vertex. One on the mesh boundary is held on the
  domain wall;
* **remove_vertex** -- the vertex goes with every face around it;
* **remove_edge** -- the two faces either side merge into one polygon. A
  boundary edge has only one face, and that face goes;
* **remove_face** -- click inside a face and it goes;
* **draw_edges** -- draw a polyline across the mesh: every edge it crosses is
  split and every face it passes through is split along it. Snap to the
  scratch points to go through existing vertices. A line cannot END inside a
  face -- a face is a closed loop of vertices -- so the ends are trimmed back to
  the last edge they crossed;
* **add_line** / **remove_line** -- the strip grammar: pick an edge, and a strip
  grows beside the whole line through it, or the whole strip through it is
  deleted. Only where the mesh is still quads: a strip that reaches a triangle,
  a pole or a polygon is refused, every other strip is not;
* **relax** -- smooth the interior, holding every boundary;
* **undo** -- one edit back; **reset** -- back to where this session started.

**Every edit happens the moment it is picked**, with no "are you sure". Each is
snapshotted first, and 'undo' is where a wrong pick goes -- the coarse editor
works the same way. Only a strip removal that would collapse a boundary asks
first, because that is the one whose consequence is not on the thing clicked.

**A successful save ENDS the command.** Nothing is written to the document until
then. A FAILED save keeps the loop open, so the edit is not lost with it.

**Saving keeps polygons as polygons.** A Rhino mesh face has four slots, so a
pentagon is stored as five triangles plus an n-gon record; the loader here reads
that record back (``helpers.mesh_from_rhino``). ``mesh_to_compas`` does not: a
command that reads this mesh with it sees the triangles and a centre vertex.

**Where the code lives.** The editing is ``compas_singular.editing.DenseMeshEditor``
and imports no ``rs``; picking, dragging and drawing are
``compas_singular.rhino.mesh_ui``. What stays here is the menu, the prompts and
the bake.

**Re-running CMD_quad_mesh regenerates from the coarse layout** and discards
everything done here. If an edit can be made on the coarse layout, make it there.
"""

# ----------------------------------------------------------------------
# imports
# ----------------------------------------------------------------------

# Development bootstrap -- delete once compas_singular is installed into Rhino's
# Python. MUST run before any compas_singular import: Rhino resets sys.path between
# runs but keeps sys.modules, so put the source on the path and drop a stale copy
# (see CMD_start for why every module, framefield included, has to go).
import sys
SINGULAR_SRC = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular\src"
if SINGULAR_SRC not in sys.path:
    sys.path.insert(0, SINGULAR_SRC)
if not getattr(sys, "compas_singular_keep_modules", False):  # set by headless tests
    for _mod in list(sys.modules):
        if _mod == "compas_singular" or _mod.startswith("compas_singular."):
            del sys.modules[_mod]

import rhinoscriptsyntax as rs
import scriptcontext as sc

from compas_singular.datastructures import QuadMesh
from compas_singular.editing import DenseMeshEditor
from compas_singular.rhino import mesh_ui
from compas_singular.rhino.helpers import bake_mesh
from compas_singular.rhino.helpers import clear_layer
from compas_singular.rhino.helpers import mesh_from_rhino
from compas_singular.rhino.helpers import read_boundary_loops
from compas_singular.rhino.project import get_settings
from compas_singular.rhino.project import ROOT, layer_path


# ----------------------------------------------------------------------
# layers and constants
# ----------------------------------------------------------------------

QUADMESH_LAYER = layer_path("QuadMesh")
EDIT_LAYER = QUADMESH_LAYER + "::Edit"
UNEDITED_LAYER = QUADMESH_LAYER + "::Unedited"
EDITED_LAYER = QUADMESH_LAYER + "::Edited"

#: Above this many edges the first drawing is enough objects to make Rhino
#: noticeably slow, so it is worth asking. Later edits only redraw what changed.
DRAW_WARN_EDGES = 4000


# ----------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------

def load_mesh():
    """The PICKED mesh, n-gons intact, with the object and layer it came from."""
    if not rs.IsLayer(QUADMESH_LAYER):
        raise RuntimeError(
            "No layer '{}' -- run CMD_quad_mesh first.".format(QUADMESH_LAYER))
    quad_layer = rs.LayerName(QUADMESH_LAYER, fullpath=True)

    def on_quad_layer(rhobj, geometry, component_index):
        layer = rs.ObjectLayer(rhobj)
        return layer == quad_layer or rs.IsLayerChildOf(layer, quad_layer)

    guid = rs.GetObject(
        message="Pick the mesh to edit, from 'QuadMesh' or a sublayer",
        filter=rs.filter.mesh, preselect=False, select=False,
        custom_filter=on_quad_layer, subobjects=False)
    if not guid:
        return None, None, None
    return (mesh_from_rhino(rs.coercemesh(guid), cls=QuadMesh),
            guid, rs.ObjectLayer(guid))


def load_walls():
    """The domain walls, for holding a moved boundary vertex on the outline.

    Optional: without the input boundaries in the document the command still
    works, the boundary simply stops being held.
    """
    try:
        spacing = get_settings()["triangulation_spacing"]
        outer_loop, inner_loops = read_boundary_loops(spacing * 0.25)
    except Exception as e:
        print("boundary walls unavailable ({}) -- a moved boundary vertex will "
              "NOT be held on the outline.".format(e))
        return []
    return [loop for loop in [outer_loop] + inner_loops if len(loop) >= 3]


def face_degrees(mesh):
    out = {}
    for fkey in mesh.faces():
        n = len(mesh.face_vertices(fkey))
        out[n] = out.get(n, 0) + 1
    return dict(sorted(out.items()))


# ----------------------------------------------------------------------
# the editor, the scene, undo
# ----------------------------------------------------------------------

def special_vertices():
    """Singularities, and any vertex where faces only touch at a corner."""
    try:
        special = set(editor.mesh.singularities())
    except Exception:
        special = set()
    special.update(editor.non_manifold_vertices())
    return special


def redraw():
    removed, added = scene.sync(editor.mesh, special=special_vertices())
    return removed, added


def edit(mutate):
    """Run one ``editor`` call with undo bookkeeping. ``(ok, notes)``.

    Snapshot first, and drop the snapshot again if the call refused -- a refusal
    never touches ``editor.mesh``, so there would be nothing to undo back to.
    """
    editor.push_undo()
    ok, notes = mutate()
    if not ok:
        editor.discard_last_undo()
    return ok, notes


def warn_non_manifold(notes):
    count = len(notes.get("non_manifold") or [])
    if count:
        print("  {} vertex/vertices now join faces that only touch at a corner "
              "(drawn in the singularity colour). add_line and remove_line are "
              "refused until that is fixed -- 'undo' takes it back.".format(count))


# ----------------------------------------------------------------------
# operations -- each loops on picks until Esc
# ----------------------------------------------------------------------

def move_vertex():
    changed = False
    while True:
        vkey = scene.pick_vertex("Select a vertex to move (Esc to return to the menu)")
        if vkey is None:
            return changed
        start = editor.mesh.vertex_coordinates(vkey)
        on_boundary = editor.is_vertex_on_boundary(vkey)
        held = on_boundary and bool(editor.walls)

        def project(point, start=start, held=held):
            # The drag is in world XY; the vertex keeps its own height, so a
            # mesh that has been given one is not flattened by a move.
            point = [point[0], point[1], start[2]]
            return editor.project_to_wall(point) if held else point

        xyz = mesh_ui.drag_point(
            "New position for vertex {}{}".format(vkey, " (held on the wall)" if held else ""),
            start,
            neighbours=[editor.mesh.vertex_coordinates(n) for n in editor.mesh.vertex_neighbors(vkey)],
            project=project)
        if xyz is None:
            return changed
        ok, _notes = edit(lambda: editor.move_vertex(vkey, xyz, project=False))
        if ok:
            changed = True
            redraw()


def remove_vertex():
    changed = False
    while True:
        vkey = scene.pick_vertex("Select a vertex to remove with its faces (Esc to return)")
        if vkey is None:
            return changed
        ok, notes = edit(lambda: editor.remove_vertex(vkey))
        if not ok:
            mesh_ui.refuse("Remove vertex", "Nothing was removed.\n\n{}.".format(notes["error"]))
            continue
        print("vertex {} removed with {} face(s); mesh now has {} face(s).".format(
            vkey, notes["faces_removed"], notes["faces"]))
        warn_non_manifold(notes)
        changed = True
        redraw()


def remove_edge():
    changed = False
    while True:
        picked = scene.pick_edge("Select an edge to remove (Esc to return to the menu)")
        if picked is None:
            return changed
        edge, _guid = picked
        ok, notes = edit(lambda: editor.remove_edge(edge))
        if not ok:
            mesh_ui.refuse("Remove edge", "Nothing was removed.\n\n{}.".format(notes["error"]))
            continue
        if notes["merged"] is None:
            print("boundary edge {}: its one face was removed; {} face(s) left.".format(
                edge, notes["faces"]))
        else:
            print("edge {} removed: faces of {} and {} sides merged into one of {}.".format(
                edge, notes["degrees"][0], notes["degrees"][1], notes["degree"]))
        warn_non_manifold(notes)
        changed = True
        redraw()


def remove_face():
    changed = False
    while True:
        point = rs.GetPoint("Click inside a face to remove it (Esc to return to the menu)")
        if point is None:
            return changed
        fkey = editor.face_at([point.X, point.Y, point.Z])
        if fkey is None:
            print("That point is not inside a face of the mesh -- click inside one.")
            continue
        ok, notes = edit(lambda: editor.remove_face(fkey))
        if not ok:
            mesh_ui.refuse("Remove face", "Nothing was removed.\n\n{}.".format(notes["error"]))
            continue
        print("face with {} sides removed; {} face(s) left.".format(notes["degree"], notes["faces"]))
        warn_non_manifold(notes)
        changed = True
        redraw()


def draw_edges():
    changed = False
    print("Draw a polyline across the mesh; Enter to finish it. Snap to the "
          "points on '{}' to go through existing vertices.".format(EDIT_LAYER))
    while True:
        points = rs.GetPolyline(
            message1="First point of the line (Esc to return to the menu)",
            message2="Next point",
            message3="Next point, Enter to cut it in")
        if not points:
            return changed
        stroke = [[p.X, p.Y, p.Z] for p in points]
        ok, notes = edit(lambda: editor.draw_edges(stroke))
        if not ok:
            mesh_ui.refuse("Draw edges", "Nothing was added.\n\n{}.".format(notes["error"]))
            continue
        print("line cut in: {} edge(s) split, {} face(s) split, {} new vertex/vertices; "
              "faces now {}.".format(notes["split_edges"], notes["split_faces"],
                                     notes["new_vertices"], face_degrees(editor.mesh)))
        left_out = []
        if notes["trimmed"]:
            left_out.append("{} loose end(s) inside a face or off the mesh".format(notes["trimmed"]))
        if notes["outside"]:
            left_out.append("{} stretch(es) outside the mesh".format(notes["outside"]))
        if notes["along_existing"]:
            left_out.append("{} stretch(es) along existing edges".format(notes["along_existing"]))
        if left_out:
            print("  left out: " + "; ".join(left_out))
        changed = True
        redraw()


def add_line():
    changed = False
    while True:
        picked = scene.pick_edge("Pick an edge; the whole line through it gains a strip (Esc to return)")
        if picked is None:
            return changed
        edge, _guid = picked
        # ``relax=True`` explicitly: the grammar creates the new row with ZERO
        # width, and the editor no longer opens it for you. Without this the
        # strip is invisible here and the mesh bakes with coincident corners.
        ok, notes = edit(lambda: editor.add_line(edge, relax=True))
        if not ok:
            mesh_ui.refuse("Add line", "No strip was added.\n\n{}".format(notes["error"]))
            continue
        print("strip added along a {} line of {} vertices: faces {} -> {}.".format(
            "closed" if notes["closed"] else "wall-to-wall", notes["polyedge_vertices"],
            notes["faces_before"], notes["faces_after"]))
        changed = True
        redraw()


def remove_line():
    changed = False
    while True:
        picked = scene.pick_edge("Pick an edge; the whole strip through it is removed (Esc to return)")
        if picked is None:
            return changed
        edge, _guid = picked
        plan = editor.plan_strip_deletion(edge)
        if not plan["ok"]:
            mesh_ui.refuse("Remove line", "Nothing was removed.\n\n{}".format(plan["reason"]))
            continue

        preserve = False
        if plan["boundaries_lost"]:
            # The one strip removal that asks: its consequence -- a hole or a
            # boundary closing up -- is not on the strip that was clicked.
            answer = mesh_ui.ask(
                "This strip COLLAPSES {} boundary/boundaries. Remove it?".format(plan["boundaries_lost"]),
                ["Yes", "PreserveBoundaries", "No"], "No")
            if not answer or answer.startswith("n"):
                print("Nothing removed.")
                continue
            preserve = answer.startswith("p")

        ok, notes = edit(lambda: editor.remove_line(edge, preserve_boundaries=preserve))
        if not ok:
            mesh_ui.refuse("Remove line", "Nothing was removed.\n\n{}".format(notes["error"]))
            continue
        print("strip removed: {} face(s){}; faces {} -> {}. The two sides welded "
              "together, so the vertices beside it moved.".format(
                  notes["faces"],
                  ", plus {} collateral strip(s)".format(len(notes["collateral"])) if notes["collateral"] else "",
                  notes["faces_before"], notes["faces_after"]))
        changed = True
        redraw()


def relax():
    value = mesh_ui.ask_integer("Smoothing passes", editor.relax_iterations, 1)
    if not value:
        return False
    editor.push_undo()
    editor.relax(iterations=value)
    redraw()
    print("relaxed: {} pass(es), every boundary held.".format(value))
    return True


def undo():
    ok, notes = editor.undo()
    if not ok:
        print("Nothing to undo.")
        return False
    redraw()
    print("Undone -- {} face(s). {} more undo(s) available.".format(notes["faces"], notes["remaining"]))
    return True


def save():
    """Bake the mesh to ``QuadMesh::Edited``. Its guids, or ``[]`` on failure.

    Whether this worked decides whether the command ends, so it is reported
    rather than assumed. ``bake_mesh`` clears the layer BEFORE adding, so a
    failed save may already have emptied it; this session's starting mesh is on
    ``QuadMesh::Unedited``, which is what the message points at.
    """
    layer = mesh_ui.ensure_layer(EDITED_LAYER, (0, 120, 200))
    try:
        bake_mesh(editor.mesh, layer)
    except Exception as e:
        mesh_ui.refuse(
            "Save mesh",
            "The mesh was NOT saved.\n\n{}\n\nThe edit is still in memory, so try "
            "again. If '{}' is now empty, this session's starting mesh is on "
            "'{}'.".format(e, layer, UNEDITED_LAYER))
        return []
    guids = rs.ObjectsByLayer(layer) or []
    if not guids:
        mesh_ui.refuse(
            "Save mesh",
            "The mesh was NOT saved -- nothing landed on '{}'.\n\nThe edit is still "
            "in memory, so try again. This session's starting mesh is on "
            "'{}'.".format(layer, UNEDITED_LAYER))
        return []
    print("Saved to '{}': {} faces {}, {} vertices.".format(
        layer, editor.mesh.number_of_faces(), face_degrees(editor.mesh),
        editor.mesh.number_of_vertices()))
    return guids


# ----------------------------------------------------------------------
# run
# ----------------------------------------------------------------------

mesh, source_guid, source_layer = load_mesh()
if mesh is None:
    print("Cancelled -- nothing was drawn or changed.")
    raise SystemExit

print("editing the mesh on '{}': {} faces {}, {} vertices, {} edges".format(
    source_layer, mesh.number_of_faces(), face_degrees(mesh),
    mesh.number_of_vertices(), mesh.number_of_edges()))

if mesh.number_of_edges() > DRAW_WARN_EDGES:
    answer = mesh_ui.ask(
        "This mesh has {} edges. Drawing them all as pickable objects will make "
        "Rhino slow. Continue?".format(mesh.number_of_edges()), ["Yes", "No"], "Yes")
    if not answer.startswith("y"):
        print("Cancelled -- nothing was drawn or changed.")
        raise SystemExit

walls = load_walls()
if not walls:
    print("  no domain walls loaded -- a moved boundary vertex will not be held "
          "on the outline.")

editor = DenseMeshEditor(mesh, walls=walls)
scene = mesh_ui.PickableMesh(EDIT_LAYER)

# This session's starting point, on its own layer, hidden. Refreshed every run:
# a backup kept from before a re-run of CMD_quad_mesh would be a mesh of a
# different layout, and stale is worse than absent.
mesh_ui.ensure_layer(UNEDITED_LAYER, (170, 170, 170))
bake_mesh(mesh, UNEDITED_LAYER)
rs.LayerVisible(UNEDITED_LAYER, False)
print("this session's starting mesh kept on '{}'".format(UNEDITED_LAYER))

# The source object is in the way while picking. Hidden rather than deleted, so
# a cancelled session leaves the document as it was found.
rs.HideObjects([source_guid])

OPERATIONS = ["move_vertex", "remove_vertex", "remove_edge", "remove_face",
              "draw_edges", "add_line", "remove_line", "relax", "undo", "reset",
              "save", "exit"]

try:
    redraw()
    lock_state = scene.unlock()
    try:
        while True:
            operation = mesh_ui.ask("next", OPERATIONS, "exit")
            if operation == "move_vertex":
                move_vertex()
            elif operation == "remove_vertex":
                remove_vertex()
            elif operation == "remove_edge":
                remove_edge()
            elif operation == "remove_face":
                remove_face()
            elif operation == "draw_edges":
                draw_edges()
            elif operation == "add_line":
                add_line()
            elif operation == "remove_line":
                remove_line()
            elif operation == "relax":
                relax()
            elif operation == "undo":
                undo()
            elif operation == "reset":
                editor.reset()
                redraw()
                print("Back to this session's starting mesh.")
            elif operation == "save":
                if save():
                    print("  CMD_quad_mesh regenerates from the coarse layout and "
                          "would discard this. Commands reading the mesh with "
                          "mesh_to_compas see polygons as triangle fans.")
                    break
                # a failed save keeps the loop, and the edit, open
            else:
                if editor.edited:
                    answer = mesh_ui.ask("Unsaved edits. Save before leaving?", ["Yes", "No"], "Yes")
                    if answer.startswith("y"):
                        if not save():
                            print("Leaving WITHOUT saving -- the save failed above.")
                    else:
                        print("Left unsaved -- '{}' is unchanged.".format(source_layer))
                break
    finally:
        scene.relock(lock_state)
finally:
    scene.clear()
    if rs.IsLayer(EDIT_LAYER):
        clear_layer(EDIT_LAYER)
    if rs.IsObject(source_guid):
        rs.ShowObjects([source_guid])
    sc.doc.Views.Redraw()
