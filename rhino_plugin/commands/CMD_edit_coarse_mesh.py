#! python3

# r: compas

"""**Step 4 -- edit the COARSE layout by hand, then commit it.**

    reads   TopologyProblem::Skeleton::Mesh          the layout
            TopologyProblem::Skeleton::Poles         its poles
            TopologyProblem::InputBoundaries::*      the domain
    writes  TopologyProblem::Skeleton::Mesh          on 'commit' only
            TopologyProblem::Skeleton::Poles
            TopologyProblem::Skeleton::Polylines     the re-traced separatrices
            TopologyProblem::Skeleton::EdgeCurves    the layout WITH its curvature
    scratch TopologyProblem::Skeleton::TempEdit      pickable corners and edges

Five operations:

* **move_vertices** -- drag a corner. A corner on the layout boundary is held
  on the domain wall, so nudging the edge of the layout does not eat the
  outline;
* **add_polyedge** -- draw a Polyline or an Arc across the layout and cut a new
  line of coarse edges into it: every edge crossed is split, every patch passed
  through is split in two. The drawn SHAPE is kept, so an arc densifies as an
  arc. The run always reaches a wall at both ends -- a run that stops inside is
  carried on along its strip rather than left to make a pole nobody drew;
* **add_strip** -- pick a run of existing corners and UNZIP it into a strip:
  each corner becomes two and a new band opens between them, redividing the
  space of the two patches either side into three. Robin's grammar rule, and it
  keeps the strip data valid, so the densities set in step 5 survive it;
* **remove_strip** -- pick a strip by its RIBBON and the whole band is deleted.
  The two sides weld together;
* **commit** -- write the edited layout back into the mesh this command read,
  in place. What that costs depends on what was edited: after moves alone the
  boundary corners are re-snapped to the walls and the strip data (and so the
  densities set in step 5) survives; after a cut the layout is welded and
  repaired, which RENUMBERS every corner and drops the strip data. Either way a
  layout that will not densify is REFUSED rather than quietly replaced.

**Nothing is written to the document until commit.** ``reset`` goes back to the
layout this session started with, and leaving without committing leaves
``Skeleton::Mesh`` exactly as it was found.

**Where the code lives, and why.** Everything in this file is Rhino: picking,
drawing, prompting, baking. The editing itself -- the cut planner, the all-quad
gate, the strip deletion, the curve map, the commit -- is
``compas_singular.editing.CoarseEditor`` and imports no ``rs`` at all, so
the same four operations can be scripted and tested with Rhino closed. The
generic half of the interaction (draw a mesh as pickable objects, pick a vertex
or an edge, drag with a live preview) is ``compas_singular.rhino.mesh_ui``,
shared with CMD_edit_quad_mesh. What stays here is what only this command needs:
drawing a curve ACROSS the layout, and the prompts that phrase this step.

**Why a line is never one edge.** A quad layout cannot gain or lose a single
edge: removing one merges two patches into a hexagon, adding one leaves a
five-sided patch. Every operation that keeps the layout all-quad runs the full
width of it. That is the one rule behind every refusal here, and the refusal
message says which half of it was hit.
"""

# ----------------------------------------------------------------------
# imports
# ----------------------------------------------------------------------

# MUST come before any compas_singular import: it fixes sys.path and purges a
# stale copy from sys.modules, and doing that after a class is bound would
# leave an existing mesh failing isinstance against the rebuilt class.
from CMD_start import import_compas_singular
import_compas_singular()


import rhinoscriptsyntax as rs
import scriptcontext as sc
import Rhino
from System.Drawing import Color


from compas_singular.editing import CoarseEditor
from compas_singular.framefield.field import CrossField
from compas_singular.rhino import mesh_ui
from compas_singular.rhino.coarse_curves import coarse_edges_to_curves
from compas_singular.rhino.helpers.helpers import bake_edge_curves
from compas_singular.rhino.helpers.helpers import bake_mesh
from compas_singular.rhino.helpers.helpers import bake_points
from compas_singular.rhino.helpers.helpers import bake_polylines
from compas_singular.rhino.helpers.helpers import clear_layer
from compas_singular.rhino.helpers.helpers import read_boundaries
from compas_singular.rhino.helpers.helpers import read_boundary_loops
from compas_singular.rhino.helpers.helpers import read_coarse
from compas_singular.rhino.helpers.helpers import read_polylines
from CMD_start import get_settings
from CMD_start import cache_path, COARSE_CACHE, FIELD_CACHE


# ----------------------------------------------------------------------
# settings and layers
# ----------------------------------------------------------------------

settings = get_settings()
SPACING = settings["triangulation_spacing"]

#: The walls the layout's boundary edges densify ALONG. Finer than the
#: background, because these points ARE the boundary from here on.
WALL_SAMPLING = SPACING * 0.25

#: Layers the helpers read and write by their short name. Left as they are:
#: ``read_coarse`` and ``bake_mesh`` both use these, and changing one without
#: the other silently splits the layout across two layers.
MESH_LAYER = "Mesh"
POLES_LAYER = "Poles"

# ``rs.AddLayer`` returns the FULL '::' path, which is what every other rs call
# needs -- a nested layer's short name is not reliably resolvable on its own,
# and handing the short name to the scratch layer used to risk creating a
# SECOND, root-level 'TempEdit' the command then never cleaned up.
EDIT_LAYER = rs.AddLayer(name="TempEdit", parent="Skeleton", color=(0, 120, 200))
POLYLINE_LAYER = rs.AddLayer(name="Polylines", parent="Skeleton")
EDGE_CURVE_LAYER = rs.AddLayer(name="EdgeCurves", parent="Skeleton", color=(0, 120, 200))


# ----------------------------------------------------------------------
# load
# ----------------------------------------------------------------------

# Read ONCE. This samples every boundary curve in the document, and the command
# used to do it twice with identical arguments and throw the first result away.
outer, inners, guides, point_features = read_boundaries(spacing=SPACING)
coarse, poles = read_coarse()
print("layout: {} patch(es), {} corner(s), {} pole(s)".format(
    coarse.number_of_faces(), coarse.number_of_vertices(), len(poles)))

# **This step no longer solves a field.** It used to call ``get_decomposition()``
# for one reason only -- ``edit_coarse`` lived on the decomposition -- and paid
# for a full solve to do it. The editor is anchored on the LAYOUT now, so all it
# needs is what step 3 already left behind:
#
#   * the traced separatrices, baked on ``Skeleton::Polylines``. They give a cut
#     across a curved edge the shape of the two halves it makes, and they are
#     republished with the user's own curves through ``editor.all_polylines``;
#   * the cross field, from the side-car ``field.json`` -- the same file step 6
#     reads. It is CARRIED, not consulted: nothing in an edit depends on it, and
#     ``default=None`` because a skeleton-route document legitimately has none.
field = CrossField.load_from_json(cache_path(FIELD_CACHE, create=False),
                                  default=None)
polylines = read_polylines(POLYLINE_LAYER)
print("separatrices: {} read from the document".format(len(polylines)))

#: What ``commit()`` puts here, and the only thing the write-back below tests.
#: It used to read ``decomposition.edit_notes``, which was a side channel through
#: an object this command no longer has.
COMMITTED = None
COMMIT_NOTES = {}

editor = CoarseEditor(coarse, field=field, loops=[outer] + inners,
                      polylines=polylines)
scene = mesh_ui.PickableMesh(EDIT_LAYER)


def redraw():
    """The layout as pickable corners and edges, poles coloured.

    An edge carrying a drawn shape is drawn as that shape -- drawing an added
    arc as the chord it is stored as would make a cut that worked look like one
    that snapped straight.
    """
    scene.draw(editor.mesh,
               edge_shape=editor.edge_shape,
               special=editor.mesh.poles() if hasattr(editor.mesh, "poles") else ())


# ----------------------------------------------------------------------
# move corners
# ----------------------------------------------------------------------

def move_vertices():
    """Drag corners until Esc. ``True`` if any moved."""
    changed = False
    while True:
        vkey = scene.pick_vertex("Select a coarse corner")
        if vkey is None:
            return changed
        on_boundary = editor.is_vertex_on_boundary(vkey)
        neighbours = [editor.mesh.vertex_coordinates(nbr)
                      for nbr in editor.mesh.vertex_neighbors(vkey)]
        xyz = mesh_ui.drag_point(
            "New position for corner {}{}".format(
                vkey, " (held on the wall)" if on_boundary else ""),
            editor.mesh.vertex_coordinates(vkey),
            neighbours=neighbours,
            # The preview shows the SNAPPED position, not the cursor, so what
            # is drawn is what the commit will get.
            project=lambda point: editor.project(point, on_boundary))
        if xyz is None:
            return changed
        moved, _notes = editor.move_vertex(vkey, xyz)
        if moved:
            changed = True
            redraw()


# ----------------------------------------------------------------------
# draw a curve across the layout -- the part only this command needs
# ----------------------------------------------------------------------

def choose_mode():
    """Polyline / Arc, or ``None`` if the user cancelled (Esc).

    **No Line.** A straight line is a two-point polyline, and having it as its
    own mode cost a decision at every cut for nothing: the polyline mode draws
    the same thing with the same two picks, and it keeps going if the run needs
    more. Polyline is offered first so Enter takes it.
    """
    go = Rhino.Input.Custom.GetOption()
    go.SetCommandPrompt("Choose path type")
    go.AddOption("Polyline")
    go.AddOption("Arc")
    if go.Get() != Rhino.Input.GetResult.Option:
        return None
    return go.Option().EnglishName


def pick_second_edge():
    """The edge to end ON. ``(guid, geometry)``, or ``(None, None)``."""
    picked = scene.pick_edge("Select the coarse edge to end ON")
    if picked is None:
        return None, None
    _edge, guid = picked
    rs.SelectObject(guid)
    return guid, rs.coercecurve(guid)


def draw_polyline(start_pt):
    """Polyline mode: free points, ending the moment one lands on the layout.

    The GetPoint is created ONCE and reused across picks. That matters: a fresh
    GetPoint per pick has an empty internal undo history, so Backspace has
    nothing to undo and the keypress is swallowed -- which is why undo did
    nothing before. Reusing one instance keeps that history alive.
    """
    points = [start_pt]
    draw_color = Color.FromArgb(40, 120, 220)

    def _draw(sender, e):
        # Reads the live list, so it always reflects the current path.
        if len(points) >= 2:
            e.Display.DrawPolyline(points, draw_color, 2)
        e.Display.DrawLine(points[-1], e.CurrentPoint, draw_color, 2)
        for point in points:
            e.Display.DrawPoint(point, draw_color)

    print("Pick points freely - Backspace removes the last point.")
    print("Placing a point ON a layout edge finishes the path (Esc to change path type).")

    gp = Rhino.Input.Custom.GetPoint()
    gp.DynamicDraw += _draw
    # The pipeline is planar, and ``mesh_from_faces`` zeroes z anyway. An
    # unconstrained free pick lands on whatever construction plane the view
    # happens to have, so the path drawn and the path cut would differ.
    gp.Constrain(Rhino.Geometry.Plane.WorldXY, False)
    try:
        while True:
            gp.SetCommandPrompt("Pick next point (Backspace = undo)")
            gp.AcceptUndo(True)   # re-armed each pick; without it, no Undo result
            result = gp.Get()

            if result == Rhino.Input.GetResult.Undo:
                if len(points) > 1:
                    points.pop()
                else:
                    print("Nothing left to undo.")
                continue

            if result != Rhino.Input.GetResult.Point:
                return None       # Esc -> back to the mode menu

            point = gp.Point()
            snapped = scene.closest_edge_point(point)
            if snapped:
                points.append(snapped)   # exact point on the layout edge
                break                    # landed on the layout -> path complete
            points.append(point)
    finally:
        gp.DynamicDraw -= _draw

    if len(points) < 2:
        return None
    return rs.AddPolyline(points)


def draw_arc(start_pt):
    """Arc mode: second edge, end point on it, then a point the arc passes
    through -- the same three picks as Rhino's native 3-point _Arc."""
    guid, geometry = pick_second_edge()
    if not guid:
        return None
    try:
        gp = Rhino.Input.Custom.GetPoint()
        gp.SetCommandPrompt("Pick end point on second curve (Esc to change path type)")
        gp.Constrain(geometry, False)
        gp.DrawLineFromPoint(start_pt, True)
        gp.Get()
        if gp.CommandResult() != Rhino.Commands.Result.Success:
            return None
        end_pt = gp.Point()

        while True:
            gp2 = Rhino.Input.Custom.GetPoint()
            gp2.SetCommandPrompt("Pick a point on the arc (Esc to change path type)")

            def _preview(sender, e):
                try:
                    arc = Rhino.Geometry.Arc(start_pt, e.CurrentPoint, end_pt)
                    if arc.IsValid:
                        e.Display.DrawCurve(Rhino.Geometry.ArcCurve(arc), Color.Blue, 2)
                except Exception:
                    pass

            gp2.DynamicDraw += _preview
            gp2.Get()
            gp2.DynamicDraw -= _preview

            if gp2.CommandResult() != Rhino.Commands.Result.Success:
                return None

            arc_id = rs.AddArc3Pt(start_pt, end_pt, gp2.Point())
            if arc_id:
                return arc_id
            print("Those three points do not define a valid arc - pick again.")
    finally:
        rs.UnselectObject(guid)


def sample_drawn_curve(guid):
    """The drawn curve as a point list, sampled finely enough to cut with.

    A line and a polyline are already point lists; an arc is not, and how
    finely it is sampled is what its coarse edge will look like after
    densification -- **the samples ARE the curve from here on**. One sample per
    tenth of a coarse edge, bounded, which is fine next to the separatrices the
    same edges are matched against.
    """
    curve = rs.coercecurve(guid)
    ok, polyline = curve.TryGetPolyline()
    if ok:
        return [[p.X, p.Y, 0.0] for p in polyline]

    length = rs.CurveLength(guid) or 0.0
    step = editor.mean_edge() * 0.1
    count = int(length / step) if step > 0 else 12
    count = max(8, min(64, count))
    return [[p.X, p.Y, 0.0] for p in rs.DivideCurve(guid, count)]


def note_extend(count):
    """Say that a short run is being carried on to the wall. Always ``True``.

    **A new line runs boundary to boundary, or it closes on itself -- there is no
    third option**, and the user is not asked to confirm a rule they cannot opt
    out of. A run that stops on an interior edge leaves the patch behind it with
    five sides, which ``solve_non_quad_faces`` then fans into a quad plus a
    triangle kept as a POLE: a singularity nobody drew. So the run is continued
    along its strip to the wall instead, and the extension is announced rather
    than offered.

    The part the user drew keeps its shape; the continuation is straight, because
    it should look like what it is.
    """
    print("  the run stopped on an interior patch edge ({} end(s)) -- carried on "
          "along the strip to the wall, because a line has to reach a boundary at "
          "both ends or close on itself.".format(count))
    return True


def add_polyedge():
    """Draw a curve across the layout and cut a new POLYEDGE into it.

    The drawn run becomes a new line of coarse edges: every edge it crosses is
    split and every patch it passes through is split in two. That is what makes
    it a polyedge rather than a strip -- it adds CORNERS at mid-edge points,
    where ``add_strip`` unzips corners that were already there.

    Repeats until Esc, so several cuts can be made in one go. The drawn object
    is scratch and is deleted either way -- what is left on screen is the
    redrawn layout, so a line that stays means the mesh has it, and a line that
    vanishes was refused, with a dialog saying why.

    Both ends must sit on the layout's OWN edges, which is why the picks are
    filtered to them: the point has to name a coarse edge to split, and a curve
    on some other layer does not name one.
    """
    changed = False
    while True:
        first_guid = None
        temp_point = None
        obj_id = None
        try:
            picked = scene.pick_edge("Select the coarse edge to start FROM")
            if picked is None:
                return changed
            _edge, first_guid = picked
            rs.SelectObject(first_guid)

            start_pt = rs.GetPointOnCurve(first_guid, "Pick start point on first curve")
            if not start_pt:
                return changed
            temp_point = rs.AddPoint(start_pt)

            # Esc in the mode menu cancels the whole command; a cancelled
            # sub-workflow returns here instead of ending it.
            while True:
                mode = choose_mode()
                if mode is None:
                    return changed
                if mode == "Arc":
                    obj_id = draw_arc(start_pt)
                else:
                    obj_id = draw_polyline(start_pt)
                if obj_id:
                    break
                print("Cancelled - choose the path type again.")
        finally:
            if first_guid:
                rs.UnselectObject(first_guid)
            if temp_point:
                rs.DeleteObject(temp_point)

        points = sample_drawn_curve(obj_id)
        rs.DeleteObject(obj_id)

        cut, _notes = editor.divide(points, extend=note_extend)
        if cut:
            note = editor.last_cut
            print("Line added: {} patch(es) cut, {} corner(s) inserted, layout "
                  "now has {} patch(es).".format(
                      note.get("faces"), note.get("corners"), note.get("faces_out")))
            changed = True
            redraw()
        else:
            mesh_ui.refuse("Add line",
                           "The line was NOT added.\n\n{}.\n\nThe layout is "
                           "unchanged.".format(editor.last_reason))
        sc.doc.Views.Redraw()


# ----------------------------------------------------------------------
# remove a line -- the strip through the picked edge
# ----------------------------------------------------------------------

def remove_strip():
    """Delete a picked strip, until Esc. ``True`` if any went.

    **The strip is picked as a SURFACE, not as an edge.** Each band is drawn as a
    ribbon down its own middle -- see ``mesh_ui.draw_strips`` for why a filled
    band cannot be picked -- so what the user clicks is the thing that goes.

    **Deleting a strip is not a local edit, and the user is told so before it
    happens.** ``plan_strip_deletion`` performs the deletion on a copy, so the
    numbers in the prompt are measured rather than estimated, and a strip that
    cannot go is refused before the question is even asked.
    """
    changed = False
    while True:
        scene.draw_strips(editor.mesh)
        sc.doc.Views.Redraw()
        skey = scene.pick_strip("Pick a strip to remove (Esc to finish)")
        if skey is None:
            scene._clear_strips()
            sc.doc.Views.Redraw()
            return changed

        info = editor.plan_strip_deletion(skey=skey)
        if not info["ok"]:
            mesh_ui.refuse(
                "Remove line",
                "Nothing was removed.\n\n{}.\n\nThe layout is unchanged."
                .format(info["reason"]))
            continue

        # The detail goes to the command history, where it can be scrolled
        # back; the PROMPT has to fit the one line Rhino shows.
        print("strip {}: {} patch(es){}. The two sides weld together, so "
              "surrounding corners move.".format(
                  info["skey"], info["faces"],
                  ", plus {} collateral strip(s) whose patches all lie inside "
                  "it".format(info["collateral"]) if info["collateral"] else ""))
        if info["boundaries_lost"]:
            # On a coarse layout a collapsed boundary means losing a hole the
            # DOMAIN still has, so the layout stops describing the problem.
            print("  WARNING: this would COLLAPSE {} boundary/boundaries -- a "
                  "hole of the domain would no longer be a hole of the layout."
                  .format(info["boundaries_lost"]))

        answer = mesh_ui.ask(
            "Remove strip: {} patch(es){}{}".format(
                info["faces"],
                ", +{} collateral".format(info["collateral"]) if info["collateral"] else "",
                ", COLLAPSES A BOUNDARY" if info["boundaries_lost"] else ""),
            ["Yes", "No"], "Yes")
        if not answer.startswith("y"):
            print("Nothing removed.")
            continue

        removed, _notes = editor.remove_strip(skey=skey)
        if removed:
            note = editor.last_deletion
            print("Line removed: strip {}, patches {} -> {}.".format(
                note.get("skey"), note.get("faces_in"), note.get("faces_out")))
            changed = True
            redraw()
            # The ribbons describe the layout as it WAS; redrawn at the top of
            # the loop, but cleared here so nothing stale is on screen meanwhile.
            scene._clear_strips()
        else:
            mesh_ui.refuse("Remove line",
                           "Nothing was removed.\n\n{}.\n\nThe layout is "
                           "unchanged.".format(editor.last_reason))

# ----------------------------------------------------------------------
# add a strip -- Robin's grammar rule, along existing corners
# ----------------------------------------------------------------------

def add_strip():
    """Pick a run of existing corners and unzip it into a strip.

    **The picks are CORNERS, not free points.** ``add_strip`` takes a polyedge of
    vertices that are already there -- that is what distinguishes it from a
    divide, which inserts new corners at mid-edge points. So the pick is filtered
    to the layout's own corner objects and each one maps straight back to a
    vertex key; the order they are picked in IS the polyedge order.

    Validity is the same rule as everywhere else here: more than two corners, and
    either closed on itself or with both ends on the layout boundary. A strip has
    to run the full width of the layout.
    """
    print("Pick corners in order along the line to unzip. Enter/Esc when done.")
    print("Both ends must be on the layout boundary, or the run must close.")
    polyedge = []
    try:
        while True:
            vkey = scene.pick_vertex(
                "Corner {} of the new line (Esc when done)".format(len(polyedge) + 1))
            if vkey is None:
                break
            if polyedge and vkey == polyedge[-1]:
                print("  that is the same corner again -- skipped.")
                continue
            # Refused at the pick, not at the end: a polyedge is a chain of
            # edges, and finding that out after the tenth corner costs all ten.
            if polyedge and vkey not in editor.mesh.halfedge[polyedge[-1]]:
                print("  corner {} is not joined to corner {} by an edge -- "
                      "skipped. Pick the next corner ALONG an edge.".format(
                          vkey, polyedge[-1]))
                continue
            polyedge.append(vkey)
            # Show the run so far. Picking corner by corner is otherwise blind:
            # the keys are on the command line, but nothing on screen says which
            # way the line is going or whether the last pick was the intended one.
            scene.show_path(polyedge, mesh=editor.mesh)
            print("  {} corner(s): {}".format(len(polyedge), polyedge))
    finally:
        scene.clear_path()
        sc.doc.Views.Redraw()

    if len(polyedge) < 3:
        if polyedge:
            mesh_ui.refuse("Add strip",
                           "A strip needs more than two corners to run along.\n\n"
                           "The layout is unchanged.")
        return False

    ok, notes = editor.add_strip(polyedge)
    if not ok:
        mesh_ui.refuse("Add strip",
                       "The strip was NOT added.\n\n{}.\n\nThe layout is "
                       "unchanged.".format(editor.last_reason))
        return False

    print("Strip added: {} corner(s) unzipped, {} opened, layout now has {} "
          "patch(es).".format(notes.get("corners"), notes.get("opened"),
                              notes.get("faces_out")))
    redraw()
    return True


def divide_strip():
    """Pick a strip by its ribbon and split it down the middle.

    The plain half of :meth:`CoarseEditor.divide` -- no drawing at all. The
    drawn-curve half is ``add_polyedge``; both do the same thing, and which strip
    gets divided is the only question either of them answers.
    """
    changed = False
    while True:
        scene.draw_strips(editor.mesh)
        sc.doc.Views.Redraw()
        skey = scene.pick_strip("Pick a strip to divide in two (Esc to finish)")
        if skey is None:
            scene._clear_strips()
            sc.doc.Views.Redraw()
            return changed

        ok, notes = editor.divide(skey=skey)
        if not ok:
            mesh_ui.refuse("Divide strip",
                           "Nothing was divided.\n\n{}.\n\nThe layout is "
                           "unchanged.".format(editor.last_reason))
            continue
        print("Strip {} divided: {} patch(es) split, layout now has {} "
              "patch(es).".format(notes.get("skey"), notes.get("faces"),
                                  notes.get("faces_out")))
        changed = True
        redraw()
        scene._clear_strips()


# ----------------------------------------------------------------------
# commit
# ----------------------------------------------------------------------

def commit():
    """Hand the layout back. ``True`` if it was adopted."""
    global COMMITTED, COMMIT_NOTES
    layout, notes = editor.commit()
    # ``is None`` rather than truthiness: commit hands back the LAYOUT, and a
    # compas mesh defines no ``__bool__``, so every mesh -- including an empty
    # one -- is truthy.
    if layout is None:
        mesh_ui.refuse(
            "Commit layout",
            "The edited layout was NOT adopted.\n\n{}\n\nThe layout is left as "
            "you edited it, so the offending corner can be moved again."
            .format(notes.get("error", editor.last_reason)))
        return False

    print("Layout adopted: {} patch(es) in, {} out, {} corner(s) snapped".format(
        notes.get("faces_in"), notes.get("faces_out"), notes.get("snapped")))
    if notes.get("off_wall"):
        print("  {} corner(s) sit on no wall -- expected only where a patch was "
              "deleted".format(notes["off_wall"]))
    # 3 is a pseudo-quad pole and is correct; anything else means a patch went
    # out with more points than it has corners, which the repair then fanned
    # into patches nobody drew.
    odd = {n: count for n, count in (notes.get("sides") or {}).items()
           if n not in (3, 4)}
    if odd:
        print("  patches with {} side(s) -- one point per CORNER, not per "
              "sample".format(sorted(odd)))
    if notes.get("lost_curves"):
        print("  {} drawn curve(s) no longer match a coarse edge and were "
              "dropped -- those edges densify straight".format(notes["lost_curves"]))
    print("  commit path: {}".format(notes.get("path")))
    COMMITTED, COMMIT_NOTES = layout, notes
    redraw()
    return True


# ----------------------------------------------------------------------
# run
# ----------------------------------------------------------------------

# The baked mesh is in the way while picking -- it still shows the corners
# where they were. Hidden rather than deleted, so a cancelled edit leaves the
# document as it was found.
mesh_guids = rs.ObjectsByLayer(MESH_LAYER)
rs.HideObjects(mesh_guids)

redraw()
lock_state = scene.unlock()

try:
    while True:
        operation = mesh_ui.ask(
            "next", ["move_vertices", "add_polyedge", "add_strip", "divide_strip",
                     "remove_strip", "reset", "commit", "finish"], "finish")

        if operation == "move_vertices":
            move_vertices()
        elif operation == "add_polyedge":
            add_polyedge()
        elif operation == "add_strip":
            add_strip()
        elif operation == "divide_strip":
            divide_strip()
        elif operation == "remove_strip":
            remove_strip()
        elif operation == "reset":
            editor.reset()
            redraw()
            print("Back to the layout this session started with.")
        elif operation == "commit":
            if commit():
                break
        else:
            break
finally:
    # Everything the command put on screen comes off, in the order that leaves
    # the document as it was found even if an operation raised.
    scene.relock(lock_state)
    scene.clear()
    clear_layer(EDIT_LAYER)
    rs.ShowObjects([guid for guid in mesh_guids if rs.IsObject(guid)])
    sc.doc.Views.Redraw()


# ----------------------------------------------------------------------
# write it back -- only if something was committed
# ----------------------------------------------------------------------

if COMMITTED is None:
    print("nothing committed -- '{}' is unchanged.".format(MESH_LAYER))
else:
    # ``commit()`` mutated the layout read at the top of this file in place, so
    # this IS ``coarse`` -- returned as well so the reference is explicit.
    committed = COMMITTED
    bake_mesh(committed, MESH_LAYER)

    # Poles come from the COMMITTED mesh, keys and coordinates both. Reading
    # keys from one mesh and coordinates from another is not a style point:
    # a rebuilding commit renumbers -- measured 2 of 8 unchanged on a disc -- so
    # the old code baked pole points at arbitrary wrong corners, or raised.
    pole_keys = committed.poles() if hasattr(committed, "poles") else []
    bake_points([committed.vertex_coordinates(vkey) for vkey in pole_keys],
                POLES_LAYER)

    # The traced separatrix network, plus any curve the user drew. These are
    # what CMD_quad_mesh matches coarse edges against, so a stale set here
    # silently costs curvature there.
    #
    # They are NOT re-traced against the edited layout -- ``_build`` is cached
    # and the field was not re-solved -- so after a corner is moved they will
    # visibly not line up with it. That is correct rather than a glitch: an
    # edge that no longer matches its old separatrix falls back to a chord, and
    # seeing the mismatch is seeing what will actually be densified.
    #
    # ``editor.all_polylines``, which is the traced network this command READ
    # from the layer plus the curves the user drew during this session. That the
    # two are published together is not cosmetic: CMD_quad_mesh rebuilds the
    # edge-curve mapping from this layer ALONE, so a drawn arc that misses this
    # bake comes out a straight chord there. The old route through the
    # decomposition had a matching trap -- ``decomposition_polylines()``
    # reassigned the list and dropped the user curves, measured 48 entries down
    # to 47 -- which is the same failure reached by a different road.
    ribs = editor.all_polylines
    _guids, skipped = bake_polylines(ribs, POLYLINE_LAYER)
    print("separatrices: {} baked".format(len(_guids)))
    if skipped:
        # Not fatal, and it must not be: this runs AFTER the mesh is written, so
        # raising here would leave the document half updated -- the layout new,
        # this layer emptied by the clear and never refilled, and the next
        # command quietly chording every edge. rs.AddPolyline RAISES on geometry
        # Rhino refuses, which is exactly how that used to happen.
        print("  {} rib(s) skipped -- Rhino refused them (coincident points at "
              "the document tolerance, or fewer than two)".format(skipped))

    # Refresh the curved layout. Left alone it would still show the shape from
    # before the edit, which is worse than not drawing it at all. Rebuilt with
    # the SAME document-only builder CMD_quad_mesh uses, so what is on screen
    # is what will actually be densified -- including a moved interior edge
    # falling back to a chord, which is the truthful state and worth seeing.
    outer_loop, inner_loops = read_boundary_loops(WALL_SAMPLING)
    curves, tally = coarse_edges_to_curves(
        committed,
        loops=[outer_loop] + inner_loops,
        polylines=read_polylines(POLYLINE_LAYER))
    _guids, skipped = bake_edge_curves(curves.values(), EDGE_CURVE_LAYER)
    print("coarse edges: {}".format(tally))
    if skipped:
        print("  {} edge curve(s) Rhino refused -- those edges densify as "
              "chords".format(skipped))
    if tally.get("wall_missed"):
        print("  {} boundary edge(s) did NOT get their wall arc -- on a hole "
              "this is what makes one circle come out round and the next a "
              "polygon".format(tally["wall_missed"]))

    notes = COMMIT_NOTES
    print("committed: {} patch(es) in, {} out, {} corner(s) snapped to a wall"
          .format(notes.get("faces_in"), notes.get("faces_out"),
                  notes.get("snapped")))
    if notes.get("off_wall"):
        print("  {} corner(s) sit on no wall -- expected only where a patch was "
              "deleted".format(notes["off_wall"]))
    # Densities stored in step 5 are keyed on edge midpoints, and this edit
    # moved some of those. Step 5 reports how many overrides survived.
    print("layout written back to '{}'".format(MESH_LAYER))

    # The side-car, refreshed to match. Steps 5 and 6 prefer it over the baked
    # mesh because it carries what a bake cannot -- strips, ``face_pole``, the
    # route -- but only while it still describes the SAME layout; ``read_layout``
    # compares the two and falls back to the document otherwise. Leaving a stale
    # one here would therefore not be wrong, just wasteful: every later step
    # would silently do the re-derivation this exists to avoid.
    committed.attributes.setdefault("route", "field")
    committed.attributes["edited"] = True
    committed.save_to_json(cache_path(COARSE_CACHE))

print("next: CMD_densities")
