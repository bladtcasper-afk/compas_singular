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

Four operations:

* **move_vertices** -- drag a corner. A corner on the layout boundary is held
  on the domain wall, so nudging the edge of the layout does not eat the
  outline;
* **add_lines** -- draw a Line, Polyline or Arc across the layout and CUT it:
  every coarse edge crossed is split, every patch passed through is split in
  two. The drawn SHAPE is kept, so an arc densifies as an arc;
* **remove_lines** -- pick an edge, and the whole STRIP through it is deleted.
  The two sides weld together;
* **commit** -- hand the layout to ``FieldDecomposition.edit_coarse``, which
  re-matches every edge to the separatrix it came from and warps that curve
  onto the moved corners. **Until this runs the edit is worth nothing**:
  densifying the moved layout directly gives straight chords and throws away
  the field alignment the whole front end exists to produce.

**Nothing is written to the document until commit.** ``reset`` goes back to the
layout this session started with, and leaving without committing leaves
``Skeleton::Mesh`` exactly as it was found.

**Where the code lives, and why.** Everything in this file is Rhino: picking,
drawing, prompting, baking. The editing itself -- the cut planner, the all-quad
gate, the strip deletion, the curve map, the commit -- is
``compas_singular.editing.CoarseLayoutEditor`` and imports no ``rs`` at all, so
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


from compas_singular.editing import CoarseLayoutEditor
from compas_singular.rhino import mesh_ui
from compas_singular.rhino.coarse_curves import coarse_edges_to_curves
from compas_singular.rhino.helpers.helpers import bake_edge_curves
from compas_singular.rhino.helpers.helpers import bake_mesh
from compas_singular.rhino.helpers.helpers import bake_points
from compas_singular.rhino.helpers.helpers import bake_polylines
from compas_singular.rhino.helpers.helpers import clean_layer
from compas_singular.rhino.helpers.helpers import read_boundaries
from compas_singular.rhino.helpers.helpers import read_boundary_loops
from compas_singular.rhino.helpers.helpers import read_coarse
from compas_singular.rhino.helpers.helpers import read_polylines
from CMD_start import get_settings
from CMD_start import get_decomposition


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

# Needed to commit: ``edit_coarse`` warps the separatrices onto the moved
# corners and snaps to the domain walls, and only the decomposition has either.
# Solved once per document and reused across steps 3, 4 and 6 -- and across
# Rhino restarts. This call used to omit ``symmetry``, which the other two
# passed, so step 4 could solve under a different group than the layout it was
# about to edit. ``CMD_start.get_decomposition`` is now the only place any of
# the three resolves anything, which is what ``resolve_symmetry`` asks for.
decomposition = get_decomposition()

# Belt and braces. The cache reconstructs a pristine object per call, so notes
# from an edit committed earlier this session cannot reach here -- if they ever
# did, a run the user CANCELS would look committed and the stale layout would
# be baked over the current one.
decomposition.edit_notes = {}

editor = CoarseLayoutEditor(decomposition, coarse=coarse, loops=[outer] + inners)
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
            # is drawn is what ``edit_coarse`` will get.
            project=lambda point: editor.project(point, on_boundary))
        if xyz is None:
            return changed
        if editor.move_vertex(vkey, xyz):
            changed = True
            redraw()


# ----------------------------------------------------------------------
# draw a curve across the layout -- the part only this command needs
# ----------------------------------------------------------------------

def choose_mode():
    """Line / Polyline / Arc, or ``None`` if the user cancelled (Esc)."""
    go = Rhino.Input.Custom.GetOption()
    go.SetCommandPrompt("Choose path type")
    go.AddOption("Line")
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


def draw_line(start_pt):
    """Line mode: pick a second edge, then a point on it."""
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
        return rs.AddLine(start_pt, gp.Point())
    finally:
        rs.UnselectObject(guid)


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


def ask_extend(count):
    """Offer to carry a cut that stopped inside the layout on to the wall."""
    answer = mesh_ui.ask(
        "The line stops on an interior patch edge ({} end(s)), which would "
        "leave a 5-sided patch. Extend the cut to the boundary?".format(count),
        ["Yes", "No"], "Yes")
    return answer.startswith("y")


def add_lines():
    """Draw a curve across the layout and CUT the layout with it.

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
                if mode == "Line":
                    obj_id = draw_line(start_pt)
                elif mode == "Polyline":
                    obj_id = draw_polyline(start_pt)
                else:
                    obj_id = draw_arc(start_pt)
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

        if editor.insert_curve(points, extend=ask_extend):
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

def remove_lines():
    """Delete the strip through a picked edge, until Esc. ``True`` if any went.

    **Deleting a strip is not a local edit, and the user is told so before it
    happens.** ``plan_strip_deletion`` performs the deletion on a copy, so the
    numbers in the prompt are measured rather than estimated, and a strip that
    cannot go is refused before the question is even asked.
    """
    changed = False
    while True:
        picked = scene.pick_edge("Pick an edge; the whole strip through it is removed")
        if picked is None:
            return changed
        edge, _guid = picked

        info = editor.plan_strip_deletion(edge)
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

        if editor.delete_strip(edge):
            note = editor.last_deletion
            print("Line removed: strip {}, patches {} -> {}.".format(
                note.get("skey"), note.get("faces_in"), note.get("faces_out")))
            changed = True
            redraw()
        else:
            mesh_ui.refuse("Remove line",
                           "Nothing was removed.\n\n{}.\n\nThe layout is "
                           "unchanged.".format(editor.last_reason))


# ----------------------------------------------------------------------
# commit
# ----------------------------------------------------------------------

def commit():
    """Hand the layout back. ``True`` if it was adopted."""
    ok, notes = editor.commit()
    if not ok:
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
            "next", ["move_vertices", "add_lines", "remove_lines", "reset",
                     "commit", "finish"], "finish")

        if operation == "move_vertices":
            move_vertices()
        elif operation == "add_lines":
            add_lines()
        elif operation == "remove_lines":
            remove_lines()
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
    clean_layer(EDIT_LAYER)
    rs.ShowObjects([guid for guid in mesh_guids if rs.IsObject(guid)])
    sc.doc.Views.Redraw()


# ----------------------------------------------------------------------
# write it back -- only if something was committed
# ----------------------------------------------------------------------

if not decomposition.edit_notes:
    print("nothing committed -- '{}' is unchanged.".format(MESH_LAYER))
else:
    committed = decomposition.mesh
    bake_mesh(committed, MESH_LAYER)

    # Poles come from the COMMITTED mesh, keys and coordinates both. Reading
    # keys from one mesh and coordinates from another is not a style point:
    # ``edit_coarse`` renumbers -- measured 2 of 8 unchanged on a disc -- so
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
    # ``decomposition.polylines`` rather than ``decomposition_polylines()``, and
    # the difference is not cosmetic: that method REASSIGNS
    # ``self.polylines = boundary + others``, discarding the user curves
    # ``commit`` has just published through ``set_user_curves`` -- measured 48
    # entries down to 47. Since CMD_quad_mesh rebuilds from this layer alone, a
    # drawn arc that is not baked here comes out a straight chord there.
    ribs = decomposition.polylines or decomposition.decomposition_polylines()
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

    notes = decomposition.edit_notes
    print("committed: {} patch(es) in, {} out, {} corner(s) snapped to a wall"
          .format(notes.get("faces_in"), notes.get("faces_out"),
                  notes.get("snapped")))
    if notes.get("off_wall"):
        print("  {} corner(s) sit on no wall -- expected only where a patch was "
              "deleted".format(notes["off_wall"]))
    # Densities stored in step 5 are keyed on edge midpoints, and this edit
    # moved some of those. Step 5 reports how many overrides survived.
    print("layout written back to '{}'".format(MESH_LAYER))

print("next: CMD_densities")
