#! python3

# r: compas

"""**Step 7 -- edit the FINAL quad mesh by hand.**

    reads   TopologyProblem::QuadMesh
    writes  TopologyProblem::QuadMesh                 (on 'save' only)
            TopologyProblem::QuadMesh::Unedited       this session's starting mesh
    scratch TopologyProblem::QuadMesh::Edit           pickable points and lines

Three operations, and they are three different things even though a user calls
two of them "adding and removing lines":

* **move_vertex** -- drag a vertex. Topology unchanged. A vertex on the mesh
  boundary is projected back onto the domain wall, so nudging the edge of the
  mesh does not eat the outline;
* **add_line** -- pick an edge, and the whole POLYEDGE through it gains a strip
  beside it (``grammar.add_strip``). One new row of elements, wall to wall;
* **remove_line** -- pick an edge, and the whole STRIP through it is deleted
  (``grammar_pattern.delete_strip``). The two sides weld together.

**A successful save ENDS the command.** Nothing is written to the document until
then, so the loop is the edit and saving is leaving it -- and the mesh on
``QuadMesh`` is what CMD_smoothen and CMD_dual read next, so handing control back
is the honest end of this step. A FAILED save keeps the loop open: dropping the
edit because the bake did not take would be the worst of both.

**Why a line is never one edge.** A quad mesh cannot gain or lose a single edge:
removing one merges two quads into a hexagon, and adding one leaves a five-sided
face. The only operations that keep every face a quad run the whole width of the
mesh -- which is what a strip is. So the pick is an edge and the unit of change
is the strip through it. Both operations are planned on a COPY and adopted only
if the copy is still all-quad, so a refusal costs nothing.

**This is the end of the line, and that is deliberate.** A hand-edited dense mesh
has no coarse layout that produces it, so re-running CMD_quad_mesh regenerates
from the layout and throws the edit away. Nothing here writes back to
``Skeleton::Mesh`` and nothing tries to replay an edit onto a new layout -- the
correspondence is not well defined once the layout changes. If an edit can be
expressed on the COARSE layout instead, do it there: CMD_edit_coarse_mesh keeps
the strip grammar, the field alignment and the boundary curvature, and a mistake
there costs one 'reset'.

**Strip operations need an all-quad mesh, and a mesh with poles is not one.**
``collect_strip`` walks ``face_opposite_edge``, and the triangle fan around a
pole has none. Rather than produce a quietly wrong strip, this command counts
non-quad faces on load and disables add_line and remove_line if there are any.
move_vertex still works, because it does not touch topology.

**Nothing is keyed by index across an edit.** Strip and polyedge keys come from
popping ``self.edges()``, so they follow vertex insertion order -- and every edit
here renumbers. So they are re-collected from the live mesh immediately before
each use and never stored between operations. That is the only safe way to use
them at all, and it is why every pick starts from an edge the user clicked rather
than from a remembered key.
"""

# ----------------------------------------------------------------------
# imports
# ----------------------------------------------------------------------

# Temporary import of compas_singular development library. MUST come before any
# compas_singular import: it fixes sys.path and purges a stale copy from
# sys.modules, and doing that after a class is bound would leave an existing mesh
# failing isinstance against the rebuilt class.
from CMD_start import import_compas_singular
import_compas_singular()

import rhinoscriptsyntax as rs
import scriptcontext as sc
import Rhino
from Rhino.Geometry import Line as RhinoLine
from Rhino.Geometry import Point3d
# ``from System.Drawing.Color import FromArgb`` is an IronPython idiom and raises
# under Rhino 8's CPython, where Color is a class and not a module.
from System.Drawing import Color

from compas.itertools import pairwise

import compas_rhino as cr
from compas_rhino.conversions import mesh_to_compas

from compas_singular.datastructures import QuadMesh
from compas_singular.datastructures.mesh_quad.grammar.add_strip import add_strip
from compas_singular.datastructures.mesh_quad.grammar.add_strip import (
    is_polyedge_valid_for_strip_addition)
from compas_singular.datastructures.mesh_quad.grammar.delete_strip import (
    strips_to_split_to_prevent_boundary_collapse)
# From grammar_pattern, NOT grammar.delete_strip: the two modules both define
# delete_strip and the package __init__ imports grammar last, so the plain name
# resolves to the other one. This is the version that merges the two sides,
# handles collateral deletions and repairs 'face_pole'.
from compas_singular.datastructures.mesh_quad.grammar_pattern import delete_strip
from compas_singular.datastructures.mesh_quad.grammar_pattern import split_strips
from compas_singular.datastructures.mesh_quad.grammar_pattern import (
    collateral_strip_deletions)
from compas_singular.datastructures.mesh_quad.grammar_pattern import (
    total_boundary_deletions)
# compas 2.x dropped mesh_smooth_centroid from compas.datastructures; the shim is
# what the rest of this codebase uses.
from compas_singular._compat import mesh_smooth_centroid
from compas_singular.rhino.coarse_curves import BoundaryLoop
from compas_singular.rhino.helpers.helpers import bake_mesh
from compas_singular.rhino.helpers.helpers import clear_layer
from compas_singular.rhino.helpers.helpers import read_boundary_loops
from compas_singular.rhino.helpers.helpers import read_mesh
from CMD_start import get_settings


# ----------------------------------------------------------------------
# layers and constants
# ----------------------------------------------------------------------

ROOT = "TopologyProblem"
QUADMESH_LAYER = ROOT + "::QuadMesh"
EDIT_LAYER = QUADMESH_LAYER + "::Edit"
UNEDITED_LAYER = QUADMESH_LAYER + "::Unedited"
EDITED_LAYER = QUADMESH_LAYER + "::Edited"

COLOR_VERTEX = (0, 0, 0)
COLOR_VERTEX_BOUNDARY = (0, 120, 200)
COLOR_VERTEX_SINGULAR = (200, 0, 120)
COLOR_EDGE = (110, 110, 110)
COLOR_PREVIEW = (255, 120, 0)

#: Above this many edges the scratch layer is enough objects to make Rhino's
#: redraw noticeably slow, so it is worth asking rather than just doing it.
DRAW_WARN_EDGES = 4000

#: Smoothing passes used to open a newly added strip. ``add_strip`` creates its
#: two vertices ON TOP of the one they replace, so without this the new strip has
#: zero width and is invisible.
RELAX_ITERATIONS = 50


def ensure_layer(path, color=None):
    """Create a ``::`` layer path, parents first, and return the full path.

    ``rs.AddLayer`` does not create intermediate parents, and ``bake_mesh``'s own
    layer creation hard-codes ``parent="TopologyProblem"``, which is wrong for a
    layer two levels down.
    """
    parts = path.split("::")
    for i in range(len(parts)):
        name = "::".join(parts[:i + 1])
        if not rs.IsLayer(name):
            rs.AddLayer(name=parts[i],
                        parent="::".join(parts[:i]) if i else None,
                        color=color if i == len(parts) - 1 else None)
    return path


# ----------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------

def load_quad_mesh():
    """The PICKED mesh as a :class:`QuadMesh`, with the object it came from.

    Returns ``(mesh, guid, layer)``. Both extras are needed downstream: the source
    object is hidden while editing and shown again on the way out, and a mesh
    picked from a sublayer must not be confused with the primal one.

    ``mesh_to_compas`` gives a plain compas ``Mesh``; the strip and polyedge
    grammar lives on ``QuadMesh``, so it is rebuilt here from vertices and faces.
    Rhino has no triangular mesh face -- it stores one as a quad with a repeated
    corner -- so the duplicate is dropped, or the face is degenerate and every
    valence through it is wrong.
    """
    if not rs.IsLayer(QUADMESH_LAYER):
        raise RuntimeError(
            "No layer '{}' -- run CMD_quad_mesh first.".format(QUADMESH_LAYER))

    quad_layer = rs.LayerName(QUADMESH_LAYER, fullpath=True)

    def quad_mesh_filter(rhobj, geometry, component_index):
        layer = rs.ObjectLayer(rhobj)
        return layer == quad_layer or rs.IsLayerChildOf(layer, quad_layer)

    guid = rs.GetObject(
        message="Pick the quad mesh to edit, from 'QuadMesh' or a sublayer",
        filter=rs.filter.mesh, preselect=False, select=False,
        custom_filter=quad_mesh_filter, subobjects=False)
    if not guid:
        print("Cancelled -- nothing was drawn or changed.")
        raise SystemExit

    # The PICKED geometry. Re-reading the layer here is what made the prompt
    # decoration: it returns the first object on QuadMesh however the user
    # answered.
    plain = mesh_to_compas(cr.objects.find_object(guid).Geometry)
    vertices, faces = plain.to_vertices_and_faces()
    faces = [[v for i, v in enumerate(f) if v != f[i - 1]] for f in faces]
    return (QuadMesh.from_vertices_and_faces(vertices, faces),
            guid, rs.ObjectLayer(guid))


def load_walls():
    """The domain walls, for holding a moved boundary vertex on the outline.

    Optional: without the input boundaries in the document the command still
    works, the boundary simply stops being held. Sampled at the resolution
    CMD_quad_mesh densified along, so a projected vertex lands where the rest of
    the dense boundary already is.
    """
    try:
        spacing = get_settings()["triangulation_spacing"]
        outer_loop, inner_loops = read_boundary_loops(spacing * 0.25)
    except Exception as e:
        print("boundary walls unavailable ({}) -- a moved boundary vertex will "
              "NOT be held on the outline.".format(e))
        return []
    return [BoundaryLoop(loop) for loop in [outer_loop] + inner_loops
            if len(loop) >= 3]


def non_quad_faces(mesh):
    return [fkey for fkey in mesh.faces() if len(mesh.face_vertices(fkey)) != 4]


def boundary_vertices(mesh):
    """Every boundary ring, not just the longest -- see :func:`relax`."""
    out = set()
    for ring in mesh.vertices_on_boundaries():
        out.update(ring)
    return out


# ----------------------------------------------------------------------
# display
# ----------------------------------------------------------------------

GUID_VERTICES = {}      # guid -> vkey
GUID_EDGES = {}         # guid -> (u, v)


def clear_scratch():
    guids = [guid for guid in list(GUID_VERTICES) + list(GUID_EDGES)
             if rs.IsObject(guid)]
    if guids:
        rs.DeleteObjects(guids)
    GUID_VERTICES.clear()
    GUID_EDGES.clear()


def draw(mesh):
    """Vertices as points and edges as lines, on the scratch layer.

    Everything is redrawn rather than updated in place: every edit renumbers
    keys, and a stale guid map is a real way to move the wrong vertex. Singular
    vertices are coloured because they are usually what the user is trying to
    move, or to avoid.
    """
    clear_scratch()
    ensure_layer(EDIT_LAYER, (0, 150, 100))

    boundary = boundary_vertices(mesh)
    try:
        singular = set(mesh.singularities())
    except Exception:
        # A convenience, not the point: a mesh mid-edit may not have the
        # consistent valences singularities() expects.
        singular = set()

    rs.EnableRedraw(False)
    try:
        for vkey in mesh.vertices():
            guid = rs.AddPoint(Point3d(*mesh.vertex_coordinates(vkey)))
            rs.ObjectLayer(guid, EDIT_LAYER)
            if vkey in singular:
                rs.ObjectColor(guid, COLOR_VERTEX_SINGULAR)
            elif vkey in boundary:
                rs.ObjectColor(guid, COLOR_VERTEX_BOUNDARY)
            else:
                rs.ObjectColor(guid, COLOR_VERTEX)
            GUID_VERTICES[guid] = vkey

        for u, v in mesh.edges():
            guid = rs.AddLine(Point3d(*mesh.vertex_coordinates(u)),
                              Point3d(*mesh.vertex_coordinates(v)))
            # AddLine returns None on a zero-length line rather than raising,
            # and ObjectLayer(None) does raise.
            if not guid:
                continue
            rs.ObjectLayer(guid, EDIT_LAYER)
            rs.ObjectColor(guid, COLOR_EDGE)
            GUID_EDGES[guid] = (u, v)
    finally:
        rs.EnableRedraw(True)
    sc.doc.Views.Redraw()


def pick_vertex():
    """A vertex of the mesh, by its scratch point. ``None`` if not one."""
    guid = rs.GetObject("Select a mesh vertex", rs.filter.point, preselect=True)
    if not guid:
        return None
    vkey = GUID_VERTICES.get(guid)
    if vkey is None:
        print("Not a vertex of the quad mesh -- pick one of the points on "
              "'{}'.".format(EDIT_LAYER))
        return None
    return vkey


def pick_edge(message):
    """An edge of the mesh, by its scratch line. ``None`` if not one."""
    guid = rs.GetObject(message, rs.filter.curve, preselect=True)
    if not guid:
        return None
    edge = GUID_EDGES.get(guid)
    if edge is None:
        print("Not an edge of the quad mesh -- pick one of the lines on "
              "'{}'.".format(EDIT_LAYER))
        return None
    return edge


# ----------------------------------------------------------------------
# move a vertex
# ----------------------------------------------------------------------

def project_to_wall(xyz, walls):
    """The nearest point of the nearest wall, or ``xyz`` if there is no wall.

    Only ever called for a vertex that is TOPOLOGICALLY on the mesh boundary.
    Projecting anything merely near a wall would drag an interior vertex that
    legitimately sits in a narrow slot onto it -- the same rule
    ``edit.snap_to_loops`` keeps, and for the same reason.
    """
    if not walls:
        return xyz
    best_d, best_p = float('inf'), None
    for wall in walls:
        d, _s, q = wall.project(xyz)
        if d < best_d:
            best_d, best_p = d, q
    return best_p if best_p is not None else xyz


def move_vertex(mesh, walls):
    """Drag one vertex, previewing the mesh as it deforms. ``True`` if moved.

    The preview shows the PROJECTED position, not the cursor, so a boundary
    vertex is drawn where it will actually end up. Constrained to world XY
    because the pipeline is planar: an unconstrained pick lands on whatever
    construction plane the view happens to have, and the vertex moved would then
    differ from the one previewed.
    """
    vkey = pick_vertex()
    if vkey is None:
        return False

    on_boundary = vkey in boundary_vertices(mesh)
    neighbours = [mesh.vertex_coordinates(nbr)
                  for nbr in mesh.vertex_neighbors(vkey)]
    start = mesh.vertex_coordinates(vkey)
    color = Color.FromArgb(*COLOR_PREVIEW)

    gp = Rhino.Input.Custom.GetPoint()
    gp.SetCommandPrompt("New position for vertex {}{}".format(
        vkey, " (held on the wall)" if on_boundary and walls else ""))
    gp.SetBasePoint(Point3d(*start), True)
    gp.Constrain(Rhino.Geometry.Plane.WorldXY, False)

    def OnDynamicDraw(sender, e):
        cp = e.CurrentPoint
        xyz = [cp.X, cp.Y, 0.0]
        if on_boundary:
            xyz = project_to_wall(xyz, walls)
        point = Point3d(*xyz)
        for nbr in neighbours:
            e.Display.DrawLine(RhinoLine(point, Point3d(*nbr)), color, 2)
        e.Display.DrawPoint(point, color)

    gp.DynamicDraw += OnDynamicDraw
    try:
        gp.Get()
        if gp.CommandResult() != Rhino.Commands.Result.Success:
            return False
        ep = gp.Point()
    finally:
        gp.DynamicDraw -= OnDynamicDraw

    xyz = [ep.X, ep.Y, 0.0]
    if on_boundary:
        xyz = project_to_wall(xyz, walls)
    mesh.vertex_attributes(vkey, 'xyz', xyz)
    return True


def move_vertices(mesh, walls):
    """Move vertices until Enter or Esc. ``True`` if any moved."""
    changed = False
    while True:
        if not move_vertex(mesh, walls):
            return changed
        changed = True
        draw(mesh)


# ----------------------------------------------------------------------
# add a line -- a strip along the polyedge through the picked edge
# ----------------------------------------------------------------------

def polyedge_through(mesh, edge):
    """``(pkey, polyedge)`` for the polyedge containing ``edge``, as a new list.

    ``collect_polyedges`` is re-run every time rather than cached: an added or
    deleted strip renumbers the keys, and a remembered pkey would then name a
    different polyedge without erroring. The list is a COPY because ``add_strip``
    pops from the list it is given, which would corrupt
    ``attributes['polyedges']`` if it were handed the stored one.
    """
    u, v = edge
    for pkey, polyedge in dict(mesh.collect_polyedges()).items():
        for a, b in pairwise(polyedge):
            if (a, b) == (u, v) or (a, b) == (v, u):
                return pkey, list(polyedge)
    return None, None


def add_line(mesh):
    """Add a strip along the polyedge through a picked edge.

    Returns ``(mesh to keep, changed)`` -- the original on a refusal, the edited
    copy on success, so a refusal leaves the mesh exactly as it was.
    """
    edge = pick_edge("Pick an edge; the whole line through it gains a strip")
    if edge is None:
        return mesh, False

    pkey, polyedge = polyedge_through(mesh, edge)
    if polyedge is None:
        print("That edge belongs to no polyedge -- nothing to add along.")
        return mesh, False

    closed = polyedge[0] == polyedge[-1]
    # >2 vertices, and closed OR both ends on a boundary. A strip has to run the
    # full width of the mesh, for the same reason a cut on the coarse layout has
    # to reach a wall: anything less leaves a face with five sides.
    if not is_polyedge_valid_for_strip_addition(mesh, polyedge):
        rs.MessageBox(
            "No strip was added.\n\nThe line through that edge has {} vertices, "
            "is not closed, and does not end on the boundary at both ends. A "
            "strip has to run the full width of the mesh -- wall to wall, or all "
            "the way round.\n\nThe mesh is unchanged.".format(len(polyedge)),
            0, "Add line")
        return mesh, False

    work = mesh.copy()
    # A PRECONDITION, not housekeeping: add_strip writes into
    # attributes['strips'] and fails without it.
    work.collect_strips()
    # Re-resolved on the copy: mesh.copy() preserves keys today, but relying on
    # that would make this break silently if it ever stopped being true.
    _pkey, work_polyedge = polyedge_through(work, edge)
    if work_polyedge is None:
        print("Lost track of that line on the working copy -- nothing changed.")
        return mesh, False

    try:
        skey, old_to_new = add_strip(work, list(work_polyedge))
    except Exception as e:
        rs.MessageBox("No strip was added.\n\nadd_strip failed: {}\n\nThe mesh is "
                      "unchanged.".format(e), 0, "Add line")
        return mesh, False

    bad = non_quad_faces(work)
    if bad:
        rs.MessageBox(
            "No strip was added.\n\nThe result would have {} face(s) that are not "
            "quads.\n\nThe mesh is unchanged.".format(len(bad)), 0, "Add line")
        return mesh, False

    # add_strip creates its two vertices ON TOP of the vertex they replace, so the
    # new strip has zero width and is invisible until the mesh is relaxed.
    # old_to_new maps old vkey -> the pair that has to separate, and those are the
    # only boundary vertices allowed to move.
    new_vertices = set()
    for pair in old_to_new.values():
        new_vertices.update(pair)
    relax(work, free=new_vertices)

    print("Line added: strip {} along polyedge {} ({} vertices, {}), "
          "faces {} -> {}".format(
              skey, pkey, len(work_polyedge),
              "closed" if closed else "wall to wall",
              mesh.number_of_faces(), work.number_of_faces()))
    return work, True


# ----------------------------------------------------------------------
# remove a line -- the strip through the picked edge
# ----------------------------------------------------------------------

def remove_line(mesh):
    """Delete the strip through a picked edge.

    **Deleting a strip is not a local edit, and the user has to be told so before
    it happens.** Three things follow from what ``delete_strip`` does, and all
    three are predicted on the untouched mesh and put in the prompt:

    * it MERGES the two sides of the strip -- for each disconnected part of the
      strip's edge network one new vertex is added at the centroid and every old
      vertex is substituted for it, so surrounding vertices move;
    * it can take other strips with it: any strip whose faces all lie inside the
      deleted one goes too;
    * it can COLLAPSE a boundary that would be left with too few splits.

    Returns ``(mesh to keep, changed)``.
    """
    edge = pick_edge("Pick an edge; the whole strip through it is removed")
    if edge is None:
        return mesh, False

    # Strips are collected on the copy so the live mesh keeps no state we did not
    # ask for, and re-collected every time because an earlier edit renumbered.
    work = mesh.copy()
    work.collect_strips()
    skey = work.edge_strip(edge)
    if skey is None:
        print("That edge belongs to no strip -- nothing to remove.")
        return mesh, False

    faces = len(work.strip_faces(skey))
    collateral = list(collateral_strip_deletions(work, [skey]))
    boundaries_lost = total_boundary_deletions(work, [skey])

    print("strip {}: {} face(s){}{}. The two sides weld together, so "
          "surrounding vertices move.".format(
              skey, faces,
              ", plus {} collateral strip(s) whose faces all lie inside "
              "it".format(len(collateral)) if collateral else "",
              ", and {} boundary/boundaries would COLLAPSE".format(
                  boundaries_lost) if boundaries_lost else ""))

    options = ["Yes", "No"]
    if boundaries_lost:
        # Pre-splitting the strips that would otherwise collapse is the
        # documented remedy. It is done here rather than through delete_strip's
        # own preserve_boundaries flag, which calls a helper --
        # boundary_strip_preserve -- that is not defined anywhere in the package
        # and raises NameError.
        options = ["Yes", "PreserveBoundaries", "No"]
    answer = (rs.GetString(
        "Remove strip: {} face(s){}{}".format(
            faces,
            ", +{} collateral".format(len(collateral)) if collateral else "",
            ", COLLAPSES A BOUNDARY" if boundaries_lost else ""),
        options[0], options) or "no").lower()

    if answer.startswith("n"):
        print("Nothing removed.")
        return mesh, False

    if answer.startswith("p"):
        to_split = strips_to_split_to_prevent_boundary_collapse(work, [skey])
        if not to_split:
            rs.MessageBox(
                "Nothing removed.\n\nThis strip cannot go without collapsing a "
                "boundary, and there is no strip left to split to prevent it."
                "\n\nThe mesh is unchanged.", 0, "Remove line")
            return mesh, False
        split_strips(work, to_split)
        print("  pre-split {} strip(s) to keep the boundaries".format(len(to_split)))
        # to_split never contains the strip being deleted, but the split rebuilt
        # the strip data, so re-resolve from the same picked edge rather than
        # trusting the skey across it.
        skey = work.edge_strip(edge)
        if skey is None:
            print("Lost track of the strip after the split -- nothing changed.")
            return mesh, False

    try:
        delete_strip(work, skey)
    except Exception as e:
        rs.MessageBox("Nothing removed.\n\ndelete_strip failed: {}\n\nThe mesh is "
                      "unchanged.".format(e), 0, "Remove line")
        return mesh, False

    if work.number_of_faces() == 0:
        rs.MessageBox("Nothing removed.\n\nThat would delete every face.\n\nThe "
                      "mesh is unchanged.", 0, "Remove line")
        return mesh, False

    bad = non_quad_faces(work)
    if bad:
        rs.MessageBox(
            "Nothing removed.\n\nThe result would have {} face(s) that are not "
            "quads.\n\nThe mesh is unchanged.".format(len(bad)), 0, "Remove line")
        return mesh, False

    print("Line removed: strip {}, faces {} -> {}, vertices {} -> {}".format(
        skey, mesh.number_of_faces(), work.number_of_faces(),
        mesh.number_of_vertices(), work.number_of_vertices()))
    return work, True


# ----------------------------------------------------------------------
# relax
# ----------------------------------------------------------------------

def relax(mesh, free=(), iterations=RELAX_ITERATIONS):
    """Smooth the interior, holding the boundary. ``free`` vertices are released.

    ``vertices_on_boundaries()`` -- PLURAL. The singular form returns only the
    longest boundary, so every hole would be left free and smoothing would round
    it off. ``free`` exists for the vertices a newly added strip has to separate:
    they sit on the boundary and must be allowed to move apart, and they are the
    only boundary vertices that may.
    """
    free = set(free)
    fixed = [w for w in boundary_vertices(mesh) if w not in free]
    mesh_smooth_centroid(mesh, kmax=iterations, fixed=fixed)
    return mesh


# ----------------------------------------------------------------------
# run
# ----------------------------------------------------------------------

mesh, source_guid, source_layer = load_quad_mesh()
print("editing the mesh on '{}'".format(source_layer))
walls = load_walls()

bad = non_quad_faces(mesh)
strips_ok = not bad

print("quad mesh: {} faces, {} vertices, {} edges".format(
    mesh.number_of_faces(), mesh.number_of_vertices(), mesh.number_of_edges()))
if bad:
    print("  {} face(s) are NOT quads -- pole fans, most likely. add_line and "
          "remove_line are disabled: collect_strip walks face_opposite_edge and a "
          "triangle has none, so a strip through a pole is not defined. "
          "move_vertex still works.".format(len(bad)))
if not walls:
    print("  no domain walls loaded -- a moved boundary vertex will not be held "
          "on the outline.")

if mesh.number_of_edges() > DRAW_WARN_EDGES:
    answer = (rs.GetString(
        "This mesh has {} edges. Drawing them all as pickable objects will make "
        "Rhino slow. Continue?".format(mesh.number_of_edges()),
        "Yes", ["Yes", "No"]) or "no").lower()
    if not answer.startswith("y"):
        print("Cancelled -- nothing was drawn or changed.")
        raise SystemExit

# This session's starting point, on its own layer, hidden. Refreshed every run
# rather than written once: bake_mesh clears only its own layer, so re-running
# CMD_quad_mesh leaves this one in place -- and a backup kept from before that
# regeneration would be a mesh of a different layout. Stale is worse than absent,
# so it always describes THIS session.
ensure_layer(UNEDITED_LAYER, (170, 170, 170))
bake_mesh(mesh, UNEDITED_LAYER)
rs.LayerVisible(UNEDITED_LAYER, False)
print("this session's starting mesh kept on '{}'".format(UNEDITED_LAYER))

snapshot = mesh.copy()
dirty = False

# The source object is in the way while picking -- it still shows the vertices
# where they were. Hidden rather than deleted, so a cancelled session leaves the
# document as it was found. Kept SEPARATE from whatever save produces: this list
# exists to be shown again, that one does not.
source_guids = [source_guid]
rs.HideObjects(source_guids)

draw(mesh)


def save(mesh):
    """Replace the mesh on ``QuadMesh``. Returns its guids, or ``[]`` on failure.

    **Whether this succeeded decides whether the command ends**, so it is
    reported rather than assumed. ``bake_mesh`` raises rather than returning
    nothing when Rhino refuses the mesh -- ``rs.AddMesh`` gives ``None`` and
    ``rs.ObjectLayer(None, ...)`` then raises -- so the failure is caught here and
    turned into an empty list.

    Note the order inside ``bake_mesh``: it clears the layer BEFORE adding, so a
    failed save has already removed the old mesh. That is recoverable rather than
    fatal only because this session's starting mesh is on ``Unedited``, which is
    what the message points at.
    """
    layer = ensure_layer(EDITED_LAYER, (0, 120, 200))
    # no clear_layer here -- bake_mesh already clears the layer before adding
    try:
        bake_mesh(mesh, layer)
    except Exception as e:
        rs.MessageBox(
            "The mesh was NOT saved.\n\n{}\n\nThe edit is still in memory, so try "
            "again. If '{}' is now empty, this session's starting mesh is on "
            "'{}'.".format(e, layer, UNEDITED_LAYER),
            0, "Save quad mesh")
        return []

    guids = rs.ObjectsByLayer(layer) or []
    if not guids:
        rs.MessageBox(
            "The mesh was NOT saved -- nothing landed on '{}'.\n\nThe edit is "
            "still in memory, so try again. This session's starting mesh is on "
            "'{}'.".format(layer, UNEDITED_LAYER),
            0, "Save quad mesh")
        return []

    print("Saved to '{}': {} faces, {} vertices.".format(
        layer, mesh.number_of_faces(), mesh.number_of_vertices()))
    return guids


try:
    while True:
        if strips_ok:
            options = ["move_vertex", "add_line", "remove_line", "relax",
                       "reset", "save", "exit"]
        else:
            options = ["move_vertex", "relax", "reset", "save", "exit"]

        operation = (rs.GetString("next", "exit", options) or "exit").lower()

        if operation == "move_vertex":
            if move_vertices(mesh, walls):
                dirty = True

        elif operation == "add_line" and strips_ok:
            mesh, changed = add_line(mesh)
            if changed:
                dirty = True
                draw(mesh)

        elif operation == "remove_line" and strips_ok:
            mesh, changed = remove_line(mesh)
            if changed:
                dirty = True
                draw(mesh)

        elif operation == "relax":
            value = rs.GetInteger("Smoothing passes", RELAX_ITERATIONS, 1)
            if value:
                relax(mesh, iterations=value)
                dirty = True
                draw(mesh)

        elif operation == "reset":
            # The whole MESH, not just the coordinates: add_line and remove_line
            # change the topology, and putting coordinates back would leave the
            # new faces in place with nothing to say where they came from.
            mesh = snapshot.copy()
            dirty = False
            draw(mesh)
            print("Back to this session's starting mesh.")

        elif operation == "save":
            guids = save(mesh)
            if guids:
                dirty = False
                print("  CMD_quad_mesh regenerates from the coarse layout and "
                      "would discard this. Edit the coarse layout instead if you "
                      "need to re-densify.")
                print("  the saved mesh is what CMD_smoothen and CMD_dual now "
                      "read.")
                break
            # A failed save leaves the edit in memory, so stay in the loop: ending
            # the command here would drop the work AND leave the layer in whatever
            # state the failed bake left it.

        else:
            if dirty:
                answer = (rs.GetString("Unsaved edits. Save before leaving?",
                                       "Yes", ["Yes", "No"]) or "yes").lower()
                if answer.startswith("y"):
                    if not save(mesh):
                        # Do not swallow it: the user asked to save on the way out
                        # and it did not happen.
                        print("Leaving WITHOUT saving -- the save failed above.")
                else:
                    print("Left unsaved -- '{}' is unchanged.".format(source_layer))
            break
finally:
    clear_scratch()
    if rs.IsLayer(EDIT_LAYER):
        clear_layer(EDIT_LAYER)
    rs.ShowObjects([guid for guid in source_guids if rs.IsObject(guid)])
    sc.doc.Views.Redraw()
