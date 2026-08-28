"""**The Rhino side of editing a mesh by hand: draw it, pick it, drag it.**

Everything here is interaction, and nothing here decides anything. It puts a
mesh on screen as pickable objects, turns a pick back into a vertex or an edge
key, drags a point with a live preview, and asks a question. What the answer
*means* belongs to :mod:`compas_singular.editing`, which is why that package
imports none of this and can run with no Rhino at all.

It lives in ``src`` rather than in a command because **two commands need the
same three hundred lines**: ``CMD_edit_coarse_mesh`` drives the coarse layout
and ``CMD_edit_quad_mesh`` the dense mesh, and both draw points and lines onto
a scratch layer, hold a guid map, pick against it and drag with a projected
preview. Anything only one command needs -- drawing an arc across the layout,
say -- stays in that command.

Three things in here are not obvious and were each a bug first:

* **redraw everything, never update in place.** Every topological edit
  renumbers keys, and a stale guid map is a real way to move the wrong vertex.
  A coarse layout is tens of patches, so there is nothing to gain by being
  clever;
* **the drag preview must show the PROJECTED point, not the cursor.** A
  boundary vertex is pulled back onto its wall, so previewing the cursor draws
  a position the commit will not produce;
* **a refusal goes in a dialog, not a print.** Rhino's command line shows one
  line and the next prompt takes it, so a printed refusal reads to the user as
  the command silently doing nothing.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


try:
    import rhinoscriptsyntax as rs
    import scriptcontext as sc
    import Rhino
    from Rhino.Geometry import Line
    from Rhino.Geometry import Point3d
    # ``from System.Drawing.Color import FromArgb`` is an IronPython idiom and
    # raises under Rhino 8's CPython, where Color is a class and not a module.
    # It raised at the END of this block, so everything above it stayed bound
    # and the failure only surfaced later as ``NameError: FromArgb``.
    from System.Drawing import Color
    RHINO = True
except ImportError:  # importable outside Rhino, so the package still imports
    RHINO = False


__all__ = [
    'DEFAULT_COLORS',
    'PickableMesh',
    'ask',
    'ask_integer',
    'drag_point',
    'ensure_layer',
    'refuse',
    'relock',
    'unlock',
]


#: Colours every hand-edit command draws with. One scheme, so a boundary corner
#: looks the same in the coarse editor and the dense one.
DEFAULT_COLORS = {
    'vertex': (0, 0, 0),
    'vertex.boundary': (0, 120, 200),
    'vertex.special': (200, 0, 120),     # a pole, or a singularity
    'edge': (110, 110, 110),
    'preview': (255, 120, 0),
}


def _require_rhino():
    if not RHINO:
        raise RuntimeError(
            'compas_singular.rhino.mesh_ui needs Rhino. The editing operations '
            'themselves do not -- see compas_singular.editing.')


# ----------------------------------------------------------------------
# layers
# ----------------------------------------------------------------------

def ensure_layer(path, color=None):
    """Create a ``::`` layer path, parents first, and return the full path.

    ``rs.AddLayer`` does not create intermediate parents, and ``bake_mesh``'s
    own layer creation hard-codes ``parent="TopologyProblem"``, which is wrong
    for a layer two levels down. Only the leaf gets the colour.
    """
    _require_rhino()
    parts = path.split('::')
    for i in range(len(parts)):
        name = '::'.join(parts[:i + 1])
        if not rs.IsLayer(name):
            rs.AddLayer(name=parts[i],
                        parent='::'.join(parts[:i]) if i else None,
                        color=color if i == len(parts) - 1 else None)
    return path


def unlock(layer, guids=()):
    """Unlock a layer, every parent of it, and the given objects.

    Returns what was locked, to hand back to :func:`relock`. Layer lock and
    object lock are separate in Rhino and either one alone is enough to make
    the points unpickable, so both are cleared. A locked PARENT locks the
    child, which is why the whole path is walked.
    """
    _require_rhino()
    parts = layer.split('::')
    state = {'layers': {}, 'objects': []}
    for i in range(len(parts)):
        name = '::'.join(parts[:i + 1])
        if rs.IsLayer(name):
            state['layers'][name] = rs.LayerLocked(name)
            if state['layers'][name]:
                rs.LayerLocked(name, False)
    for guid in guids:
        if rs.IsObject(guid) and rs.IsObjectLocked(guid):
            state['objects'].append(guid)
            rs.UnlockObject(guid)
    return state


def relock(state):
    """Put back exactly what :func:`unlock` took off, and nothing else."""
    _require_rhino()
    for guid in state.get('objects', []):
        if rs.IsObject(guid):
            rs.LockObject(guid)
    for layer, was_locked in state.get('layers', {}).items():
        if was_locked and rs.IsLayer(layer):
            rs.LayerLocked(layer, True)


# ----------------------------------------------------------------------
# asking
# ----------------------------------------------------------------------

def ask(message, options, default=None):
    """A command-line menu. Returns the answer lowercased, or ``''`` on Esc.

    Lowercased because every caller compares against a lowercase name, and
    ``rs.GetString`` returns whatever case the option was declared in.
    """
    _require_rhino()
    answer = rs.GetString(message, default if default is not None else options[0],
                          list(options))
    return (answer or '').lower()


def ask_integer(message, default, minimum=None, maximum=None):
    """``rs.GetInteger``, returning ``None`` on Esc rather than a surprise."""
    _require_rhino()
    value = rs.GetInteger(message, default, minimum, maximum)
    return None if value is None else int(value)


def refuse(title, message):
    """Tell the user why nothing happened, in a DIALOG.

    Rhino's command line shows one line and the next prompt overwrites it, so a
    printed refusal reads as the command silently doing nothing -- which is the
    worst possible reading of a command that deliberately changed nothing.
    """
    _require_rhino()
    rs.MessageBox(message, 0, title)


# ----------------------------------------------------------------------
# dragging
# ----------------------------------------------------------------------

def drag_point(prompt, start, neighbours=(), project=None, colors=None):
    """Drag one point with a live preview. ``[x, y, 0.0]``, or ``None`` on Esc.

    Parameters
    ----------
    start : list
        Where the point is now. Used as the drag base point.
    neighbours : list, optional
        Points to rubber-band to, so the mesh is seen deforming rather than one
        dot moving.
    project : callable, optional
        ``xyz -> xyz``, applied to the cursor before anything is drawn. This is
        what makes the preview honest: a boundary vertex is shown where it will
        actually end up, not where the cursor is.

    Constrained to world XY because the pipeline is planar. An unconstrained
    pick lands on whatever construction plane the view happens to have, and the
    point moved would then differ from the one previewed.
    """
    _require_rhino()
    colors = colors or DEFAULT_COLORS
    color = Color.FromArgb(*colors['preview'])
    neighbours = [list(point) for point in neighbours]

    def _at(x, y):
        xyz = [x, y, 0.0]
        return project(xyz) if project is not None else xyz

    gp = Rhino.Input.Custom.GetPoint()
    gp.SetCommandPrompt(prompt)
    gp.SetBasePoint(Point3d(*start), True)
    gp.Constrain(Rhino.Geometry.Plane.WorldXY, False)

    def on_dynamic_draw(sender, e):
        cp = e.CurrentPoint
        point = Point3d(*_at(cp.X, cp.Y))
        for nbr in neighbours:
            e.Display.DrawLine(Line(point, Point3d(*nbr)), color, 2)
        e.Display.DrawPoint(point, color)

    gp.DynamicDraw += on_dynamic_draw
    try:
        gp.Get()
        if gp.CommandResult() != Rhino.Commands.Result.Success:
            return None
        ep = gp.Point()
    finally:
        # Detached explicitly: the handler closes over ``neighbours``, and a
        # GetPoint that outlives the call would keep drawing a stale mesh.
        gp.DynamicDraw -= on_dynamic_draw
    return _at(ep.X, ep.Y)


# ----------------------------------------------------------------------
# the mesh on screen
# ----------------------------------------------------------------------

class PickableMesh(object):
    """A mesh drawn as points and lines the user can actually click.

    Owns the two guid maps that are the entire pick model -- there is no Rhino
    user data on the objects, because a redraw throws them all away anyway.

    Parameters
    ----------
    layer : str
        Full ``::`` path of the scratch layer. Created if missing.
    colors : dict, optional
        Overrides for :data:`DEFAULT_COLORS`.

    Examples
    --------
    >>> scene = PickableMesh('TopologyProblem::Skeleton::TempEdit')
    >>> scene.draw(mesh, special=mesh.poles())
    >>> edge = scene.pick_edge('Pick an edge')
    >>> scene.clear()
    """

    def __init__(self, layer, colors=None):
        _require_rhino()
        self.layer = layer
        self.colors = dict(DEFAULT_COLORS, **(colors or {}))
        self.guid_vertices = {}
        self.guid_edges = {}

    # -- drawing ---------------------------------------------------------

    def draw(self, mesh, edge_shape=None, special=()):
        """Redraw the whole mesh. Everything, every time -- see the module doc.

        Parameters
        ----------
        edge_shape : callable, optional
            ``(u, v) -> point list or None``. An edge that carries a drawn
            SHAPE is drawn as that shape rather than as its chord, or a cut
            that worked looks like one that snapped straight. It is still one
            object and still picked the same way: a polyline is a curve to
            ``rs.GetObject``.
        special : iterable, optional
            Vertices to colour as poles or singularities -- usually what the
            user is trying to reach, or to avoid.
        """
        self.clear()
        ensure_layer(self.layer)

        boundary = set()
        for ring in mesh.vertices_on_boundaries():
            boundary.update(ring)
        special = set(special)

        rs.EnableRedraw(False)
        try:
            for vkey in mesh.vertices():
                guid = rs.AddPoint(Point3d(*mesh.vertex_coordinates(vkey)))
                if not guid:
                    continue
                rs.ObjectLayer(guid, self.layer)
                if vkey in special:
                    key = 'vertex.special'
                elif vkey in boundary:
                    key = 'vertex.boundary'
                else:
                    key = 'vertex'
                rs.ObjectColor(guid, self.colors[key])
                self.guid_vertices[guid] = vkey

            for u, v in mesh.edges():
                shape = edge_shape(u, v) if edge_shape is not None else None
                if shape is not None and len(shape) > 2:
                    guid = rs.AddPolyline([Point3d(*point) for point in shape])
                else:
                    guid = rs.AddLine(Point3d(*mesh.vertex_coordinates(u)),
                                      Point3d(*mesh.vertex_coordinates(v)))
                # AddLine/AddPolyline return None on geometry Rhino will not
                # take rather than raising, and ObjectLayer(None) does raise.
                if not guid:
                    continue
                rs.ObjectLayer(guid, self.layer)
                rs.ObjectColor(guid, self.colors['edge'])
                self.guid_edges[guid] = (u, v)
        finally:
            rs.EnableRedraw(True)
        sc.doc.Views.Redraw()

    def clear(self):
        """Delete everything this object drew. Safe to call twice."""
        guids = [guid for guid in list(self.guid_vertices) + list(self.guid_edges)
                 if rs.IsObject(guid)]
        if guids:
            rs.DeleteObjects(guids)
        self.guid_vertices = {}
        self.guid_edges = {}

    # -- locking ---------------------------------------------------------

    def unlock(self):
        """Unlock the scratch layer and its objects. Pass the result to :meth:`relock`."""
        return unlock(self.layer,
                      list(self.guid_vertices) + list(self.guid_edges))

    def relock(self, state):
        relock(state)

    # -- picking ---------------------------------------------------------

    def pick_vertex(self, message='Select a vertex', preselect=True):
        """A vertex key, or ``None`` if the user pressed Esc.

        Re-prompts on anything that is not ours, for the same reason
        :meth:`pick_edge` does: clicking slightly off is a miss, not a decision
        to leave, and ending the loop on one would be maddening.
        """
        while True:
            guid = rs.GetObject(message, rs.filter.point, preselect=preselect)
            if not guid:
                return None
            vkey = self.guid_vertices.get(guid)
            if vkey is not None:
                return vkey
            print("Not a vertex of this mesh -- pick one of the points on "
                  "{!r}.".format(self.layer))

    def pick_edge(self, message='Select an edge', preselect=True):
        """``((u, v), guid)``, or ``None`` if the user pressed Esc.

        The guid comes back as well as the key because picking a point ALONG
        the edge -- ``rs.GetPointOnCurve``, or constraining a drag to it --
        needs the Rhino object, and looking it up again by key would mean
        keeping a second map.

        Re-prompts rather than giving up: picking a curve that is not one of
        ours is a miss, not a decision to stop, and the alternative is a
        command that quietly ends when the user clicks slightly off.
        """
        while True:
            guid = rs.GetObject(message, rs.filter.curve, preselect=preselect)
            if not guid:
                return None
            edge = self.guid_edges.get(guid)
            if edge is not None:
                return edge, guid
            print("Not an edge of this mesh -- pick one of the lines on "
                  "{!r}.".format(self.layer))

    def closest_edge_point(self, point, tol=None):
        """``point`` snapped onto one of OUR edges if it lies on one, else ``None``.

        Only this mesh's edges are searched. Any other curve in the document --
        an input boundary, a baked separatrix, a construction line -- would
        anchor a path somewhere the mesh has no edge, and the operation would
        then refuse something that looked finished.
        """
        if tol is None:
            tol = sc.doc.ModelAbsoluteTolerance * 2.0
        for guid in self.guid_edges:
            curve = rs.coercecurve(guid)
            if not curve:
                continue
            success, t = curve.ClosestPoint(point)
            if not success:
                continue
            candidate = curve.PointAt(t)
            if candidate.DistanceTo(point) <= tol:
                return candidate
        return None
