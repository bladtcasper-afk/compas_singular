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

Four things in here are not obvious and were each a bug first:

* **redraw everything, never update in place.** Every topological edit
  renumbers keys, and a stale guid map is a real way to move the wrong vertex.
  A coarse layout is tens of patches, so there is nothing to gain by being
  clever;
* **the drag preview must show the PROJECTED point, not the cursor.** A
  boundary vertex is pulled back onto its wall, so previewing the cursor draws
  a position the commit will not produce;
* **a refusal goes in a dialog, not a print.** Rhino's command line shows one
  line and the next prompt takes it, so a printed refusal reads to the user as
  the command silently doing nothing;
* **every ``rs.Add*`` here RAISES on geometry Rhino will not take** -- none of
  them returns ``None`` -- and one refused FACE invalidates a whole mesh, so a
  band of good quads goes with it. That is what hid every strip with a pole:
  the triangle in it was padded into a quad with two corners on the same point.
  What Rhino refuses is measured in :func:`_append_face`; and a swallowed
  refusal is worse than a visible one, because an object that was never drawn
  cannot be picked and nothing on screen says why.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import struct


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
    'density_colors',
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
    'strip': (0, 160, 190),
    'strip.highlight': (255, 120, 0),
    'path': (255, 120, 0),
    'ortho': (255,100,100),
    'diagonal': (50,255,50),
    'fan': (50,50,255),
}

#: The two ends of the density colour scale: the lowest density in the layout is
#: drawn LIGHT, the highest DARK. The ends of matplotlib's 'Blues'.
DENSITY_COLOR_LOW = (198, 219, 239)
DENSITY_COLOR_HIGH = (8, 48, 107)


def density_colors(strip_densities, low=DENSITY_COLOR_LOW, high=DENSITY_COLOR_HIGH):
    """``{skey: (r, g, b)}``, a blue darkening with density, RELATIVE to this layout.

    The lightest shade goes to the smallest density present and the darkest to
    the largest, so the scale always uses its full range. When every strip has
    the same density there is no range to spread, and all are drawn mid-scale
    rather than all-light, which would read as "all at the minimum".
    """
    if not strip_densities:
        return {}
    lo = min(strip_densities.values())
    hi = max(strip_densities.values())
    out = {}
    for skey, d in strip_densities.items():
        t = 0.5 if hi == lo else (d - lo) / float(hi - lo)
        out[skey] = tuple(int(round(low[i] + (high[i] - low[i]) * t)) for i in range(3))
    return out

#: How much of a pseudo-quad face is drawn, toward its own centroid. Big enough
#: to click on, small enough not to swallow the faces beside it.
POLE_FACE_SCALE = 0.6


def _shrunk_face(points, scale):
    """``points`` pulled toward their own centroid. As many as came in."""
    n = float(len(points))
    cx = sum(p[0] for p in points) / n
    cy = sum(p[1] for p in points) / n
    cz = sum(p[2] for p in points) / n
    return [[cx + (p[0] - cx) * scale,
             cy + (p[1] - cy) * scale,
             cz + (p[2] - cz) * scale] for p in points]


def _f32(point):
    """``point`` rounded the way a Rhino mesh vertex stores it.

    A mesh vertex is a ``Point3f``, so two corners 1e-9 apart are the SAME
    corner to Rhino -- and coincidence is what it refuses. Measured: 1e-9 apart
    is refused, 1e-7 apart is not.
    """
    return struct.unpack('<3f', struct.pack('<3f', point[0], point[1], point[2]))


def _append_face(vertices, points):
    """Append ``points`` as fresh corners of ``vertices``; return their faces.

    **A face may not have two coincident corners, and one that does invalidates
    the WHOLE mesh.** ``doc.Objects.AddMesh`` then returns an empty guid and
    ``rs.AddMesh`` RAISES. Measured against openNURBS 8: a quad whose corners 2
    and 3 sit at the same point is refused, so is one pinched across either
    diagonal, and a single such face carries a mesh of good quads down with it.
    A mesh whose VERTEX LIST merely repeats a point is fine -- which is why a
    ribbon spanning several faces, where consecutive faces put two corners on
    the same rung, has always drawn.

    So the corners are collapsed at float32 first, and what is left is spelled
    the way Rhino spells it: a triangle is a quad with its last index REPEATED
    (``rs.AddMesh``'s own docstring documents exactly that), never a fourth
    corner placed on top of the third. Returns ``[]`` for what is no longer a
    face at all.
    """
    kept = []
    for point in points:
        key = _f32(point)
        if any(key == other for other, _corner in kept):
            continue
        kept.append((key, list(point)))
    if len(kept) < 3:
        return []
    i = len(vertices)
    vertices += [corner for _key, corner in kept]
    if len(kept) == 3:
        return [[i, i + 1, i + 2, i + 2]]
    if len(kept) == 4:
        return [[i, i + 1, i + 2, i + 3]]
    # An n-gon has no business in a quad layout, but keeping only its first four
    # corners is how ``bake_mesh`` used to lose a fifth of the dual. Fan it.
    return [[i, i + k, i + k + 1, i + k + 1] for k in range(1, len(kept) - 1)]


def _lerp(a, b, t):
    """The point ``t`` of the way from ``a`` to ``b``."""
    return [a[0] + (b[0] - a[0]) * t,
            a[1] + (b[1] - a[1]) * t,
            a[2] + (b[2] - a[2]) * t]


class _Ribbon(object):
    """One strip's ribbon as ONE welded mesh: neighbouring faces share corners.

    Built face by face, but a corner is never appended twice. That is the whole
    difference from :func:`_append_face`, which gives every face fresh corners --
    correct, but Rhino then shades each quad on its own and the band reads as a
    row of tiles with seams between them.

    Two things make the sharing exact rather than approximate:

    * **a rung point is computed the same way from both sides.** The next face
      walks the shared rung the other way round, and ``lerp(a, b, t)`` and
      ``lerp(b, a, 1 - t)`` are the same point but not always the same bits. So
      it is always computed from the lower vertex key;
    * **corners are merged at float32, across the whole ribbon.** Rhino stores
      ``Point3f``, and two corners equal at float32 are the same corner to it --
      one face with two of them invalidates the entire mesh (see
      :func:`_append_face`). Merging per ribbon instead of per face is what lets a
      collapsed rung merge with the face on its other side too.
    """

    def __init__(self, mesh, half):
        self.mesh = mesh
        self.half = half
        self.vertices = []
        self.faces = []
        self._index = {}                # float32 point -> vertex index

    def _add(self, point):
        key = _f32(point)
        if key not in self._index:
            self._index[key] = len(self.vertices)
            self.vertices.append(list(point))
        return self._index[key]

    def rung_point(self, a, b, t):
        """The point ``t`` of the way from vertex ``a`` to ``b``, as an index."""
        if b < a:
            a, b, t = b, a, 1.0 - t
        return self._add(_lerp(self.mesh.vertex_coordinates(a),
                               self.mesh.vertex_coordinates(b), t))

    def corner(self, vkey):
        """A mesh vertex itself, as an index -- the pole a ribbon tapers into."""
        return self._add(self.mesh.vertex_coordinates(vkey))

    def add_face(self, indices):
        """Add a face; dropped if fewer than three distinct corners are left.

        Repeats are removed wherever they are, not only when consecutive: a quad
        pinched across a diagonal is refused by Rhino just as a collapsed side is.
        A triangle is spelled with its last index repeated, as ``rs.AddMesh``
        documents.
        """
        kept = []
        for index in indices:
            if index not in kept:
                kept.append(index)
        if len(kept) == 3:
            self.faces.append(kept + [kept[2]])
        elif len(kept) == 4:
            self.faces.append(kept)


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
        self.guid_faces = {}
        self.guid_strips = {}
        self.guid_path = []
        self._path_colors = {}          # guid -> colour before show_path
        self.guid_ftext = {}
        self.guid_stext = {}            # strip label dot -> skey

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

        chorded, lost = [], []

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
                guid = None
                if shape is not None and len(shape) > 2:
                    # ``rs.AddPolyline`` RAISES on points Rhino will not take --
                    # it does NOT return None -- and it judges coincidence at the
                    # DOCUMENT tolerance, not at 1e-9, so a drawn arc that
                    # doubles back on itself by half a tolerance is geometrically
                    # fine and still unbakeable (see helpers.clean_polyline_points,
                    # where the same refusal once cost a whole commit). The edge
                    # falls back to its chord rather than going missing: that is
                    # what the densification does with it anyway.
                    try:
                        guid = rs.AddPolyline([Point3d(*point) for point in shape])
                    except Exception:                             # noqa: BLE001
                        chorded.append((u, v))
                if not guid:
                    try:
                        guid = rs.AddLine(Point3d(*mesh.vertex_coordinates(u)),
                                          Point3d(*mesh.vertex_coordinates(v)))
                    except Exception:                             # noqa: BLE001
                        # A zero-length chord. Nothing can be drawn for it, and
                        # ``ObjectLayer(None)`` would raise out of the redraw.
                        guid = None
                if not guid:
                    lost.append((u, v))
                    continue
                rs.ObjectLayer(guid, self.layer)
                rs.ObjectColor(guid, self.colors['edge'])
                self.guid_edges[guid] = (u, v)
        finally:
            rs.EnableRedraw(True)
        sc.doc.Views.Redraw()
        if chorded:
            print("{} edge(s) drawn as a chord -- Rhino refused the curve."
                  .format(len(chorded)))
        if lost:
            print("{} edge(s) could not be drawn, and cannot be picked: {}."
                  .format(len(lost), ", ".join("{}-{}".format(*e) for e in lost)))

    def update_face(self, mesh, fkey):
        _require_rhino()

        guid = None
        for d_guid, d_fkey in self.guid_faces.items():
            if d_fkey == fkey:
                guid = d_guid
                break

        # In theory the order of this one is always the same in the current flows
        # Thus, finding the text_guid could be simplified here.
        text_guid = None
        for d_guid, d_fkey in self.guid_ftext.items():
            if d_fkey == fkey:
                text_guid = d_guid
                break
        print(text_guid)
        if guid:
            pattern = mesh.get_face_pattern(fkey)
            rs.ObjectColor(guid, self.colors[pattern])

            if text_guid:
                rs.TextObjectText(text_guid, text=pattern)

    def draw_faces(self, mesh, highlight=(), highlight_pattern=False, clear=True):
        _require_rhino()
        if clear:
            self._clear_faces()
        ensure_layer(self.layer)

        work = mesh.copy()
        keys = work.faces()
        highlight = set(highlight)
        missing = []

        rs.EnableRedraw(False)
        try:
            for fkey in keys:
                vkeys = work.face_vertices(fkey)
                vcoordinates = [work.vertex_attributes(vkey, "xyz") for vkey in vkeys]

                #check for triangular faces to dupe point for rhino
                if len(vcoordinates) == 3:
                    face = [0,1,2,2]
                else:
                    face = [0,1,2,3]

                try:
                    guid = rs.AddMesh(vcoordinates, [face])
                except Exception as error:                        # noqa: BLE001
                    print("strip {}: Rhino refused the band -- {}".format(fkey, error))
                    guid = None
                if not guid:
                    missing.append(fkey)
                    continue
                rs.ObjectLayer(guid, self.layer)
                if highlight_pattern:
                    pattern = mesh.get_face_pattern(fkey)
                    rs.ObjectColor(guid, self.colors[pattern])
                    centroid = rs.MeshAreaCentroid(guid)
                    text_height = rs.MeshArea(guid)[1]/100
                    text_guid = rs.AddText(pattern, point_or_plane=centroid, height=text_height, justification=131074)
                    rs.ObjectLayer(text_guid, self.layer)
                    rs.ObjectColor(text_guid, (0,0,0))
                    self.guid_ftext[text_guid] = fkey
                else:
                    rs.ObjectColor(guid, self.colors[
                        'strip.highlight' if fkey in highlight else 'strip'])
                self.guid_faces[guid] = fkey
        finally:
            rs.EnableRedraw(True)
        sc.doc.Views.Redraw()
        if missing:
            print("{} of {} strips could not be drawn, and cannot be picked: {}."
                  .format(len(missing), len(keys),
                          ", ".join(str(skey) for skey in missing)))
        return dict(self.guid_strips)

    def draw_strips(self, mesh, strips=None, width=0.5, highlight=(), clear=True,
                    collect=True, strip_colors=None, labels=None):
        """**Draw strips as pickable ribbons.** Returns ``{guid: skey}``.

        A strip is a band of faces, and the obvious way to draw one -- fill its
        faces -- cannot be picked: **every face belongs to TWO strips**, so the
        filled bands overlap everywhere and whichever was drawn last takes every
        click. Instead each band is drawn as a RIBBON down its own middle: within
        each face, the quad between its two rungs, narrowed to ``width`` of them.
        Two crossing strips then overlap only in a small square at the centre of
        the face they share, and the rest of each ribbon is uniquely clickable.

        Parameters
        ----------
        mesh : :class:`QuadMesh`
        strips : iterable, optional
            Which strips to draw. Default is all of them. Pass a short list to
            show only the candidates a command is offering.
        width : float, optional
            The ribbon width as a fraction of the rung it spans. 0.5 by default
            -- wide enough to hit, narrow enough to leave the crossing strip
            visible underneath it.
        highlight : iterable, optional
            Strips to draw in the highlight colour, for showing which band a
            drawn curve has resolved to before the user commits to it.
        clear : bool, optional
            Remove previously drawn strips first. Off when overlaying.
        collect : bool, optional
            Re-collect the strips on a copy before drawing -- the default, and
            right for an editor whose mesh just changed. Pass False to draw the
            strips the mesh ALREADY carries, under the keys it carries them: a
            caller holding data keyed by those strips (``strips_density``) needs
            the picked key to be its key, and a re-collection is free to number
            the same bands differently. Collected anyway if the mesh has none.
        strip_colors : dict, optional
            ``{skey: (r, g, b)}``, overriding the strip colour per strip. A
            highlighted strip still takes the highlight colour.
        labels : dict, optional
            ``{skey: text}``, drawn as a text dot on the middle rung of each
            strip. A rung belongs to exactly one strip, so unlike a face centre
            it is never shared with the crossing band's label.

        Notes
        -----
        **Each ribbon is ONE welded mesh** (:class:`_Ribbon`): consecutive faces
        share the corners on their common rung, so the band shades as one
        continuous patch rather than a row of tiles, and a closed strip closes.

        **A pseudo-quad is part of the ribbon, not skipped.** Where the strip
        passes through it between the two pole-edges it is an ordinary ribbon
        face with one side at the pole; where the strip ENDS at the pole, the
        ribbon tapers into it as a triangle. Only a strip no ribbon face can be
        built for falls back to its faces drawn shrunk. Skipping pseudo-quads
        was the obvious thing and exactly wrong: measured on a plate with a circular
        hole and four poles, of the 12 strips that DO get a ribbon **none could
        be deleted**, while all 8 single-pseudo-quad strips could. Drawing only
        the ribbons showed the user every strip they could not act on and hid the
        only ones they could.

        A pseudo-quad is a TRIANGLE, and every face here goes through
        :func:`_append_face` for that reason: a triangle padded with a fourth
        corner on top of its third is a face Rhino refuses, and it takes the
        whole strip mesh -- ribbons included -- down with it. That is why the
        strips with a pole used to draw nothing at all.
        """
        _require_rhino()
        if clear:
            self._clear_strips()
        ensure_layer(self.layer)

        if collect or not mesh.attributes.get('strips'):
            work = mesh.copy()
            work.collect_strips()
        else:
            work = mesh            # read only below
        keys = list(work.attributes['strips']) if strips is None else list(strips)
        highlight = set(highlight)
        strip_colors = strip_colors or {}
        labels = labels or {}
        half = max(0.05, min(0.49, float(width) / 2.0))
        missing = []

        rs.EnableRedraw(False)
        try:
            for skey in keys:
                if skey not in work.attributes['strips']:
                    continue
                rungs = set(frozenset(edge) for edge in work.strip_edges(skey))
                ribbon = _Ribbon(work, half)
                for fkey in work.strip_faces(skey):
                    here = [(u, v) for u, v in work.face_halfedges(fkey)
                            if u != v and frozenset((u, v)) in rungs]
                    opposite = work.face_opposite_edge(*here[0]) if here else None
                    if opposite is None:
                        continue
                    u, v = here[0]
                    w, x = opposite
                    if w == x:
                        # The far side of a pseudo-quad collapsed to its pole:
                        # the ribbon tapers into the pole instead of stopping.
                        ribbon.add_face([ribbon.rung_point(u, v, 0.5 - half),
                                         ribbon.rung_point(u, v, 0.5 + half),
                                         ribbon.corner(w)])
                    else:
                        # (u, v, w, x) go round the face, so u pairs with x and
                        # v with w -- the ribbon runs u-side, v-side, w-side,
                        # x-side.
                        ribbon.add_face([ribbon.rung_point(u, v, 0.5 - half),
                                         ribbon.rung_point(u, v, 0.5 + half),
                                         ribbon.rung_point(x, w, 0.5 + half),
                                         ribbon.rung_point(x, w, 0.5 - half)])
                vertices, faces = ribbon.vertices, ribbon.faces
                if not faces:
                    # Nothing to run a ribbon across -- a strip no face could be
                    # walked for. Draw its faces shrunk instead, so it can at
                    # least still be picked (see Notes).
                    for fkey in work.strip_faces(skey):
                        faces += _append_face(vertices, _shrunk_face(
                            [work.vertex_coordinates(c) for c in work.face_vertices(fkey)],
                            POLE_FACE_SCALE))
                if not faces:
                    missing.append(skey)
                    continue
                # ``rs.AddMesh`` is safe HERE and nowhere else in this
                # workflow: it keeps only the first four corners of a face, and
                # every face built above has exactly four indices.
                #
                # It RAISES rather than returning None when Rhino refuses the
                # mesh, so the call is guarded -- but a strip that goes missing
                # is a strip the user can see no way to select, so the guard
                # SAYS SO instead of dropping it in silence.
                try:
                    guid = rs.AddMesh(vertices, faces)
                except Exception as error:                        # noqa: BLE001
                    print("strip {}: Rhino refused the band -- {}".format(skey, error))
                    guid = None
                if not guid:
                    missing.append(skey)
                    continue
                rs.ObjectLayer(guid, self.layer)
                if skey in highlight:
                    color = self.colors['strip.highlight']
                else:
                    color = strip_colors.get(skey, self.colors['strip'])
                rs.ObjectColor(guid, color)
                self.guid_strips[guid] = skey

                if skey in labels:
                    rungs = [(u, v) for u, v in work.strip_edges(skey) if u != v]
                    if rungs:
                        u, v = rungs[len(rungs) // 2]
                        mid = _lerp(work.vertex_coordinates(u),
                                    work.vertex_coordinates(v), 0.5)
                        dot = rs.AddTextDot(str(labels[skey]), Point3d(*mid))
                        if dot:
                            rs.ObjectLayer(dot, self.layer)
                            self.guid_stext[dot] = skey
        finally:
            rs.EnableRedraw(True)
        sc.doc.Views.Redraw()
        if missing:
            print("{} of {} strips could not be drawn, and cannot be picked: {}."
                  .format(len(missing), len(keys),
                          ", ".join(str(skey) for skey in missing)))
        return dict(self.guid_strips)

    def show_path(self, vkeys, mesh=None, color=None):
        """**Show the run picked so far, so the user can see what they are building.**

        Picking a polyedge corner by corner is otherwise blind: the command line
        lists the keys, but nothing on screen says which way the line is going or
        whether the last pick was the one intended. Call it after every pick; it
        replaces whatever it showed before.

        **It RECOLOURS this scene's own edge and corner objects; it does not draw
        a line on top of them.** It used to add a polyline through the picks, and
        in Rhino nothing appeared: consecutive corners of a polyedge are joined by
        a coarse edge, so that polyline lay exactly on top of the grey edge
        objects at the same width, and Rhino does not promise which of two
        coincident curves it draws last. Recolouring the edge itself cannot lose
        that race, and it lights up the FIRST pick too, where a line needs two.

        **Not selection either.** Selecting the objects would highlight them just
        as well, but ``pick_vertex`` asks with ``preselect=True``, so the next
        prompt would take the highlighted corner as its answer before the user
        clicked anything.

        Parameters
        ----------
        vkeys : list[int]
            The corners picked so far, in order.
        mesh : Mesh, optional
            Only needed for a pair of corners with no edge object between them
            (not joined, or the edge could not be drawn); that pair gets a chord
            polyline so the run still reads as one line. Without it the pair is
            skipped.
        """
        _require_rhino()
        self.clear_path()
        color = color or self.colors['path']
        vkeys = list(vkeys)

        guid_of_vertex = {vkey: guid for guid, vkey in self.guid_vertices.items()}
        guid_of_edge = {frozenset(edge): guid for guid, edge in self.guid_edges.items()}

        targets = [guid_of_vertex[vkey] for vkey in vkeys if vkey in guid_of_vertex]
        chords = []
        for u, v in zip(vkeys[:-1], vkeys[1:]):
            guid = guid_of_edge.get(frozenset((u, v)))
            if guid is not None:
                targets.append(guid)
            elif mesh is not None:
                chords.append((mesh.vertex_coordinates(u), mesh.vertex_coordinates(v)))

        for guid in targets:
            if guid in self._path_colors or not rs.IsObject(guid):
                continue
            self._path_colors[guid] = rs.ObjectColor(guid)
            rs.ObjectColor(guid, color)

        if chords:
            ensure_layer(self.layer)
            for a, b in chords:
                try:
                    guid = rs.AddPolyline([Point3d(*a), Point3d(*b)])
                except Exception:                                 # noqa: BLE001
                    guid = None
                if not guid:
                    continue
                rs.ObjectLayer(guid, self.layer)
                rs.ObjectColor(guid, color)
                self.guid_path.append(guid)

        sc.doc.Views.Redraw()
        return list(self._path_colors) + list(self.guid_path)

    def clear_path(self):
        """Put back the colours :meth:`show_path` changed. Safe to call twice."""
        for guid, original in self._path_colors.items():
            if rs.IsObject(guid):
                rs.ObjectColor(guid, original)
        self._path_colors = {}
        guids = [guid for guid in self.guid_path if rs.IsObject(guid)]
        if guids:
            rs.DeleteObjects(guids)
        self.guid_path = []

    def _clear_faces(self):
        guids = [guid for guid in self.guid_faces if rs.IsObject(guid)]
        if guids:
            rs.DeleteObjects(guids)
        self.guid_faces = {}
    
    def _clear_strips(self):
        guids = [guid for guid in list(self.guid_strips) + list(self.guid_stext)
                 if rs.IsObject(guid)]
        if guids:
            rs.DeleteObjects(guids)
        self.guid_strips = {}
        self.guid_stext = {}

    def _clear_text(self):
        guids = [guid for guid in self.guid_ftext if rs.IsObject(guid)]
        if guids:
            rs.DeleteObjects(guids)
        self.guid_ftext = {}

    def clear(self):
        """Delete everything this object drew. Safe to call twice."""
        guids = [guid for guid in (list(self.guid_vertices)
                                   + list(self.guid_edges)
                                   + list(self.guid_strips)
                                   + list(self.guid_stext)
                                   + list(self.guid_path)
                                   + list(self.guid_faces)
                                   + list(self.guid_ftext))
                 if rs.IsObject(guid)]
        if guids:
            rs.DeleteObjects(guids)
        self.guid_vertices = {}
        self.guid_edges = {}
        self.guid_strips = {}
        self.guid_stext = {}
        self.guid_path = {}
        self._path_colors = {}

    # -- locking ---------------------------------------------------------

    def unlock(self):
        """Unlock the scratch layer and its objects. Pass the result to :meth:`relock`."""
        return unlock(self.layer,
                      list(self.guid_vertices) + list(self.guid_edges)
                      + list(self.guid_strips) + list(self.guid_path))

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

    def pick_face(self, message='Select a face', preselect=True):
        def face_filter(rh_object, geometry, component_index):
            fkey = self.guid_faces.get(rh_object.id)
            return fkey is not None
   
        while True:
            guid = rs.GetObject(message, rs.filter.mesh, preselect=preselect, custom_filter=face_filter)
            if not guid:
                return None
            fkey = self.guid_faces.get(guid)
            if fkey is not None:
                return fkey
            print("Not a face of this layout -- pick one of the bands on "
                  "{!r}.".format(self.layer))

    def pick_faces(self, message='Select faces', preselect=True):
        def face_filter(rh_object, geometry, component_index):
            fkey = self.guid_faces.get(rh_object.id)
            return fkey is not None
   
        while True:
            guids = rs.GetObjects(message, rs.filter.mesh, preselect=preselect, custom_filter=face_filter)
            if not guids:
                return None

            fkeys = []
            for guid in guids:
                fkey = self.guid_faces.get(guid)
                fkeys.append(fkey)
            if fkeys != []:
                return fkeys
            print("Not a face of this layout -- pick one of the bands on "
                  "{!r}.".format(self.layer))        

    def pick_strip(self, message='Select a strip', preselect=True):
        """A strip key, or ``None`` if the user pressed Esc.

        Re-prompts on anything that is not one of ours, for the same reason
        :meth:`pick_edge` does: clicking slightly off is a miss, not a decision
        to stop.
        """
        while True:
            guid = rs.GetObject(message, rs.filter.mesh, preselect=preselect)
            if not guid:
                return None
            skey = self.guid_strips.get(guid)
            if skey is not None:
                return skey
            print("Not a strip of this layout -- pick one of the bands on "
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
