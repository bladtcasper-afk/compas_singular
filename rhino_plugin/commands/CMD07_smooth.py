#! python3
# r: compas
# r: pydantic
# r: compas_fd
"""Smooth a dense quad mesh -- the whole of it, or a region -- from one option line.

One command for what ``CMD_smoothen`` (whole / region) and ``CMD_smoothen_guide`` (guide
curves) do separately. Both are left as they are. Pick the mesh, answer Whole or Region,
set the options and press Enter to smooth::

    Whole   Algorithm  Boundary  FixedVertices  Guides  Iterations  Damping
    Region  Region  Blend  Boundary  FixedVertices  Guides  Iterations  Damping

* **Algorithm** -- Area, Centroid, CenterOfMass or ForceDensity. A region is always
  smoothed by area: ``smoothing_region`` has no other rule, because centroid equalises
  edge lengths and fights the grading a mesh is meant to have.
* **Boundary** -- Sliding runs every boundary vertex along its own outline with the
  corners pinned; Fixed pins the whole boundary; Free holds none of it, corners included,
  so only FixedVertices and guides anchor the mesh. Free is what ForceDensity form finding
  wants: fix the supports and let the outline find its own shape.
* **FixedVertices** -- vertices that do not move at all: pick them, or take the vertex at
  every point on the PointFeatures layer. A fixed vertex wins over everything else,
  including a guide chain running through it.
* **Guides** -- pick a guide curve and the chain along it is proposed: the longest run of
  ONE polyedge that follows it (see ``compas_singular.editing.guide_chain`` for why it is
  chosen rather than built). Edit it with Add / Remove / Clear and set how it holds.
  AllGuides proposes a chain for every curve on the Guides layer at once. Picking a guide
  that already has a chain edits that chain. A guide may be any curve: the chain is chosen
  on its samples, and a curved guide's vertices land ON the curve, not on a chord of them.

Nothing is attached until Enter. Every chain is chosen on the mesh as it was picked, and
all of them are moved onto their guides together, so the result does not depend on the
order the guides were picked in.

**ForceDensity cannot hold a vertex on a curve.** ``relaxation`` takes a fixed set and
nothing else, so under it Boundary=Sliding holds only the corners -- the rest of the
outline is held by stiff edges alone and can drift off the wall -- and a guide vertex is
moved onto its guide and pinned there, whatever its hold says. The command prints what it
could not honour.

The smoothing itself is in ``compas_singular.datastructures.mesh.smoothing`` and the chain
selection in ``compas_singular.editing.guide_chain``. This file is the picks, the option
lines, the preview and the bake -- plus :func:`smooth_whole` / :func:`smooth_region`, which
put the options together and take no Rhino input, so they can be run without Rhino.
"""

from compas_singular.rhino.project import LAYER_DATA

import Rhino
import System
import rhinoscriptsyntax as rs
import scriptcontext as sc
import compas_rhino as cr
from compas.geometry import Point
from compas_rhino.conversions import curve_to_compas
from compas_singular.rhino.helpers import mesh_from_rhino

from compas_singular.datastructures import QuadMesh
from compas_singular.datastructures.mesh.smoothing import automated_boundary_constraints
from compas_singular.datastructures.mesh.smoothing import constrained_smoothing
from compas_singular.datastructures.mesh.smoothing import mesh_boundary_corners
from compas_singular.datastructures.mesh.smoothing import relaxation
from compas_singular.datastructures.mesh.smoothing import smoothing_region
from compas_singular.editing import GuideCurve
from compas_singular.editing import attach_chain
from compas_singular.editing import chain_quality
from compas_singular.editing import collect_polyedges
from compas_singular.editing import guide_chain
from compas_singular.editing import mean_edge_length
from compas_singular.editing.guide_chain import DEFAULT_MAX_ANGLE
from compas_singular.editing.guide_chain import DEFAULT_TOLERANCE_FACTOR
from compas_singular.rhino.helpers import bake_mesh, clear_layer, curve_points

ALGORITHMS = ["Area", "Centroid", "CenterOfMass", "ForceDensity"]
BOUNDARY_MODES = ["Sliding", "Fixed", "Free"]

#: Where a smoothed region is baked, next to the per-algorithm layers of a whole-mesh
#: smooth. The same name CMD_smoothen bakes a region to.
REGION_OUTPUT_LAYER = "Relaxed"

PREVIEW_LAYER = "SmoothPreview"
CORE_COLOUR = (255, 0, 0)
BLEND_COLOUR = (255, 160, 0)
GUIDE_COLOUR = (200, 0, 255)
FIXED_COLOUR = (0, 110, 255)

#: The full paths, from the one place the layer map is written down. rhinoscriptsyntax
#: resolves a bare name with FindName, which returns the FIRST layer of that name anywhere
#: in the document.
GUIDE_LAYER = LAYER_DATA.get("Guides", ("Guides", None))[0]
POINT_LAYER = LAYER_DATA.get("PointFeatures", ("PointFeatures", None))[0]

#: Rings of vertices outside a region over which the damping falls to zero. A region
#: smoothed at full strength up to a hard edge leaves a CREASE there -- a fixed ring next
#: to a fully relaxed one is a kink, the sort of defect smoothing is meant to remove.
DEFAULT_BLEND = 3

#: How much heavier a boundary edge is than an interior one under ForceDensity. The value
#: CMD_smoothen runs it with.
DEFAULT_Q_FACTOR = 10.0

#: A point feature further than this from every vertex, in mean edge lengths, is not a
#: vertex of this mesh, and fixing whatever vertex happens to be nearest would be wrong.
POINT_FEATURE_REACH = 0.5

#: What option_line returns when the user presses Enter.
ENTER = "<enter>"


# ==============================================================================
# Putting the options together -- no Rhino input from here to the UI section
# ==============================================================================

def boundary_vertices(mesh):
    # NOT vertices_on_boundary(): in COMPAS 2 that returns the LONGEST boundary only, so on
    # a mesh with a hole the hole is left out -- and left free to blow up.
    return set(vertex for loop in mesh.vertices_on_boundaries() for vertex in loop)


def grow(mesh, core, rings=1):
    """The core plus ``rings`` of neighbours around it."""
    out = set(core)
    frontier = set(core)
    for _ in range(rings):
        frontier = set(n for vertex in frontier for n in mesh.vertex_neighbors(vertex)) - out
        out |= frontier
    return out


def taper(mesh, core, blend):
    """Per-vertex damping weight: 1 in the core, falling to 0 outside the blend.

    The same taper ``smoothing_region`` builds for itself. It is built here instead, so the
    region that is previewed is exactly the region that moves, and so fixed vertices can be
    given a weight of 0 inside it.
    """
    weights = dict((vertex, 1.0) for vertex in core)
    frontier = set(core)
    for ring in range(1, blend + 1):
        frontier = set(n for vertex in frontier for n in mesh.vertex_neighbors(vertex)) - set(weights)
        if not frontier:
            break
        weight = 1.0 - ring / float(blend + 1)
        for vertex in frontier:
            weights[vertex] = weight
    return weights


def propose_chain(mesh, guide, tolerance, max_angle, polyedges):
    """The chain for one guide, why it is empty if it is, and how it measures."""
    selected, info = guide_chain(mesh, guide, tolerance=tolerance, max_angle=max_angle,
                                 polyedges=polyedges)
    quality = chain_quality(mesh, selected, guide) if selected else None
    return selected, info, quality


def attach_guides(mesh, guides, fixed=(), allowed=None):
    """Move every guide's chain onto its guide, and return what holds each vertex there.

    Every move is worked out on the mesh as it is BEFORE any of them is applied: a guide
    attached second must not choose where to go on a mesh the first one has already pulled
    about, or the result would depend on the order the guides were picked in. A vertex
    claimed by two guides is held by the one listed last.

    A BOUNDARY vertex is never moved onto a guide; ``attach_chain`` holds it on its own
    outline, so at most it slides. Moving it would take the outline with it.

    Parameters
    ----------
    guides : list[dict]
        ``{'guide': GuideCurve, 'chain': [vertex, ...], 'hold': 'fixed' | 'sliding'}``.
    fixed : iterable[int]
        Vertices left out of every chain. A fixed vertex does not move, not even onto a
        guide.
    allowed : set[int], optional
        If given, every vertex not in it is left out as well.
    """
    fixed = set(fixed)
    moves, constraints, overlap = {}, {}, set()
    skipped_fixed, skipped_outside, attached_guides = 0, 0, 0

    for entry in guides:
        chain = []
        for vertex in entry["chain"]:
            if vertex in fixed:
                skipped_fixed += 1
            elif allowed is not None and vertex not in allowed:
                skipped_outside += 1
            else:
                chain.append(vertex)
        if not chain:
            continue
        attached_guides += 1
        entry_moves, entry_constraints = attach_chain(mesh, chain, entry["guide"],
                                                      hold=entry["hold"])
        overlap.update(vertex for vertex in entry_constraints if vertex in constraints)
        moves.update(entry_moves)
        constraints.update(entry_constraints)

    for vertex, xyz in moves.items():
        mesh.vertex_attributes(vertex, "xyz", xyz)

    return {
        "constraints": constraints,
        "moved": set(moves),
        # interior vertices held ON the guide rather than pinned where they landed
        "sliding": set(vertex for vertex in moves if not isinstance(constraints[vertex], Point)),
        # boundary vertices of a chain: not moved, held on their own outline
        "anchored": set(constraints) - set(moves),
        "guides": attached_guides,
        "overlap": len(overlap),
        "skipped_fixed": skipped_fixed,
        "skipped_outside": skipped_outside,
    }


def _report(attached, pinned, notes, **extra):
    report = {
        "guides": attached["guides"],
        "attached": len(attached["constraints"]),
        "moved": len(attached["moved"]),
        "anchored": len(attached["anchored"]),
        "overlap": attached["overlap"],
        "skipped_fixed": attached["skipped_fixed"],
        "skipped_outside": attached["skipped_outside"],
        "pinned": len(pinned),
        "notes": notes,
    }
    report.update(extra)
    return report


def boundary_holds(mesh, boundary):
    """``(constraints, pinned)`` for a Boundary setting, read before anything moves.

    ``sliding`` -- every boundary vertex slides along its own outline, corners pinned.
    ``fixed``   -- every boundary vertex pinned.
    ``free``    -- nothing: the boundary, corners included, is smoothed like the interior,
                   and only fixed vertices and guides hold the mesh.
    """
    if boundary == "sliding":
        return automated_boundary_constraints(mesh), set()
    if boundary == "fixed":
        return {}, boundary_vertices(mesh)
    if boundary == "free":
        return {}, set()
    raise ValueError("boundary must be sliding, fixed or free, got {!r}".format(boundary))


def _release_anchors(constraints, attached, boundary):
    """On a Free boundary a chain's boundary ends are not held on the outline either.

    attach_chain constrains them to their wall as it was; with the rest of that wall free to
    move, holding just the ends to it would pin the outline at two points.
    """
    if boundary == "free":
        for vertex in attached["anchored"]:
            constraints.pop(vertex, None)


def smooth_whole(mesh, algorithm="area", boundary="sliding", fixed=(), guides=(),
                 kmax=100, damping=0.5, q_factor=DEFAULT_Q_FACTOR):
    """Smooth the whole mesh in place with the options of the command. Returns a report.

    Raises
    ------
    ValueError
        If ForceDensity would run with nothing fixed.
    """
    algorithm = algorithm.lower()
    boundary = boundary.lower()
    fixed = set(fixed)
    walls = boundary_vertices(mesh)
    notes = []

    if algorithm == "forcedensity":
        if boundary == "fixed":
            pinned = set(walls)
        elif boundary == "free":
            # the form is found from the fixed vertices and the guides alone; the boundary
            # edges are still BoundaryStiffness times stiffer, which is what keeps the
            # outline from collapsing onto the supports
            pinned = set()
        else:
            pinned = set(mesh_boundary_corners(mesh))
            notes.append(
                "ForceDensity cannot slide a boundary: only the {} corner(s) are held, and "
                "the rest of the outline is held by stiff edges alone, so it can drift off "
                "the wall. Use Boundary=Fixed, or another algorithm, to keep it on the "
                "wall.".format(len(pinned)))
        attached = attach_guides(mesh, guides, fixed)
        # an interior vertex on a guide is pinned where it landed -- the only hold there is
        pinned |= attached["moved"]
        if attached["sliding"]:
            notes.append(
                "ForceDensity cannot hold a vertex on a guide: {} vertices with a Sliding "
                "hold were pinned where they landed instead.".format(len(attached["sliding"])))
        pinned |= fixed
        if not pinned:
            raise ValueError(
                "ForceDensity needs at least one fixed vertex and there is none: no fixed "
                "boundary, no corners, no fixed vertices and no guides. Fix some vertices, "
                "or set Boundary=Fixed or Sliding.")
        relaxation(mesh, fixed="manual", fixed_vertices=sorted(pinned), q_factor=q_factor)
        return _report(attached, pinned, notes)

    # the constraints read the outline before any guide moves anything
    constraints, pinned = boundary_holds(mesh, boundary)

    attached = attach_guides(mesh, guides, fixed)
    constraints.update(attached["constraints"])
    _release_anchors(constraints, attached, boundary)
    if boundary == "free" and not fixed and not constraints:
        notes.append(
            "the boundary is free and nothing is fixed or on a guide, so nothing holds the "
            "mesh: it shrinks a little with every iteration.")
    # An attached BOUNDARY vertex is held on its own outline by attach_chain -- the same
    # constraint a sliding boundary gives it. On a Fixed boundary that is what lets a
    # cable's anchored end follow the chain along the wall, and it still cannot leave it.
    pinned -= set(attached["constraints"])
    # and a fixed vertex is fixed, whatever else was going to hold it
    pinned |= fixed
    for vertex in fixed:
        constraints.pop(vertex, None)

    if algorithm == "centerofmass":
        # compas' centre-of-mass rule takes the polygon of a vertex's neighbours and RAISES
        # ("At least three points required") on a vertex with two -- which is what every
        # corner of a quad mesh has. A Fixed boundary never reaches it; Sliding and Free do.
        lonely = set(vertex for vertex in mesh.vertices()
                     if vertex not in pinned and len(mesh.vertex_neighbors(vertex)) < 3)
        # a corner already pinned by a Point constraint loses nothing by being fixed instead
        mobile = [vertex for vertex in lonely if not isinstance(constraints.get(vertex), Point)]
        if mobile:
            notes.append(
                "CenterOfMass cannot move a vertex with only two neighbours: {} of them, which "
                "would have moved, were pinned instead.".format(len(mobile)))
        pinned |= lonely
        for vertex in lonely:
            constraints.pop(vertex, None)

    constrained_smoothing(mesh, kmax=kmax, damping=damping, constraints=constraints,
                          algorithm=algorithm, fixed=sorted(pinned))
    return _report(attached, pinned, notes)


def smooth_region(mesh, core, blend=DEFAULT_BLEND, boundary="sliding", fixed=(), guides=(),
                  kmax=50, damping=0.5):
    """Smooth a region of the mesh in place with the options of the command. Returns a report.

    Guide chains are clipped to the CORE. A chain vertex outside it would be moved onto its
    guide where the smoothing is too weak to blend the move in -- or, outside the blend,
    where there is no smoothing at all -- and that is a crease.
    """
    boundary = boundary.lower()
    core = set(core)
    fixed = set(fixed)
    weights = taper(mesh, core, blend)
    constraints, pinned = boundary_holds(mesh, boundary)

    attached = attach_guides(mesh, guides, fixed, allowed=core)
    constraints.update(attached["constraints"])
    _release_anchors(constraints, attached, boundary)
    pinned -= set(attached["constraints"])
    pinned |= fixed

    # smoothing_region fixes a vertex whose weight is not positive, and only projects the
    # constraints of the vertices it moves
    for vertex in pinned:
        if vertex in weights:
            weights[vertex] = 0.0
        constraints.pop(vertex, None)

    # a dict, never None: None would make smoothing_region build sliding boundary
    # constraints of its own, and a Fixed boundary would slide
    smoothing_region(mesh, weights, kmax=kmax, damping=damping, constraints=constraints)
    moving = len([vertex for vertex, weight in weights.items() if weight > 0.0])
    return _report(attached, set(vertex for vertex in pinned if vertex in weights), [],
                   core=len(core), moving=moving)


# ==============================================================================
# Rhino
# ==============================================================================

def quad_mesh_filter(rhobj, geometry, component_index):
    layer = rs.ObjectLayer(rhobj)
    quad_layer = rs.LayerName("QuadMesh", fullpath=True)
    return layer == quad_layer or rs.IsLayerChildOf(quad_layer, layer)


def _xyz(point):
    if hasattr(point, "X"):
        return [point.X, point.Y, point.Z]
    return [point[0], point[1], point[2]]


def _index(added):
    """The index an ``AddOption*`` call returns.

    The overloads that take their option by ``ref`` -- Toggle, Integer, Double -- return
    ``(index, option)`` under Rhino 8's CPython, not the index (checked in Rhino). 0 means
    Rhino refused the option, e.g. because the name was already used on this line.
    """
    return added[0] if isinstance(added, tuple) else added


def option_line(prompt, options):
    """Show one command-line option line and wait.

    ``options`` is a list of ``(key, option)``: ``option`` is a plain option name, or a
    function that adds the option to the GetOption and returns what that call returned.

    Returns ``(key, option)`` for a click, ``(ENTER, None)`` for Enter and ``(None, None)``
    for Esc. A typed Integer or Double option has already stored its new value when this
    returns.
    """
    go = Rhino.Input.Custom.GetOption()
    go.SetCommandPrompt(prompt)
    go.AcceptNothing(True)
    keys = {}
    for key, option in options:
        index = go.AddOption(option) if isinstance(option, str) else _index(option(go))
        if index > 0:
            keys[index] = key
    result = go.Get()
    if result == Rhino.Input.GetResult.Nothing:
        return ENTER, None
    if result != Rhino.Input.GetResult.Option:
        return None, None
    return keys.get(go.OptionIndex()), go.Option()


def guide_curves():
    """Every curve on the Guides layer and its sublayers.

    These are the curves CMD_boundary_selection collected and the layout was generated
    from, so holding the mesh to them holds it to what it came from.
    """
    if not rs.IsLayer(GUIDE_LAYER):
        return []
    found, seen = [], set()
    for layer in [GUIDE_LAYER] + (rs.LayerChildren(GUIDE_LAYER) or []):
        for guid in rs.ObjectsByLayer(layer) or []:
            if rs.IsCurve(guid) and str(guid) not in seen:
                seen.add(str(guid))
                found.append(guid)
    return found


def build_guide(curve_id):
    """A Rhino curve as a GuideCurve: samples to choose the chain on, the curve to land on."""
    # max_edge only sets how finely a curved guide is sampled; a polyline guide is taken
    # at its own vertices either way
    points = curve_points(curve_id, 0.25)
    # curve_points DROPS the closing point -- a loop is closed by convention there. A ring
    # handed over like that is an arc with a gap in it, and the arc-length parameter
    # cannot wrap across the seam.
    if rs.IsCurveClosed(curve_id) and len(points) > 2:
        points = points + points[:1]
    return GuideCurve(points, curve=guide_geometry(curve_id))


def guide_geometry(curve_id):
    """The guide as a compas curve on World XY, for its vertices to land ON -- or None.

    None for a polyline: its samples are its own vertices, so landing on them already is
    landing on the guide. An arc, a circle or a NURBS curve is sampled at 0.25, and a vertex
    landing on those samples lands on a chord, up to a sagitta off what was drawn.

    Flattened onto World XY because ``curve_points`` puts every sample at z = 0, and the
    curve has to lie where the samples do: a guide drawn at another height would otherwise
    lift the vertices attached to it off the plan. None too if Rhino cannot flatten it,
    which leaves the samples to land on, as before.
    """
    curve = rs.coercecurve(curve_id)
    if curve is None or curve.TryGetPolyline()[0]:
        return None
    flat = Rhino.Geometry.Curve.ProjectToPlane(curve, Rhino.Geometry.Plane.WorldXY)
    if flat is None:
        return None
    return curve_to_compas(flat)


class Preview(object):
    """Coloured points on their own layer, showing what the smoothing will hold.

    **Drawn with RhinoCommon, not rhinoscriptsyntax, and only when it has changed.**
    ``rs.ObjectLayer`` and ``rs.ObjectColor`` each end in a full ``doc.Views.Redraw()``, as
    do ``rs.AddPoints`` and ``rs.DeleteObjects``. Colouring N points one call at a time was
    2N redraws, measured at 11 ms each in a 3 300-object document -- seconds per click, on
    every pass round the menu, even for a click on Algorithm that changes nothing on screen.
    Here every point is added with its attributes already set, and the viewport is redrawn
    once.
    """

    def __init__(self, doc, layer=PREVIEW_LAYER):
        self.doc = doc
        self.layer = layer
        self.guids = []
        self.shown = None

    def _layer_index(self):
        index = self.doc.Layers.FindByFullPath(self.layer, -1)
        if index < 0:
            index = self.doc.Layers.Add(self.layer, System.Drawing.Color.Red)
        elif self.shown is None:
            # left behind by a run that did not get to purge it
            for rhino_object in self.doc.Objects.FindByLayer(self.doc.Layers[index]) or []:
                self.doc.Objects.Delete(rhino_object.Id, True)
        return index

    def draw(self, colours, coordinates):
        """``colours``: {vertex: (r, g, b)}; ``coordinates``: vertex -> [x, y, z]."""
        signature = frozenset(colours.items())
        if signature == self.shown:
            return False
        index = self._layer_index()
        self.clear(redraw=False)
        groups = {}
        for vertex, colour in colours.items():
            groups.setdefault(colour, []).append(vertex)
        for colour, vertices in groups.items():
            attributes = Rhino.DocObjects.ObjectAttributes()
            attributes.LayerIndex = index
            attributes.ColorSource = Rhino.DocObjects.ObjectColorSource.ColorFromObject
            attributes.ObjectColor = System.Drawing.Color.FromArgb(colour[0], colour[1], colour[2])
            for vertex in vertices:
                x, y, z = coordinates(vertex)
                self.guids.append(self.doc.Objects.AddPoint(Rhino.Geometry.Point3d(x, y, z), attributes))
        self.shown = signature
        self.doc.Views.Redraw()
        return True

    def clear(self, redraw=True):
        for guid in self.guids:
            self.doc.Objects.Delete(guid, True)
        self.guids = []
        if redraw:
            self.shown = None
            self.doc.Views.Redraw()


class SmoothCommand(object):
    """The state of the option lines, the picks that change it, and the preview."""

    def __init__(self, mesh_id, region):
        self.mesh_id = mesh_id
        self.region = region
        geometry = cr.objects.find_object(mesh_id).Geometry
        # QuadMesh, not a plain Mesh: the guide chains are chosen from its polyedges
        self.mesh = QuadMesh.from_vertices_and_faces(
            *mesh_from_rhino(geometry).to_vertices_and_faces())
        self.rhino_vertices = rs.MeshVertices(mesh_id)
        # nothing moves until the options are settled, so this stays true for every pick
        self.coordinates = [(vertex, self.mesh.vertex_coordinates(vertex))
                            for vertex in self.mesh.vertices()]
        self.average = mean_edge_length(self.mesh)
        self._polyedges = None

        self.preview = Preview(sc.doc)

        self.algorithm = 0
        self.boundary = "Sliding"
        self.fixed = []
        self.guides = []
        self.core = []
        self.blend = Rhino.Input.Custom.OptionInteger(DEFAULT_BLEND, 0, 1000)
        self.iterations = Rhino.Input.Custom.OptionInteger(50 if region else 100, 1, 100000)
        self.damping = Rhino.Input.Custom.OptionDouble(0.5, 0.0, 1.0)
        self.q_factor = Rhino.Input.Custom.OptionDouble(DEFAULT_Q_FACTOR, 0.001, 1.0e6)
        self.tolerance = Rhino.Input.Custom.OptionDouble(float(DEFAULT_TOLERANCE_FACTOR), 0.05, 5.0)
        self.max_angle = Rhino.Input.Custom.OptionDouble(float(DEFAULT_MAX_ANGLE), 1.0, 90.0)

    # ------------------------------------------------------------------
    # picks
    # ------------------------------------------------------------------

    def nearest_vertex(self, xyz):
        """The compas vertex closest to a point, and how far it is.

        Matched by POSITION, not by index. mesh_to_compas and the bake round-trip both go
        through float conversions, so an index that looks like it should line up is not
        something to bet a constraint on.
        """
        best, best_d = None, None
        for vertex, point in self.coordinates:
            d = (point[0] - xyz[0]) ** 2 + (point[1] - xyz[1]) ** 2 + (point[2] - xyz[2]) ** 2
            if best_d is None or d < best_d:
                best, best_d = vertex, d
        return best, (best_d or 0.0) ** 0.5

    def pick_vertices(self, message):
        picked = rs.GetMeshVertices(self.mesh_id, message)
        if not picked:
            return []
        return [self.nearest_vertex(_xyz(self.rhino_vertices[i]))[0] for i in picked]

    def vertices_in_curves(self):
        curve_ids = rs.GetObjects(message="Pick one or more CLOSED curves around the zone",
                                  filter=rs.filter.curve, preselect=False, select=False,
                                  minimum_count=1) or []
        closed = [curve_id for curve_id in curve_ids if rs.IsCurveClosed(curve_id)]
        if len(closed) < len(curve_ids):
            print("{} of the picked curves not closed -- skipped.".format(
                len(curve_ids) - len(closed)))
        inside = []
        for vertex, xyz in self.coordinates:
            # 1 = inside, 2 = on the curve
            if any(rs.PointInPlanarClosedCurve(xyz, curve_id) in (1, 2) for curve_id in closed):
                inside.append(vertex)
        return inside

    def point_feature_vertices(self):
        """The vertex at every point on the PointFeatures layer."""
        if not rs.IsLayer(POINT_LAYER):
            print("No {} layer -- no point features to fix.".format(POINT_LAYER))
            return []
        found, far = [], 0
        for guid in rs.ObjectsByLayer(POINT_LAYER) or []:
            if not rs.IsPoint(guid):
                continue
            vertex, distance = self.nearest_vertex(_xyz(rs.PointCoordinates(guid)))
            if distance > POINT_FEATURE_REACH * self.average:
                far += 1
                continue
            found.append(vertex)
        print("{} point feature(s) matched to a vertex.".format(len(found)))
        if far:
            print("  {} skipped: no vertex within {} edge length(s) of them, so they are not "
                  "vertices of this mesh.".format(far, POINT_FEATURE_REACH))
        return found

    def polyedges(self):
        if self._polyedges is None:
            nonquad = [f for f in self.mesh.faces() if len(self.mesh.face_vertices(f)) != 4]
            if nonquad:
                print("note: {} non-quad face(s) -- polyedges stop at those.".format(len(nonquad)))
            # collected ONCE: it costs O(edges squared) and does not depend on the guide
            self._polyedges = collect_polyedges(self.mesh)
        return self._polyedges

    @staticmethod
    def _add(target, vertices):
        for vertex in vertices:
            if vertex not in target:
                target.append(vertex)

    # ------------------------------------------------------------------
    # preview
    # ------------------------------------------------------------------

    def show(self, pending=()):
        """Draw what will happen: fixed blue, guide chains purple, region core red, blend orange.

        One point per vertex, by priority, so the preview shows what holds each vertex.
        ``pending`` are guide entries being edited, drawn in place of their saved versions.
        Costs nothing when the picture has not changed -- see :class:`Preview`.
        """
        colours = {}
        if self.region:
            for vertex, weight in taper(self.mesh, self.core, self.blend.CurrentValue).items():
                colours[vertex] = CORE_COLOUR if weight >= 1.0 else BLEND_COLOUR
        allowed = set(self.core) if self.region else None
        for entry in self.guide_entries(pending):
            for vertex in entry["chain"]:
                if allowed is None or vertex in allowed:
                    colours[vertex] = GUIDE_COLOUR
        for vertex in self.fixed:
            colours[vertex] = FIXED_COLOUR
        # the mesh does not move until the options are settled, so a vertex key stands for
        # its position here
        self.preview.draw(colours, self.mesh.vertex_coordinates)

    def guide_entries(self, pending=()):
        keys = set(entry["key"] for entry in pending)
        return [entry for entry in self.guides if entry["key"] not in keys] + list(pending)

    # ------------------------------------------------------------------
    # the main option line
    # ------------------------------------------------------------------

    def menu(self):
        """The main option line. True to smooth, False when cancelled."""
        while True:
            self.show()
            force_density = not self.region and ALGORITHMS[self.algorithm] == "ForceDensity"
            options = []
            if self.region:
                options.append(("region", lambda go: go.AddOption("Region", str(len(self.core)))))
                options.append(("blend", lambda go: go.AddOptionInteger("Blend", self.blend)))
            else:
                options.append(("algorithm", lambda go: go.AddOptionList(
                    "Algorithm", ALGORITHMS, self.algorithm)))
            options.append(("boundary", lambda go: go.AddOptionList(
                "Boundary", BOUNDARY_MODES, BOUNDARY_MODES.index(self.boundary))))
            options.append(("fixed", lambda go: go.AddOption("FixedVertices", str(len(self.fixed)))))
            options.append(("guides", lambda go: go.AddOption("Guides", str(len(self.guides)))))
            if force_density:
                # no iterations and no damping to set: relaxation solves to equilibrium
                options.append(("stiffness", lambda go: go.AddOptionDouble(
                    "BoundaryStiffness", self.q_factor)))
            else:
                options.append(("iterations", lambda go: go.AddOptionInteger(
                    "Iterations", self.iterations)))
                options.append(("damping", lambda go: go.AddOptionDouble("Damping", self.damping)))

            what = "a region" if self.region else "the whole mesh"
            key, option = option_line("Smooth {}. Press Enter to run".format(what), options)
            if key is None:
                return False
            if key == ENTER:
                if self.region and not self.core:
                    print("The region is empty -- add vertices to it first.")
                    continue
                return True
            if key == "algorithm":
                self.algorithm = option.CurrentListOptionIndex
            elif key == "boundary":
                self.boundary = BOUNDARY_MODES[option.CurrentListOptionIndex]
            elif key == "region":
                self.edit_region()
            elif key == "fixed":
                self.edit_fixed()
            elif key == "guides":
                self.edit_guides()
            # Blend, Iterations, Damping, BoundaryStiffness: Rhino has stored the value

    # ------------------------------------------------------------------
    # region
    # ------------------------------------------------------------------

    def pick_region(self):
        """The first pick of the region. True if it selected anything."""
        key, _ = option_line("Select the zone with closed curves or vertices. Press Enter for Curve",
                             [("curve", "Curve"), ("vertices", "Vertices")])
        if key is None:
            return False
        if key == "vertices":
            self._add(self.core, self.pick_vertices("Pick vertices in the zone"))
        else:
            self._add(self.core, self.vertices_in_curves())
        print("region: {} core vertices.".format(len(self.core)))
        return bool(self.core)

    def edit_region(self):
        while True:
            self.show()
            prompt = "Region: {} core vertices, {} with the blend. Press Enter when done".format(
                len(self.core), len(taper(self.mesh, self.core, self.blend.CurrentValue)))
            key, _ = option_line(prompt, [("curve", "Curve"), ("add", "Add"), ("remove", "Remove"),
                                          ("grow", "Grow"), ("clear", "Clear")])
            if key in (None, ENTER):
                return
            if key == "curve":
                self._add(self.core, self.vertices_in_curves())
            elif key == "add":
                self._add(self.core, self.pick_vertices("Pick vertices to ADD to the region"))
            elif key == "remove":
                drop = set(self.pick_vertices("Pick vertices to REMOVE from the region"))
                self.core = [vertex for vertex in self.core if vertex not in drop]
            elif key == "grow":
                self._add(self.core, grow(self.mesh, self.core, 1))
            elif key == "clear":
                self.core = []

    # ------------------------------------------------------------------
    # fixed vertices
    # ------------------------------------------------------------------

    def edit_fixed(self):
        while True:
            self.show()
            prompt = "Fixed vertices: {}. Press Enter when done".format(len(self.fixed))
            key, _ = option_line(prompt, [("add", "Add"), ("points", "PointFeatures"),
                                          ("remove", "Remove"), ("clear", "Clear")])
            if key in (None, ENTER):
                return
            if key == "add":
                self._add(self.fixed, self.pick_vertices("Pick vertices to FIX"))
            elif key == "points":
                self._add(self.fixed, self.point_feature_vertices())
            elif key == "remove":
                drop = set(self.pick_vertices("Pick fixed vertices to FREE"))
                self.fixed = [vertex for vertex in self.fixed if vertex not in drop]
            elif key == "clear":
                self.fixed = []

    # ------------------------------------------------------------------
    # guides
    # ------------------------------------------------------------------

    def edit_guides(self):
        while True:
            self.show()
            count = len(set(vertex for entry in self.guides for vertex in entry["chain"]))
            prompt = ("Guides: {} with {} vertices (Distance is in edge lengths). Press Enter "
                      "when done").format(len(self.guides), count)
            # Distance and Angle are the two gates the proposal is chosen by, and they do
            # different jobs. DISTANCE is how far off the guide a vertex may sit -- and it is
            # pulled onto the guide by exactly that far, so it sets both coverage and damage.
            # ANGLE is whether the polyedge still runs the guide's way, which is what cuts a
            # chain where the line veers off near its ends. 90 turns the angle gate off.
            # They apply to the next proposal; a chain already made keeps its vertices.
            key, _ = option_line(prompt, [
                ("all", "AllGuides"), ("pick", "Pick"), ("remove", "Remove"), ("clear", "Clear"),
                ("distance", lambda go: go.AddOptionDouble("Distance", self.tolerance)),
                ("angle", lambda go: go.AddOptionDouble("Angle", self.max_angle)),
            ])
            if key in (None, ENTER):
                return
            if key == "pick":
                self.pick_guide()
            elif key == "all":
                self.all_guides()
            elif key == "remove":
                self.remove_guide()
            elif key == "clear":
                self.guides = []

    def propose(self, guide):
        """Propose a chain for one guide and say what it is. Returns (chain, info)."""
        selected, info, quality = propose_chain(
            self.mesh, guide, self.tolerance.CurrentValue * self.average,
            self.max_angle.CurrentValue, self.polyedges())
        if not selected:
            return selected, info
        # coverage says how much of the guide it got; alignment says whether the mesh has a
        # course along this guide at all -- low alignment means it has not, and no
        # selection will fix that
        print("  {} vertices: {:.0f}% of the guide, alignment {:.2f}, worst face angle if "
              "attached {:.0f}.".format(len(selected), 100 * quality["coverage"],
                                        quality["alignment"], quality["min_face_angle"] or 0.0))
        if quality["alignment"] < 0.8:
            print("  note: this guide does not follow a course of the mesh.")
        if quality["folded_faces"]:
            print("  WARNING: attaching this one turns {} face(s) inside out. It is reaching "
                  "{:.1f} edge lengths to do it -- lower Distance to stop it.".format(
                      quality["folded_faces"], quality["off_max"]))
        if self.region:
            core = set(self.core)
            outside = len([vertex for vertex in selected if vertex not in core])
            if outside:
                print("  {} of them are outside the region and will be left alone.".format(outside))
        return selected, info

    def pick_guide(self):
        curve_id = rs.GetObject(message="Pick a guide curve", filter=rs.filter.curve,
                                preselect=False, select=False)
        if not curve_id:
            return
        key = str(curve_id)
        existing = [entry for entry in self.guides if entry["key"] == key]
        if existing:
            print("editing the chain this guide already has.")
            entry = existing[0]
        else:
            print("guide:")
            guide = build_guide(curve_id)
            selected, info = self.propose(guide)
            if not selected:
                print("  {} -- use Add to pick the chain by hand.".format(info["reason"]))
            entry = {"key": key, "guide": guide, "chain": selected, "hold": "fixed"}

        edited = self.edit_chain(entry)
        if edited is None:
            print("  chain left as it was." if existing else "  discarded.")
            return
        self.guides = [other for other in self.guides if other["key"] != key]
        if not edited["chain"]:
            print("  no vertices -- this guide holds nothing.")
            return
        self.guides.append(edited)

    def edit_chain(self, entry):
        """Add / Remove / Clear on one chain. The edited entry, or None when discarded."""
        entry = dict(entry, chain=list(entry["chain"]))
        while True:
            self.show(pending=[entry])
            prompt = "Chain: {} vertices. Press Enter to keep it, Esc to discard".format(
                len(entry["chain"]))
            # Fixed pins an interior vertex where it lands on the guide; Sliding re-projects
            # it every iteration, so it may travel ALONG the guide but never leave it. There
            # is no ordering guard on Sliding -- two vertices may pass each other.
            key, _ = option_line(prompt, [
                ("add", "Add"), ("remove", "Remove"), ("clear", "Clear"),
                ("hold", lambda go: go.AddOption("Hold", entry["hold"].capitalize())),
            ])
            if key is None:
                return None
            if key == ENTER:
                return entry
            if key == "add":
                self._add(entry["chain"], self.pick_vertices("Pick vertices to ADD to the chain"))
            elif key == "remove":
                drop = set(self.pick_vertices("Pick vertices to REMOVE from the chain"))
                entry["chain"] = [vertex for vertex in entry["chain"] if vertex not in drop]
            elif key == "clear":
                entry["chain"] = []
            elif key == "hold":
                entry["hold"] = "sliding" if entry["hold"] == "fixed" else "fixed"

    def all_guides(self):
        curve_ids = guide_curves()
        if not curve_ids:
            print("No curves on the {} layer -- pick the guides one at a time instead.".format(
                GUIDE_LAYER))
            return
        present = set(entry["key"] for entry in self.guides)
        print("proposing a chain for every guide on the {} layer: {} curve(s), distance {:.2f} "
              "edge lengths, angle {:.0f} degrees.".format(
                  GUIDE_LAYER, len(curve_ids), self.tolerance.CurrentValue,
                  self.max_angle.CurrentValue))

        new, refused = [], []
        for number, curve_id in enumerate(curve_ids, 1):
            if str(curve_id) in present:
                print("guide {}: already has a chain -- kept as it is.".format(number))
                continue
            print("guide {}:".format(number))
            guide = build_guide(curve_id)
            selected, info = self.propose(guide)
            if not selected:
                print("  skipped -- {}.".format(info["reason"]))
                refused.append(number)
                continue
            new.append({"key": str(curve_id), "guide": guide, "chain": selected, "hold": "fixed"})

        if refused:
            print("{} guide(s) got nothing: {}. Widen Distance or Angle, or Pick them to choose "
                  "the chain by hand.".format(len(refused), ", ".join(str(i) for i in refused)))
        if not new:
            return

        # everything is on screen before the one question is asked
        hold = "fixed"
        count = len(set(vertex for entry in new for vertex in entry["chain"]))
        while True:
            self.show(pending=new)
            prompt = "{} new chain(s), {} vertices. Press Enter to add them, Esc to discard".format(
                len(new), count)
            key, _ = option_line(prompt, [("hold", lambda go: go.AddOption("Hold", hold.capitalize()))])
            if key is None:
                print("discarded.")
                return
            if key == ENTER:
                break
            hold = "sliding" if hold == "fixed" else "fixed"
        for entry in new:
            entry["hold"] = hold
        self.guides.extend(new)

    def remove_guide(self):
        curve_id = rs.GetObject(message="Pick the guide curve to remove", filter=rs.filter.curve,
                                preselect=False, select=False)
        if not curve_id:
            return
        before = len(self.guides)
        self.guides = [entry for entry in self.guides if entry["key"] != str(curve_id)]
        if len(self.guides) == before:
            print("That curve has no chain.")

    # ------------------------------------------------------------------
    # run
    # ------------------------------------------------------------------

    def run(self):
        """Smooth, and return (report, the layer to bake to, a line saying what was done)."""
        kmax = int(self.iterations.CurrentValue)
        damping = float(self.damping.CurrentValue)
        if self.region:
            blend = int(self.blend.CurrentValue)
            report = smooth_region(self.mesh, self.core, blend=blend, boundary=self.boundary,
                                   fixed=self.fixed, guides=self.guides, kmax=kmax, damping=damping)
            summary = ("smoothed a region: {} core vertices, {} moving with a blend of {} rings, "
                       "boundary {}, {} iterations at damping {}.").format(
                report["core"], report["moving"], blend, self.boundary.lower(), kmax, damping)
            return report, REGION_OUTPUT_LAYER, summary

        algorithm = ALGORITHMS[self.algorithm]
        report = smooth_whole(self.mesh, algorithm, boundary=self.boundary, fixed=self.fixed,
                              guides=self.guides, kmax=kmax, damping=damping,
                              q_factor=float(self.q_factor.CurrentValue))
        if algorithm == "ForceDensity":
            how = "boundary stiffness {}".format(self.q_factor.CurrentValue)
        else:
            how = "{} iterations at damping {}".format(kmax, damping)
        summary = "smoothed the whole mesh: {}, boundary {}, {}.".format(
            algorithm, self.boundary.lower(), how)
        return report, algorithm, summary


def print_report(report, fixed_count, boundary):
    anchored = "held nowhere, like the rest of the free boundary" if boundary.lower() == "free" \
        else "they slide along it"
    if fixed_count:
        print("  {} fixed vertices.".format(fixed_count))
    if report["attached"] or report["skipped_fixed"] or report["skipped_outside"]:
        print("  {} guide(s): {} vertices attached, {} moved onto their guide, {} anchored "
              "on the boundary ({}).".format(
                  report["guides"], report["attached"], report["moved"], report["anchored"],
                  anchored))
    if report["overlap"]:
        print("  note: {} vertices were in two chains and are held by the guide picked "
              "last.".format(report["overlap"]))
    if report["skipped_fixed"]:
        print("  note: {} chain vertices are fixed vertices, and stayed where they "
              "are.".format(report["skipped_fixed"]))
    if report["skipped_outside"]:
        print("  note: {} chain vertices are outside the region and were left "
              "alone.".format(report["skipped_outside"]))
    for note in report["notes"]:
        print("  note: " + note)


def main():
    if not rs.IsLayer("QuadMesh"):
        print("No dense quad mesh yet -- the QuadMesh layer does not exist.")
        return

    mesh_id = rs.GetObject(
        message="Pick a dense quad mesh to smooth from the layer QuadMesh and its sublayers",
        filter=rs.filter.mesh, preselect=False, select=False, custom_filter=quad_mesh_filter)
    if not mesh_id:
        print("No mesh picked.")
        return

    key, _ = option_line("Smooth what? Press Enter for Whole",
                         [("whole", "Whole"), ("region", "Region")])
    if key is None:
        print("Cancelled -- nothing smoothed.")
        return

    command = SmoothCommand(mesh_id, region=(key == "region"))
    try:
        if command.region and not command.pick_region():
            print("No region selected -- nothing smoothed.")
            return
        if not command.menu():
            print("Cancelled -- nothing smoothed.")
            return
    finally:
        if rs.IsLayer(PREVIEW_LAYER):
            rs.PurgeLayer(PREVIEW_LAYER)

    try:
        report, layer, summary = command.run()
    except Exception as exc:
        # nothing is baked, so the mesh in the document is untouched
        print("Smoothing failed, nothing baked: {}: {}".format(type(exc).__name__, exc))
        return

    rs.AddLayer("Smoothened", parent="QuadMesh")
    rs.AddLayer(layer, parent="Smoothened")
    clear_layer(layer, clean_sublayers=True)
    bake_mesh(command.mesh, layer)
    print(summary)
    print_report(report, len(command.fixed), command.boundary)
    print("  baked to Smoothened::{}.".format(layer))


if __name__ == "__main__":
    main()
