#! python3

# r: compas

"""Shared plumbing for the ``ff_0*`` commands. Not a command -- import only.

Everything here is conversion, layer handling and session state. No geometry is
computed in this file; every result comes from ``compas_singular`` or
``framefield``.

**The document holds the state, not this module.** Each step bakes its result
on its own layer and the next step reads it back, so the workflow survives a
save, a reopen, and running the steps out of order. The two things geometry
cannot express -- the solver settings and the per-strip densities -- go in
document user text as small JSON blobs. The one thing that is neither, the
solved field, is never stored: it is 92 KB at 0.5 spacing and the solve is
DETERMINISTIC, so it is rebuilt on demand and cached in ``scriptcontext.sticky``
for the rest of the session.

Densities are keyed GEOMETRICALLY, and that is not a style choice. Strip indices
come from ``collect_strips`` popping ``self.edges()``, so they follow vertex
insertion order -- and a layout welded back from Rhino has whatever order the
geometry came back in. Measured: bake a 4-strip layout, read it back, 0 of 4
indices still point at the same strip. A density map keyed by index does not
error, it silently applies the wrong numbers. So a density is stored against the
rounded midpoint of one remembered edge, and the strip that owns that edge is
looked up again on the other side.
"""
import json
import sys


# ----------------------------------------------------------------------
# where the code is -- EDIT THESE if the repository moves
# ----------------------------------------------------------------------

SINGULAR_SRC = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular\src"
NEW_APPROACH_PATH = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular\examples\New approach"
RHINO_PLUGIN_PATH = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular"


def ensure_paths():
    """Put compas_singular, framefield and rhino_plugin on ``sys.path``.

    Idempotent, and safe to call at ANY time -- unlike :func:`setup_paths`.
    **Call this immediately before every late import**, because this module's
    body runs only once. Rhino's script engine keeps one interpreter alive
    across runs, so the second time a script does ``from ff_common import ...``
    the module is served from ``sys.modules`` with its body skipped, while
    ``sys.path`` has been re-initialised by the script engine in between. The
    module-level imports below survive that -- they are already bound -- but a
    late import inside a function does not, and fails with a bare
    ``ModuleNotFoundError`` that looks nothing like a path problem.
    """
    for path in (SINGULAR_SRC, NEW_APPROACH_PATH, RHINO_PLUGIN_PATH):
        if path not in sys.path:
            sys.path.insert(0, path)


def setup_paths():
    """:func:`ensure_paths`, plus a one-time purge of a stale ``compas_singular``.

    A ``compas_singular`` imported before the path was fixed stays in
    ``sys.modules`` and would be served again, so drop it first.

    **Module-body use only.** Purging ``sys.modules`` after this module has
    bound ``CoarsePseudoQuadMesh`` would let a later import rebuild the class
    object, and a mesh made by one would not be ``isinstance`` of the other --
    a far more confusing failure than the import error this exists to prevent.
    """
    for module in list(sys.modules):
        if module == "compas_singular" or module.startswith("compas_singular."):
            del sys.modules[module]
    ensure_paths()


setup_paths()

import rhinoscriptsyntax as rs          # noqa: E402
import scriptcontext as sc              # noqa: E402

from compas.tolerance import TOL        # noqa: E402

from compas_singular.datastructures import CoarsePseudoQuadMesh  # noqa: E402


# ----------------------------------------------------------------------
# layers
# ----------------------------------------------------------------------

ROOT = "TopologyProblem"

OUTER_LAYER = ROOT + "::InputBoundaries::Outer"
INNER_LAYER = ROOT + "::InputBoundaries::Inner"
GUIDES_LAYER = ROOT + "::InputBoundaries::Guides"

CROSSES_LAYER = ROOT + "::Field::Crosses"
SINGULARITIES_LAYER = ROOT + "::Field::Singularities"
SEPARATRICES_LAYER = ROOT + "::Field::Separatrices"

MESH_LAYER = ROOT + "::Skeleton::Mesh"
POLES_LAYER = ROOT + "::Skeleton::Poles"
EDIT_LAYER = ROOT + "::Skeleton::Edit"
DENSITY_LAYER = ROOT + "::Skeleton::Densities"

QUADMESH_LAYER = ROOT + "::QuadMesh"


def ensure_layer(path, color=None):
    """Create a ``::`` layer path, parents first, and return it."""
    parts = path.split("::")
    for i in range(len(parts)):
        name = "::".join(parts[:i + 1])
        if not rs.IsLayer(name):
            rs.AddLayer(name=parts[i],
                        parent="::".join(parts[:i]) if i else None,
                        color=color if i == len(parts) - 1 else None)
    return path


def clear_layer(path):
    if rs.IsLayer(path):
        guids = rs.ObjectsByLayer(path)
        if guids:
            rs.DeleteObjects(guids)


def add_to_layer(guid, path):
    rs.ObjectLayer(guid, layer=path)
    return guid


# ----------------------------------------------------------------------
# settings, in document user text
# ----------------------------------------------------------------------

SETTINGS_KEY = "ff.settings"

DEFAULT_SETTINGS = {
    #: Background triangulation spacing. NOT the quad size -- the field is
    #: solved on this. Finer means a better field and a slower solve; 0.3 on a
    #: 20 x 14 m plate is 2 438 vertices and 2.1 s.
    "spacing": 0.5,
    #: 'tangent' makes elements run ALONG the guides, 'perpendicular' across.
    "mode": "tangent",
    #: Target quad edge length for strips with no explicit density.
    "target_length": 0.5,
    #: Whether patch INTERIORS are integrated from the field. Real on the field
    #: route; measured to COST quality on a skeleton layout, whose patch edges
    #: do not follow the field -- aspect 1.99 -> 3.20 on a square with a cable.
    "field_aware": True,
}


def get_settings():
    blob = rs.GetDocumentUserText(SETTINGS_KEY)
    settings = dict(DEFAULT_SETTINGS)
    if blob:
        settings.update(json.loads(blob))
    return settings


def set_settings(settings):
    rs.SetDocumentUserText(SETTINGS_KEY, json.dumps(settings))
    return settings


# ----------------------------------------------------------------------
# curves -> point lists
# ----------------------------------------------------------------------

def curve_points(guid, max_edge):
    """A curve as a list of points: corners kept, curvature sampled.

    A polyline is taken at its own vertices, because sampling it would round the
    corners and a corner is exactly what the field has to see as a corner.
    Anything else -- an arc, a circle, a NURBS curve -- is divided by length,
    because every consumer treats the loop as straight segments between the
    points given, so an arc handed over as two points IS a chord.

    ``max_edge`` should be the background spacing or finer: sampled coarser than
    the background, the field never sees the wall it is meant to meet at a right
    angle.
    """
    curve = rs.coercecurve(guid)
    ok, polyline = curve.TryGetPolyline()
    if ok:
        points = [[p.X, p.Y, 0.0] for p in polyline]
    else:
        count = max(8, int(round(curve.GetLength() / max(max_edge, 1e-6))))
        points = []
        for t in curve.DivideByCount(count, True):
            p = curve.PointAt(t)
            points.append([p.X, p.Y, 0.0])
    if len(points) > 1 and _distance(points[0], points[-1]) < 1e-9:
        points = points[:-1]
    return points


def _distance(a, b):
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5


def polygon_area(points):
    """Shoelace, unsigned. Only used to tell the outer loop from the holes."""
    total = 0.0
    for i in range(len(points)):
        a, b = points[i], points[(i + 1) % len(points)]
        total += a[0] * b[1] - b[0] * a[1]
    return abs(total) * 0.5


def read_boundaries(spacing=None):
    """``(outer, inners, guides)`` as point lists, from the input layers."""
    if spacing is None:
        spacing = get_settings()["spacing"]

    outer_ids = rs.ObjectsByLayer(OUTER_LAYER)
    if not outer_ids:
        raise RuntimeError(
            "No outer boundary on '{}' -- run ff_01_boundaries first.".format(
                OUTER_LAYER))

    outer = curve_points(outer_ids[0], spacing)
    inners = [curve_points(guid, spacing)
              for guid in rs.ObjectsByLayer(INNER_LAYER) or []]
    guides = [curve_points(guid, spacing)
              for guid in rs.ObjectsByLayer(GUIDES_LAYER) or []]
    return outer, inners, guides


# ----------------------------------------------------------------------
# the coarse layout, on the document
# ----------------------------------------------------------------------

def read_coarse():
    """``(CoarsePseudoQuadMesh, poles)`` from ``Skeleton::Mesh`` + ``::Poles``.

    Rhino has no triangular mesh face -- it stores one as a quad with a repeated
    corner -- so the duplicate is dropped. A face that keeps it is degenerate
    and gets thrown out further down without explanation.
    """
    from compas_rhino.conversions import mesh_to_compas
    from compas_rhino.conversions import point_to_compas
    import compas_rhino as cr

    guids = rs.ObjectsByLayer(MESH_LAYER)
    if not guids:
        raise RuntimeError(
            "No coarse layout on '{}' -- run ff_03_coarse first.".format(
                MESH_LAYER))

    mesh = mesh_to_compas(cr.objects.find_object(guids[0]).Geometry)
    vertices, faces = mesh.to_vertices_and_faces()
    faces = [[v for i, v in enumerate(f) if v != f[i - 1]] for f in faces]

    poles = [list(point_to_compas(cr.objects.find_object(guid).Geometry))
             for guid in rs.ObjectsByLayer(POLES_LAYER) or []]

    coarse = CoarsePseudoQuadMesh.from_vertices_and_faces_with_poles(
        vertices, faces, poles)
    return coarse, poles


def bake_coarse(coarse):
    """Replace the layout on ``Skeleton::Mesh``, and its poles.

    The geometry IS the state the next step reads, so the old objects go. Poles
    are rewritten too: a stale pole point would make a later step register a
    collapsed corner that is no longer there.
    """
    from compas_rhino.conversions import point_to_rhino

    ensure_layer(MESH_LAYER, (0, 0, 0))
    ensure_layer(POLES_LAYER, (200, 0, 120))
    clear_layer(MESH_LAYER)
    clear_layer(POLES_LAYER)

    index = {v: i for i, v in enumerate(coarse.vertices())}
    vertices = [coarse.vertex_attributes(v, "xyz") for v in coarse.vertices()]
    faces = [[index[v] for v in coarse.face_vertices(f)] for f in coarse.faces()]
    guid = add_to_layer(rs.AddMesh(vertices, faces), MESH_LAYER)

    poles = coarse.poles() if hasattr(coarse, "poles") else []
    for vkey in poles:
        add_to_layer(rs.AddPoint(point_to_rhino(coarse.vertex_coordinates(vkey))),
                     POLES_LAYER)
    return guid


# ----------------------------------------------------------------------
# the field, solved once per session
# ----------------------------------------------------------------------

STICKY_KEY = "ff.decomposition"


def get_decomposition(force=False, verbose=True):
    """The live ``FieldDecomposition`` for the current inputs.

    Cached in ``sticky`` under a hash of the inputs, so a moved boundary curve
    misses the cache and re-solves rather than silently meshing the old outline.
    The solve is deterministic -- two builds give identical singularities and
    identical coarse vertices -- so a rebuild reproduces the field exactly and
    there is no drift to defend against.
    """
    ensure_paths()          # this module's body ran once; sys.path may not have
    from framefield.decomposition import FieldDecomposition
    import hashlib
    import time

    settings = get_settings()
    outer, inners, guides = read_boundaries(settings["spacing"])

    blob = json.dumps([[[round(c, 6) for c in p] for p in loop]
                       for loop in [outer] + inners + guides]
                      + [settings["spacing"], settings["mode"]])
    key = STICKY_KEY + "." + hashlib.md5(blob.encode("utf-8")).hexdigest()

    if not force:
        cached = sc.sticky.get(key)
        if cached is not None:
            if verbose:
                print("field: reusing this session's solve")
            sc.sticky[STICKY_KEY] = cached
            return cached

    start = time.time()
    decomposition = FieldDecomposition.from_boundary(
        outer, inners or None, guides=guides or None, mode=settings["mode"],
        target_length=settings["spacing"])
    if verbose:
        print("field: solved and traced in {:.1f} s".format(time.time() - start))

    sc.sticky[key] = decomposition
    sc.sticky[STICKY_KEY] = decomposition
    return decomposition


def print_warnings(warnings):
    """Print each distinct warning once.

    ``repair_notes`` accumulates for the life of the decomposition, and the
    decomposition outlives a command -- so a layout adopted in step 4 and again
    in step 6 reports its repair twice, which reads as two problems.
    """
    seen = set()
    for warning in warnings:
        if warning not in seen:
            seen.add(warning)
            print("  warn: " + warning)


# ----------------------------------------------------------------------
# densities, keyed geometrically
# ----------------------------------------------------------------------

DENSITY_KEY = "ff.densities"


def edge_key(a, b):
    """The identity a density is stored against: the edge's rounded midpoint.

    Order-independent by construction, and rounded the way the rest of the
    codebase rounds -- ``TOL.geometric_key``, 3 decimals, the same resolution
    ``from_polylines`` matches endpoints at.
    """
    mid = [(a[i] + b[i]) / 2.0 for i in range(3)]
    return TOL.geometric_key(mid)


def get_densities():
    blob = rs.GetDocumentUserText(DENSITY_KEY)
    return json.loads(blob) if blob else {}


def set_densities(densities):
    rs.SetDocumentUserText(DENSITY_KEY, json.dumps(densities))
    return densities


def strip_edge_map(coarse):
    """``{edge midpoint key: skey}`` for every edge of every strip.

    Built once and shared: this is the lookup that replaces the strip index
    everywhere it would otherwise be used across a bake.
    """
    out = {}
    for skey in coarse.strips():
        for u, v in coarse.strip_edges(skey):
            if u == v:
                continue
            out[edge_key(coarse.vertex_coordinates(u),
                         coarse.vertex_coordinates(v))] = skey
    return out


def strip_of_edge(coarse, a, b, lookup=None):
    """Which strip owns the edge between ``a`` and ``b``. ``None`` if no match."""
    lookup = strip_edge_map(coarse) if lookup is None else lookup
    return lookup.get(edge_key(a, b))


def apply_densities(coarse, target_length, densities=None, verbose=True):
    """Set every strip's density: the target length, then the stored overrides.

    Returns ``(applied, lost)`` -- how many overrides found their strip, and how
    many did not because the layout changed under them. A lost override is worth
    saying out loud: quietly meshing at the default is how a density the user set
    goes missing.
    """
    if densities is None:
        densities = get_densities()

    coarse.collect_strips()
    coarse.set_strips_density_target(target_length)

    lookup = strip_edge_map(coarse)
    applied = lost = 0
    for key, density in densities.items():
        skey = lookup.get(key)
        if skey is None:
            lost += 1
            continue
        coarse.set_strip_density(skey, int(density))
        applied += 1

    if verbose and (applied or lost):
        print("densities: {} override(s) applied, {} lost to a changed layout"
              .format(applied, lost))
    return applied, lost
