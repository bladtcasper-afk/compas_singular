"""The Rhino project the ``CMD_`` commands share: settings, layers and side-cars.

Everything a command used to import from another command lives here, so that no
command imports another. Three things:

* **Settings** -- ``DEFAULT_SETTINGS``, stored on the document as user text under
  ``SETTINGS_KEY`` so they travel with the ``.3dm``. :func:`get_settings` /
  :func:`set_settings` are the only way to read and write them, and the
  ``resolve_*`` functions turn a stored setting into what a solve needs.
* **Layers** -- ``ROOT`` and ``LAYER_DATA``, the one place the layer tree is
  written down. :func:`layer_path` gives a full path by its short name.
* **Side-cars** -- the JSON files beside the ``.3dm`` that carry what a bake
  cannot: :func:`cache_path`, the ``*_CACHE`` names and :func:`read_layout`.

**Importing this module does nothing.** It reads no document, writes no settings
and creates no layers. ``rhinoscriptsyntax`` is optional so the Rhino-free parts
(``resolve_relax``, ``resolve_symmetry``, ``resolve_densities``, the layer names)
import headless; anything that touches the document raises there.

**What is NOT here: the ``sys.path`` bootstrap.** A module inside
``compas_singular`` cannot put ``compas_singular`` on the path, and the purge of a
stale copy would delete this module. Each command carries that block itself until
``compas_singular`` is installed into Rhino's Python.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import json
import os
import tempfile

try:
    import rhinoscriptsyntax as rs
except ImportError:
    rs = None


__all__ = [
    'SETTINGS_KEY',
    'DEFAULT_SETTINGS',
    'get_settings',
    'set_settings',
    'resolve_relax',
    'resolve_symmetry',
    'resolve_densities',
    'has_closed_guide',
    'ROOT',
    'LAYER_DATA',
    'layer_path',
    'ensure_layers',
    'COARSE_CACHE',
    'FIELD_CACHE',
    'DENSE_CACHE',
    'cache_directory',
    'cache_path',
    'read_layout',
]


def _require_rhino():
    if rs is None:
        raise RuntimeError("compas_singular.rhino.project: this function needs Rhino "
                           "(rhinoscriptsyntax is not importable here).")


# ==============================================================================
# settings
# ==============================================================================

SETTINGS_KEY = "settings"
DEFAULT_SETTINGS = {
    #: Background triangulation spacing. NOT the quad size -- the field is
    #: solved on this. Finer means a better field and a slower solve; 0.3 on a
    #: 20 x 14 m plate is 2 438 vertices and 2.1 s.
    "triangulation_spacing": 0.5,
    #: 'tangent' makes elements run ALONG the guides, 'perpendicular' across.
    "guide_allignment": "tangent",
    #: Target quad edge length for strips with no explicit density.
    "target_length": 0.5,
    #: Elements across every strip with no explicit density, when
    #: ``density_mode`` is ``"density"``.
    "target_density": 5,
    #: Which global rule step 5 last applied: ``"length"`` or ``"density"``.
    #: Only consulted when a layout has no densities saved on it.
    "density_mode": "length",
    #: Whether patch INTERIORS are integrated from the field. Real on the field
    #: route; measured to COST quality on a skeleton layout, whose patch edges
    #: do not follow the field -- aspect 1.99 -> 3.20 on a square with a cable.
    "field_aware": True,
    "density_key": "density_key",
    "symmetry": "auto",
    "relax": "auto",
}


def get_settings():
    """The document's settings, over the defaults. A key the document lacks gets its default."""
    _require_rhino()
    data = rs.GetDocumentUserText(SETTINGS_KEY)
    settings = dict(DEFAULT_SETTINGS)
    if data:
        settings.update(json.loads(data))
    return settings


def set_settings(settings):
    """Store ``settings`` on the document. Returns them."""
    _require_rhino()
    rs.SetDocumentUserText(SETTINGS_KEY, json.dumps(settings))
    return settings


def resolve_relax(settings, guides):
    """Whether to solve with diffusion + normalisation, for THIS input."""
    value = settings.get("relax", "auto")
    if isinstance(value, str):
        value = value.strip().lower()
        if value == "auto":
            return bool(guides)
        return value in ("on", "yes", "true", "1")
    return bool(value)


def has_closed_guide():
    """Is any curve on the Guides layer closed? A closed guide is the one case
    measured where relaxation makes the layout WORSE (13 patches -> 1)."""
    _require_rhino()
    for guid in rs.ObjectsByLayer(layer_path("Guides")) or []:
        if rs.IsCurveClosed(guid):
            return True
    return False


def resolve_symmetry(settings):
    """The symmetry group to solve under, for THIS document.

    ``'auto'`` detects it from outer + holes + guides; JSON ``null`` disables it
    and takes the pre-2026-08-14 route. Kept next to ``resolve_relax`` because
    the three commands that build a decomposition MUST agree -- steps 4 and 6
    re-solve, and a different group there is a different layout under the
    user's edits.
    """
    value = settings.get("symmetry", "auto")
    if isinstance(value, str):
        value = value.strip().lower()
        if value in ("", "none", "off", "no", "false", "0"):
            return None
        return "auto"
    return None if not value else "auto"


def resolve_densities(coarse, settings, verbose=True):
    """The layout's saved densities if every strip has one; otherwise the global rule.

    Shared by step 5 (``CMD_densities``) and step 6 (``CMD_quad_mesh``). A layout
    from the side-car that step 5 saved keeps its densities as they are. One read
    from the baked mesh, saved before step 5 ran, or whose strips changed in step
    4 gets every strip set by ``settings["density_mode"]``: ``"length"`` sizes
    each strip from ``target_length``, ``"density"`` gives them all
    ``target_density``. Densities picked per strip are not recovered then -- they
    were keyed by strips that may no longer exist.

    Returns
    -------
    bool
        True when the densities had to be re-derived.
    """
    if coarse.has_densities():
        return False
    coarse.attributes['strips_density'] = {}
    if settings.get("density_mode") == "density":
        coarse.set_strips_density(int(settings["target_density"]))
    else:
        coarse.set_strips_density_target(settings["target_length"])
    if verbose:
        print("densities: none saved on the layout -- every strip set by the "
              "target {}".format(settings.get("density_mode", "length")))
    return True


# ==============================================================================
# layers
# ==============================================================================

ROOT = "TopologyProblem"

#: Short name -> (full layer path, colour or None). ``CMD_start`` creates every
#: one of these, parents first, in this order.
LAYER_DATA = {
    "Problem": (ROOT, None),
    "Input": (ROOT + "::InputBoundaries", None),
    "Outer": (ROOT + "::InputBoundaries::Outer", (255, 0, 0)),
    "Inner": (ROOT + "::InputBoundaries::Inner", (0, 255, 0)),
    "Guides": (ROOT + "::InputBoundaries::Guides", (255, 127, 0)),
    "PointFeatures": (ROOT + "::PointFeatures", (0, 255, 0)),
    "Skeleton": (ROOT + "::Skeleton", None),
    "Poles": (ROOT + "::Skeleton::Poles", None),
    "Polylines": (ROOT + "::Skeleton::Polylines", None),
    "Mesh": (ROOT + "::Skeleton::Mesh", None),
    "EdgeCurves": (ROOT + "::Skeleton::EdgeCurves", (0, 120, 200)),
    "TempEdit": (ROOT + "::Skeleton::TempEdit", (0, 120, 200)),
    "Densities": (ROOT + "::Attributes::Densities", None),
    "Patterns": (ROOT + "::Attributes::Patterns", None),
    "QuadMesh": (ROOT + "::QuadMesh", None),
    "Dual": (ROOT + "::QuadMesh::Dual", None),
    "Smoothened": (ROOT + "::QuadMesh::Smoothened", None),
    "Area": (ROOT + "::QuadMesh::Smoothened::Area", None),
}


def layer_path(name):
    """The full ``::`` path of a project layer, by its short name in ``LAYER_DATA``.

    Use the full path, not the short name, with ``rhinoscriptsyntax``: a bare name
    resolves with FindName, which returns the FIRST layer of that name anywhere in
    the document.
    """
    return LAYER_DATA[name][0]


def ensure_layers():
    """Create every project layer that does not exist yet. Deletes nothing."""
    _require_rhino()
    for name, color in LAYER_DATA.values():
        if not rs.IsLayer(name):
            rs.AddLayer(name=name, color=color)


# ==============================================================================
# side-cars: where this document keeps what geometry cannot carry
# ==============================================================================
#
# A Rhino mesh is vertices and faces. Everything a ``CoarsePseudoQuadMesh``
# KNOWS -- ``strips``, ``strips_density``, ``polyedges``, ``face_pole``, and
# whatever a step chose to record in ``attributes`` -- is dropped by the bake,
# and a cross field cannot be baked at all. Both classes serialise themselves
# (``save_to_json`` / ``load_from_json``); what belongs HERE is the only
# document-shaped part of the question, which folder.

#: The side-cars steps 3 to 6 pass between them.
COARSE_CACHE = "coarse"
FIELD_CACHE = "field"
DENSE_CACHE = "dense"


def cache_directory(create=True):
    """The folder this document's side-cars live in.

    **Per DOCUMENT, not per script.** A ``cache`` folder next to the command
    files would be shared by every model Rhino opens, so two open documents
    would silently overwrite each other's layout. This is a hidden folder beside
    the ``.3dm`` instead, so it travels when the model is copied with its folder
    and is obviously absent when it is not.

    An UNSAVED document has no path and falls back to the system temp directory.
    Nothing written there survives a machine restart; save the model first if
    that matters.
    """
    _require_rhino()
    folder = None
    name = rs.DocumentName()
    path = rs.DocumentPath()
    if name and path:
        folder = os.path.join(path, "." + os.path.splitext(name)[0] + ".compas_singular")
    if folder is None:
        folder = os.path.join(tempfile.gettempdir(), "compas_singular", "unsaved")
    if create and not os.path.isdir(folder):
        os.makedirs(folder)
    return folder


def cache_path(name, create=True):
    """``<this document's cache>/<name>.json``. Nothing is read or written."""
    return os.path.join(cache_directory(create=create), name + ".json")


def read_layout(verbose=True):
    """**The coarse layout, with everything it knows.** ``(coarse, poles, source)``.

    Two copies of the layout exist and they answer different questions. The
    BAKED mesh on ``Skeleton::Mesh`` is what the user sees and edits, and it is
    vertices and faces at single precision. The side-car JSON written beside it
    carries the same layout WITH its ``attributes`` -- strips, ``strips_density``,
    ``polyedges``, ``face_pole``, the route that built it -- none of which
    survives a bake.

    So: prefer the side-car, but only after checking it still describes what is
    on the layer. They diverge whenever something outside these commands moved
    the mesh, and meshing a layout the user cannot see is worse than losing the
    attributes.

    The check is the ROUNDED VERTEX SET, not vertex order: a bake renumbers, so
    order proves nothing. Rounded because ``rs.AddMesh`` stores ``Point3f`` and
    a corner read back is a single-precision copy of the one that was written --
    about 2e-6 out on a 20-unit plate, which no exact comparison survives.

    Returns
    -------
    (CoarsePseudoQuadMesh, list, str)
        ``source`` is ``'cache'`` or ``'document'``, so a caller can say which
        it got.
    """
    from compas.tolerance import TOL
    from compas_singular.datastructures import CoarsePseudoQuadMesh
    from compas_singular.rhino.helpers import read_coarse

    baked, poles = read_coarse()
    cached = CoarsePseudoQuadMesh.load_from_json(
        cache_path(COARSE_CACHE, create=False), default=None)

    if cached is None:
        if verbose:
            print("layout: from the document -- no side-car yet, so strips and "
                  "densities are re-derived")
        return baked, poles, "document"

    def gkeys(mesh):
        return set(TOL.geometric_key(mesh.vertex_coordinates(v)) for v in mesh.vertices())

    if gkeys(cached) != gkeys(baked):
        if verbose:
            print("layout: from the document -- the side-car describes a "
                  "DIFFERENT layout ({} vs {} corners). Something moved the "
                  "baked mesh outside these commands; re-run step 3 or 4 to "
                  "bring them back together."
                  .format(cached.number_of_vertices(), baked.number_of_vertices()))
        return baked, poles, "document"

    if verbose:
        print("layout: from the side-car -- {} patch(es), route {!r}".format(
            cached.number_of_faces(), cached.attributes.get("route", "unknown")))
    return cached, poles, "cache"
