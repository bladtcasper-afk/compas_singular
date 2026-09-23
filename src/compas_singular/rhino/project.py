"""The Rhino project the ``CMD_`` commands share: settings, layers and the layout.

Everything a command used to import from another command lives here, so that no
command imports another. Three things:

* **Settings** -- :func:`get_settings` / :func:`set_settings`, a plain dict view
  of the session's :class:`~compas_singular.settings.Settings`. The ``resolve_*``
  functions turn a stored setting into what a solve needs.
* **Layers** -- ``ROOT`` and ``LAYER_DATA``, the one place the layer tree is
  written down. :func:`layer_path` gives a full path by its short name.
* **The layout** -- :func:`read_layout`, a copy of the session's coarse layout.

Everything a command keeps between runs is in the session
(:class:`~compas_singular.rhino.session.RhinoSession`), stored inside the
``.3dm``. The JSON side-cars that used to sit beside it are only read once, to
move an old document over (:mod:`compas_singular.rhino.legacy`).

**Importing this module does nothing.** It reads no document, writes no settings
and creates no layers. ``rhinoscriptsyntax`` is optional so the Rhino-free parts
(``resolve_relax``, ``resolve_symmetry``, ``resolve_densities``, the layer names)
import headless; anything that touches the document raises there.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas_singular.settings import Settings

try:
    import rhinoscriptsyntax as rs
except ImportError:
    rs = None


__all__ = [
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
    'read_layout',
    'layout_polylines',
]


def _require_rhino():
    if rs is None:
        raise RuntimeError("compas_singular.rhino.project: this function needs Rhino "
                           "(rhinoscriptsyntax is not importable here).")


# ==============================================================================
# settings
# ==============================================================================

#: Every setting at its default, as the dict :func:`get_settings` returns.
DEFAULT_SETTINGS = Settings().model_dump()


def get_settings():
    """The session's settings, as a dict a command can change and hand to :func:`set_settings`."""
    _require_rhino()
    from compas_singular.rhino.session import RhinoSession
    return RhinoSession.current().settings.model_dump()


def set_settings(settings):
    """Store ``settings`` (a dict) in the session, and record it. Returns them as stored."""
    _require_rhino()
    from compas_singular.rhino.session import RhinoSession
    session = RhinoSession.current()
    session.settings = Settings.model_validate(settings)
    session.record('Settings')
    return session.settings.model_dump()


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
    step 5 gave densities keeps them as they are. One that has none -- step 5
    never ran, its strips changed in step 4, or it was imported from an old
    drawing -- gets every strip set by ``settings["density_mode"]``: ``"length"`` sizes
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
# the layout
# ==============================================================================

def read_layout(verbose=True):
    """**The session's coarse layout, as a COPY.**

    A copy because every caller changes it -- densities, patterns, snapped
    corners -- and may then be cancelled. It goes back into the session only
    when the caller assigns it and records (see
    :class:`~compas_singular.rhino.session.RhinoSession`).

    Raises
    ------
    RuntimeError
        If the session holds no layout. The mesh drawn on ``Skeleton::Mesh`` is
        NOT read instead: it is vertices and faces only, a view of the layout. A
        document from before sessions has its drawn layout imported once
        (:mod:`compas_singular.rhino.legacy`).
    """
    _require_rhino()
    from compas_singular.rhino.session import RhinoSession

    layout = RhinoSession.current().coarse
    if layout is None:
        raise RuntimeError("This document holds no coarse layout. Run CMD_coarse_mesh "
                           "or CMD_read_coarse_mesh first.")
    coarse = layout.copy()
    if verbose:
        print("layout: {} patch(es), route {!r}".format(
            coarse.number_of_faces(), coarse.attributes.get("route", "unknown")))
    return coarse


def layout_polylines(coarse):
    """The polylines ``coarse``'s edges take their shape from (``shape_polylines``).

    A layout recorded before 2026-09-18 carries none; its polylines are still
    drawn on ``Skeleton::Polylines``, so they are read from there instead.
    """
    polylines = coarse.shape_polylines()
    if polylines:
        return polylines
    _require_rhino()
    from compas_singular.rhino.helpers import read_polylines
    layer = layer_path("Polylines")
    return read_polylines(layer) if rs.IsLayer(layer) else []
