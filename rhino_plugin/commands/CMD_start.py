#! python3

# r: compas


import rhinoscriptsyntax as rs
import json
import os
import tempfile

#Rhino settings
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
    for guid in rs.ObjectsByLayer("Guides") or []:
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


def get_settings():
    data = rs.GetDocumentUserText(SETTINGS_KEY)
    settings = dict(DEFAULT_SETTINGS)
    if data:
        settings.update(json.loads(data))
    return settings


def set_settings(settings):
    rs.SetDocumentUserText(SETTINGS_KEY, json.dumps(settings))
    return settings


# ----------------------------------------------------------------------
# where this document keeps what geometry cannot carry
# ----------------------------------------------------------------------
#
# A Rhino mesh is vertices and faces. Everything a ``CoarsePseudoQuadMesh``
# KNOWS -- ``strips``, ``strips_density``, ``polyedges``, ``face_pole``, and
# whatever a step chose to record in ``attributes`` -- is dropped by the bake,
# and a cross field cannot be baked at all. Both classes serialise themselves
# (``save_to_json`` / ``load_from_json``); what belongs HERE is the only
# document-shaped part of the question, which folder.
#
# This is plugin policy, not a library capability, which is why it lives beside
# the settings rather than in ``compas_singular``.


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


#: The two side-cars steps 3 and 6 pass between them.
COARSE_CACHE = "coarse"
FIELD_CACHE = "field"
DENSE_CACHE = "dense"


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
    from compas_singular.rhino.helpers.helpers import read_coarse

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

# NOTE do the same for layer names

settings = get_settings()
set_settings(settings)
print(settings)

#Rhino layer names

LAYERS_KEY = "layers"
LAYER_DATA = {
    "Problem": ("TopologyProblem", None),
    "Input": ("TopologyProblem::InputBoundaries", None),
    "Outer": ("TopologyProblem::InputBoundaries::Outer", (255, 0, 0)),
    "Inner": ("TopologyProblem::InputBoundaries::Inner", (0, 255, 0)),
    "Guides": ("TopologyProblem::InputBoundaries::Guides", (255, 127, 0)),
    "PointFeatures": ("TopologyProblem::PointFeatures", (0, 255, 0)),
    "Skeleton": ("TopologyProblem::Skeleton", None),
    "Poles": ("TopologyProblem::Skeleton::Poles", None),
    "Polylines": ("TopologyProblem::Skeleton::Polylines", None),
    "Mesh": ("TopologyProblem::Skeleton::Mesh", None),
    "EdgeCurves": ("TopologyProblem::Skeleton::EdgeCurves", (0, 120, 200)),
    "TempEdit": ("TopologyProblem::Skeleton::TempEdit", (0, 120, 200)),
    "Densities": ("TopologyProblem::Attributes::Densities", None),
    "Patterns": ("TopologyProblem::Attributes::Patterns", None),
    "QuadMesh": ("TopologyProblem::QuadMesh", None),
    "Dual": ("TopologyProblem::QuadMesh::Dual", None),
    "Smoothened": ("TopologyProblem::QuadMesh::Smoothened", None),
    "Area": ("TopologyProblem::QuadMesh::Smoothened::Area", None),
}

#Rhino file setup
SINGULAR_SRC = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular\src"


def ensure_paths():
    """Put compas_singular on ``sys.path``. Purges NOTHING.

    Safe to call at any time, unlike :func:`import_compas_singular` -- see the
    warning there. Call this immediately before a late import inside a function:
    this module's body runs once, while Rhino's script engine re-initialises
    ``sys.path`` between runs, so a late import can otherwise fail with a bare
    ``ModuleNotFoundError`` that looks nothing like a path problem.
    """
    import sys
    if SINGULAR_SRC not in sys.path:
        sys.path.insert(0, SINGULAR_SRC)


def import_compas_singular():
    """:func:`ensure_paths`, plus a one-time purge of a stale compas_singular.

    **Module-body use only, before any compas_singular import.** Rhino's script
    engine keeps its interpreter alive between runs, so a compas_singular
    imported before the ``sys.path`` fix stays in ``sys.modules`` and would
    otherwise be served again.

    **``framefield`` must be purged WITH everything else, and exempting it was
    a real bug** (2026-08-28). It imports ``Mesh`` and ``CoarsePseudoQuadMesh``
    at module level, so a ``framefield`` kept across a purge goes on holding the
    PREVIOUS run's class objects while the command that just imported
    ``compas_singular`` holds new ones. Nothing looks wrong until two of them
    meet: ``isinstance`` is silently false, and ``CMD_coarse_mesh`` failed with

        PicklingError: Can't pickle <class '...mesh.Mesh'>:
        it's not the same object as compas_singular.datastructures.mesh.mesh.Mesh

    which was a pickle detecting the mismatch rather than causing it.
    Re-importing ``framefield`` costs about a fifth of a second, and every Rhino
    command runs once.
    """
    import sys
    ensure_paths()
    for _mod in list(sys.modules):
        if _mod == "compas_singular" or _mod.startswith("compas_singular."):
            del sys.modules[_mod]
import_compas_singular()


# ----------------------------------------------------------------------
# running CMD_start AS a command
# ----------------------------------------------------------------------

def reset_project():
    """Empty the project layers and recreate them. This command's actual job.

    **Called only under the ``__name__`` guard below, and that guard is
    load-bearing rather than tidiness.** Every other ``CMD_`` file opens with
    ``from CMD_start import ...``, which executes this module's body -- so
    before the guard existed, the first ``CMD_`` command of a Rhino session
    DELETED the user's whole ``TopologyProblem`` layer as a side effect of an
    import, before they had asked for anything. It only ever looked survivable
    because Rhino keeps one interpreter alive across runs, so it fired once per
    session rather than once per command.

    The body below is unchanged from when it ran at import; only its
    reachability is.
    """
    #Setup of project folder
    project_folder = "TopologyProblem"
    if rs.IsLayer(project_folder):
        if rs.IsLayer("Default"):
            rs.CurrentLayer("Default")
        else:
            print("Will not be able to clear problem folder. Set another mayer as active.")

        sublayers = rs.LayerChildren(project_folder)

        for layer in sublayers:
            layer_objects = rs.ObjectsByLayer(layer)
            rs.PurgeLayer(layer)
        objects = rs.ObjectsByLayer(project_folder)
        rs.DeleteObjects(objects)
    else:
        rs.AddLayer(name=project_folder)

    #Setup of subfolders
    for name, color in LAYER_DATA.values():
        rs.AddLayer(name=name, color=color)

    # Say so. A reset that silently stops happening -- if a future Rhino ran
    # this file under a name other than "__main__" -- would look like a command
    # that did nothing, and the absence of this line is how you would spot it.
    print("CMD_start: project '{}' reset.".format(project_folder))
    print("note: this deletes every input, layout and mesh under '{}' -- nothing "
          "from an earlier session survives.".format(project_folder))
    print("next: CMD_boundary_selection to define a new domain.")
    return project_folder


if __name__ == "__main__":
    reset_project()

#Some general informational links
#Filter numbers: https://developer.rhino3d.com/api/rhinoscript/selection_methods/filterobjects.htm




