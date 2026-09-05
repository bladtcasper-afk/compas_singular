#! python3

# r: compas

"""**Step 1 -- the domain.** Pick the outer boundary, the holes, and the guides.

    writes  TopologyProblem::InputBoundaries::Outer     one closed curve
            TopologyProblem::InputBoundaries::Inner     holes
            TopologyProblem::InputBoundaries::Guides    cables / force lines
            document user text 'ff.settings'            spacing, mode

Selected curves are COPIED onto these layers, so the drawing you picked from is
left alone and re-running the step cannot slowly eat it.

Guides are what make the field do something a mesher would not do by itself: a
cable across a plain square forces no topology at all, so without a guide there
is nothing for the field to align to and every route gives the same grid. Their
``mode`` is a design decision and never a default -- ``tangent`` makes elements
run ALONG the guide, ``perpendicular`` across it.

The spacing set here is the BACKGROUND TRIANGULATION spacing, which is what the
field is solved on -- not the quad size, which is step 6. Finer gives a better
field and a slower solve: on a 20 x 14 m plate, 1.0 is 258 vertices and 0.2 s,
0.3 is 2 438 vertices and 2.1 s.
"""
from ff_common import GUIDES_LAYER
from ff_common import INNER_LAYER
from ff_common import OUTER_LAYER
from ff_common import clear_layer
from ff_common import curve_points
from ff_common import ensure_layer
from ff_common import get_settings
from ff_common import polygon_area
from ff_common import set_settings

import rhinoscriptsyntax as rs

from compas.geometry import Polygon
from compas.geometry import is_polygon_in_polygon_xy


def closed_filter(rhobj, geometry, component_index):
    return bool(rs.IsCurveClosed(rhobj))


def pick(message, multiple=False, closed=True):
    kwargs = dict(filter=rs.filter.curve, preselect=False, select=True)
    if closed:
        kwargs["custom_filter"] = closed_filter
    if multiple:
        return rs.GetObjects(message=message, **kwargs) or []
    guid = rs.GetObject(message=message, **kwargs)
    return [guid] if guid else []


def place(guids, layer):
    """Copy onto ``layer`` -- never move: the source drawing is not ours."""
    copies = rs.CopyObjects(guids) or []
    for guid in copies:
        rs.ObjectLayer(guid, layer=layer)
    return copies


def set_outer():
    guids = pick("Pick ONE closed curve as the outer boundary")
    if guids:
        clear_layer(OUTER_LAYER)
        place(guids, OUTER_LAYER)


def set_inner(mode):
    if mode == "new":
        clear_layer(INNER_LAYER)
    if mode == "delete":
        def on_inner(rhobj, geometry, component_index):
            return rs.ObjectLayer(rhobj) == INNER_LAYER
        guids = rs.GetObjects(message="Select holes to remove",
                              filter=rs.filter.curve, custom_filter=on_inner)
        if guids:
            rs.DeleteObjects(guids)
        return
    place(pick("Pick closed curves as holes", multiple=True), INNER_LAYER)


def set_guides():
    # Open curves are fine and usual here -- a cable runs across the domain.
    place(pick("Pick guide curves -- cables, force lines", multiple=True,
               closed=False), GUIDES_LAYER)


def validate():
    """Drop holes that are not inside the outer boundary.

    Cheap to check here and expensive to find later: a hole outside the domain
    survives the triangulation and comes back as an unexplained missing patch.
    """
    outer_ids = rs.ObjectsByLayer(OUTER_LAYER)
    if not outer_ids:
        print("No outer boundary set.")
        return

    spacing = get_settings()["spacing"]
    outer = Polygon([p for p in curve_points(outer_ids[0], spacing)])
    for guid in rs.ObjectsByLayer(INNER_LAYER) or []:
        inner = Polygon([p for p in curve_points(guid, spacing)])
        if not is_polygon_in_polygon_xy(outer, inner):
            rs.DeleteObjects(guid)
            print("A hole was not inside the outer boundary -- removed.")


# ----------------------------------------------------------------------
# run
# ----------------------------------------------------------------------

ensure_layer(OUTER_LAYER, (255, 0, 0))
ensure_layer(INNER_LAYER, (0, 255, 0))
ensure_layer(GUIDES_LAYER, (0, 120, 255))

settings = get_settings()

while True:
    # rs.GetString returns the option lowercased, whatever case it was offered in.
    option = (rs.GetString(message="Domain",
                           defaultString="Continue",
                           strings=["Outer", "Inner", "Guides", "Spacing",
                                    "Mode", "Continue"]) or "continue").lower()

    if option == "outer":
        set_outer()
    elif option == "inner":
        how = (rs.GetString(message="Holes", defaultString="Add",
                            strings=["Add", "New", "Delete"]) or "add").lower()
        set_inner(how)
    elif option == "guides":
        set_guides()
    elif option == "spacing":
        value = rs.GetReal("Background spacing (NOT the quad size)",
                           settings["spacing"], 1e-3)
        if value:
            settings["spacing"] = value
            set_settings(settings)
    elif option == "mode":
        answer = rs.GetString(message="Elements run ALONG or ACROSS the guides?",
                              defaultString=settings["mode"],
                              strings=["tangent", "perpendicular"])
        if answer:
            settings["mode"] = answer.lower()
            set_settings(settings)
    else:
        break

validate()

outer_ids = rs.ObjectsByLayer(OUTER_LAYER)
if not outer_ids:
    raise RuntimeError("No outer boundary selected.")

outer = curve_points(outer_ids[0], settings["spacing"])
inners = [curve_points(g, settings["spacing"])
          for g in rs.ObjectsByLayer(INNER_LAYER) or []]
guides = rs.ObjectsByLayer(GUIDES_LAYER) or []

print("domain: outer {} points, area {:.1f}; {} hole(s); {} guide(s)".format(
    len(outer), polygon_area(outer), len(inners), len(guides)))
print("spacing {}, mode '{}'".format(settings["spacing"], settings["mode"]))
print("next: ff_02_field")
