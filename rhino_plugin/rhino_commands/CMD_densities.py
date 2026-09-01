#! python3

# r: compas

"""**Step 5 -- how many elements across each strip.**

    reads   TopologyProblem::Skeleton::Mesh
    writes  document user text 'ff.densities'
    scratch TopologyProblem::Skeleton::Densities   pickable edges + labels

``pick`` an edge and give it a number, or set a ``target`` length that sizes
every strip that has no explicit number. A strip is a band of quads running
across the layout: setting one edge sets the whole band, because a quad mesh
cannot subdivide one patch without subdividing everything in line with it.

**The densities are stored against geometry, never against the strip index, and
that is the one trap in this whole design.** ``collect_strips`` builds its list
by popping ``self.edges()``, so indices follow vertex insertion order -- and a
layout that has been baked to Rhino and read back has whatever order the
geometry came back in. Measured: bake a 4-strip layout, read it back, **0 of 4
indices still point at the same strip**, permuted 0->2 1->3 2->0 3->1. A density
map keyed by index does not raise; it silently applies the wrong numbers to the
wrong bands. So each density is stored against the rounded midpoint of the edge
that was picked, and the strip owning that edge is looked up again on the other
side. Editing the layout in step 4 can move that midpoint, and an override that
no longer finds its strip is reported as lost rather than dropped in silence.
"""

#Temporary import of compas_singular development library
from CMD_start import import_compas_singular
import_compas_singular()

import json

import rhinoscriptsyntax as rs

from compas_rhino.conversions import point_to_rhino

from compas_singular.rhino.helpers.helpers import read_mesh, bake_mesh, clean_layer, read_coarse
clear_layer = clean_layer
from compas_singular.rhino.helpers.helpers import curve_to_polyline, read_boundary_loops, bake_edge_curves, read_boundaries
from compas_singular.rhino.coarse_curves import coarse_edges_to_curves, snap_corners_to_walls
from CMD_start import set_settings, get_settings
from compas.tolerance import TOL   

DENSITY_LAYER = "Densities"
settings = get_settings()

coarse, poles = read_coarse()

def get_densities():
    density_key = settings["density_key"]
    blob = rs.GetDocumentUserText(density_key)
    return json.loads(blob) if blob else {}


def set_densities(densities):
    density_key = settings["density_key"]
    rs.SetDocumentUserText(density_key, json.dumps(densities))
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


def edge_key(a, b):
    """The identity a density is stored against: the edge's rounded midpoint.

    Order-independent by construction, and rounded the way the rest of the
    codebase rounds -- ``TOL.geometric_key``, 3 decimals, the same resolution
    ``from_polylines`` matches endpoints at.
    """
    mid = [(a[i] + b[i]) / 2.0 for i in range(3)]
    return TOL.geometric_key(mid)



GUIDS = {}          # guid -> (u, v)
def draw(coarse):
    """Every coarse edge as a pickable line, plus one label per strip.

    The label sits at the midpoint of the strip's first edge and carries the
    density that strip will actually be meshed at -- the override if it has one,
    otherwise what the target length works out to. Showing the resolved number
    rather than the override is the point: a strip with no override is not
    unset, it is set by the target.
    """
    rs.AddLayer(DENSITY_LAYER, parent="Skeleton")
    clear_layer(DENSITY_LAYER)
    GUIDS.clear()

    overrides = get_densities()
    for u, v in coarse.edges():
        a = coarse.vertex_coordinates(u)
        b = coarse.vertex_coordinates(v)
        guid = rs.AddLine(point_to_rhino(a), point_to_rhino(b))
        rs.ObjectLayer(guid, layer=DENSITY_LAYER)
        GUIDS[guid] = (u, v)

    for skey in coarse.strips():
        edges = [(u, v) for u, v in coarse.strip_edges(skey) if u != v]
        if not edges:
            continue
        u, v = edges[len(edges) // 2]
        a = coarse.vertex_coordinates(u)
        b = coarse.vertex_coordinates(v)
        mid = [(a[i] + b[i]) / 2.0 for i in range(3)]
        density = coarse.get_strip_density(skey)
        marked = "*" if any(edge_key(coarse.vertex_coordinates(x),
                                     coarse.vertex_coordinates(y)) in overrides
                            for x, y in edges) else ""
        guid = rs.AddTextDot("{}{}".format(density, marked), point_to_rhino(mid))
        rs.ObjectLayer(guid, layer=DENSITY_LAYER)


def pick_strip(coarse, lookup):
    """``(skey, edge_key)`` for an edge the user picks. ``(None, None)`` if not."""
    guid = rs.GetObject(message="Pick an edge across the strip to set",
                        filter=rs.filter.curve, preselect=True)
    if not guid or guid not in GUIDS:
        return None, None
    u, v = GUIDS[guid]
    a = coarse.vertex_coordinates(u)
    b = coarse.vertex_coordinates(v)
    key = edge_key(a, b)
    return lookup.get(key), key


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


densities = get_densities()

applied, lost = apply_densities(coarse, settings["target_length"], densities)
if lost:
    print("  a lost override is one whose edge no longer exists -- usually "
          "because step 4 moved that corner. Re-pick those strips.")

lookup = strip_edge_map(coarse)
print("layout: {} patch(es), {} strip(s); target length {}".format(
    coarse.number_of_faces(), len(list(coarse.strips())),
    settings["target_length"]))

draw(coarse)

try:
    while True:
        option = (rs.GetString(message="Densities", defaultString="Finish",
                               strings=["Pick", "Target", "Clear", "Finish"])
                  or "finish").lower()
        if option == "pick":
            skey, key = pick_strip(coarse, lookup)
            if skey is None:
                print("Not an edge of the layout.")
                continue
            current = coarse.get_strip_density(skey)
            value = rs.GetInteger("Elements across this strip", current, 1)
            if value is None:
                continue
            densities[key] = int(value)
            set_densities(densities)
            coarse.set_strip_density(skey, int(value))
            draw(coarse)

        elif option == "target":
            value = rs.GetReal("Target quad edge length",
                               settings["target_length"], 1e-3)
            if value:
                settings["target_length"] = value
                set_settings(settings)
                # Re-resolve: the target sets every strip that has no override,
                # and the labels have to show the new numbers.
                apply_densities(coarse, value, densities, verbose=False)
                draw(coarse)

        elif option == "clear":
            densities = set_densities({})
            apply_densities(coarse, settings["target_length"], densities,
                            verbose=False)
            draw(coarse)
            print("all overrides cleared -- every strip follows the target.")

        else:
            break
finally:
    clear_layer(DENSITY_LAYER)

print("densities: {} strip(s) with an explicit number, the rest at target {}"
      .format(len(densities), settings["target_length"]))
