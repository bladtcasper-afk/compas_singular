#! python3

# r: compas

"""**Step 2 -- solve the cross field and look at it.**

    reads   TopologyProblem::InputBoundaries::*
    writes  TopologyProblem::Field::Crosses        display only
            TopologyProblem::Field::Singularities  display only
            TopologyProblem::Field::Separatrices   display only

Nothing on these layers is ever read back. They exist because the field is the
one part of the pipeline with no visible output of its own, and every later
surprise -- a patch in a silly place, a layout that would not densify -- is
already visible here if you know what to look at.

**The singularities are the thing to look at.** Each is where a quad mesh must
put an irregular vertex: red for index +1, which becomes a valence-3 vertex,
blue for -1, which becomes valence-5. Their count and placement decide the whole
patch layout, so a spray of them along a curved wall means the layout is about
to be a mess -- and the usual cause is a sampled arc turning more than 45
degrees between points, which a cross field cannot tell apart from a corner.
Re-sample the boundary finer and they go away.

**The field is never serialized.** It is 92 KB at 0.5 spacing and 240 KB at 0.3,
against a solve that is deterministic and costs 0.4-2.1 s, so it is rebuilt on
demand and cached in ``scriptcontext.sticky`` for the session.
"""
from ff_common import CROSSES_LAYER
from ff_common import SEPARATRICES_LAYER
from ff_common import SINGULARITIES_LAYER
from ff_common import clear_layer
from ff_common import ensure_layer
from ff_common import get_decomposition
from ff_common import get_settings
from ff_common import print_warnings

import rhinoscriptsyntax as rs

from compas_rhino.conversions import point_to_rhino


#: One cross per this many background vertices. A 0.5-spacing background is
#: ~940 vertices, and 940 crosses is 1 880 lines of visual noise that also
#: makes the viewport crawl.
DECIMATE = 4


def draw_crosses(decomposition, decimate=DECIMATE):
    """Two crossed lines per sampled background vertex.

    Drawn at 40% of the background spacing so neighbouring crosses stay apart:
    at full spacing they touch and the field reads as a solid hatch.
    """
    ensure_layer(CROSSES_LAYER, (120, 120, 120))
    clear_layer(CROSSES_LAYER)

    field = decomposition.field
    mesh = field.background.mesh
    scale = decomposition.background.target_length * 0.4

    count = 0
    for i, vkey in enumerate(mesh.vertices()):
        if i % decimate:
            continue
        origin = mesh.vertex_coordinates(vkey)
        # Two of the four arms are enough -- the other two are their negatives,
        # and a cross drawn as four half-arms is the same picture at twice the
        # object count.
        for direction in field.directions(vkey)[:2]:
            a = [origin[k] - direction[k] * scale for k in range(3)]
            b = [origin[k] + direction[k] * scale for k in range(3)]
            guid = rs.AddLine(point_to_rhino(a), point_to_rhino(b))
            rs.ObjectLayer(guid, layer=CROSSES_LAYER)
            count += 1
    return count


def draw_singularities(decomposition):
    ensure_layer(SINGULARITIES_LAYER, (255, 0, 0))
    clear_layer(SINGULARITIES_LAYER)

    mesh = decomposition.field.background.mesh
    tally = {}
    for fkey, index in decomposition.field.singularities():
        guid = rs.AddPoint(point_to_rhino(mesh.face_centroid(fkey)))
        rs.ObjectLayer(guid, layer=SINGULARITIES_LAYER)
        rs.ObjectColor(guid, (255, 0, 0) if index > 0 else (0, 90, 255))
        tally[index] = tally.get(index, 0) + 1
    return tally


def draw_separatrices(decomposition):
    """The curves the coarse layout's edges are cut from."""
    ensure_layer(SEPARATRICES_LAYER, (160, 160, 160))
    clear_layer(SEPARATRICES_LAYER)

    count = 0
    for polyline in decomposition.decomposition_polylines():
        points = [point_to_rhino([p[0], p[1], p[2] if len(p) > 2 else 0.0])
                  for p in polyline]
        if len(points) > 1:
            guid = rs.AddPolyline(points)
            rs.ObjectLayer(guid, layer=SEPARATRICES_LAYER)
            count += 1
    return count


# ----------------------------------------------------------------------
# run
# ----------------------------------------------------------------------

settings = get_settings()
print("spacing {}, mode '{}'".format(settings["spacing"], settings["mode"]))

resolve = (rs.GetString(message="Re-solve from scratch?", defaultString="No",
                        strings=["No", "Yes"]) or "no").lower() == "yes"

decomposition = get_decomposition(force=resolve)

crosses = draw_crosses(decomposition)
tally = draw_singularities(decomposition)
ribs = draw_separatrices(decomposition)

print("background: {} vertices".format(
    decomposition.background.mesh.number_of_vertices()))
print("field: {} cross lines drawn (1 in {} vertices)".format(crosses, DECIMATE))
print("singularities: {}".format(
    ", ".join("index {:+d}: {}".format(k, v) for k, v in sorted(tally.items()))
    or "none"))
print("separatrices: {} polyline(s)".format(ribs))

# ``arm_mismatch`` is the warning that matters here. If the probe circles around
# a singularity do not produce exactly 4 - index launch directions, the tracer
# falls back to a partial set and the layout comes out with the wrong number of
# patches -- and nothing about the coarse mesh afterwards says that happened.
print_warnings(decomposition.warnings())

print("next: ff_03_coarse")
