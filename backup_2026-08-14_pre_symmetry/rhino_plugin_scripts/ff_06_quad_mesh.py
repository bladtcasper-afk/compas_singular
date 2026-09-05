#! python3

# r: compas

"""**Step 6 -- densify the layout into the quad mesh.**

    reads   TopologyProblem::Skeleton::Mesh, ::Poles, InputBoundaries::*
            document user text 'ff.settings', 'ff.densities'
    writes  TopologyProblem::QuadMesh

**Why the layout goes back through ``edit_coarse`` before being densified.**
The densifier maps each coarse EDGE to the curve it should follow, and that map
is keyed by the vertex keys of the decomposition's own layout. A layout read
back from Rhino has different keys, so handing it over directly would silently
fall back to straight chords for every edge -- and on a curved domain the mesh
then loses the area between chord and wall outright: 65% coverage on an ellipse,
74% on a disc. Adopting it first puts the read-back layout IN the decomposition,
so every edge finds its separatrix, its boundary arc, or its warp.

**``field_aware`` is a real feature on the field route and a marginal one on a
skeleton layout.** It integrates patch INTERIORS from the field instead of
filling them with a bilinear blend of their own boundaries, which is the only
way a cable reaches the inside of a patch it does not split. But a skeleton
layout's edges do not follow the field, so interior relaxation fights the patch
boundaries it was handed: measured on a square with a hole and a cable, aspect
ratio 1.99 -> 3.20 and min angle 35.5 -> 30.3 degrees, for an alignment gain of
0.9 degrees. Turn it off for skeleton layouts.

Densities come from step 5 and are re-resolved here against the current layout,
so an override whose edge no longer exists is reported rather than lost.
"""
from ff_common import QUADMESH_LAYER
from ff_common import apply_densities
from ff_common import clear_layer
from ff_common import ensure_layer
from ff_common import get_decomposition
from ff_common import get_settings
from ff_common import print_warnings
from ff_common import read_coarse
from ff_common import set_settings

import rhinoscriptsyntax as rs
import scriptcontext as sc

from compas_rhino.conversions import mesh_to_rhino


settings = get_settings()
coarse, poles = read_coarse()
decomposition = get_decomposition()

answer = (rs.GetString(message="Integrate patch interiors from the field?",
                       defaultString="Yes" if settings["field_aware"] else "No",
                       strings=["Yes", "No"]) or "yes").lower()
settings["field_aware"] = answer == "yes"
set_settings(settings)

# Adopt the read-back layout. ``strict=False`` because this step should still
# produce a mesh from a layout that is merely imperfect -- step 4 is where a bad
# edit is supposed to be refused, loudly, while it can still be fixed.
decomposition.edit_coarse(coarse, strict=False)
coarse = decomposition.mesh

applied, lost = apply_densities(coarse, settings["target_length"])
if lost:
    print("  re-pick those strips in ff_05 -- they are meshing at the target.")

decomposition.field_aware = settings["field_aware"]
dense = decomposition.densify()

ensure_layer(QUADMESH_LAYER, (0, 0, 0))
clear_layer(QUADMESH_LAYER)
guid = sc.doc.Objects.AddMesh(mesh_to_rhino(dense))
rs.ObjectLayer(guid, layer=QUADMESH_LAYER)

print("mesh: {} faces from {} patch(es), field_aware={}".format(
    dense.number_of_faces(), coarse.number_of_faces(), settings["field_aware"]))

quality = decomposition.quality(mesh=dense)
print("quality: min angle {:.1f}, max angle {:.1f}, aspect {:.2f}, "
      "coverage {:.2f}".format(quality["min_angle"], quality["max_angle"],
                               quality["aspect_max"], quality["coverage"]))
# Coverage is the one to read first: it is the meshed area over the domain area,
# and anything below 1 means the mesh does not fill the shape -- a chorded
# boundary, or a patch that went missing.
if quality["coverage"] < 0.99:
    print("  coverage below 1: the mesh does not fill the domain. Usually a "
          "boundary edge densified as a chord -- check the layout's corners "
          "actually sit on the walls.")

edges = decomposition.edit_notes.get("edges")
if edges:
    # How each coarse edge found its curve. 'chord' is the one that costs area.
    print("edges: {}".format(edges))

print_warnings(decomposition.warnings())

print("baked to '{}'".format(QUADMESH_LAYER))
