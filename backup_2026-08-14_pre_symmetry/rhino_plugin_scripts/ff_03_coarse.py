#! python3

# r: compas

"""**Step 3 -- cut the coarse patch layout.**

    reads   TopologyProblem::InputBoundaries::*
    writes  TopologyProblem::Skeleton::Mesh    the layout, as a Rhino mesh
            TopologyProblem::Skeleton::Poles   collapsed corners, if any

Two routes, and they are not interchangeable:

* **field** -- patches are cut by the separatrices of the cross field, so their
  edges follow it and a guide curve actually steers the layout. This is the
  route the rest of the pipeline is built around.
* **skeleton** -- the original ``compas_singular`` medial-axis decomposition.
  It ignores the field entirely. Worth having because it is robust on domains
  where the separatrix network fails to close, and because it is the baseline
  the field route is measured against: over 18 domains, field won 15 and
  skeleton 1.

**The skeleton route needs a densely sampled outline and the field route does
not.** ``SkeletonDecomposition`` builds its medial axis from the circumcentres
of a triangulation whose only vertices are the boundary points, so a
corners-only square raises ``KeyError: popitem(): dictionary is empty``. The
field route builds its own interior instead. Both are fed the outline sampled at
the background spacing here, which satisfies the first and does no harm to the
second.

A layout is baked as a MESH, not as polylines: polylines duplicate every shared
corner, and step 4 has to be able to drag one corner and move all the patches
that meet there.
"""
from ff_common import MESH_LAYER
from ff_common import bake_coarse
from ff_common import get_decomposition
from ff_common import get_settings
from ff_common import print_warnings
from ff_common import read_boundaries

import rhinoscriptsyntax as rs

from compas_singular.algorithms import SkeletonDecomposition
from compas_singular.algorithms import boundary_triangulation


def from_field():
    decomposition = get_decomposition()
    coarse = decomposition.decomposition_mesh()

    # ``quad_mesh`` ALWAYS returns a mesh, and so does this: on a domain whose
    # separatrix network did not close, what comes back is a fallback that
    # covers the shape and ignores the field. Nothing about the mesh says which
    # happened, so the warnings are the only place it is visible.
    discarded = any("DISCARDED" in w for w in decomposition.warnings())
    return coarse, decomposition.warnings(), discarded


def from_skeleton():
    settings = get_settings()
    outer, inners, guides = read_boundaries(settings["spacing"])
    if guides:
        print("note: the skeleton route ignores guide curves entirely.")

    trimesh = boundary_triangulation(outer_boundary=outer,
                                     inner_boundaries=inners,
                                     point_features=[])
    decomposition = SkeletonDecomposition.from_mesh(trimesh)
    return decomposition.decomposition_mesh([]), [], False


# ----------------------------------------------------------------------
# run
# ----------------------------------------------------------------------

route = (rs.GetString(message="Layout route", defaultString="Field",
                      strings=["Field", "Skeleton"]) or "field").lower()

if route == "skeleton":
    coarse, warnings, discarded = from_skeleton()
else:
    coarse, warnings, discarded = from_field()

sides = {}
for fkey in coarse.faces():
    n = len(coarse.face_vertices(fkey))
    sides[n] = sides.get(n, 0) + 1

bake_coarse(coarse)

print("layout ({} route): {} patch(es), {} corner(s), {} pole(s)".format(
    route, coarse.number_of_faces(), coarse.number_of_vertices(),
    len(coarse.poles()) if hasattr(coarse, "poles") else 0))
# 4 is a patch and 3 is a pseudo-quad pole; anything else means the layout was
# repaired into shapes nobody drew.
print("sides: {}".format(sides))
print_warnings(warnings)

if discarded:
    print("The field layout was DISCARDED and a fallback is on the layer in its "
          "place. It covers the domain but ignores the field, and it is not a "
          "layout worth hand-editing -- try a finer spacing in ff_01, or the "
          "skeleton route.")

print("baked to '{}'".format(MESH_LAYER))
print("next: ff_04_edit_coarse (optional), then ff_05_densities")
