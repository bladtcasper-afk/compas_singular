#! python3

# r: compas

"""**Densify the coarse layout into the quad mesh.**

    reads   TopologyProblem::Skeleton::{Mesh, Poles, Polylines}
            TopologyProblem::InputBoundaries::{Outer, Inner}
    writes  TopologyProblem::QuadMesh
            TopologyProblem::Skeleton::EdgeCurves

**A coarse edge is a straight chord, and the coarse layout has to stay that
way.** Strips, densities, poles and ``add_strip`` are all defined on the
topological quad graph, so the SHAPE of an edge cannot live in the layout -- it
is handed to ``densification`` separately, as one polyline per edge. Without it
every edge densifies as the chord between its corners and the mesh loses the
area between chord and wall outright: measured on the field route, 65% coverage
on an ellipse and 74% on a disc, and on a radius-5 disc the dense boundary sits
1.46 units inside the circle.

That mapping used to be dropped here. The layout is baked as a Rhino mesh, and a
Rhino mesh is vertices and faces -- it cannot carry a per-edge polyline -- so
this step read the layout back and called a bare ``densification()``.
``coarse_edges_to_curves`` rebuilds the mapping from the document instead: the
walls come from the input curves on ``InputBoundaries``, the interior branches
from ``Skeleton::Polylines``. Nothing is cached and nothing is keyed by vertex
index, so it survives a save, a reopen, an edit in CMD_edit_coarse_mesh, and the
single-precision round trip ``rs.AddMesh`` puts every corner through.

Read the ``coarse edges`` tally it prints. ``chord`` is the count that costs
area, and on a curved domain it should be interior edges only.
"""

#Temporary import of compas_singular development library
from CMD_start import import_compas_singular
import_compas_singular()

import rhinoscriptsyntax as rs
import Rhino

import compas_rhino as cr
from compas_rhino.conversions import curve_to_compas_polyline, curve_to_compas_circle, mesh_to_compas, point_to_rhino, point_to_compas

from math import pi

import json

import compas
from compas.tolerance import TOL  
from compas.datastructures.mesh.mesh import Mesh
from compas.datastructures.mesh.duality import mesh_dual
from compas.geometry import Point, Polyline, Line, Polygon, is_polygon_in_polygon_xy
from compas.scene import Scene
from compas_singular.datastructures import CoarseQuadMesh, QuadMesh, PseudoQuadMesh
from compas_singular.algorithms import boundary_triangulation, SkeletonDecomposition
from compas_singular.datastructures import CoarsePseudoQuadMesh
from compas_singular.rhino.helpers.helpers import bake_mesh, bake_edge_curves, read_boundary_loops, read_polylines, read_coarse, read_boundaries, clean_layer
from compas_singular.rhino.coarse_curves import coarse_edges_to_curves, snap_corners_to_walls
from CMD_start import get_settings, set_settings
from CMD_start import get_decomposition

settings = get_settings()
coarse, poles = read_coarse()


outer, inners, guides, point_features = read_boundaries(spacing=settings["triangulation_spacing"])
# Solved once per document and reused across steps 3, 4 and 6 -- and across
# Rhino restarts. ``CMD_start.get_decomposition`` owns the settings, the
# resolvers and the staleness rules; see its docstring. The object is
# reconstructed per call, so ``edit_coarse`` below cannot leak into step 4.
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


def get_densities():
    density_key = settings["density_key"]
    blob = rs.GetDocumentUserText(density_key)
    return json.loads(blob) if blob else {}


def set_densities(densities):
    density_key = settings["density_key"]
    rs.SetDocumentUserText(density_key, json.dumps(densities))
    return densities


def edge_key(a, b):
    """The identity a density is stored against: the edge's rounded midpoint.

    Order-independent by construction, and rounded the way the rest of the
    codebase rounds -- ``TOL.geometric_key``, 3 decimals, the same resolution
    ``from_polylines`` matches endpoints at.
    """
    mid = [(a[i] + b[i]) / 2.0 for i in range(3)]
    return TOL.geometric_key(mid)


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

applied, lost = apply_densities(coarse, settings["target_length"])
if lost:
    print("  re-pick those strips -- they are meshing at the target.")

decomposition.field_aware = settings["field_aware"]
dense = decomposition.densify()

layer = rs.AddLayer("QuadMesh", parent="TopologyProblem")
QUADMESH_LAYER = layer

clean_layer(QUADMESH_LAYER, clean_sublayers=True)
bake_mesh(dense, QUADMESH_LAYER)


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

for warning in decomposition.warnings():
    print("warning: " + warning)

print("baked to '{}'".format(QUADMESH_LAYER))
