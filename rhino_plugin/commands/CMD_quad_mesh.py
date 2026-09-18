#! python3

# r: compas

"""**Densify the coarse layout into the quad mesh.**

    reads   TopologyProblem::Skeleton::{Mesh, Poles, Polylines}
            TopologyProblem::InputBoundaries::{Outer, Inner}
            <document cache>/coarse.json   the layout WITH its attributes
            <document cache>/field.json    the cross field, if step 3 solved one
    writes  TopologyProblem::QuadMesh

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

# Development bootstrap -- delete once compas_singular is installed into Rhino's
# Python. MUST run before any compas_singular import: Rhino resets sys.path between
# runs but keeps sys.modules, so put the source on the path and drop a stale copy
# (see CMD_start for why every module, framefield included, has to go).
import sys
SINGULAR_SRC = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular\src"
if SINGULAR_SRC not in sys.path:
    sys.path.insert(0, SINGULAR_SRC)
if not getattr(sys, "compas_singular_keep_modules", False):  # set by headless tests
    for _mod in list(sys.modules):
        if _mod == "compas_singular" or _mod.startswith("compas_singular."):
            del sys.modules[_mod]

import rhinoscriptsyntax as rs

from compas_singular.datastructures import coarse_edges_to_curves, snap_corners_to_walls
from compas_singular.framefield.field import CrossField
from compas_singular.rhino.helpers import bake_mesh, clear_layer
from compas_singular.rhino.helpers import read_boundaries, read_boundary_loops, read_polylines
from compas_singular.framefield.quality import mesh_quality

from compas_singular.rhino.project import get_settings, set_settings
from compas_singular.rhino.project import cache_path, read_layout, FIELD_CACHE, DENSE_CACHE
from compas_singular.rhino.project import resolve_relax, resolve_symmetry
# One implementation, shared with CMD_densities.
from compas_singular.rhino.project import resolve_densities
from compas_singular.rhino.project import ROOT

#The walls the layout's boundary edges densify ALONG, as a multiple of the
#background spacing. Finer than the background, because these points are the
#boundary from here on: the final mesh follows the real curve only as closely as
#these follow it.
WALL_SAMPLING_FACTOR = 0.25


def main():
    settings = get_settings()

    # ------------------------------------------------------------------
    # what the previous steps left behind
    # ------------------------------------------------------------------
    # The side-car when it still matches the baked layout, the document
    # otherwise -- see ``read_layout``. Trusting the side-car unconditionally
    # used to mean a stale coarse.json (one CMD_edit_coarse_mesh never
    # finished writing, say) silently densified whatever it had, with the
    # baked ``Mesh`` layer showing the real, edited layout and nobody told.
    coarse, _poles, _source = read_layout()

    # The field is a side-car because it cannot be baked: it is 443 background
    # vertices and a complex number per vertex, and no Rhino object holds that.
    # ``default=None`` because a skeleton-route document legitimately has none.
    field = CrossField.load_from_json(cache_path(FIELD_CACHE, create=False),
                                      default=None)

    # ------------------------------------------------------------------
    # is the field still this document's field?
    # ------------------------------------------------------------------
    # A field is a function of its domain, and densifying with one solved for a
    # DIFFERENT outline fails nowhere: every patch interior still integrates,
    # the mesh still welds, the quality gate still passes. It is simply aligned
    # to a shape that is no longer there. ``mismatch`` is the only thing in the
    # pipeline that can catch it.
    outer, inners, guides, _point_features = read_boundaries(
        spacing=settings["triangulation_spacing"])
    if field is not None:
        # ``resolve_relax`` / ``resolve_symmetry``, not the defaults: they are
        # what step 3 SOLVED with, and ``mismatch`` compares solver settings as
        # well as geometry. Checking against ``relax=False`` would report
        # "solver settings changed" on every guided document and throw away a
        # perfectly good field.
        why = field.mismatch(outer, inners, guides=guides,
                             mode=settings["guide_allignment"],
                             target_length=settings["triangulation_spacing"],
                             relax=resolve_relax(settings, guides),
                             symmetry=resolve_symmetry(settings))
        if why:
            print("field ignored: {} since it was solved. Re-run step 3.".format(why))
            field = None


    options = (
        ("Field", "Off", "On"),
        ("Boundary_curvature", "False", "True"),
        ("Skeleton_curvature", "False", "True")
    )

    option_defaults = [False, True, True]

    options = rs.GetBoolean("Densify the coarse layout.", options, option_defaults)
    if options is None:
        print("Command cancelled. Nothing changed.")
        return
    field_aware, boundary_curvature, skeleton_curvature = options

    settings["field_aware"] = field_aware
    set_settings(settings)

    if field is None and settings["field_aware"]:
        print("no field available -- patch interiors will be plain Coons. This "
              "is expected on the skeleton route.")

    # ------------------------------------------------------------------
    # the shape of every coarse edge
    # ------------------------------------------------------------------
    # A coarse edge is a straight chord and the layout has to keep it that way,
    # so the SHAPE of each edge is handed to ``densification`` separately. Built
    # from the DOCUMENT -- walls plus the branches on ``Skeleton::Polylines`` --
    # so it survives a save, a reopen and the single-precision round trip
    # ``rs.AddMesh`` puts every corner through.
    wall_sampling = settings["triangulation_spacing"] * WALL_SAMPLING_FACTOR
    outer_loop, inner_loops = read_boundary_loops(wall_sampling)
    loops = [outer_loop] + inner_loops

    # Before the mapping, never after: the arcs are anchored on corner positions.
    snapped, worst = snap_corners_to_walls(coarse, loops=loops)
    if snapped:
        print("snapped {} boundary corner(s) onto their wall, worst {:.4f}".format(
            snapped, worst))

    # ``AddLayer`` on an existing layer returns its FULL '::' path without
    # recreating it, and the full path is what ``ObjectsByLayer`` needs -- a
    # nested layer's short name is not resolvable on its own. Same idiom as
    # CMD_edit_coarse_mesh.
    polyline_layer = rs.AddLayer(name="Polylines", parent="Skeleton")
    edges_to_curves, tally = coarse_edges_to_curves(
        coarse, loops=loops, polylines=read_polylines(polyline_layer))

    # ------------------------------------------------------------------
    # densities, then the mesh
    # ------------------------------------------------------------------
    # The densities step 5 saved on the layout. Re-derived only when the layout
    # has none -- read from the baked mesh, or step 5 never ran -- because
    # ``densification`` raises KeyError on the first strip without one.
    resolve_densities(coarse, settings)

    if settings["field_aware"] and field is not None:
        dense = coarse.densification(boundary_curvature=boundary_curvature, skeleton_curvature=skeleton_curvature, field=field)
    else:
        print(coarse.edges_to_curves())
        dense = coarse.quad_mesh(boundary_curvature=boundary_curvature, skeleton_curvature=skeleton_curvature)

    layer = rs.AddLayer("QuadMesh", parent=ROOT)
    clear_layer(layer, clean_sublayers=True)
    bake_mesh(dense, layer)
    dense.save_to_json(cache_path(DENSE_CACHE))

    # ------------------------------------------------------------------
    # what to read in the output
    # ------------------------------------------------------------------
    print("mesh: {} faces from {} patch(es), field_aware={}".format(
        dense.number_of_faces(), coarse.number_of_faces(),
        settings["field_aware"] and field is not None))

    quality = mesh_quality(dense)
    print("quality: min angle {:.1f}, max angle {:.1f}, aspect {:.2f}, "
          "{:.1%} of corners below {:.0f} deg".format(
              quality["min_angle"], quality["max_angle"], quality["aspect_max"],
              quality["share_below"], quality["low_angle"]))

    # ``chord`` is the count that costs area. On a curved domain it should be
    # interior edges only -- a boundary edge that chords means a corner is not
    # on its wall. An interior corner DRAGGED in step 4 also chords here: the
    # warp that would keep its curvature needs the live traced network, which
    # this step no longer builds.
    print("coarse edges: {}".format(dict(tally)))
    if tally.get("wall_missed"):
        print("  {} boundary edge(s) did NOT get their wall arc -- on a hole "
              "this is what makes one circle come out round and the next a "
              "polygon".format(tally["wall_missed"]))

    print("baked to '{}'".format(layer))
    print("note: re-running this step regenerates from the coarse layout and "
          "overwrites any hand edits made in CMD_edit_quad_mesh.")
    print("next: CMD_smoothen / CMD_smoothen_guide to relax the mesh, CMD_dual for "
          "the dual mesh, or CMD_edit_quad_mesh to hand-edit it.")


if __name__ == "__main__":
    main()
