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

#Temporary import of compas_singular development library
from CMD_start import import_compas_singular
import_compas_singular()

import rhinoscriptsyntax as rs

from compas_singular.datastructures import coarse_edges_to_curves, snap_corners_to_walls
from compas_singular.framefield.field import CrossField
from compas_singular.rhino.helpers.helpers import bake_mesh, clear_layer
from compas_singular.rhino.helpers.helpers import read_boundaries, read_boundary_loops, read_polylines
from compas_singular.framefield.quality import mesh_quality

from CMD_start import get_settings, set_settings
from CMD_start import cache_path, read_layout, FIELD_CACHE
from CMD_start import resolve_relax, resolve_symmetry
# One implementation, shared with CMD_densities -- see density_common.
from density_common import apply_densities

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

    answer = (rs.GetString(message="Integrate patch interiors from the field?",
                           defaultString="Yes" if settings["field_aware"] else "No",
                           strings=["Yes", "No"]) or "yes").lower()
    settings["field_aware"] = answer == "yes"
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
    applied, lost = apply_densities(coarse, settings["target_length"])
    if lost:
        print("  re-pick those strips -- they are meshing at the target.")

    if settings["field_aware"] and field is not None:
        dense = coarse.densification(edges_to_curves=edges_to_curves, field=field)
    else:
        dense = coarse.densification(edges_to_curves=edges_to_curves)

    layer = rs.AddLayer("QuadMesh", parent="TopologyProblem")
    clear_layer(layer, clean_sublayers=True)
    bake_mesh(dense, layer)

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


if __name__ == "__main__":
    main()
