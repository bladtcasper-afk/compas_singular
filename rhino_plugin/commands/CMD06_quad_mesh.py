#! python3

# r: compas
# r: pydantic

"""**Densify the coarse layout into the quad mesh.**

    reads   the session's layout WITH its attributes, and its field if step 3 solved one
            TopologyProblem::InputBoundaries::{Outer, Inner}
    writes  the session's dense mesh, drawn on TopologyProblem::QuadMesh::Dense

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
``coarse_edges_to_curves`` rebuilds the mapping: the walls come from the input
curves on ``InputBoundaries``, the interior branches from the layout's own
``shape_polylines``. Nothing is keyed by vertex index, so it survives an edit in
CMD_edit_coarse_mesh.

Read the ``coarse edges`` tally it prints. ``chord`` is the count that costs
area, and on a curved domain it should be interior edges only.
"""


import rhinoscriptsyntax as rs

from compas_singular.datastructures import coarse_edges_to_curves, snap_corners_to_walls
from compas_singular.rhino.helpers import clear_layer
from compas_singular.rhino.helpers import read_boundaries, read_boundary_loops
from compas_singular.framefield.quality import mesh_quality

from compas_singular.rhino.project import get_settings, set_settings
from compas_singular.rhino.project import read_layout
from compas_singular.rhino.project import resolve_relax, resolve_field_symmetry, resolve_spacing
from compas_singular.rhino.session import RhinoSession
# One implementation, shared with CMD_densities.
from compas_singular.rhino.project import resolve_densities
from compas_singular.rhino.project import ROOT, layer_path

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
    # A COPY of the session's layout: snapping and densities below change it,
    # and it goes back into the session only when this step records.
    coarse = read_layout()

    # The field lives in the session because it cannot be baked: it is 443
    # background vertices and a complex number per vertex, and no Rhino object
    # holds that. ``None`` on the skeleton route, legitimately.
    field = RhinoSession.current().field

    # ------------------------------------------------------------------
    # is the field still this document's field?
    # ------------------------------------------------------------------
    # A field is a function of its domain, and densifying with one solved for a
    # DIFFERENT outline fails nowhere: every patch interior still integrates,
    # the mesh still welds, the quality gate still passes. It is simply aligned
    # to a shape that is no longer there. ``mismatch`` is the only thing in the
    # pipeline that can catch it.
    spacing = resolve_spacing(settings)
    outer, inners, guides, _point_features = read_boundaries(spacing=spacing)
    if field is not None:
        # ``resolve_relax`` / ``resolve_field_symmetry``, not the defaults: they are
        # what step 3 SOLVED with, and ``mismatch`` compares solver settings as
        # well as geometry. Checking against ``relax=False`` would report
        # "solver settings changed" on every guided document and throw away a
        # perfectly good field.
        why = field.mismatch(outer, inners, guides=guides,
                             mode=settings["guide_alignment"],
                             target_length=settings["triangulation_spacing"],
                             relax=resolve_relax(settings, guides),
                             symmetry=resolve_field_symmetry(settings))
        if why:
            print("field ignored: {} since it was solved. Re-run step 3.".format(why))
            field = None


    options = (
        ("Field", "Off", "On"),
        ("Boundary_curvature", "False", "True"),
        ("Skeleton_curvature", "False", "True")
    )

    option_defaults = [False, True, True]

    if coarse.attributes['decomposition_type'] == 'field':
        option_defaults[0] = True

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
    # from the walls in the DOCUMENT plus the layout's own shape polylines.
    wall_sampling = spacing * WALL_SAMPLING_FACTOR
    outer_loop, inner_loops = read_boundary_loops(wall_sampling)
    loops = [outer_loop] + inner_loops

    # Before the mapping, never after: the arcs are anchored on corner positions.
    snapped, worst = snap_corners_to_walls(coarse, loops=loops)
    if snapped:
        print("snapped {} boundary corner(s) onto their wall, worst {:.4f}".format(
            snapped, worst))

    edges_to_curves, tally = coarse_edges_to_curves(
        coarse, loops=loops, polylines=coarse.shape_polylines())

    # ------------------------------------------------------------------
    # densities, then the mesh
    # ------------------------------------------------------------------
    # The densities step 5 saved on the layout. Re-derived only when the layout
    # has none -- read from the baked mesh, or step 5 never ran -- because
    # ``densification`` raises KeyError on the first strip without one.
    resolve_densities(coarse, settings)

    if settings["field_aware"] and field is not None:
        dense = coarse.quad_mesh(boundary_curvature=boundary_curvature, skeleton_curvature=skeleton_curvature, field=field)
    else:
        print(coarse.edges_to_curves())
        dense = coarse.quad_mesh(boundary_curvature=boundary_curvature, skeleton_curvature=skeleton_curvature)

    # What was under QuadMesh belonged to the previous mesh: an edited copy, a
    # smoothed one, a dual. Recording draws the new one on QuadMesh::Dense.
    clear_layer(rs.AddLayer("QuadMesh", parent=ROOT), clean_sublayers=True)
    layer = layer_path("Dense")
    session = RhinoSession.current()
    session.coarse = coarse          # with the densities it was meshed at
    session.dense = dense

    session.record("Quad mesh")

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

    print("drawn on '{}'".format(layer))
    print("note: re-running this step regenerates from the coarse layout and "
          "overwrites any hand edits made in CMD_edit_quad_mesh.")
    print("next: CMD_smoothen / CMD_smoothen_guide to relax the mesh, CMD_dual for "
          "the dual mesh, or CMD_edit_quad_mesh to hand-edit it.")


if __name__ == "__main__":
    main()
