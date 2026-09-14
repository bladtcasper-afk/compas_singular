#! python3

# r: compas

"""**Step 5 -- how many elements across each strip.**

    reads   TopologyProblem::Skeleton::Mesh
    writes  <document cache>/coarse.json   the layout WITH strips_density
    scratch TopologyProblem::Attributes::Densities   pickable strips + labels

``pick`` a strip and give it a number, or set a global ``target`` -- a length
or a density -- that sets every strip at once. A strip is a band of quads
running across the layout, and it is drawn as one: a ribbon shaded light to
dark blue from the smallest density in the layout to the largest, with its
density on it. A quad mesh cannot subdivide one patch without subdividing
everything in line with it, which is why the band and not the edge is the unit.

The edge-picking version of this command is kept in
``backup_2026-09-14_density_edge_pick/`` (with the ``mesh_ui`` it ran against).

**The densities are stored by strip key, on the layout, and the layout is saved
to the side-car** -- like the face patterns ``CMD_dense_pattern`` sets. Strip
keys are NOT stable across a bake (``collect_strips`` numbers strips in edge
order, and a layout read back from Rhino comes back in whatever order the
geometry did), so they are only trusted together with the strips they were
saved with: ``read_layout`` uses the side-car only while it still describes the
baked layout, and the strips are never re-collected on it. A layout that has
lost its densities -- read from the document, or rebuilt in step 3 or 4 -- is
re-set from the global target by ``resolve_densities``; per-strip picks are gone
at that point, because the strips they belonged to may be too.
"""

#Temporary import of compas_singular development library
from CMD_start import import_compas_singular
import_compas_singular()

import rhinoscriptsyntax as rs

from compas_singular.rhino import mesh_ui
from compas_singular.rhino.helpers.helpers import clear_layer
# This step reads its layout from the document (read_coarse) and touches no
# boundary geometry, so the wall helpers that used to be imported here --
# read_boundary_loops, coarse_edges_to_curves, snap_corners_to_walls,
# read_boundaries, curve_to_polyline, bake_edge_curves -- were never called.
from CMD_start import set_settings, get_settings, read_layout
from CMD_start import cache_path, COARSE_CACHE
# One implementation, shared with CMD_quad_mesh.
from CMD_start import resolve_densities

DENSITY_LAYER = "TopologyProblem::Attributes::Densities"

def draw(scene, coarse):
    """Every strip as a pickable ribbon, shaded by density, with its number on it.

    Drawn the way ``CMD_edit_coarse_mesh`` draws strips for removing and
    dividing -- ``mesh_ui.draw_strips`` -- but with ``collect=False``: what is
    picked has to come back as the key ``strips_density`` is stored under, and a
    re-collection may number the same bands differently.

    The label is the density the strip will actually be meshed at. The shade is
    relative: lightest blue for the smallest density in this layout, darkest for
    the largest.
    """
    strip_densities = coarse.get_strip_densities()
    scene.draw_strips(coarse, collect=False,
                      strip_colors=mesh_ui.density_colors(strip_densities),
                      labels={skey: d for skey, d in strip_densities.items()})


def main():
    settings = get_settings()
    # The side-car when it still matches the baked layout, the document
    # otherwise. Either way the strips below are this layout's own.
    coarse, poles, _source = read_layout()

    # The densities saved on the layout by a previous run, when it has a full
    # set; every strip set by the last global target otherwise.
    resolve_densities(coarse, settings)

    print("layout: {} patch(es), {} strip(s); target {} ({})".format(
        coarse.number_of_faces(), len(list(coarse.strips())),
        settings[target_setting(settings)], settings["density_mode"]))

    scene = mesh_ui.PickableMesh(DENSITY_LAYER)

    # Hidden while editing, so the ribbons are not drawn over the baked layout
    # and the old quad mesh. ``LayerVisible`` raises on a missing layer, and
    # QuadMesh does not exist until step 6 has run once.
    if rs.IsLayer("TopologyProblem::Skeleton::Mesh"):
        rs.LayerVisible("TopologyProblem::Skeleton::Mesh", False)
    if rs.IsLayer("TopologyProblem::QuadMesh"):
        rs.LayerVisible("TopologyProblem::QuadMesh", False)

    draw(scene, coarse)

    try:
        while True:
            option = (rs.GetString(message="Densities", defaultString="Finish",
                                   strings=["Pick", "Target_length", "Target_density", "Clear", "Finish"])
                      or "finish").lower()
            if option == "pick":
                # Strip after strip until Esc, redrawn after each so the shade
                # and number of the one just set are visible before the next.
                while True:
                    skey = scene.pick_strip("Pick a strip to set its density (Esc to finish)")
                    if skey is None:
                        break
                    current = coarse.get_strip_density(skey)
                    value = rs.GetInteger("Elements across strip {}".format(skey),
                                          current, 1)
                    if value is None:
                        continue
                    coarse.set_strip_density(skey, int(value))
                    draw(scene, coarse)

            elif option == "target_length":
                value = rs.GetReal("Target quad edge length",
                                   settings["target_length"], 1e-3)
                if value:
                    settings["target_length"] = value
                    settings["density_mode"] = "length"
                    set_settings(settings)
                    # Every strip, picked ones included: the target replaces
                    # all densities, and the labels show the new numbers.
                    coarse.set_strips_density_target(value)
                    draw(scene, coarse)
            
            elif option == "target_density":
                value = rs.GetInteger("Elements across every strip",
                                      settings["target_density"], 1)
                if value:
                    settings["target_density"] = value
                    settings["density_mode"] = "density"
                    set_settings(settings)
                    coarse.set_strips_density(int(value))
                    draw(scene, coarse)

            elif option == "clear":
                # Emptied first, so resolve_densities re-sets every strip from
                # the current target instead of keeping the picked numbers.
                coarse.attributes['strips_density'] = {}
                resolve_densities(coarse, settings, verbose=False)
                draw(scene, coarse)
                print("picks cleared -- every strip follows the target.")

            else:
                break
    finally:
        scene.clear()
        clear_layer(DENSITY_LAYER)
        # Back on however the command ends, Esc and errors included.
        if rs.IsLayer("TopologyProblem::Skeleton::Mesh"):
            rs.LayerVisible("TopologyProblem::Skeleton::Mesh", True)
        if rs.IsLayer("TopologyProblem::QuadMesh"):
            rs.LayerVisible("TopologyProblem::QuadMesh", True)

    # Onto the side-car, like the patterns in CMD_dense_pattern: the densities
    # are an attribute of the layout, and a bake cannot carry attributes. Not in
    # the ``finally`` -- a command that raised must not overwrite a good cache.
    coarse.save_to_json(cache_path(COARSE_CACHE))
    print("densities: saved on the layout ({} strip(s))".format(
        len(coarse.get_strip_densities())))


def target_setting(settings):
    """The settings key of the global rule currently in force."""
    return "target_density" if settings["density_mode"] == "density" else "target_length"



if __name__ == "__main__":
    main()
