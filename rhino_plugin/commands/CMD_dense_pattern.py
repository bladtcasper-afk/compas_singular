#! python3

# r: compas

#Temporary import of compas_singular development library
from CMD_start import import_compas_singular
import_compas_singular()

import rhinoscriptsyntax as rs

from compas_singular.datastructures import coarse_edges_to_curves, snap_corners_to_walls
from compas_singular.framefield.field import CrossField
from compas_singular.rhino.helpers.helpers import bake_mesh, clear_layer
from compas_singular.rhino.helpers.helpers import read_boundaries, read_boundary_loops, read_polylines
from compas_singular.framefield.quality import mesh_quality
from compas_singular.rhino import mesh_ui

from CMD_start import get_settings, set_settings
from CMD_start import cache_path, read_layout, FIELD_CACHE, COARSE_CACHE
from CMD_start import resolve_relax, resolve_symmetry
from compas_singular.datastructures.mesh_quad_coarse.patterns import PATTERNS


coarse, _poles, _source = read_layout()

scene = mesh_ui.PickableMesh("TopologyProblem::Attributes::TempPattern")

rs.LayerVisible(layer="Skeleton", visible=False)
rs.LayerVisible(layer="QuadMesh", visible=False)

while True:
    scene.draw_faces(coarse, highlight_pattern=True)

    pattern = rs.GetString("Apply a dense mesh patterns to the coarse mesh.", defaultString="Finish", strings=PATTERNS+["Finish"])
    pattern = (pattern or "Finish").lower()

    if pattern=="finish":
        break

    mode = rs.GetString("Pick the faces or assign the pattern to all faces.", defaultString="Pick", strings=["Pick", "Pick_multiple", "All"])
    mode = (mode or "Pick").lower()
    if mode=="all":
        coarse.set_global_face_pattern(pattern)
        scene._clear_faces()
        scene._clear_text()
        continue
    elif mode=="pick_multiple":
        fkeys = scene.pick_faces()
        if fkeys is not None:
            rs.EnableRedraw(False)
            for fkey in fkeys:
                if coarse.get_face_pattern(fkey) != pattern:
                    coarse.set_face_pattern(fkey, pattern)
                    print(f"Successfully added {pattern} pattern to face {fkey}")
                scene.update_face(coarse, fkey)
            rs.EnableRedraw(True)
        else:
            continue
    elif mode=="pick":
        while True:
            fkey = scene.pick_face()
            if fkey is not None:
                if coarse.get_face_pattern(fkey) != pattern:
                    coarse.set_face_pattern(fkey, pattern)
                    print(f"Successfully added {pattern} pattern to face {fkey}")
                scene.update_face(coarse, fkey)
            else:
                break
    scene._clear_faces()
    scene._clear_text()

rs.LayerVisible(layer="Skeleton", visible=True)
rs.LayerVisible(layer="QuadMesh", visible=True)
scene.clear()

coarse.save_to_json(cache_path(COARSE_CACHE))




