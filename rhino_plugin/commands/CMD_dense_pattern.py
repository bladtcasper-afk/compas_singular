#! python3

# r: compas

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

from compas_singular.datastructures import coarse_edges_to_curves, snap_corners_to_walls, CoarsePseudoQuadMesh
from compas_singular.framefield.field import CrossField
from compas_singular.rhino.helpers import bake_mesh, clear_layer
from compas_singular.rhino.helpers import read_boundaries, read_boundary_loops, read_polylines
from compas_singular.framefield.quality import mesh_quality
from compas_singular.rhino import mesh_ui

from compas_singular.rhino.project import get_settings, set_settings
from compas_singular.rhino.project import cache_path, read_layout, FIELD_CACHE, COARSE_CACHE
from compas_singular.rhino.project import resolve_relax, resolve_symmetry
from compas_singular.rhino.project import ROOT
from compas_singular.datastructures.mesh_quad_coarse.patterns import PATTERNS


coarse, _poles, _source = read_layout()
guids = rs.ObjectsByLayer("Mesh")
if not guids:
    raise RuntimeError(
        "No coarse layout on '{}' ".format(
            "Mesh"))
coarse = CoarsePseudoQuadMesh.load_from_json(
        cache_path(COARSE_CACHE, create=False), default=None)


scene = mesh_ui.PickableMesh(ROOT + "::Attributes::TempPattern")

rs.LayerVisible(layer="Skeleton", visible=False)
rs.LayerVisible(layer="QuadMesh", visible=False)

while True:
    scene._clear_faces()
    scene._clear_text()
    scene.draw_faces(coarse, highlight_pattern=True)

    pattern = rs.GetString("Apply a dense mesh patterns to the coarse mesh.", defaultString="Finish", strings=PATTERNS+["Finish"])
    pattern = (pattern or "Finish").lower()

    if pattern=="finish":
        break

    mode = rs.GetString("Pick the faces or assign the pattern to all faces.", defaultString="Pick", strings=["Pick", "Pick_multiple", "All"])
    mode = (mode or "Exit").lower()
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
    else:
        continue
    scene._clear_faces()
    scene._clear_text()

rs.LayerVisible(layer="Skeleton", visible=True)
rs.LayerVisible(layer="QuadMesh", visible=True)
scene.clear()

coarse.save_to_json(cache_path(COARSE_CACHE))
custom = sum(1 for p in coarse.dense_patterns().values() if p != "ortho")
print("patterns: {} patch(es) set away from the default 'ortho'".format(custom))
print("note: re-running CMD_coarse_mesh, or a topology-changing edit in "
      "CMD_edit_coarse_mesh, replaces the layout and drops these patterns.")
print("next: CMD_quad_mesh to densify.")




