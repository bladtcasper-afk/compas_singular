#! python3

# r: compas

#Temporary import of compas_singular development library
from CMD_start import import_compas_singular
import_compas_singular()


import rhinoscriptsyntax as rs
import scriptcontext as sc
import Rhino

import compas_rhino as cr
from compas_rhino.conversions import polyline_to_rhino, point_to_rhino

import compas
from compas.geometry import Point, Line, Polyline, Polygon, is_polygon_in_polygon_xy
from compas.scene import Scene
from compas_singular.algorithms import boundary_triangulation, SkeletonDecomposition
from compas_singular.datastructures import CoarsePseudoQuadMesh
from compas_singular.rhino.helpers.helpers import clean_layer, curve_to_polyline, read_boundary_loops, bake_edge_curves, bake_polylines, read_boundaries
from compas_singular.rhino.coarse_curves import coarse_edges_to_curves, snap_corners_to_walls
from CMD_start import get_settings
from CMD_start import get_decomposition

settings = get_settings()

#How finely a CURVED input is sampled. A polyline input is still taken at its
#own vertices, so nothing changes for the polyline boundaries this workflow has
#always used -- but an arc or a NURBS curve is now usable at all, where
#``curve_to_compas_polyline`` raised ConversionError and the command died here.
CURVE_SAMPLING = settings["triangulation_spacing"]
#The walls the layout's boundary edges densify ALONG. Finer than the background,
#because these points are the boundary from here on: the final mesh follows the
#real curve only as closely as these follow it.
WALL_SAMPLING = settings["triangulation_spacing"] * 0.25

#Finding an checking of boundaries
outer, inners, guides, point_features = read_boundaries(spacing=settings["triangulation_spacing"])

#Coarse mesh generation

#New approach: frame-field decomposition -> all-quad mesh
#See examples/New approach/README.md -- "the two lines that matter":
#  d = FieldDecomposition.from_boundary(outline, guides=cables)
#  mesh = d.quad_mesh(target_length=0.5)

def coarse_from_field():
    QUAD_TARGET_LENGTH = 1

    # Solved once per document and reused across steps 3, 4 and 6 -- and
    # across Rhino restarts. ``CMD_start.get_decomposition`` owns the settings,
    # the resolvers and the staleness rules; see its docstring.
    decomposition = get_decomposition()
    print(decomposition.field.report())
    # Print the group: a smaller group than the drawing suggests means the
    # layout is ALLOWED to be less symmetric, and nothing downstream complains.
    group = decomposition.symmetry
    if group is None or group.trivial:
        print("symmetry: none detected -- the layout is under no obligation to have any")
    else:
        print("symmetry: order {} about ({:.3f}, {:.3f}) -- {}".format(
            len(group), group.centre[0], group.centre[1], ", ".join(group.names())))

    skeleton = decomposition.decomposition_polylines()
    coarse_mesh = decomposition.decomposition_mesh()

    if True:
        pass

    #dense_mesh = decomposition.quad_mesh(target_length=QUAD_TARGET_LENGTH)
    #scene=Scene()
    #scene.add(dense_mesh)
    #scene.draw()

    return coarse_mesh, skeleton

#From Skeleton
def coarse_from_skeleton():
    BACKGROUND_TARGET_LENGTH = settings["triangulation_spacing"]

    def densify(curve, target_length):
        line_pts = curve.points if isinstance(curve, Polyline) else list(curve)
        densified_pts = []
        for a, b in zip(line_pts, line_pts[1:]):
            line = Line(a, b)
            resolution = max(1, int(round(line.length / target_length)))
            # drop the closing point -- the next segment starts there
            densified_pts.extend(line.to_polyline(n=resolution).points[:-1])
        return densified_pts

    dense_outer = densify(outer,BACKGROUND_TARGET_LENGTH)

    curves = inners.copy()
    dense_inners = []
    for curve in curves:
        dense_inners.append(densify(curve,BACKGROUND_TARGET_LENGTH))

    def to_coordinates(pts: list[Point]):
        return [[pt.x, pt.y, pt.z] for pt in pts]

    trimesh = boundary_triangulation(outer_boundary=to_coordinates(dense_outer), inner_boundaries=[to_coordinates(dense_inner) for dense_inner in dense_inners], point_features=point_features)
    decomposition = SkeletonDecomposition.from_mesh(trimesh) #plotting this returns basic trimesh, but it contains the info and methods
    coarse_mesh: CoarsePseudoQuadMesh = decomposition.decomposition_mesh(point_features)
    skeleton = decomposition.decomposition_polylines()
    
    return coarse_mesh, skeleton

coarse_mesh, skeleton = None, None
while True:
    mode = rs.GetString(message="Choose coarse mesh generation", defaultString="Skeleton", strings=["FrameField", "Skeleton", "Exit"])
    mode = (mode or "Exit").lower()

    if mode == "Skeleton".lower():
        coarse_mesh, skeleton = coarse_from_skeleton()
        break
    elif mode== "FrameField".lower():
        coarse_mesh, skeleton = coarse_from_field()
        break
    elif mode == "Exit".lower():
        break
    else:
        print(f"Unrecognised mode: {mode}")

    
#Bake results
if coarse_mesh and skeleton:
    clean_layer("Skeleton", clean_sublayers = True) #Clean Layer before adding objects

    #BEFORE the bake, or the baked mesh and the edge curves disagree about where
    #the corners are. A boundary corner of the layout IS a point of the domain
    #boundary, but nothing upstream puts it there -- the background
    #triangulation places it, and on a curved wall that means on a CHORD of the
    #wall. An arc is pinned to its two corners, so those corners would be the
    #only points of the dense boundary not on the wall. Measured on a plate with
    #three r=0.7 holes near the wall: worst corner 0.106 off its own hole, 15.1%
    #of the radius. Snapping costs nothing -- same patches, same face count,
    #quality unchanged to a few tenths of a degree.
    outer_loop, inner_loops = read_boundary_loops(WALL_SAMPLING)
    snapped, worst = snap_corners_to_walls(
        coarse_mesh, loops=[outer_loop] + inner_loops)
    if snapped:
        print("snapped {} boundary corner(s) onto their wall, worst {:.4f}".format(
            snapped, worst))

    rs.AddLayer(name="Skeleton", parent="TopologyProblem")

    layer = rs.AddLayer(name="Poles", parent="Skeleton")
    vkeys = coarse_mesh.poles()
    pts = [coarse_mesh.vertex_coordinates(v) for v in vkeys]
    points = [point_to_rhino(pt) for pt in pts]
    guids = rs.AddPoints(points)
    for guid in guids:
        rs.ObjectLayer(guid, layer=layer)

    #``rs.AddPolyline`` RAISES -- 'Unable to add polyline to document' -- on a
    #rib Rhino will not take, and it judges that at the DOCUMENT tolerance, not
    #at 1e-9: a separatrix that doubles back on itself by half a tolerance is
    #geometrically fine and still unbakeable. That used to kill the command
    #here, which cost the MESH as well, because the mesh is baked below this.
    #``bake_polylines`` cleans each rib and guards the add.
    layer = rs.AddLayer(name="Polylines", parent="Skeleton")
    _guids, skipped = bake_polylines(skeleton, layer)
    print("separatrices: {} baked".format(len(_guids)))
    if skipped:
        print("  {} rib(s) skipped -- Rhino refused them (coincident points at "
              "the document tolerance, or fewer than two)".format(skipped))

    layer = rs.AddLayer(name="Mesh", parent="Skeleton")
    # face_vertices gives vertex KEYS, the vertex list is positional, and the
    # two only agree while the keys are 0..n-1. Repair and weld leave gaps, and
    # the baked faces would then reference the wrong corners -- silently.
    index = {v: i for i, v in enumerate(coarse_mesh.vertices())}
    vertices = [coarse_mesh.vertex_attributes(v, "xyz") for v in coarse_mesh.vertices()]
    faces = [[index[v] for v in coarse_mesh.face_vertices(f)] for f in coarse_mesh.faces()]
    guid = rs.AddMesh(vertices, faces)
    rs.ObjectLayer(guid, layer=layer)

    #The layout WITH its curvature. A Rhino mesh has straight edges, so on a
    #curved domain the mesh above is drawn cutting the corner off its own
    #boundary -- which reads as the workflow having lost the curve when only
    #the display has. These polylines are what CMD_quad_mesh densifies along,
    #rebuilt there from the same document data rather than read back from here.
    curves, tally = coarse_edges_to_curves(
        coarse_mesh, loops=[outer_loop] + inner_loops, polylines=skeleton)
    # AddLayer returns the FULL '::' path, which is what every rs call below
    # needs -- a nested layer's short name is not resolvable on its own.
    edge_layer = rs.AddLayer(name="EdgeCurves", parent="Skeleton", color=(0, 120, 200))
    _guids, skipped = bake_edge_curves(curves.values(), edge_layer)
    print("coarse edges: {}".format(tally))
    if skipped:
        print("  {} edge curve(s) Rhino refused -- those edges densify as "
              "chords".format(skipped))
    if tally["chord"]:
        print("  {} edge(s) densify as a straight chord -- interior edges with "
              "no traced branch, which is expected, or a boundary corner that "
              "is not on a wall, which is not".format(tally["chord"]))
    if tally.get("wall_missed"):
        print("  {} boundary edge(s) did NOT get their wall arc -- on a hole "
              "this is what makes one circle come out round and the next a "
              "polygon".format(tally["wall_missed"]))

