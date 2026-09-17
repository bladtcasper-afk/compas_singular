#! python3

# r: compas

#Temporary import of compas_singular development library
from CMD_start import import_compas_singular
import_compas_singular()

import rhinoscriptsyntax as rs
import Rhino

import compas_rhino as cr
from compas_rhino.conversions import curve_to_compas_polyline, curve_to_compas_circle

from math import pi

import compas
from compas.datastructures.mesh.mesh import Mesh
from compas.datastructures.mesh.duality import mesh_dual
from compas.geometry import Point, Polyline, Line, Polygon, is_polygon_in_polygon_xy
from compas.scene import Scene
from compas_singular.datastructures import CoarseQuadMesh, QuadMesh
from compas_singular.algorithms import boundary_triangulation, SkeletonDecomposition
from compas_singular.datastructures import CoarsePseudoQuadMesh
from compas_singular.rhino.helpers.helpers import clear_layer, read_boundaries
from CMD_start import get_settings, set_settings

settings = get_settings()

# NOTE ADD FILTER TO PREVENT DOUBLE BOUNDARY SELECTION
def closed_filter(rhobj, geometry, component_index):
    if isinstance(geometry, Rhino.Geometry.Curve):
        return rs.IsCurveClosed(rhobj)
    elif isinstance(geometry, Rhino.Geometry.Surface):
        return rs.IsSurfaceClosed(rhobj)
    elif isinstance(geometry, Rhino.Geometry.Mesh):
        return rs.IsMeshClosed(rhobj)

    return False

def select_boundary():
    curve = rs.GetObject(message="Pick a closed polyline as boundary", filter=rs.filter.curve, preselect=False, select=True, subobjects=False, custom_filter=closed_filter)
    if curve:
        curve = rs.CopyObject(curve)
    return curve

def select_boundaries():
    curves = rs.GetObjects(message="Pick closed polylines as boundaries", filter=rs.filter.curve, group=True, preselect=False, select=True, minimum_count=0, maximum_count=10, custom_filter=closed_filter)
    if not curves:
        return []
    return rs.CopyObjects(curves) or []

def edit_outer():
    existing = rs.ObjectsByLayer("Outer")
    rs.ObjectLayer(existing, "TrashBin")
    outer_boundary = select_boundary()
    if outer_boundary:
        rs.ObjectLayer(outer_boundary, layer="Outer")

def edit_inner():
    mode = rs.GetString(message="Inner boundary selection mode", defaultString="Add", strings=["Add", "Delete", "New"])
    if mode is None:
        return
    mode = mode.lower()

    if mode == "new":
        new_inner_boundaries()
    elif mode == "delete":
        delete_inner_boundaries()
    elif mode == "add":
        add_inner_boundaries()
    else:
        print("Unrecognised mode '{}', defaulting to Add".format(mode))
        add_inner_boundaries()

def new_inner_boundaries():
    existing = rs.ObjectsByLayer("Inner")
    rs.ObjectLayer(existing, "TrashBin")
    for boundary in select_boundaries():
        rs.ObjectLayer(boundary, layer="Inner")

def add_inner_boundaries():
    for boundary in select_boundaries():
        rs.ObjectLayer(boundary, layer="Inner")

def delete_inner_boundaries():
    def inner_filter(rhobj, geometry, component_index):
        return rs.ObjectLayer(rhobj) == "TopologyProblem::InputBoundaries::Inner"

    to_delete = rs.GetObjects(message="Select inner boundaries to delete", filter=rs.filter.curve, group=True, preselect=False, select=True, minimum_count=0, custom_filter=inner_filter)
    if to_delete:
        rs.ObjectLayer(to_delete, "TrashBin")
    else:
        print("No boundaries selected for deletion.")

def edit_guides():
    mode = rs.GetString(message="Inner boundary selection mode", defaultString="Add", strings=["Add", "Delete", "New"])
    if mode is None:
        return
    mode = mode.lower()

    if mode == "new":
        new_guides()
    elif mode == "delete":
        delete_guides()
    elif mode == "add":
        add_guides()
    else:
        print("Unrecognised mode '{}', defaulting to Add".format(mode))
        add_guides()

def new_guides():
    existing = rs.ObjectsByLayer("Guides")
    rs.ObjectLayer(existing, "TrashBin")

    curves = rs.GetObjects(message="Pick polylines as guides.", filter=rs.filter.curve, group=True, preselect=False, select=True, minimum_count=0)
    curves = rs.CopyObjects(curves) if curves else []
    for guide in curves:
        rs.ObjectLayer(guide, layer="Guides")

def add_guides():
    curves = rs.GetObjects(message="Pick polylines as guides.", filter=rs.filter.curve, group=True, preselect=False, select=True, minimum_count=0)
    curves = rs.CopyObjects(curves) if curves else []
    for guide in curves:
        rs.ObjectLayer(guide, layer="Guides")

def delete_guides():
    def guide_filter(rhobj, geometry, component_index):
        return rs.ObjectLayer(rhobj) == "TopologyProblem::InputBoundaries::Guides"

    to_delete = rs.GetObjects(message="Select guides to delete", filter=rs.filter.curve, group=True, preselect=False, select=True, minimum_count=0, custom_filter=guide_filter)
    if to_delete:
        rs.ObjectLayer(to_delete, "TrashBin")
    else:
        print("No guides selected for deletion.")    

def edit_point_features():
    mode = rs.GetString(message="Point feature selection mode", defaultString="Add", strings=["Add", "Delete", "New"])
    if mode is None:
        return
    mode = mode.lower()

    if mode == "new":
        new_pts()
    elif mode == "delete":
        delete_pts()
    elif mode == "add":
        add_pts()
    else:
        print("Unrecognised mode '{}', defaulting to Add".format(mode))
        add_pts()

def new_pts():
    existing = rs.ObjectsByLayer("PointFeatures")
    rs.ObjectLayer(existing, "TrashBin")

    curves = rs.GetObjects(message="Pick points as point features.", filter=rs.filter.point, group=True, preselect=False, select=True, minimum_count=0)
    curves = rs.CopyObjects(curves) if curves else []
    for guide in curves:
        rs.ObjectLayer(guide, layer="PointFeatures")

def add_pts():
    curves = rs.GetObjects(message="Pick points as point features.", filter=rs.filter.point, group=True, preselect=False, select=True, minimum_count=0)
    curves = rs.CopyObjects(curves) if curves else []
    for guide in curves:
        rs.ObjectLayer(guide, layer="PointFeatures")

def delete_pts():
    def guide_filter(rhobj, geometry, component_index):
        return rs.ObjectLayer(rhobj) == "TopologyProblem::InputBoundaries::PointFeatures"

    to_delete = rs.GetObjects(message="Select point features to delete", filter=rs.filter.point, group=True, preselect=False, select=True, minimum_count=0, custom_filter=guide_filter)
    if to_delete:
        rs.ObjectLayer(to_delete, "TrashBin")
    else:
        print("No point features selected for deletion.")    



#Boundary Selection
rs.AddLayer(name="Outer", parent="InputBoundaries", color=(255, 0, 0))
rs.AddLayer(name="Inner", parent="InputBoundaries", color=(0, 255, 0))
rs.AddLayer(name="Guides", parent="InputBoundaries", color=(255, 127, 0))
rs.AddLayer(name="PointFeatures", parent="InputBoundaries", color=(0, 255, 0))
rs.AddLayer(name="TrashBin", parent="TopologyProblem", color=(90, 90, 90), visible=False)

while True:
    section = rs.GetString(message="Edit boundary", defaultString="Continue",
                           strings=["Outer", "Inner", "Background_Triangulation",
                                    "Point_Features", "Guides", "Guide_Allignment",
                                    "Field_Solver", "Clear", "Continue"])
    section = (section or "Continue").lower()
    # NOTE the output always contains NO capitals!
    if section == "Outer".lower():
        edit_outer()
    elif section == "Inner".lower():
        edit_inner()
    elif section == "guides":
        edit_guides()
    elif section=="Guide_Allignment".lower():
        answer = rs.GetString(message="Elements run ALONG or ACROSS the guides?",
                              defaultString=settings["guide_allignment"],
                              strings=["tangent", "perpendicular"])
        if answer:
            settings["guide_allignment"] = answer.lower()
            set_settings(settings)
    elif section=="Field_Solver".lower():
        current = settings.get("relax", "auto")
        default = "Auto" if str(current).lower() == "auto" else ("On" if current else "Off")
        answer = rs.GetString(message="Relax the field solve? Auto = on when guides exist",
                              defaultString=default, strings=["Auto", "On", "Off"])
        if answer:
            answer = answer.lower()
            settings["relax"] = "auto" if answer == "auto" else (answer == "on")
            set_settings(settings)
    elif section=="Background_Triangulation".lower():
        value = rs.GetReal("Background spacing (NOT the quad size)",
                           settings["triangulation_spacing"], 1e-3)
        if value:
            settings["triangulation_spacing"] = value
            set_settings(settings)
    elif section=="Point_Features".lower():
        edit_point_features()
    elif section=="Clear".lower():
        clear_layer("InputBoundaries", True)
    elif section == "Continue".lower():
        break
    else:
        print("Unrecognised option '{}'".format(section))

outer, inner, guides, point_features = read_boundaries(spacing=settings["triangulation_spacing"])

print("Boundary selection completed.")
print("1 outer boundary curve")
print(f"{len(inner)} inner boundary curves")
print(f"{len(guides)} guide curves")
print(f"{len(point_features)} point features")
print("next: CMD_coarse_mesh to generate a coarse layout, or CMD_read_coarse_mesh to read a hand-drawn one.")

"""
def densify(curve, resolution: int) -> list[Point]:
    #Important to discretise the sides and never the entire polyline to make sure that hard corners are kept as a point
    if isinstance(curve, list) and isinstance(curve[0], Point):
        polyline = Polyline(curve)
        line_pts = curve
    elif isinstance(curve, Polyline):
        polyline = curve
        line_pts = curve.points
    else:
        raise TypeError(f"Wrong input type, expected Polyline or list[Point], got {type(curve)}")

    lines = []
    for indx in range(len(line_pts)-1):
        lines.append(Line(line_pts[indx], line_pts[indx+1]))
    
    densified_pts = []
    for line in lines:
        #As divide_by_count is not implemented yet, this is a small work around using Polyline. Alternative to doing a range division
        # densified_pts.extend(pts = line.divide_by_count(resolution))
        polyline: Polyline = line.to_polyline(n=resolution)
        densified_pts.extend(polyline.points)
    return densified_pts

#note how for resolution one (keep same number of points) the trimesh is not complete. Higher resolution is required for decomposition in step2
outer = densify(outer,20)

curves = list(inner)
inner = []
for curve in curves:
    inner.append(densify(curve,10))


scene = Scene()
scene.clear()
for curve in inner:
    for pt in curve:
        scene.add(pt)
for pt in outer:
    scene.add(pt)

#-----------------------------------------
#Start with boundary triangulation
#triangulation.py
def to_coordinates(pts: list[Point]):
    return [[pt.x, pt.y, pt.z] for pt in pts]

point_features = []
polyline_features = [[2.5,3.5,0]]

#Note that the outer boundary should be singular, but the inenr boundary can be a list.
trimesh = boundary_triangulation(outer_boundary=to_coordinates(outer), inner_boundaries=[inner], point_features=point_features)


#-----------------------------------------
#Skeleton decomposition: turns trimesh into a coarse quad-patch layout -> CoarsePseudoQuadMesh
#Uses same feature values as in trimesh
#decomposition.py

decomposition = SkeletonDecomposition.from_mesh(trimesh) #plotting this returns basic trimesh, but it contains the info and methods

coarse_mesh: CoarsePseudoQuadMesh = decomposition.decomposition_mesh(point_features) #decomposed CoarsePseudoQuadMesh  derived from SkeletonDecomposition
#Each path is a Coons-patch placeholder. Coons patch helps divide complicated area geometries

#-----------------------------------------
#density -> CoarseQuadMesh

#this methods collects the strips and returns skey and vertix indices. Afterwards skey can be collected using .strips
coarse_mesh.collect_strips()



#using skey the density of each strip can also be set individually
# coarse_mesh.set_strips_density_target(t=5)
# coarse_mesh.set_strips_density_func(func=, func_args=)
coarse_mesh.set_strips_density_target(t=1)

coarse_mesh.densification() #the density is stored in the object, but not yet visible.
dense_mesh: QuadMesh = coarse_mesh.get_quad_mesh() #Move from pseudo quad to quad


#Viewer setup


scene.add(dense_mesh)
scene.add(coarse_mesh)
scene.draw()"""
"""
step1 = viewer.scene.add_group(name="step1: boundary")
for pt in outer:
    step1.add(pt)
step1.add(outer_curve)
step1.add(inner_curve)
step1.add(trimesh, name="trimesh")

step2 = viewer.scene.add_group(name="step2: trimesh")
step2.add(decomposition)
step2.add(coarse_mesh)

step3 = viewer.scene.add_group(name="step3: densification")
step3.add(dense_mesh)

#-----------------------------------------
#strip visualisation
step4 = viewer.scene.add_group(name="step4: strips")

#.strips() gives skeys -> .strip_faces(skey) gives the faces in a strip using fkeys -> face_vertices(fkey) gives the vertices of a face
dense_mesh.collect_strips()
for skey in coarse_mesh.strips():
    strip = viewer.scene.add_group(name="strip", parent=step4)
    fkeys = coarse_mesh.strip_faces(skey)
    print(fkeys)
    for fkey in fkeys:
        # fkey -> ordered vertex keys -> vertex coordinates -> Polygon
        vkeys = coarse_mesh.face_vertices(fkey)
        face_pts = [Point(*coarse_mesh.vertex_coordinates(vkey)) for vkey in vkeys]
        face_polygon = Polygon(face_pts)
        strip.add(face_polygon)


viewer.show()"""
