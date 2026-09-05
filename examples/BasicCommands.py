
import compas
from compas_singular.datastructures import CoarseQuadMesh, QuadMesh
from compas_viewer import Viewer
from compas.datastructures.mesh.mesh import Mesh
from compas.datastructures.mesh.duality import mesh_dual
from compas.geometry import Point, Polyline, Line, Polygon

from compas_singular.algorithms import boundary_triangulation, SkeletonDecomposition
from compas_singular.datastructures import CoarsePseudoQuadMesh

from math import pi

#-----------------------------------------
#Define boundary curves
outer = [Point(5,5,0), Point(5,-5,0), Point(-5,-5,0), Point(-5,5,0), Point(0,7,0)]
inner = [Point(1,1,0), Point(1,-1,0), Point(-1,-1,0), Point(-1,1,0)]

outer += [outer[0]]
inner += [inner[0]]

outer_curve = Polyline(outer + [outer[0]])
inner_curve = Polyline(inner + [inner[0]])

def densify(curve: list[Point] | Polyline, resolution: int) -> list[Point]:
    #Important to discretise the sides and never the entire polyline to make sure that hard corners are kept as a point
    if isinstance(curve, list) and isinstance(curve[0], Point):
        polyline = Polyline(curve)
        line_pts = curve
    elif isinstance(curve, Polyline):
        polyline = curve
        line_pts = curve.split_at_corners(pi/10)
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
outer = densify(outer,4)
inner = densify(inner,2)

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
coarse_mesh.set_strips_density_target(t=0.5)

coarse_mesh.densification() #the density is stored in the object, but not yet visible.
dense_mesh: QuadMesh = coarse_mesh.get_quad_mesh() #Move from pseudo quad to quad


#Viewer setup
viewer = Viewer()


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


viewer.show()

