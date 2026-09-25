
import compas
from compas_singular.datastructures import CoarseQuadMesh, QuadMesh
from compas_viewer import Viewer
from compas_viewer.scene import Tag
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
outer = densify(outer, 20)
inner = densify(inner, 20)

#-----------------------------------------
#Start with boundary triangulation
#triangulation.py
def to_coordinates(pts: list[Point]):
    return [[pt.x, pt.y, pt.z] for pt in pts]

point_features = []
polyline_features = [[2.5,3.5,0]]

#Note that the outer boundary should be singular, but the inenr boundary can be a list.
trimesh = boundary_triangulation(outer_boundary=to_coordinates(outer), inner_boundaries=[to_coordinates(inner)], point_features=point_features)


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
coarse_mesh.set_strips_density(d=5)

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

#-----------------------------------------
#editing the quad mesh dense_mesh: QuadMesh


#find the appropriate polyedge

from compas_singular.datastructures.mesh_quad.grammar.add_strip import (
    add_strip,
    is_polyedge_valid_for_strip_addition,
)
#compas 2.x dropped mesh_smooth_centroid from compas.datastructures; compas_singular ships a shim
from compas.datastructures.mesh.smoothing import mesh_smooth_centroid

#NOTE both .edges() and collect_polyedges() can be used to find the mesh edges, but polyedges are not compatible with edge_line
viewer_polyedges = viewer.scene.add_group(name="polyedges")
viewer_test = viewer.scene.add_group(name="test edges")

"""
for polyedge in coarse_mesh.edges():
    print(polyedge)
    viewer_polyedges.add(coarse_mesh.edge_line(polyedge), name=str(polyedge))
    polyedge_line = Polyline([coarse_mesh.vertex_point(v) for v in polyedge])
    viewer_test.add(Polyline([coarse_mesh.vertex_point(v) for v in polyedge]))
print(len(list(coarse_mesh.edges())))
"""

#collect_polyedges() yields (pkey, polyedge) pairs, where polyedge is a flat list of vertex keys
#of arbitrary length -- NOT an edge. So edge_line() would silently keep only the first two.
for pkey, polyedge in coarse_mesh.collect_polyedges():
    viewer_polyedges.add(Polyline([coarse_mesh.vertex_point(v) for v in polyedge]), name=str(pkey))
print(len(dict(coarse_mesh.collect_polyedges())))


#-----------------------------------------
#adding a strip to the dense mesh
#
#add_strip(mesh, polyedge) takes ONE polyedge = flat list of vertex keys, e.g. [12, 13, 19, 25].
#For several polyedges use add_strips(mesh, [pe1, pe2, ...]).
#The vertex keys must belong to the mesh being edited: densification() rebuilds the mesh from
#Coons patches, so coarse_mesh keys (0..24) do NOT map onto dense_mesh keys (0..424).

quad_mesh = dense_mesh.copy()
quad_mesh.collect_strips()      #add_strip updates attributes['strips'], so it must be populated
quad_mesh_polyedges = dict(quad_mesh.collect_polyedges())

#a polyedge is only valid if it has >2 vertices and is either closed or ends on the boundary twice
valid = [pkey for pkey, pe in quad_mesh_polyedges.items()
         if is_polyedge_valid_for_strip_addition(quad_mesh, pe)]
print("valid polyedges for strip addition:", valid)

for pkey, polyedge in quad_mesh.collect_polyedges():
    pl = Polyline([quad_mesh.vertex_point(v) for v in polyedge])
    viewer_polyedges.add(pl, name=str(pkey))
    viewer_polyedges.add(Tag(text=str(pkey), position=pl.points[0]))
    

# for pkey in valid[:1]:
#     #pass a copy: add_strip pops from the list in place and would corrupt attributes['polyedges']
#     skey, old_to_new = add_strip(quad_mesh, list(quad_mesh_polyedges[pkey]))
#     print("added strip", skey, "along polyedge", pkey)

c = quad_mesh.copy()

#NOTE quad_mesh_polyedges was collected BEFORE any edit. add_strip deletes the vertices it splits,
#so keys stay valid for a second call only if the two polyedges are disjoint (41 and 38 are --
#they are parallel). For crossing polyedges, re-run collect_polyedges() between calls or remap
#through old_to_new.
skey_a, old_to_new_a = add_strip(c, list(quad_mesh_polyedges[41]))
skey_b, old_to_new_b = add_strip(c, list(quad_mesh_polyedges[38]))
print("added strips", skey_a, "and", skey_b, "| faces:", quad_mesh.number_of_faces(), "->", c.number_of_faces())

#add_strip creates the two new vertices ON TOP of the vertex they replace
#(add_vertex(attr_dict=mesh.vertex[v])), so a new strip has zero width and is invisible until the
#mesh is relaxed. old_to_new maps old vkey -> (v1, v2): the pair that has to separate.
new_vertices = {v for mapping in (old_to_new_a, old_to_new_b)
                for pair in mapping.values() for v in pair}

#NOTE vertices_on_boundary() returns only the LONGEST boundary, so the inner hole would be left
#free and smoothing would round it off. vertices_on_boundaries() (plural) gives every ring.
#Pinning every original boundary vertex keeps the outline and the hole; only the newly created
#boundary vertices stay free so the strips can open up.
fixed = [v for ring in c.vertices_on_boundaries() for v in ring if v not in new_vertices]
mesh_smooth_centroid(c, kmax=50, fixed=fixed)

viewer.scene.add(c, name=f"quad_mesh + strips {skey_a}, {skey_b}")



viewer.show()
