
import compas
from compas_singular.datastructures import CoarseQuadMesh, QuadMesh
from compas_viewer import Viewer
from compas_viewer.scene import Tag
from compas.datastructures.mesh.mesh import Mesh
from compas.datastructures.mesh.duality import mesh_dual
from compas.geometry import Point, Polyline, Line, Polygon

from compas_singular.algorithms import boundary_triangulation, SkeletonDecomposition
from compas_singular.datastructures import CoarsePseudoQuadMesh
from compas_singular.datastructures.mesh_quad import delete_strips

from math import pi
import re

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

point_features = [[2.5,3.5,0]]
polyline_features = []

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

from math import ceil
from compas.geometry import vector_average

def density_by_target(skey, target_length):
    lengths = [coarse_mesh.edge_length(u, v)
               for u, v in coarse_mesh.strip_edges(skey) if u != v]
    return ceil(vector_average(lengths) / target_length)

coarse_mesh.set_strips_density_func(func=density_by_target, func_args=0.5)
coarse_mesh.densification()
dense_mesh: QuadMesh = coarse_mesh.get_quad_mesh()
step3.add(dense_mesh)

def density_clamped(skey, args):
    mesh, target, d_min, d_max = args
    lengths = [mesh.edge_length(u, v) for u, v in mesh.strip_edges(skey) if u != v]
    d = ceil(vector_average(lengths) / target)
    return max(d_min, min(d_max, d))

coarse_mesh.set_strips_density_func(
    func=density_clamped,
    func_args=(coarse_mesh, 0.3, 1, 8)
)
coarse_mesh.densification()
dense_mesh: QuadMesh = coarse_mesh.get_quad_mesh()
step3.add(dense_mesh)

coarse_mesh.set_strips_density_target(t=0.5)
coarse_mesh.densification()
dense_mesh: QuadMesh = coarse_mesh.get_quad_mesh()
step3.add(dense_mesh, name="target length")

#-----------------------------------------
#strip visualisation
step4 = viewer.scene.add_group(name="step4: strips")

#.strips() gives skeys -> .strip_faces(skey) gives the faces in a strip using fkeys -> face_vertices(fkey) gives the vertices of a face
dense_mesh.collect_strips()
for skey in coarse_mesh.strips():
    strip = viewer.scene.add_group(name=f"strip skey: {int(skey)}", parent=step4)
    fkeys = coarse_mesh.strip_faces(skey)
    print(fkeys)
    union = None
    for fkey in fkeys:
        # fkey -> ordered vertex keys -> vertex coordinates -> Polygon
        vkeys = coarse_mesh.face_vertices(fkey)
        face_pts = [Point(*coarse_mesh.vertex_coordinates(vkey)) for vkey in vkeys]
        face_polygon = Polygon(face_pts)
        strip.add(face_polygon)
        if union is None:
            union = face_polygon
        else:
            union = union.boolean_union(face_polygon)
    strip.add(union)
    strip.add(Tag(text=str(skey), position=face_polygon.centroid))



viewer.show()


#-----------------------------------------
#Interactive strip-based mesh editing

def _get_strip_indices(args: str) -> list[int]:
    """Extract the strip indices from an argument string like '-[0, 2, 4]'."""
    match = re.search(r"\[(.*?)\]", args)
    if not match:
        raise ValueError("expected strip indices as a list, e.g. -[0, 2, 4]")
    items = [item.strip() for item in match.group(1).split(",")]
    return [int(item) for item in items if item]


def _show_mesh_with_strips(viewer: Viewer, mesh) -> None:
    """Clear the viewer and redraw the mesh with a tag on every current strip."""
    group = viewer.scene.add_group(name="mesh editor")
    group.add(mesh)

    strips_group = viewer.scene.add_group(name="strips", parent=group)
    for skey in mesh.strips():
        strip = viewer.scene.add_group(name=f"strip skey: {int(skey)}", parent=strips_group)
        union = None
        for fkey in mesh.strip_faces(skey):
            vkeys = mesh.face_vertices(fkey)
            face_pts = [Point(*mesh.vertex_coordinates(vkey)) for vkey in vkeys]
            face_polygon = Polygon(face_pts)
            strip.add(face_polygon)
            union = face_polygon if union is None else union.boolean_union(face_polygon)
        if union is not None:
            strip.add(Tag(text=str(skey), position=union.centroid))


def mesh_editor(mesh, viewer: Viewer = None) -> Viewer:
    """Interactively edit a quad mesh strip-by-strip.

    Shows the mesh labelled with its current strip indices in the viewer, then
    prompts in the terminal for an edit. Repeats, showing the updated mesh
    after every edit, until the user types 'f' to finish.

    Parameters
    ----------
    mesh : CoarseQuadMesh | CoarsePseudoQuadMesh
        The mesh to edit. Modified in place.
    viewer : Viewer, optional
        Viewer instance to reuse. A new one is created if not given.

    Returns
    -------
    Viewer
        The viewer used, left showing the final state of the mesh.

    """


    if not list(mesh.strips()):
        mesh.collect_strips()

    while True:
        _show_mesh_with_strips(viewer, mesh)

        viewer.show()
        viewer.scene.remove(step1)
        viewer.scene.remove(step2)
        viewer.scene.remove(step3)
        print(f"Current strip indices: {sorted(mesh.strips())}")
        print("Options:")
        print("  d -[strip_index, strip_index, ...]   delete strips")
        print("  f                                    finished editing")
        raw = input("Edit option: ").strip()

        if raw.lower() == "f":
            break

        command, _, args = raw.partition(" ")
        command = command.strip().lower()

        if command == "d":
            try:
                skeys = _get_strip_indices(args)
            except ValueError as error:
                print(f"Could not parse strip indices: {error}")
                continue

            valid_skeys = list(mesh.strips())
            unknown = [skey for skey in skeys if skey not in valid_skeys]
            if unknown:
                print(f"Unknown strip indices, ignoring: {unknown}")
            skeys = [skey for skey in skeys if skey in valid_skeys]

            delete_strips(mesh=mesh, skeys=skeys)
        else:
            print(f"Unknown option: {command!r}")

    return viewer


mesh_editor(coarse_mesh, viewer=viewer)

