
import os
import re
from math import atan2, pi, sin

import compas
from compas.colors import Color
from compas.datastructures import Mesh
from compas.datastructures.mesh.duality import mesh_dual
from compas.datastructures.mesh.mesh import Mesh
from compas.files import OBJ
from compas.geometry import (
    Frame,
    Line,
    NurbsCurve,
    Plane,
    Point,
    Polygon,
    Polyline,
    Translation,
    intersection_segment_polyline_xy,
    normalize_vector,
    subtract_vectors,
)
from compas.geometry.transformation import Transformation
from compas.tolerance import TOL
from compas_singular.algorithms import SkeletonDecomposition, boundary_triangulation
from compas_singular.datastructures import (
    CoarsePseudoQuadMesh,
    CoarseQuadMesh,
    QuadMesh,
)
from compas_singular.datastructures.mesh_quad import delete_strips
from compas_singular.guide_lines import guide_band_mesh, guide_line_mesh
from compas_viewer import Viewer
from compas_viewer.scene import Tag

viewer = Viewer()


def boundary_edges_to_curves(coarse, decomposition):
    """Map every edge of a coarse (pseudo) quad mesh to the actual boundary/skeleton
    polyline it was collapsed from during SkeletonDecomposition.decomposition_mesh().

    `decomposition_mesh()` builds the coarse mesh via `Mesh.from_polylines`, which
    keeps only each polyline's two endpoints as coarse vertices -- so by default
    `densification()` draws a straight 2-point chord per coarse edge and any curve
    the boundary had between two coarse (corner) vertices is lost. The original
    points are still recoverable from `decomposition.polylines`
    (`decomposition.decomposition_polyline`), so we look them back up here and feed
    them to `densification(edges_to_curves=...)` instead.

    Every coarse edge must be present or `densification` KeyErrors on the ones
    missing, so interior (non-boundary) edges fall back to their own straight chord.
    """
    edges_to_curves = {}
    for u, v in coarse.edges():
        gk_u = TOL.geometric_key(coarse.vertex_coordinates(u))
        gk_v = TOL.geometric_key(coarse.vertex_coordinates(v))
        curve = decomposition.decomposition_polyline(gk_u, gk_v)
        if curve is None:
            curve = [coarse.vertex_coordinates(u), coarse.vertex_coordinates(v)]
        elif TOL.geometric_key(curve[0]) != gk_u:
            # decomposition_polyline() returns the stored polyline regardless of
            # which requested key matched which end -- reverse it if it runs v -> u.
            curve = list(reversed(curve))
        edges_to_curves[u, v] = curve
    return edges_to_curves


def quad_dual(mesh: Mesh) -> QuadMesh:
    """Dual of an all-quad mesh, closed on the boundary with QUADS.

    Every primal VERTEX becomes a dual FACE: a regular valence-4 vertex gives a
    quad, so the dual is still a quad mesh everywhere the primal was regular, and
    a singular valence-n vertex gives an n-gon. That is the point of dualising
    here -- a singularity stops being a joint where n mesh lines gather (which a
    block system cannot build) and becomes a single n-sided BLOCK. The index is
    conserved, only which element carries it changes.

    `compas.datastructures.Mesh.dual(include_boundary=True)` always inserts the
    boundary vertex itself into its closure face, which turns every REGULAR
    (valence-3) boundary vertex into a pentagon -- 44 of them on the pentagon
    plate. Here the vertex is inserted only at a domain CORNER (valence 2), which
    is the one place it is needed to reach four sides, so a regular boundary ring
    comes out all quads. Both closures put their boundary on primal edge midpoints
    and corners, which lie exactly ON the primal boundary, so the outline is kept.
    """
    index: dict[str, int] = {}
    points: list[list[float]] = []

    def add(xyz) -> int:
        gk = TOL.geometric_key(xyz)
        if gk not in index:
            index[gk] = len(points)
            points.append([float(c) for c in xyz])
        return index[gk]

    faces = []
    for v in mesh.vertices():
        #ordered=True gives the incident faces as a cycle, so the dual face comes
        #out wound consistently; for a boundary vertex it starts and ends at the wall
        fkeys = [f for f in mesh.vertex_faces(v, ordered=True) if f is not None]
        if not mesh.is_vertex_on_boundary(v):
            faces.append([add(mesh.face_centroid(f)) for f in fkeys])
            continue
        nbrs = mesh.vertex_neighbors(v, ordered=True)
        ring = [add(mesh.edge_midpoint((v, nbrs[0])))]
        ring += [add(mesh.face_centroid(f)) for f in fkeys]
        ring.append(add(mesh.edge_midpoint((v, nbrs[-1]))))
        if len(nbrs) == 2:
            ring.append(add(mesh.vertex_coordinates(v)))
        faces.append(ring)

    dual = QuadMesh.from_vertices_and_faces(points, faces)
    dual.unify_cycles()
    return dual


#-----------------------------------------
#Define boundary curves
outer = [Point(6,5,0), Point(5,-5,0), Point(-5,-5,0), Point(-6,5,0), Point(0,7,0)]
outer = [Point(6,10,0), Point(6,-10,0), Point(-6,-10,0), Point(-6,10,0)]
inner = []


HERE = os.path.dirname(__file__)
FILE = os.path.join(HERE, 'curve_outline.obj')
obj = OBJ(FILE)
obj.read()

outer = []
for v in obj.vertices:
    outer.append(Point(*v))


ref = Frame.from_plane(Plane.from_points(outer))
for pt in outer:
    T = Transformation.from_frame_to_frame(ref, Frame.worldXY())
    pt.transform(T)

#The obj has no edge/connectivity info (just "p" point elements), and Rhino
#exported them out of curve order, so recover the loop by angle around the
#centroid -- valid because the flattened points are star-shaped around it
#(verified: this gives near-uniform segment lengths, nearest-neighbour chaining
#does not).
def order_closed_points(points: list[Point]) -> list[Point]:
    cx = sum(pt.x for pt in points) / len(points)
    cy = sum(pt.y for pt in points) / len(points)
    return sorted(points, key=lambda pt: atan2(pt.y - cy, pt.x - cx))

outer = order_closed_points(outer)

outer += [outer[0]]

outer_curve = Polyline(outer + [outer[0]])
def to_coordinates(pts: list[Point]):
    return [[pt.x, pt.y, pt.z] for pt in pts]

from compas_singular.algorithms import boundary_triangulation
from compas_singular.algorithms import SkeletonDecomposition

trimesh = boundary_triangulation(outer_boundary=outer, inner_boundaries=[])
decomposition = SkeletonDecomposition.from_mesh(trimesh)
coarse = decomposition.decomposition_mesh([])
coarse.collect_strips()
coarse.set_strips_density_target(t=0.5)
coarse.densification() #the density is stored in the object, but not yet visible.
dense_mesh: QuadMesh = coarse.get_quad_mesh() #Move from pseudo quad to quad
viewer.scene.add(dense_mesh)
viewer.show()

def bow(start, end, bulge, samples=48):
    """A half-sine bow off the chord -- flat in curvature where it meets the
    boundary rather than kinking into it. bulge=0 gives the straight chord."""
    chord = subtract_vectors(end, start)
    normal = normalize_vector([-chord[1], chord[0], 0.0])
    return [[start[i] + chord[i] * (k / samples)
                + normal[i] * bulge * sin(pi * k / samples) for i in range(3)]
            for k in range(samples + 1)]

guides = {
    'straight': Polyline([Point(0.22, 4, 0), Point(-0.22, -4, 0)]),
    'oblique':  Polyline([Point(-1, -3.8, 0), Point(3, -3.1, 0)]),   # ~10 deg
    'curved':   Polyline([Point(*p) for p in bow([-3.5, -1.5, 0], [3.5, -1.5, 0], 1.1)]),
}

cable_curves = [guides['straight'], guides["curved"]]

for polyline in list(guides.values()):
    viewer.scene.add(polyline)

intersection = None
for segment in cable_curves[0].lines:
    viewer.scene.add(segment)
    # NOTE using the 3D intersection_segment_polyline function return (None, None) for some reason.
    intersection =intersection_segment_polyline_xy(segment, cable_curves[1])
    print(intersection)
    if intersection != (None, None):
        print("intersection found")
        viewer.scene.add(Point(*intersection))
        break
point_features = []
point_features.append(intersection)




#Note that the outer boundary should be singular, but the inenr boundary can be a list.
trimesh = boundary_triangulation(outer_boundary=to_coordinates(outer), inner_boundaries=[], point_features=point_features)


#-----------------------------------------
#Skeleton decomposition: turns trimesh into a coarse quad-patch layout -> CoarsePseudoQuadMesh
#Uses same feature values as in trimesh
#decomposition.py

decomposition = SkeletonDecomposition.from_mesh(trimesh) #plotting this returns basic trimesh, but it contains the info and methods
viewer.scene.add(trimesh)

coarse_mesh: CoarsePseudoQuadMesh = decomposition.decomposition_mesh(point_features) #decomposed CoarsePseudoQuadMesh  derived from SkeletonDecomposition
#Each path is a Coons-patch placeholder. Coons patch helps divide complicated area geometries
viewer.scene.add(coarse_mesh)
#-----------------------------------------
#density -> CoarseQuadMesh

#this methods collects the strips and returns skey and vertix indices. Afterwards skey can be collected using .strips
coarse_mesh.collect_strips()


#using skey the density of each strip can also be set individually
# coarse_mesh.set_strips_density_target(t=5)
# coarse_mesh.set_strips_density_func(func=, func_args=)
coarse_mesh.set_strips_density_target(t=5)

#Without edges_to_curves, densification() draws a straight chord between each pair
#of coarse vertices -- that's why the coarse/dense boundary looks faceted even
#though the source trimesh boundary was smooth. Feed back the real boundary curve
#per coarse edge (recovered from the decomposition) so the dense boundary follows it.
edges_to_curves = boundary_edges_to_curves(coarse_mesh, decomposition)
coarse_mesh.densification(edges_to_curves=edges_to_curves) #the density is stored in the object, but not yet visible.
dense_mesh: QuadMesh = coarse_mesh.get_quad_mesh() #Move from pseudo quad to quad


viewer.scene.add(dense_mesh)



#-----------------------------------------
#Cable placement: guide_line_mesh puts the cable ON THE MESH EDGES (the joint
#between two blocks); guide_band_mesh with block_rows=1 (its default) puts it
#down the MIDDLE of a row of blocks instead, which is where a cable threaded
#through blocks actually belongs. Two dense meshes, shown side by side, with the
#row of blocks each cable threads through highlighted so the centring is visible
#rather than assumed.



cable_curves = [guides['straight'], guides["curved"]]

#`boundary` defaults to coarse_mesh's own straight-chorded polygon, which is what
#snap/relax/swing measure fold-checks and boundary reprojection against -- so even
#with edges_to_curves keeping the new mesh's OWN boundary curved, any vertex that
#relax/swing touches gets pulled back onto the straight polygon. Pass the real fine
#boundary (still held in `outer`, pre-decomposition) so that reprojection follows
#the true curve too.
fine_boundary = to_coordinates(outer)
dense_edge, _ = guide_line_mesh(coarse_mesh, cable_curves, 0.5, edges_to_curves=edges_to_curves, boundary=fine_boundary)

#The band cannot seat BOTH of these cables on this plate -- one band blocks every
#slot the other could use (add_strip deletes and duplicates the vertices of the
#polyedge it inflates, so bands are not allowed to cross or touch). The line seats
#both, and the dual section below is what then gets them mid-block anyway. Fail
#soft rather than killing the script, so the comparison still shows what it can.
try:
    dense_block, block_info = guide_band_mesh(coarse_mesh, cable_curves, 0.5, edges_to_curves=edges_to_curves, boundary=fine_boundary)
except ValueError as exc:
    print(f'guide_band_mesh refused these cables: {exc}')
    dense_block, block_info = None, None

def face_centroid(mesh, fkey):
    pts = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
    return [sum(p[i] for p in pts) / len(pts) for i in range(3)]

def closest_on_polyline(point, poly_pts):
    """(distance, foot) of the closest point on a chain of points."""
    best_d, best_q = float('inf'), None
    for a, b in zip(poly_pts[:-1], poly_pts[1:]):
        ab = subtract_vectors(b, a)
        length2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2]
        if length2 == 0.0:
            continue
        ap = subtract_vectors(point, a)
        t = max(0.0, min(1.0, (ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / length2))
        q = [a[i] + ab[i] * t for i in range(3)]
        d = sum((point[i] - q[i]) ** 2 for i in range(3)) ** 0.5
        if d < best_d:
            best_d, best_q = d, q
    return best_d, best_q

def cable_block_faces(mesh, centreline_pts):
    """Faces the cable actually threads through -- its centroid closer to the
    cable than the face's own half-diagonal, which is what 'runs through this
    block' means."""
    faces = []
    for fkey in mesh.faces():
        pts = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
        c = face_centroid(mesh, fkey)
        d, _ = closest_on_polyline(c, centreline_pts)
        radius = max(((c[0] - p[0]) ** 2 + (c[1] - p[1]) ** 2 + (c[2] - p[2]) ** 2) ** 0.5
                     for p in pts)
        if d < radius:
            faces.append(fkey)
    return faces

cable_compare = viewer.scene.add_group(name="cable placement: edge vs mid-block")

g_edge = viewer.scene.add_group(name="guide_line_mesh -- cable ON THE EDGE",
                                parent=cable_compare)
g_edge.add(dense_edge, show_points=False, show_lines=True, show_faces=True,
          name="dense mesh")
for cable_name, curve in (('curved', guides['curved']), ('oblique', guides['oblique'])):
    g_edge.add(curve, linecolor=Color.red(), linewidth=4, name=f"cable: {cable_name}")


shift = Translation.from_vector([0.0, 26.0, 0.0])
if dense_block is not None:
    g_block = viewer.scene.add_group(name="guide_band_mesh -- cable MID-BLOCK (block_rows=1)",
                                     parent=cable_compare)
    g_block.add(dense_block.transformed(shift), show_points=False, show_lines=True,
               show_faces=True, name="dense mesh")
    for cable_name, pts in zip(('curved', 'oblique'), block_info['centrelines']):
        curve = Polyline([Point(*p) for p in pts]).transformed(shift)
        g_block.add(curve, linecolor=Color.red(), linewidth=4, name=f"cable: {cable_name}")
        for fkey in cable_block_faces(dense_block, pts):
            face_pts = [Point(*dense_block.vertex_coordinates(v))
                       for v in dense_block.face_vertices(fkey)]
            block = Polygon(face_pts).transformed(shift)
            g_block.add(block, surfacecolor=Color.orange(), name=f"block threaded by {cable_name}")


#-----------------------------------------
#Dual: move the singularity from the JOINTS to the BLOCKS.
#
#A singularity in the quad mesh is a vertex where 3 or 5 mesh lines gather. The
#block system needs a BLOCK there, not a joint. Dualising gives exactly that: the
#valence-n vertex becomes an n-gon face, while every regular valence-4 vertex
#becomes a quad -- so nothing else about the mesh changes except a half-block
#shift of the whole grid. There is no need to dualise only locally.
#
#Duality also SWAPS which guide function centres the cable. A primal edge-chain
#dualises to a ROW OF BLOCKS running along it, so guide_line_mesh -- which puts
#the cable on the joint in the primal -- puts it MID-BLOCK in the dual. That
#matters because guide_line_mesh is the only one that supports CROSSING cables;
#guide_band_mesh refuses them outright. So line + dual is how a crossing cable
#gets to be mid-block, which neither function can manage on its own.
#
#CAVEAT, and the printout below shows it: a POLE is already the dual situation --
#it is a triangular face, not a gathering joint -- so dualising trades it the
#WRONG way, turning each pole into a valence-3 joint. `dense_mesh` comes from a
#CoarsePseudoQuadMesh with a point_feature, so it has poles and its dual keeps
#some irregular joints. `dense_edge` is genuinely all-quad, and its dual reaches
#zero: every one of its singularities becomes a block. Dualise all-quad meshes.

dual_plain = quad_dual(dense_mesh)
dual_cable = quad_dual(dense_edge)


def block_degrees(mesh):
    counts = {}
    for fkey in mesh.faces():
        n = len(mesh.face_vertices(fkey))
        counts[n] = counts.get(n, 0) + 1
    return dict(sorted(counts.items()))


def irregular_joints(mesh):
    return [v for v in mesh.vertices()
            if not mesh.is_vertex_on_boundary(v) and len(mesh.vertex_neighbors(v)) != 4]


for label, primal, dual in (('plain', dense_mesh, dual_plain), ('cable', dense_edge, dual_cable)):
    print(f"{label}: primal {len(irregular_joints(primal))} irregular JOINTS, blocks {block_degrees(primal)}"
          f"  ->  dual {len(irregular_joints(dual))} irregular JOINTS, blocks {block_degrees(dual)}")

dual_group = viewer.scene.add_group(name="dual: singularity as a BLOCK, cable mid-block")

g_dual = viewer.scene.add_group(name="dual of the plain mesh", parent=dual_group)
g_dual.add(dual_plain, show_points=False, show_lines=True, show_faces=True, name="dual mesh")
#the non-quad faces are exactly the singularities, now buildable as single blocks
for fkey in dual_plain.faces():
    if len(dual_plain.face_vertices(fkey)) != 4:
        pts = [Point(*dual_plain.vertex_coordinates(v)) for v in dual_plain.face_vertices(fkey)]
        g_dual.add(Polygon(pts), surfacecolor=Color.cyan(),
                   name=f"singular block ({len(pts)} sides)")

dual_shift = Translation.from_vector([0.0, -26.0, 0.0])
g_dual_cable = viewer.scene.add_group(name="dual of guide_line_mesh -- cable now MID-BLOCK",
                                      parent=dual_group)
g_dual_cable.add(dual_cable.transformed(dual_shift), show_points=False, show_lines=True,
                 show_faces=True, name="dual mesh")
for cable_name, curve in zip(('straight', 'curved'), cable_curves):
    pts = to_coordinates(curve.points)
    g_dual_cable.add(Polyline([Point(*p) for p in pts]).transformed(dual_shift),
                     linecolor=Color.red(), linewidth=4, name=f"cable: {cable_name}")
    for fkey in cable_block_faces(dual_cable, pts):
        face_pts = [Point(*dual_cable.vertex_coordinates(v))
                    for v in dual_cable.face_vertices(fkey)]
        g_dual_cable.add(Polygon(face_pts).transformed(dual_shift),
                         surfacecolor=Color.orange(), name=f"block threaded by {cable_name}")


coarse_mesh_improved = coarse_mesh.copy()
coarse_mesh_improved.set_strip_density_target(skey=3, t=0.6) #Both strips are responsible for one direction of lines in the quads.
coarse_mesh_improved.set_strip_density_target(skey=4, t=0.6)

coarse_mesh.densification(edges_to_curves=edges_to_curves)
coarse_mesh_improved.densification(edges_to_curves=boundary_edges_to_curves(coarse_mesh_improved, decomposition))

quad_mesh = coarse_mesh.get_quad_mesh()
quad_mesh_improved = coarse_mesh_improved.get_quad_mesh()

viewer.scene.add(quad_mesh)
viewer.scene.add(quad_mesh_improved)

step4 = viewer.scene.add_group(name="step4: strips")

#.strips() gives skeys -> .strip_faces(skey) gives the faces in a strip using fkeys -> face_vertices(fkey) gives the vertices of a face
coarse_mesh.collect_strips()
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