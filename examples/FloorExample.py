
import re
from math import pi

import compas
from compas.datastructures.mesh.duality import mesh_dual
from compas.datastructures.mesh.mesh import Mesh
from compas.geometry import Line, Point, Polygon, Polyline
from compas_singular.algorithms import SkeletonDecomposition, boundary_triangulation
from compas_singular.datastructures import (
    CoarsePseudoQuadMesh,
    CoarseQuadMesh,
    QuadMesh,
)
from compas_singular.datastructures.mesh_quad import delete_strips
from compas.colors import Color
from compas.geometry import Translation
from compas_viewer import Viewer
from compas_viewer.scene import Tag
from compas_singular.guide_lines import guide_line_mesh, guide_band_mesh

viewer = Viewer()

#-----------------------------------------
#Define boundary curves
outer = [Point(6,5,0), Point(5,-5,0), Point(-5,-5,0), Point(-6,5,0), Point(0,7,0)]
outer = [Point(6,10,0), Point(6,-10,0), Point(-6,-10,0), Point(-6,10,0)]
inner = []

outer += [outer[0]]

outer_curve = Polyline(outer + [outer[0]])

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
#NOTE while the number is not important anymore when it starts generating correct meshes, it can influence symmetry. Use even values.
outer = densify(outer,20)

#-----------------------------------------
#Start with boundary triangulation
#triangulation.py
def to_coordinates(pts: list[Point]):
    return [[pt.x, pt.y, pt.z] for pt in pts]

trimesh = boundary_triangulation(outer_boundary=outer, inner_boundaries=[])
viewer.scene.add(trimesh)

decomposition = SkeletonDecomposition.from_mesh(trimesh)
coarse_quad: CoarsePseudoQuadMesh= decomposition.decomposition_mesh([])

coarse_quad.collect_strips()
coarse_quad.set_strips_density_target(t=0.5)

from compas.geometry import subtract_vectors, normalize_vector
from math import sin


#-----------------------------------------
#Cable placement: guide_line_mesh puts the cable ON THE MESH EDGES (the joint
#between two blocks); guide_band_mesh with block_rows=1 (its default) puts it
#down the MIDDLE of a row of blocks instead, which is where a cable threaded
#through blocks actually belongs. Two dense meshes, shown side by side, with the
#row of blocks each cable threads through highlighted so the centring is visible
#rather than assumed.

def bow(start, end, bulge, samples=48):
    """A half-sine bow off the chord -- flat in curvature where it meets the
    boundary rather than kinking into it. bulge=0 gives the straight chord."""
    chord = subtract_vectors(end, start)
    normal = normalize_vector([-chord[1], chord[0], 0.0])
    return [[start[i] + chord[i] * (k / samples)
                + normal[i] * bulge * sin(pi * k / samples) for i in range(3)]
            for k in range(samples + 1)]

guides = {
    'straight': Polyline([Point(-8, -1.5, 0), Point(8, -1.5, 0)]),
    'oblique':  Polyline([Point(-8, -4.7, 0), Point(8, -3.3, 0)]),   # ~8.5 deg
    'curved':   Polyline([Point(*p) for p in bow([-5.5, -1.5, 0], [5.5, -1.5, 0], 1.1)]),
}

cable_curves = [guides['curved'], guides['oblique']]

dense_edge, _ = guide_line_mesh(coarse_quad, cable_curves, 0.5)
dense_block, block_info = guide_band_mesh(coarse_quad, cable_curves, 0.5)

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

coarse_quad_improved = coarse_quad.copy()
coarse_quad_improved.set_strip_density_target(skey=3, t=0.6) #Both strips are responsible for one direction of lines in the quads.
coarse_quad_improved.set_strip_density_target(skey=4, t=0.6)

coarse_quad.densification()
coarse_quad_improved.densification()

quad_mesh = coarse_quad.get_quad_mesh()
quad_mesh_improved = coarse_quad_improved.get_quad_mesh()

viewer.scene.add(quad_mesh)
viewer.scene.add(quad_mesh_improved)

step4 = viewer.scene.add_group(name="step4: strips")

#.strips() gives skeys -> .strip_faces(skey) gives the faces in a strip using fkeys -> face_vertices(fkey) gives the vertices of a face
coarse_quad.collect_strips()
for skey in coarse_quad.strips():
    strip = viewer.scene.add_group(name=f"strip skey: {int(skey)}", parent=step4)
    fkeys = coarse_quad.strip_faces(skey)
    print(fkeys)
    union = None
    for fkey in fkeys:
        # fkey -> ordered vertex keys -> vertex coordinates -> Polygon
        vkeys = coarse_quad.face_vertices(fkey)
        face_pts = [Point(*coarse_quad.vertex_coordinates(vkey)) for vkey in vkeys]
        face_polygon = Polygon(face_pts)
        strip.add(face_polygon)
        if union is None:
            union = face_polygon
        else:
            union = union.boolean_union(face_polygon)
    strip.add(union)
    strip.add(Tag(text=str(skey), position=face_polygon.centroid))




viewer.show()