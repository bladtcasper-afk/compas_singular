import compas

from compas.geometry import is_polygon_in_polygon_xy, is_point_in_polygon_xy, Polygon, Polyline
from compas_rhino.conversions import point_to_rhino, polyline_to_rhino
from compas_rhino.conversions import vertices_and_faces_to_rhino
from compas_rhino.conversions import mesh_to_compas
from compas_rhino.conversions import point_to_compas
import compas_rhino as cr
from compas_singular.datastructures import CoarsePseudoQuadMesh

if compas.RHINO:
	import rhinoscriptsyntax as rs
	import scriptcontext as sc
	import System


def clear_layer(layer, clean_sublayers=False):
	if not rs.IsLayer(layer):
		print(f"No layer named {layer}. Deleted nothing.")
		return 0
	
	layer = rs.LayerName(layer, fullpath = True)
	count = 0 # Number of deleted objects
	objs = rs.ObjectsByLayer(layer)
	if objs:
		rs.DeleteObjects(objs)
		count += len(objs)

	if clean_sublayers:
		all_layers = rs.LayerNames()
		if all_layers:
			prefix = layer + '::'
			for other in all_layers:
				if other != layer and other.startswith(prefix):
					objs = rs.ObjectsByLayer(other)
					if objs:
						rs.DeleteObjects(objs)
						count += len(objs)

	return count


def bake_points(points, layer, color=None, clear_existing=True):
	if not rs.IsLayer(layer):
		rs.AddLayer(layer, color)
	if clear_existing:
		clear_layer(layer)

	guids = rs.AddPoints([point_to_rhino(pt) for pt in points])
	for guid in guids:
		rs.ObjectLayer(guid, layer=layer)
	return guids


def clean_polyline_points(points, tol=None):
	"""``points`` with consecutive duplicates dropped, at RHINO's tolerance.

	Rhino rejects a polyline whose consecutive points coincide, and it judges
	that at the document's absolute tolerance -- not at 1e-9. A separatrix that
	doubles back on itself by half a document tolerance is geometrically fine
	and still unbakeable, so the cleaning has to use the same yardstick Rhino
	will.

	Returns ``[]`` when nothing bakeable is left, which the caller must treat as
	"skip this one" rather than as an error.
	"""
	if tol is None:
		try:
			tol = sc.doc.ModelAbsoluteTolerance
		except Exception:
			tol = 1e-6
	out = []
	for point in points:
		point = list(point)[:3]
		if out and _distance(out[-1], point) <= tol:
			continue
		out.append(point)
	return out if len(out) >= 2 else []


def bake_polylines(polylines, layer, color=None, clear_existing=True):
	"""Bake polylines, skipping any Rhino will not take. ``(guids, skipped)``.

	**``rs.AddPolyline`` RAISES** -- ``Unable to add polyline to document`` --
	when Rhino refuses the geometry; it does not return ``None``. That matters
	more than it sounds: these bakes run AFTER the mesh has been written, so an
	exception here leaves the document half updated, with this layer cleared and
	not repopulated, and the next command then silently finds no polylines and
	densifies every edge as a chord.

	So each one is cleaned first and the add is guarded anyway, and the number
	skipped is RETURNED rather than swallowed -- a caller that bakes fewer
	curves than it was given should say so.
	"""
	if not rs.IsLayer(layer):
		rs.AddLayer(layer, color)
	if clear_existing:
		clear_layer(layer)

	guids = []
	skipped = 0
	for polyline in polylines:
		points = clean_polyline_points(polyline)
		if not points:
			skipped += 1
			continue
		try:
			guid = rs.AddPolyline(polyline_to_rhino(Polyline(points)))
		except Exception:
			guid = None
		if not guid:
			skipped += 1
			continue
		rs.ObjectLayer(guid, layer=layer)
		guids.append(guid)
	return guids, skipped


def bake_mesh(mesh, layer, color=None, clear_existing=True):
	"""Bake a compas mesh, n-gons included. Returns the guid.

	**Not ``rs.AddMesh``.** That function takes a face of any length and does
	``AddFace(face[0], face[1], face[2], face[3])`` -- so a pentagon is baked as
	the QUAD of its first four corners and the remaining triangle is a hole in
	the middle of the mesh. Nothing reports it: the result is a valid Rhino mesh,
	just not the one that was handed in. The dual of any quad mesh with a
	singularity has such faces (a valence-5 joint dualises to a pentagon, a
	boundary singularity likewise), which is where it shows.

	``vertices_and_faces_to_rhino`` is what ``scene.draw`` already uses, and it
	fans a face of more than four corners into triangles around its centroid and
	then groups them as a ``MeshNgon``, so Rhino displays and snaps to the
	polygon while the geometry underneath is complete. ``disjoint=False`` keeps
	the vertices WELDED -- the default ``True`` gives every face its own copy of
	its corners, and this workflow reads its meshes back out of the document,
	where an unwelded mesh has a boundary around every face.
	"""
	if not rs.IsLayer(layer):
		rs.AddLayer(layer, color, parent="TopologyProblem")
		print("Layer did not exist. Added the layer.")
	if clear_existing:
		clear_layer(layer)

	# ``face_vertices`` gives vertex KEYS while ``vertices`` is positional, and
	# the two only agree while the keys happen to be 0..n-1. A repair, a weld or
	# a deleted patch leaves gaps, and then the baked faces reference the wrong
	# corners -- silently, as a valid mesh of the wrong shape. Map explicitly.
	index = {v: i for i, v in enumerate(mesh.vertices())}
	vertices = [mesh.vertex_attributes(v, "xyz") for v in mesh.vertices()]
	faces = [[index[v] for v in mesh.face_vertices(f)] for f in mesh.faces()]

	geometry = vertices_and_faces_to_rhino(vertices, faces, disjoint=False)
	guid = sc.doc.Objects.AddMesh(geometry)
	if guid == System.Guid.Empty:
		raise Exception("Unable to add the mesh to the document.")
	sc.doc.Views.Redraw()
	rs.ObjectLayer(guid, layer=layer)
	return guid

#: Fewest points a CLOSED curved input is divided into, whatever ``max_edge``
#: says. A curve is divided by CHORD LENGTH here, so the count a loop gets is
#: its own length over the background spacing -- and a small feature is short.
#: At the default spacing of 0.5 a round hole of radius 0.64 or less came out as
#: a regular OCTAGON: 45 degrees of turn per vertex, which is exactly the angle
#: a cross field cannot tell from a corner (see ``13_limits.py``). Rhino's own
#: Convert samples by DEVIATION instead and never does this, which is why the
#: same boundary worked once it had been converted by hand first.
#:
#: 28 caps the turn at 12.9 degrees, comfortably clear of that limit, and is
#: measured: it puts every curved domain in the test set on the field route,
#: while 32 and above start costing background triangle quality -- the loop's
#: chords fall well below ``target_length`` and ``background._densify_loop``
#: only ever subdivides, so it cannot coarsen them back. Measured minimum
#: background angle on a rounded plate with a radius 0.5 hole: 20.8 deg at the
#: old floor of 8, 11.0 at 28, 9.3 at 32, 6.3 at 48.
MIN_CLOSED_CURVE_POINTS = 28

#: A CLOSED loop's point count is rounded UP to a multiple of this.
#:
#: ``framefield.symmetry.Symmetry.detect`` matches by POINT SET at a 1e-6
#: tolerance, so a symmetric drawing is only detected as symmetric if its
#: SAMPLES are symmetric too. Dividing a circle into ``n`` equal parts gives a
#: set closed under a quarter turn only when ``n % 4 == 0``; ``n`` even gives
#: order 4 and ``n`` odd order 2. The count came straight from
#: ``length / max_edge``, so it was effectively arbitrary: measured on a square
#: plate with a centred hole at spacing 0.5, radii 0.5-2.0 detected the full D4
#: (count 28, from the floor above, which happens to be a multiple of 4) while
#: r=2.5 gave count 31 and collapsed to order 2, and r=3.0 gave count 38 and
#: order 4. Rounding up costs at most three points.
#:
#: This fixes the COUNT, not the PHASE. ``DivideByCount`` starts at the curve's
#: seam, so the same 28-point circle rotated 5 degrees still measures order 4:
#: the mirror axes no longer land on samples. Draw holes with the seam on a
#: symmetry axis (Rhino's own circles have theirs at angle 0, which is why this
#: is rarely hit), or pass an explicit group instead of ``'auto'``.
SYMMETRY_SAMPLE_MULTIPLE = 4


def curve_points(guid, max_edge):
	"""A Rhino curve as a list of points: corners kept, curvature sampled.

	A polyline is taken at its OWN vertices, because sampling it would round the
	corners and a corner is exactly what the decomposition has to see as a
	corner. Anything else -- an arc, a circle, a NURBS curve -- is divided by
	length, because every consumer downstream treats a loop as straight segments
	between the points given, so an arc handed over as two points IS a chord.

	This is what ``curve_to_compas_polyline`` cannot do: it calls
	``TryGetPolyline`` and RAISES ``ConversionError`` when the curve is not one,
	which is why a genuinely curved outer boundary could not be used at all.

	A CLOSED curve is never divided into fewer than
	:data:`MIN_CLOSED_CURVE_POINTS` points, however coarse ``max_edge`` is --
	see that constant for the measurement. Division by chord length ties a
	loop's point count to its own LENGTH, so a small hole came out with eight
	points no matter what, and the field front end cannot read a shape that
	coarse as curved.

	A closed loop's count is then rounded up to a multiple of
	:data:`SYMMETRY_SAMPLE_MULTIPLE` so that symmetry detection can see the
	symmetry the user drew -- see that constant.

	The closing point is dropped -- the loop is closed by convention, and a
	repeated point puts a zero-length segment in every arc that crosses it.
	"""
	curve = rs.coercecurve(guid)
	ok, polyline = curve.TryGetPolyline()
	if ok:
		points = [[p.X, p.Y, 0.0] for p in polyline]
	else:
		count = int(round(curve.GetLength() / max(max_edge, 1e-6)))
		# A closed loop gets a floor of its own; an open guide does not need
		# one, having no winding for the field to read.
		count = max(MIN_CLOSED_CURVE_POINTS if curve.IsClosed else 8, count, 8)
		if curve.IsClosed:
			count = -(-count // SYMMETRY_SAMPLE_MULTIPLE) * SYMMETRY_SAMPLE_MULTIPLE
		points = []
		for t in curve.DivideByCount(count, True):
			p = curve.PointAt(t)
			points.append([p.X, p.Y, 0.0])
	while len(points) > 1 and _distance(points[0], points[-1]) < 1e-9:
		points = points[:-1]
	return points


def _distance(a, b):
	return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5


def curve_to_polyline(guid, max_edge):
	""":func:`curve_points` as a CLOSED compas :class:`Polyline`.

	Closed in the ``first point repeated at the end`` sense, which is what
	``curve_to_compas_polyline`` returns for a closed Rhino polyline and what
	``Polyline.is_closed`` and the existing callers expect. For an input that
	really is a polyline this returns exactly what ``curve_to_compas_polyline``
	did, point for point -- the only thing that changes is that a curved input
	now works instead of raising.
	"""
	points = curve_points(guid, max_edge)
	return Polyline(points + points[:1])


def read_boundary_loops(max_edge):
	"""``(outer, inners)`` as point lists, sampled from the Rhino curves.

	The loops for :func:`compas_singular.rhino.coarse_curves.coarse_edges_to_curves`.
	Sample FINER than the background triangulation: these points are the wall
	from here on, and ``densification`` walks them with ``Polyline.point_at``,
	so the final mesh's boundary is only as smooth as what is handed in. The
	bound is one segment's sagitta -- measured 1.5e-4 on a radius-5 disc sampled
	at 400 points, against 1.46 units for the straight-chord densification it
	replaces.
	"""
	outer_ids = rs.ObjectsByLayer("Outer")
	if not outer_ids:
		raise RuntimeError("No outer boundary selected.")

	outer = curve_points(outer_ids[0], max_edge)
	inners = [curve_points(guid, max_edge)
			  for guid in rs.ObjectsByLayer("Inner") or []]
	return outer, inners


def read_polylines(layer):
	"""Every polyline on a layer as a point list. ``[]`` if the layer is empty.

	Used for the decomposition's branches on ``Skeleton::Polylines``. Unlike the
	coarse MESH, these are baked as curves and so keep full double precision.
	"""
	out = []
	for guid in rs.ObjectsByLayer(layer) or []:
		curve = rs.coercecurve(guid)
		if curve is None:
			continue
		ok, polyline = curve.TryGetPolyline()
		if not ok:
			continue
		out.append([[p.X, p.Y, p.Z] for p in polyline])
	return out


def bake_edge_curves(curves, layer, color=None):
	"""Bake one polyline per coarse edge -- the layout WITH its curvature.

	The coarse layout is baked as a Rhino mesh, and a Rhino mesh has straight
	edges: on a curved domain it is drawn cutting the corner off its own
	boundary, which reads as the workflow having lost the curve when in fact
	only the display has. These are what the layout actually densifies along.

	Returns ``(guids, skipped)`` -- see :func:`bake_polylines`. An edge whose
	curve Rhino will not take densifies as a chord, which is worth printing
	rather than losing.
	"""
	polylines = [Polyline([list(p) for p in curve]) for curve in curves
				 if len(curve) >= 2]
	return bake_polylines(polylines, layer, color=color)


def read_boundaries(spacing=None):
	"""``(outer, inners, guides, poles)``: polylines closed except the guides.

	``spacing`` is the length a CURVED input is divided by. A polyline input is
	still taken at its own vertices, so for the polyline boundaries this
	workflow has always used, the result is what
	``curve_to_compas_polyline`` gave, point for point.
	"""
	max_edge = 0.125 if spacing is None else spacing
	outer_ids = rs.ObjectsByLayer("Outer")
	inner_ids = rs.ObjectsByLayer("Inner") or []
	guide_ids = rs.ObjectsByLayer("Guides") or []
	pole_ids = rs.ObjectsByLayer("PointFeatures") or []

	if not outer_ids:
		raise RuntimeError("No outer boundary selected.")

	outer = curve_to_polyline(outer_ids[0], max_edge)

	inners = []
	for obj_id in inner_ids:
		inner_compas = curve_to_polyline(obj_id, max_edge)
		if is_polygon_in_polygon_xy(Polygon(outer), Polygon(inner_compas)):
			inners.append(inner_compas)
		else:
			rs.DeleteObjects(obj_id)
			print("This inner boundary did not lie inside the outer boundary. It has been removed from the selection.")

	for boundary in inners + [outer]:
		if not boundary.is_closed:
			raise RuntimeError("Not all boundaries are closed.")

	guides = []
	for guide_id in guide_ids:
		# A guide is an OPEN curve -- do not close it the way a loop is closed.
		guide_compas = Polyline(curve_points(guide_id, max_edge))
		if is_polygon_in_polygon_xy(Polygon(outer), Polygon(guide_compas)):
			guides.append(guide_compas)
		else:
			rs.DeleteObjects(guide_id)
			print("This guide did not lie inside the outer boundary. It has been removed from the selection.")

	poles = []
	for pole_id in pole_ids:
		point = cr.objects.find_object(pole_id)
		point_compas = point_to_compas(point.Geometry)
		if is_point_in_polygon_xy(point_compas, Polygon(outer)):
			if all(not is_point_in_polygon_xy(point_compas, Polygon(inner)) for inner in inners):
				poles.append(point_compas)
			else:
				rs.DeleteObjects(pole_id)
				print("This pole lies inside an inner boundary. It has been removed from the selection.")
		else:
			rs.DeleteObjects(pole_id)
			print("This pole did not lie inside the outer boundary. It has been removed from the selection.")

	return outer, inners, guides, poles


def read_coarse():
	guids = rs.ObjectsByLayer("Mesh")
	if not guids:
		raise RuntimeError(
			"No coarse layout on '{}' ".format(
				"Mesh"))
	print(cr.objects.find_object(guids[0]))
	mesh = mesh_to_compas(cr.objects.find_object(guids[0]).Geometry)
	vertices, faces = mesh.to_vertices_and_faces()
	faces = [[v for i, v in enumerate(f) if v != f[i - 1]] for f in faces]

	poles = [list(point_to_compas(cr.objects.find_object(guid).Geometry))
			for guid in rs.ObjectsByLayer("Poles") or []]

	coarse_mesh = CoarsePseudoQuadMesh.from_vertices_and_faces_with_poles(
		vertices, faces, poles)
	return coarse_mesh, poles

def read_mesh(layer):
	guids = rs.ObjectsByLayer(layer)
	if not guids:
		raise RuntimeError(
			"No mesh on the layer '{}' ".format(
				"Mesh"))

	mesh = mesh_to_compas(cr.objects.find_object(guids[0]).Geometry)
	return mesh