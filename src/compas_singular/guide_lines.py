"""Guide lines for a coarse quad mesh -- two self-contained implementations.

Both functions take a COARSE QUAD MESH and ONE OR MORE GUIDE CURVES and return
a dense quad mesh whose lines follow them (a course along each guide, interfaces
across it at close to 90 degrees).

    guide_line_mesh(coarse, guides)   GUIDE F -- snap existing polyedges onto the
                                      curves. No topology change.
    guide_band_mesh(coarse, guides)   GUIDE B -- inflate polyedges into bands with
                                      add_strip, snap each band's two new rails
                                      onto parallel offsets of its curve.

``guides`` is one curve or a list of curves. Each is any of: a
``compas.geometry.Polyline``, a list of ``Point``, or a list of ``[x, y, z]``.
Curves are used in XY only (z is dropped) -- the fold guard, the slot scoring
and the boundary intersection are all planar.

Each function is written with every helper nested inside it, so either can be
lifted out of this file on its own. The duplication between them is deliberate.

Which to use
------------
The first question is not curved-versus-oblique, it is WHERE ON THE QUADS the
guide should land.

* ON THE EDGES between quads -- ``guide_line_mesh``. The guide becomes a chain of
  mesh edges.
* DOWN THE MIDDLE of a row of quads -- ``guide_band_mesh``, which is what its
  default ``block_rows=1`` gives. The band densifies to a single row of faces
  with the guide as its centreline.

If the quads are BLOCKS and the guide is a cable threaded through them, the
second is the one you want: a cable on an edge is a cable in the joint between
two blocks, which is not where a cable goes. Measured on the 12x20 plate, the
cable's distance off the centre of its block, in units of half the block width
(0 = dead centre, 1 = on the block edge):

    block_rows=1 (default)   median 0.03, worst 0.07
    block_rows=None (old)    median 0.98, worst 1.13     -- i.e. on the edge

The narrow bands this implies also place far more readily than the wide ones the
old ``half_width=1.2`` default gave: two guides on that plate are REJECTED at
1.2 and seat cleanly at 0.35.

Otherwise the two are now close to interchangeable, which they were not before
``swing`` (see point 4 below); with it, line and band agree to a few tenths of a
degree on both oblique and curved guides. What is left:

* CROSSING guides: only the line supports them, so a crossing cable cannot be
  put mid-block. Bands cannot cross -- ``add_strip`` would tear an earlier
  band's rails off their curve -- and the line has nowhere to put the shared
  block. See below.
* Strongly oblique (past ~20-24 deg on a small plate): the band's outer rail
  runs off the boundary corner before the line does and the fold guard refuses
  it. The line still returns something, degraded. Shrinking ``half_width`` buys
  the range back, and with ``block_rows=1`` it is already small.
* The band changes topology, one added strip per guide, and the line changes
  none -- but NEITHER adds a singularity. Both match the unguided mesh exactly
  (0 on the rectangle, 1 on the pentagon), because a strip added along a
  wall-to-wall polyedge leaves every vertex four-valent.

Multiple guides
---------------
Every guide needs its own slot -- an interior wall-to-wall polyedge -- and a
clean decomposition has very few. The pentagon floor plate has exactly ONE.

Refinement is the slot supply. A topological quad split adds no singularities
(on an all-quad mesh every new vertex comes out valence 4) and multiplies slots
fast; measured on that plate:

    refine   faces   slots   rail lengths
      0         7      1     [3]
      1        28      7     [5, 5, 5, 5, 5, 7, 7]
      2       112     19     [9 x 13, 13]
      3       448     43     [17 x 14, ...]

so ``refine`` is chosen automatically to supply enough slots (and, for curved
guides, long enough rails). Note that TWO straight guides already force
``refine=1``: at refine 0 there is only one slot in the whole mesh.

But refinement trades one resource for the other, and this is the thing to
understand before placing guides. Snapping moves a slot while its neighbouring
rows stay put, so a guide can only travel about as far as the local row spacing
before the faces in between invert -- and refining PACKS THE ROWS CLOSER. More
slots, less reach each. Two parallel guides on this plate fail outright that
way: the slot at y = 0.62 has its neighbouring row at y = 1.05, so a guide at
y = 1.5 folds four faces.

Two things fix it, and together they are what make several guides workable:

* ``relax`` -- centroid-relax everything that is NOT on a guide once a
  placement would otherwise fold, so the other rows redistribute instead of
  being run over. Boundary vertices slide ALONG the outline (domain corners
  pinned), which matters as much as the interior: without it the slot's end
  vertex cannot get past its neighbour on the wall and the guide is stuck in a
  window a fraction of a metre wide. This is not the dense-mesh smoothing a
  guide exists to avoid -- rails stay exactly on their curves and only coarse
  patch corners move. It never runs when a placement already succeeds, so
  single-guide results are unchanged.
* backtracking -- candidate slots are tried in score order and one that still
  folds is undone and the next tried, because the best-ALIGNED slot is often
  not a REACHABLE one.

Measured on the pentagon plate, interfaces across the guides:

    guides                 line              band
    two parallel           88.1, 100%        88.9, 100%
    three parallel         88.5,  87%        90.0,  91%
    two crossing           88.0,  96%        refused (bands cannot cross)

which is the single-guide 88.1 essentially unchanged -- adding guides costs
nothing as long as each can find a slot it can reach.

What still fails is a guide with no near-parallel line anywhere near it. This
plate's decomposition has no vertical line within x < 1.2 of the centre, so a
guide at x = 0.5 is rejected however much you refine; at x = 2.0 it places
cleanly. Refinement subdivides the lines that exist, it does not invent lines
in new directions -- only ``add_strip`` does that.

The two functions assign guides to slots differently, because one of them
changes topology as it goes:

* ``guide_line_mesh`` snaps only, so all slots are known upfront and it does a
  GLOBAL greedy assignment -- repeatedly take the best remaining (guide, slot)
  pair. The order you pass the guides in does not matter.
* ``guide_band_mesh`` runs ``add_strip`` per guide, which invalidates polyedge
  keys, so it works SEQUENTIALLY in the order you pass them and re-picks from
  the updated mesh each time. Order matters; put the guide you care most about
  first.

Crossing guides (line only)
---------------------------
Two guides running across each other are assigned two polyedges that share a
vertex. That vertex has to sit on BOTH curves, so it is placed at the curves'
intersection rather than at either one's own parameter -- which is the correct
answer, and it keeps both guides satisfied because the intersection lies on
each. If the two curves do not actually intersect near the shared vertex you get
a clear error instead of a silent compromise.

``guide_band_mesh`` refuses crossings. ``add_strip`` deletes and duplicates the
vertices of the polyedge it inflates, so a band crossing an earlier band would
tear that band's rails off their curve. It excludes any slot touching an earlier
guide's rails and tells you to use ``guide_line_mesh``.

The constraints you have to respect
-----------------------------------
1. Only an INTERIOR WALL-TO-WALL polyedge can carry a guide: open, both ends on
   the boundary, not itself part of the boundary.

2. Snapping slides the slot's end vertices ALONG the boundary. Push a guide too
   far and the coarse outline folds back on itself, which densifies into a
   silently self-overlapping slab. Both functions guard on the signed area of
   every coarse face in XY and raise instead.

3. A guide can only bend where the coarse rail HAS vertices -- densification
   fills each patch with a discrete Coons patch, whose sides are straight
   between coarse vertices. A 3-vertex rail snapped to a curve gives a CHEVRON.
   Hence the automatic refinement.

4. SNAPPING ALONE DOES NOT TURN THE COURSES THAT CROSS THE GUIDE -- it bends the
   row the guide sits on and nothing else. Each rail vertex moves along the chord
   normal only (that is what minimises movement), so it keeps the x of the coarse
   vertex above it, the patch side joining them stays exactly as it was, and
   every dense course in that patch leaves the guide at 90 deg MINUS the guide's
   local tilt. Measured against a no-guide control, the crossing family is
   0.0 deg off its original direction at every tilt band.

   The ``swing`` parameter fixes it: the coarse courses leaving each rail are put
   on the guide's local normal, decaying over ``swing`` rows. In the most tilted
   band, rectangle+curved 73.3 -> 89.5 deg, pentagon+curved 74.2 -> 90.0,
   pentagon+oblique 81.3 -> 89.3. It costs quad quality -- turning the cross
   field has to skew the quads somewhere, worst-corner-angle 5th percentile
   79.4 -> 73.7 on the rectangle -- but without it the whole distortion budget is
   spent on the guide's own angle, the one place it is not wanted.

   Even with ``swing`` there is no exact guarantee: a Coons patch interpolates
   between OPPOSITE sides, so a course leaves at exactly 90 deg only where the
   patch across the guide has a parallel side.

   When measuring this yourself, BIN BY THE GUIDE'S LOCAL TILT and compare
   against a no-guide control. A plain median hides the failure, because for a
   bowed guide the tilt is worst at the ends -- exactly what ``measure()`` drops
   as boundary fan. Same mesh, varying only that margin: 0 -> 76.4 deg / 62%
   within 15; 1.5 -> 78.7 / 79%; 2.0 -> 80.6 / 88%.

5. The guide need not reach the boundary -- it is extended along its end
   tangents, or trimmed if it already pokes out. It must run roughly wall to
   wall and stay monotone along its own chord; a guide that doubles back cannot
   be carried by a single polyedge crossing the plate, and raises.

Run this file directly for a worked example on the pentagon floor plate:

    python guide_lines.py            # numbers only
    python guide_lines.py --view     # opens the compas_viewer scene
"""

from compas.geometry import add_vectors
from compas.geometry import distance_point_point
from compas.geometry import dot_vectors
from compas.geometry import intersection_segment_segment_xy
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas.itertools import pairwise

from compas_singular.datastructures import CoarseQuadMesh
from compas_singular.datastructures.mesh_quad.grammar.add_strip import add_strip


__all__ = ['guide_line_mesh', 'guide_band_mesh', 'guide_feature_mesh']


# =============================================================================
# GUIDE F -- snap existing polyedges onto the guides
# =============================================================================

def guide_line_mesh(coarse, guides, target_length=0.6, refine=None, boundary=None,
                    min_rail_vertices=9, max_refine=3, relax=50, swing=2,
                    edges_to_curves=None):
    """Densify ``coarse`` with one of its polyedges snapped onto each guide.

    Reuses polyedges the decomposition already contains, so there is no topology
    change at all -- the patch layout is untouched and only the guide regions
    move. This is also the only one of the two that supports guides which CROSS.

    Each guide ends up ON THE MESH EDGES. If the quads are blocks and the guide
    is a cable threaded through them, that puts the cable in the joint between
    two blocks; use ``guide_band_mesh`` instead, whose default ``block_rows=1``
    runs each guide down the middle of its own row of blocks. Measured, distance
    off the centre of the block in units of half its width: 1.00 here, 0.00
    there. The exception is crossing guides, which only this function can place
    at all -- a crossing cable cannot be put mid-block.

    Parameters
    ----------
    coarse : CoarseQuadMesh
        The coarse mesh to guide. Not modified -- it is copied first.
    guides : curve | list[curve]
        One guide curve or several. A curve is a ``Polyline``, a list of
        ``Point``, or a list of ``[x, y, z]``; used in XY. Each is extended
        along its end tangents to the mesh boundary, or trimmed back to it.
        Guides are assigned to slots by a global greedy match, so the order you
        pass them in does not matter.
    target_length : float, optional
        Target quad edge length for densification.
    refine : int, optional
        Topological quad splits to apply before snapping. These supply the slots
        the guides need and the rail vertices a curved guide bends at. ``None``
        (default) picks the smallest number that gives every guide a slot and
        every curved guide a rail of at least ``min_rail_vertices``, capped at
        ``max_refine``. Note that two guides already need ``refine >= 1`` on a
        typical clean decomposition.
    boundary : list[[x, y, z]], optional
        The domain outline the guides' ends must land on. Defaults to the coarse
        mesh's own boundary polygon, which is what you want unless the real
        domain boundary is finer than the coarse mesh (a curved outline chorded
        by coarse edges).
    min_rail_vertices, max_refine : int, optional
        Bounds for the automatic ``refine``.
    relax : int, optional
        Centroid-relaxation iterations applied to the coarse interior when a
        placement would otherwise fold, with the guides and the boundary held
        fixed. 0 disables it. This is what makes several guides placeable at
        all -- see the note on reach below. It never runs when a placement
        already succeeds, so single-guide results are untouched.
    swing : int, optional
        Coarse rows either side of each rail whose courses are turned onto the
        guide's local normal, decaying linearly with distance. This is what makes
        the quads perpendicular TO THE GUIDE rather than merely perpendicular to
        each other -- snapping on its own leaves the crossing courses exactly
        where they were. 2 (default) is the measured sweet spot; 3 reaches no
        further and costs more quad quality; 0 disables it and reproduces the
        pre-``swing`` behaviour. Skipped, with a note, if it would fold the
        outline.
    edges_to_curves : dict, optional
        ``{(u, v): [[x, y, z], ...]}`` keyed by ``coarse``'s OWN vertex indices
        (e.g. from a helper like ``boundary_edges_to_curves(coarse,
        decomposition)``), mapping a coarse edge to the real curve it was
        chorded from -- see ``CoarsePseudoQuadMesh.densification``. Without it,
        the final densification chords every coarse edge as a straight line, so
        a curved domain boundary comes out faceted even though the source mesh
        was smooth. Used only for edges of the FINAL coarse mesh whose two
        endpoints still sit exactly where ``coarse`` had them -- an edge inside
        a ``refine`` quad split has no original counterpart at all (the split
        replaces it with several shorter ones), and one a guide's
        snap/relax/swing moved has no known curve any more (the vertex slid
        along the real boundary, not along the stale one); both fall back to a
        straight chord. On a small coarse mesh ``swing``'s reach (default 2)
        alone can cover most of the boundary, so most of the benefit shows up
        with ``refine=0`` and ``swing=0``, or on a larger/more refined mesh
        where untouched boundary edges remain.

    Returns
    -------
    (QuadMesh, dict)
        The dense quad mesh, and a dict with ``coarse`` (the edited coarse
        mesh), ``guides`` (the clipped guide polylines, for drawing or
        measuring), ``rails`` (the snapped vertex chain per guide, in the order
        the guides were given), ``refine`` and ``note``.

    Raises
    ------
    ValueError
        Not enough interior wall-to-wall polyedges to carry every guide; a guide
        does not cross the domain; a guide doubles back on its own chord; two
        crossing guides do not actually intersect; or snapping folded the coarse
        outline.
    """

    # -- input -----------------------------------------------------------------

    def as_curve_list(obj):
        """One curve or several -> list of curves, without guessing wrongly.

        A ``Polyline`` has ``.points``; a bare list whose first item is a number
        pair/triple is a single curve; anything else is a list of curves.
        """
        if hasattr(obj, 'points'):
            return [obj]
        seq = list(obj)
        if not seq:
            raise ValueError('no guide curves given')
        if hasattr(seq[0], 'points'):
            return seq
        try:
            if isinstance(seq[0][0], (int, float)):
                return [seq]                            # a single curve of points
        except (TypeError, IndexError):
            pass
        return seq

    def as_points(curve):
        """Polyline | list[Point] | list[[x, y, z]]  ->  list of [x, y, 0.0]."""
        raw = getattr(curve, 'points', curve)
        pts = [[float(p[0]), float(p[1]), 0.0] for p in raw]
        clean = [pts[0]]
        for p in pts[1:]:
            if distance_point_point(p, clean[-1]) > 1e-9:
                clean.append(p)
        if len(clean) < 2:
            raise ValueError('a guide curve degenerates to a single point')
        return clean

    # -- geometry helpers ------------------------------------------------------

    def closest(p, poly):
        """(distance, unit tangent, foot) of the closest point on a polyline."""
        best = (float('inf'), None, None)
        for a, b in pairwise(poly):
            ab = subtract_vectors(b, a)
            length2 = dot_vectors(ab, ab)
            if length2 == 0.0:
                continue
            t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / length2))
            q = [a[i] + ab[i] * t for i in range(3)]
            d = distance_point_point(p, q)
            if d < best[0]:
                best = (d, normalize_vector(ab), q)
        return best

    def curve_intersection(poly_a, poly_b, near):
        """Where two guide polylines cross, taking the crossing closest to ``near``."""
        best, best_d = None, float('inf')
        for a, b in pairwise(poly_a):
            if distance_point_point(a, b) == 0.0:
                continue
            for p, q in pairwise(poly_b):
                if distance_point_point(p, q) == 0.0:
                    continue
                x = intersection_segment_segment_xy((a, b), (p, q))
                if x is None:
                    continue
                x = [x[0], x[1], 0.0]
                d = distance_point_point(x, near)
                if d < best_d:
                    best, best_d = x, d
        return best

    def boundary_polygon(mesh):
        """The mesh's outer boundary as an ordered closed list of points."""
        loops = mesh.boundaries()
        if not loops:
            raise ValueError('the coarse mesh has no boundary')
        loop = max(loops, key=len)                      # outer, if there are holes
        return [mesh.vertex_coordinates(v) for v in loop]

    def clip_to_boundary(pts, poly):
        """Extend the guide along its end tangents to the boundary, or trim it.

        The crossings taken are the ones BRACKETING the middle of the guide, not
        the outermost ones, so a ray that re-enters a non-convex domain further
        out does not stretch the guide across the gap.
        """
        far = 1.0e4
        head = add_vectors(pts[0], scale_vector(normalize_vector(
            subtract_vectors(pts[0], pts[1])), far))
        tail = add_vectors(pts[-1], scale_vector(normalize_vector(
            subtract_vectors(pts[-1], pts[-2])), far))
        ext = [head] + pts + [tail]

        hits = []
        for i, (a, b) in enumerate(pairwise(ext)):
            ab = subtract_vectors(b, a)
            length2 = dot_vectors(ab, ab)
            if length2 == 0.0:
                continue
            for p, q in pairwise(list(poly) + [poly[0]]):
                if distance_point_point(p, q) == 0.0:
                    continue                            # duplicated corner sample
                x = intersection_segment_segment_xy((a, b), (p, q))
                if x is None:
                    continue
                x = [x[0], x[1], 0.0]
                t = dot_vectors(subtract_vectors(x, a), ab) / length2
                hits.append((i + t, x))

        anchor = 1.0 + 0.5 * (len(pts) - 1)             # ext-index of the guide's middle
        before = [h for h in hits if h[0] <= anchor]
        after = [h for h in hits if h[0] >= anchor]
        if not before or not after:
            raise ValueError('a guide does not cross the domain -- check that its '
                             'middle lies inside the coarse mesh boundary')
        s0, first = max(before, key=lambda h: h[0])
        s1, last = min(after, key=lambda h: h[0])

        out = [first] + [ext[j] for j in range(len(ext)) if s0 < j < s1] + [last]
        clean = [out[0]]
        for p in out[1:]:
            if distance_point_point(p, clean[-1]) > 1e-9:
                clean.append(p)
        if len(clean) < 2:
            raise ValueError('a guide clips to a degenerate segment')
        return clean

    def chord_parametrise(pts):
        """(start, end, point_at) with t = the CHORD PROJECTION, not arc length.

        This is what makes snapping a minimal move: a vertex keeps its position
        along the span and only travels sideways onto the curve. It requires the
        guide to be monotone along its own chord, which is exactly the class of
        guides a single wall-to-wall polyedge can carry.
        """
        start, end = pts[0], pts[-1]
        chord = subtract_vectors(end, start)
        chord2 = dot_vectors(chord, chord)
        if chord2 == 0.0:
            raise ValueError('a guide starts and ends at the same point')
        us = [dot_vectors(subtract_vectors(p, start), chord) / chord2 for p in pts]
        for a, b in pairwise(us):
            if b <= a + 1e-12:
                raise ValueError('a guide doubles back on its own chord; a single '
                                 'polyedge cannot carry it')

        def point_at(t):
            t = min(1.0, max(0.0, t))
            for i in range(len(us) - 1):
                if us[i] <= t <= us[i + 1]:
                    span = us[i + 1] - us[i]
                    f = 0.0 if span == 0.0 else (t - us[i]) / span
                    return [pts[i][k] + (pts[i + 1][k] - pts[i][k]) * f for k in range(3)]
            return list(pts[-1])

        return start, end, point_at

    def is_curved(pts):
        """Whether the guide bows off its chord enough to need rail vertices."""
        start, end = pts[0], pts[-1]
        length = distance_point_point(start, end)
        if length == 0.0:
            return False
        d = scale_vector(subtract_vectors(end, start), 1.0 / length)
        normal = [-d[1], d[0], 0.0]
        bow = max(abs(dot_vectors(subtract_vectors(p, start), normal)) for p in pts)
        return bow > 1e-6 * length

    # -- coarse-mesh helpers ---------------------------------------------------

    def split_curve_at_half(curve):
        """Arc-length midpoint of a polyline, and the two halves either side of
        it (both include the split point, so rejoining reproduces the input)."""
        seglens = [distance_point_point(curve[i], curve[i + 1]) for i in range(len(curve) - 1)]
        total = sum(seglens)
        if total == 0.0:
            i = len(curve) // 2
            return curve[: i + 1], curve[i:]
        target = total / 2.0
        acc = 0.0
        for i, seglen in enumerate(seglens):
            if acc + seglen >= target or i == len(seglens) - 1:
                t = 0.0 if seglen == 0.0 else max(0.0, min(1.0, (target - acc) / seglen))
                a, b = curve[i], curve[i + 1]
                point = [a[k] + (b[k] - a[k]) * t for k in range(3)]
                return curve[: i + 1] + [point], [point] + curve[i + 1:]
            acc += seglen
        i = len(curve) // 2
        return curve[: i + 1], curve[i:]

    def quad_split(vertices, faces, edges_to_curves=None):
        """Catmull-Clark TOPOLOGICAL split: every n-gon becomes n quads.

        Original vertices keep their position -- no smoothing. On an all-quad
        mesh every new vertex (face centre and edge midpoint alike) comes out
        valence 4, so this buys slots and rail resolution without adding a
        singularity.

        edges_to_curves : dict, optional
            ``{(u, v): [[x, y, z], ...]}`` keyed by the SAME positional indices
            as ``vertices``/``faces``. When a split edge has a curve, its new
            midpoint is placed at the curve's own arc-length midpoint (not the
            plain average of the two endpoints) and the curve is chopped in two
            for the resulting sub-edges. Without this, refining a curved
            boundary before densification throws the curve away one split
            early -- the new vertices would sit on straight chords no matter
            what ``edges_to_curves`` densification is later given. Returns a
            third value: the same mapping, re-keyed to the output vertices,
            covering only edges that had a curve to begin with.
        """
        new_vertices = [list(p) for p in vertices]
        edges_to_curves = edges_to_curves or {}
        edge_mid = {}
        new_curves = {}

        def curve_for(u, v):
            if (u, v) in edges_to_curves:
                return edges_to_curves[u, v]
            if (v, u) in edges_to_curves:
                return list(reversed(edges_to_curves[v, u]))
            return None

        def mid(u, v):
            key = (min(u, v), max(u, v))
            if key not in edge_mid:
                idx = len(new_vertices)
                edge_mid[key] = idx
                curve = curve_for(u, v)
                if curve and len(curve) >= 2:
                    first, second = split_curve_at_half(curve)
                    new_vertices.append(list(first[-1]))
                    new_curves[u, idx] = first
                    new_curves[idx, v] = second
                else:
                    a, b = vertices[u], vertices[v]
                    new_vertices.append([(a[i] + b[i]) / 2.0 for i in range(3)])
            return edge_mid[key]

        new_faces = []
        for fv in faces:
            n = len(fv)
            centre = len(new_vertices)
            new_vertices.append([sum(vertices[v][i] for v in fv) / n for i in range(3)])
            mids = [mid(fv[i], fv[(i + 1) % n]) for i in range(n)]
            for i in range(n):
                new_faces.append([fv[i], mids[i], centre, mids[i - 1]])
        return new_vertices, new_faces, new_curves

    def prepared(times, edges_to_curves=None):
        """A quad-split copy of the input mesh, with strips collected, plus
        ``edges_to_curves`` carried through the splits and re-keyed to the
        output mesh's own vertex indices (see ``quad_split``).

        NOTE: rebuilding drops mesh attributes, including pseudo-quad pole data
        -- do not refine a pole-bearing mesh.
        """
        mesh = coarse.copy()
        curves = dict(edges_to_curves) if edges_to_curves else {}
        for _ in range(times):
            index = {vkey: i for i, vkey in enumerate(mesh.vertices())}
            vertices = [mesh.vertex_coordinates(v) for v in mesh.vertices()]
            faces = [[index[v] for v in mesh.face_vertices(f)] for f in mesh.faces()]
            curves = {(index[u], index[v]): c for (u, v), c in curves.items()
                      if u in index and v in index}
            vertices, faces, curves = quad_split(vertices, faces, curves)
            mesh = CoarseQuadMesh.from_vertices_and_faces(vertices, faces)
        mesh.collect_strips()
        return mesh, curves

    def all_slots(mesh):
        """Every polyedge that can carry a guide: open, interior, wall to wall.

        Interior AND wall-to-wall: both ends on the boundary but not the whole
        chain, or snapping would drag the outline of the slab around. Closed
        rings are excluded -- they cannot carry an open guide.
        """
        out = []
        for pkey, polyedge in mesh.collect_polyedges():
            polyedge = list(polyedge)
            if len(polyedge) < 3 or polyedge[0] == polyedge[-1]:
                continue
            if not (mesh.is_vertex_on_boundary(polyedge[0])
                    and mesh.is_vertex_on_boundary(polyedge[-1])):
                continue
            if all(mesh.is_vertex_on_boundary(v) for v in polyedge):
                continue                                        # the boundary itself
            out.append((pkey, polyedge))
        return out

    def score(mesh, polyedge, pts):
        """(alignment, score) of a slot against a guide.

        Alignment is |cos| between the two chords, 1 == parallel; it is returned
        separately so a slot pointing the wrong way can be rejected outright
        rather than merely ranked low. The score subtracts the slot's mean
        distance from the curve, normalised by guide length so it is scale-free.
        """
        chord = normalize_vector(subtract_vectors(pts[-1], pts[0]))
        length = distance_point_point(pts[0], pts[-1])
        cps = [mesh.vertex_coordinates(v) for v in polyedge]
        align = abs(dot_vectors(normalize_vector(subtract_vectors(cps[-1], cps[0])), chord))
        offset = sum(closest(p, pts)[0] for p in cps) / len(cps)
        return align, align - 1.5 * offset / length

    def snap_rail(mesh, rail, pts, start, end, point_at, pinned):
        """Move a chain of coarse vertices onto a guide.

        Each vertex keeps its projection parameter on the straight chord and
        only moves across onto the curve, so a straight guide moves everything
        the minimum distance and a curved one just adds the bow. The two ends
        are boundary vertices and must stay on the boundary, so they go to the
        guide's endpoints -- which are its boundary crossings.

        A vertex already claimed by an earlier guide is where two guides cross.
        It has to lie on both curves, so it goes to their intersection instead.
        That still satisfies the earlier guide -- the intersection is on its
        curve too, the vertex only slides along it.
        """
        chord = subtract_vectors(end, start)
        chord2 = dot_vectors(chord, chord)
        last = len(rail) - 1
        for i, vkey in enumerate(rail):
            p = mesh.vertex_coordinates(vkey)
            if vkey in pinned:
                q = curve_intersection(pinned[vkey], pts, p)
                if q is None:
                    raise ValueError(
                        'two guides are assigned polyedges that share a vertex but '
                        'the curves never cross -- move them apart, or pass them '
                        'one at a time')
            elif i == 0 or i == last:
                t = 0.0 if distance_point_point(p, start) <= distance_point_point(p, end) else 1.0
                q = point_at(t)
            else:
                t = min(1.0, max(0.0, dot_vectors(subtract_vectors(p, start), chord) / chord2))
                q = point_at(t)
            mesh.vertex[vkey]['x'], mesh.vertex[vkey]['y'], mesh.vertex[vkey]['z'] = q
            pinned[vkey] = pts

    def walk_column(mesh, source, first, steps):
        """The quad polyedge leaving ``source`` through ``first``, ``steps`` long.

        Continues straight across each four-valent vertex -- of the neighbours of
        the vertex reached, the continuation is the one sharing no face with the
        vertex arrived from. Stops early at the boundary or at a singularity,
        which is what should happen: neither has a column to turn.
        """
        out, previous, current = [first], source, first
        for _ in range(steps - 1):
            if mesh.is_vertex_on_boundary(current):
                break
            neighbours = mesh.vertex_neighbors(current)
            if len(neighbours) != 4:
                break
            behind = set(mesh.vertex_faces(previous))
            ahead = [v for v in neighbours
                     if v != previous and not (set(mesh.vertex_faces(v)) & behind)]
            if len(ahead) != 1:
                break
            previous, current = current, ahead[0]
            out.append(current)
        return out

    def swing_columns(mesh, chains, steps, poly, protected):
        """Turn the coarse courses leaving each rail onto the guide's normal.

        Snapping bends the row the guide sits on and NOTHING else. It moves each
        rail vertex along the chord normal only, so the vertex keeps the x of the
        coarse vertex above it, so the patch side joining them is untouched, so
        every dense course in that patch still runs the way it always did. The
        across-guide angle is then just 90 deg minus the guide's local tilt --
        measured against a no-guide control, the crossing family sits 0.0 deg off
        its original direction at every tilt band. The guide bends its own course
        and does not turn the ones that cross it.

        So turn them here: walk each coarse polyedge leaving a rail and put its
        k-th vertex on the guide's normal ray at the same arc distance from the
        rail, blended back to where it was with weight ``1 - k / steps``. A
        vertex reached from two rails takes the average, which is what two
        crossing guides should give it. Boundary vertices are projected back onto
        the outline -- a wall has to stay a wall.

        Measured, across-guide angle in the most tilted band: rectangle+curved
        73.3 -> 89.5 deg, pentagon+curved 74.2 -> 90.0, pentagon+oblique
        81.3 -> 89.3, no folds. The cost is real and unavoidable, because turning
        the cross field has to skew the quads somewhere: worst-corner-angle 5th
        percentile 79.4 -> 73.7 on the rectangle. Without it the entire budget is
        spent on the guide's own angle.
        """
        ring = list(poly) + [poly[0]]

        def project(p):
            """Back onto the domain outline, so a wall vertex stays on the wall."""
            best, best_d = p, float('inf')
            for a, b in pairwise(ring):
                ab = subtract_vectors(b, a)
                length2 = dot_vectors(ab, ab)
                if length2 == 0.0:
                    continue
                t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / length2))
                q = [a[i] + ab[i] * t for i in range(3)]
                d = distance_point_point(p, q)
                if d < best_d:
                    best, best_d = q, d
            return best

        proposed = {}
        for rail, pts in chains:
            on_rail = set(rail)
            for vkey in rail:
                anchor = mesh.vertex_coordinates(vkey)
                _, tangent, _ = closest(anchor, pts)
                normal = [-tangent[1], tangent[0], 0.0]
                for first in mesh.vertex_neighbors(vkey):
                    if first in on_rail or first in protected:
                        continue                # along the guide, or on another rail
                    side = 1.0 if dot_vectors(subtract_vectors(
                        mesh.vertex_coordinates(first), anchor), normal) > 0.0 else -1.0
                    arc, previous = 0.0, anchor
                    for k, v in enumerate(walk_column(mesh, vkey, first, steps)):
                        p = mesh.vertex_coordinates(v)
                        arc += distance_point_point(previous, p)
                        previous = p
                        if v in protected:
                            continue            # never move a vertex that is on a guide
                        weight = 1.0 - k / float(steps)
                        target = [anchor[i] + side * arc * normal[i] for i in range(3)]
                        moved = [p[i] * (1.0 - weight) + target[i] * weight for i in range(3)]
                        if mesh.is_vertex_on_boundary(v):
                            moved = project(moved)
                        proposed.setdefault(v, []).append(moved)

        for vkey, candidates in proposed.items():
            p = [sum(c[i] for c in candidates) / len(candidates) for i in range(3)]
            mesh.vertex[vkey]['x'], mesh.vertex[vkey]['y'], mesh.vertex[vkey]['z'] = p

    def folded_faces(mesh):
        """Faces with non-positive signed area in XY -- a folded coarse outline.
        Densifying one produces a self-overlapping slab, silently."""
        bad = []
        for fkey in mesh.faces():
            pts = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
            area2 = sum(p[0] * q[1] - q[0] * p[1] for p, q in zip(pts, pts[1:] + pts[:1]))
            if area2 <= 1e-9:
                bad.append(fkey)
        return bad

    def relaxed(mesh, fixed, poly, iterations):
        """Centroid-relax the coarse interior, guides and boundary held fixed.

        Snapping moves a slot but leaves its neighbouring rows where they were,
        so a guide travelling much more than the local row spacing folds the
        faces in between -- and refining makes that WORSE, because it packs the
        rows closer. Letting the rows that are not on a guide redistribute is
        what buys the reach back.

        This is not the dense-mesh smoothing a guide exists to avoid. The rails
        stay exactly on their curves; only coarse patch corners away from the
        guides move, and densification still builds its Coons patches between
        them. Measured on two parallel guides: 4 folded faces before, 0 after,
        and interfaces at 88.3 deg -- the same as a single guide.
        """
        ring = list(poly) + [poly[0]]

        def project(p):
            """Back onto the domain outline, so a sliding vertex stays on it."""
            best, best_d = p, float('inf')
            for a, b in pairwise(ring):
                ab = subtract_vectors(b, a)
                length2 = dot_vectors(ab, ab)
                if length2 == 0.0:
                    continue
                t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / length2))
                q = [a[i] + ab[i] * t for i in range(3)]
                d = distance_point_point(p, q)
                if d < best_d:
                    best, best_d = q, d
            return best

        # domain corners must not slide, or the outline gets rounded off
        corners = []
        for i in range(len(ring) - 1):
            u = normalize_vector(subtract_vectors(ring[i], ring[i - 2 if i == 0 else i - 1]))
            v = normalize_vector(subtract_vectors(ring[i + 1], ring[i]))
            if dot_vectors(u, v) < 0.985:                       # over ~10 deg of turn
                corners.append(ring[i])

        interior, sliding = [], set()
        for v in mesh.vertices():
            if v in fixed:
                continue
            if mesh.is_vertex_on_boundary(v):
                p = mesh.vertex_coordinates(v)
                if not any(distance_point_point(p, c) < 1e-6 for c in corners):
                    sliding.add(v)
            else:
                interior.append(v)

        free = interior + sorted(sliding)
        for _ in range(iterations):
            moved = {}
            for v in free:
                nbrs = mesh.vertex_neighbors(v)
                if not nbrs:
                    continue
                ps = [mesh.vertex_coordinates(u) for u in nbrs]
                c = [sum(p[i] for p in ps) / len(ps) for i in range(3)]
                moved[v] = project(c) if v in sliding else c
            for v, p in moved.items():
                mesh.vertex[v]['x'], mesh.vertex[v]['y'], mesh.vertex[v]['z'] = p

    def boundary_scrambled(mesh, poly):
        """True if snapping has slid a boundary vertex past its neighbour.

        The fold guard cannot catch this once the interior has been relaxed:
        relaxation pulls the inverted faces back to positive area while the
        OUTLINE itself stays crossed over, which densifies into a slab that
        overlaps along its own edge. The invariant that survives relaxation is
        that the boundary vertices keep their cyclic order along the domain
        outline, so check that directly.
        """
        ring = list(poly) + [poly[0]]
        cum = [0.0]
        for a, b in pairwise(ring):
            cum.append(cum[-1] + distance_point_point(a, b))

        def arc(p):
            best, best_d = 0.0, float('inf')
            for i, (a, b) in enumerate(pairwise(ring)):
                ab = subtract_vectors(b, a)
                length2 = dot_vectors(ab, ab)
                if length2 == 0.0:
                    continue
                t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / length2))
                q = [a[j] + ab[j] * t for j in range(3)]
                d = distance_point_point(p, q)
                if d < best_d:
                    best, best_d = cum[i] + t * distance_point_point(a, b), d
            return best

        try:
            loops = mesh.boundaries()
        except (KeyError, TypeError):
            return True                         # the loop no longer even walks
        if not loops:
            return True
        loop = max(loops, key=len)
        s = [arc(mesh.vertex_coordinates(v)) for v in loop]
        # one descent going round if the loop runs with the outline, all but one
        # if it runs against it; anything else means two vertices have swapped
        descents = sum(1 for a, b in pairwise(s + s[:1]) if b < a - 1e-9)
        return descents not in (1, len(s) - 1)

    # -- run -------------------------------------------------------------------

    curves = as_curve_list(guides)
    poly = boundary_polygon(coarse) if boundary is None else [
        [float(p[0]), float(p[1]), 0.0] for p in boundary]

    pts_list = [clip_to_boundary(as_points(c), poly) for c in curves]
    params = [chord_parametrise(p) for p in pts_list]
    curved = [is_curved(p) for p in pts_list]
    n = len(pts_list)

    def place(k):
        """Seat every guide on a slot at refinement level ``k``.

        Candidates are tried in score order and a placement that folds the
        outline is UNDONE and the next one tried. That backtracking is what
        makes several guides workable, because refinement trades one resource
        for the other: it buys more slots, but packs the rows closer together so
        each slot can travel less far before it runs into its neighbour. The
        best-aligned slot is often not a reachable one.

        Returns ((mesh, mesh_curves, chosen), None) or (None, reason).
        """
        mesh, mesh_curves = prepared(k, edges_to_curves)
        slots = all_slots(mesh)
        if len(slots) < n:
            return None, '{} guide(s) but only {} interior wall-to-wall polyedge(s)'.format(
                n, len(slots))

        ranked = {}
        for gi in range(n):
            cands = []
            for si, (pkey, polyedge) in enumerate(slots):
                if curved[gi] and len(polyedge) < min_rail_vertices:
                    continue                    # too few vertices: would snap to a chevron
                align, s = score(mesh, polyedge, pts_list[gi])
                if align < 0.5:
                    continue                    # over 60 deg off; not this guide's slot
                cands.append((s, si))
            if not cands:
                return None, 'guide {} has no slot aligned with it{}'.format(
                    gi, ' and long enough to bend' if curved[gi] else '')
            cands.sort(reverse=True)
            ranked[gi] = cands

        # most confident guide first, so a guide with only one plausible slot is
        # not left holding whatever a less fussy one did not take
        order = sorted(range(n), key=lambda gi: -ranked[gi][0][0])
        used, pinned, chosen = set(), {}, {}
        held = set()
        for gi in order:
            seated = False
            for _, si in ranked[gi]:
                if si in used:
                    continue
                before = {v: mesh.vertex_coordinates(v) for v in mesh.vertices()}
                keep = dict(pinned)
                try:
                    snap_rail(mesh, slots[si][1], pts_list[gi], *params[gi], pinned=pinned)
                    bad = folded_faces(mesh)
                    if bad and relax:
                        relaxed(mesh, held | set(slots[si][1]), poly, relax)
                        bad = folded_faces(mesh)
                    if not bad and boundary_scrambled(mesh, poly):
                        bad = ['boundary']      # relaxation can hide this one
                except ValueError:
                    bad = [None]                # crossing guides that never meet
                if not bad:
                    used.add(si)
                    chosen[gi] = slots[si]
                    held |= set(slots[si][1])
                    seated = True
                    break
                for v, p in before.items():
                    mesh.vertex[v]['x'], mesh.vertex[v]['y'], mesh.vertex[v]['z'] = p
                pinned.clear()
                pinned.update(keep)
            if not seated:
                return None, ('guide {} cannot be seated: every aligned slot either folds '
                              'the outline or is taken'.format(gi))
        return (mesh, mesh_curves, chosen), None

    reason = None
    if refine is None:
        result = None
        for k in range(max_refine + 1):
            result, reason = place(k)
            if result:
                refine = k
                break
    else:
        result, reason = place(refine)

    if not result:
        raise ValueError(
            'could not place the guides ({}). Move them towards positions the mesh '
            'already has lines at, raise max_refine, or lower min_rail_vertices -- note '
            'that refining supplies more slots but shortens how far each can '
            'travel.'.format(reason))

    work, work_curves, chosen = result
    rails = [chosen[gi][1] for gi in range(n)]

    swung = ''
    if swing:
        before = {v: work.vertex_coordinates(v) for v in work.vertices()}
        swing_columns(work, list(zip(rails, pts_list)), swing, poly,
                      {v for rail in rails for v in rail})
        if folded_faces(work) or boundary_scrambled(work, poly):
            for v, p in before.items():
                work.vertex[v]['x'], work.vertex[v]['y'], work.vertex[v]['z'] = p
            swung = '; swing SKIPPED (would fold the outline)'
        else:
            swung = '; swung x{}'.format(swing)

    def boundary_arc_between(pa, pb, poly):
        """Arc of the closed boundary polygon `poly` between the points on it
        closest to `pa` and `pb`, oriented pa -> pb. `poly` need not already be
        closed (both `boundary_polygon`'s open loop and a caller-supplied
        explicitly-closed curve are accepted).

        This is what lets an edge whose endpoints snap/relax/swing MOVED still
        densify along the true curve: both endpoints are still ON the
        boundary (that is what boundary reprojection guarantees), just not
        where they started, so re-deriving the arc from their CURRENT
        positions works where reusing the ORIGINAL curve (which no longer
        meets them) cannot.
        """
        pts = list(poly)
        if distance_point_point(pts[0], pts[-1]) > 1e-9:
            pts = pts + [pts[0]]
        n = len(pts) - 1                            # unique vertices
        if n < 2:
            return None
        doubled = pts[:-1] + pts[:-1]                # two laps, no wraparound needed
        segs = list(pairwise(doubled))
        cum = [0.0]
        for a, b in segs:
            cum.append(cum[-1] + distance_point_point(a, b))
        total = cum[n]                               # length of exactly one lap

        def project(p):
            best_d, best_s = float('inf'), None
            for (a, b), s0 in zip(segs[:n], cum[:n]):
                L = distance_point_point(a, b)
                if L == 0:
                    continue
                ab = subtract_vectors(b, a)
                t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / (L * L)))
                q = [a[k] + ab[k] * t for k in range(3)]
                d = distance_point_point(p, q)
                if d < best_d:
                    best_d, best_s = d, s0 + t * L
            return best_d, best_s

        _, sa = project(pa)
        _, sb = project(pb)
        if sa is None or sb is None:
            return None

        fwd = (sb - sa) % total
        bwd = (sa - sb) % total
        reverse = fwd > bwd
        s0, length = (sb, bwd) if reverse else (sa, fwd)
        s1 = s0 + length                             # <= s0 + total < 2 * total: safe, no wrap

        def point_at(s):
            for i, (c0, c1) in enumerate(zip(cum, cum[1:])):
                if c0 - 1e-9 <= s <= c1 + 1e-9:
                    a, b = segs[i]
                    L = c1 - c0
                    t = 0.0 if L == 0 else (s - c0) / L
                    return [a[k] + (b[k] - a[k]) * t for k in range(3)]
            return list(doubled[-1])

        arc = [point_at(s0)]
        for i, c in enumerate(cum):
            if s0 < c < s1:
                arc.append(list(doubled[i]))
        arc.append(point_at(s1))
        if reverse:
            arc.reverse()
        return arc

    def curves_for_densification(mesh, given):
        """For a boundary edge, prefer re-deriving the curve from `poly`
        between the endpoints' CURRENT positions (see `boundary_arc_between`)
        -- this is what makes the boundary densify along the true curve even
        where a guide moved it. Falls back to `given` (curves carried through
        refine by `prepared`, keyed by `mesh`'s OWN vertex indices where a
        split happened -- see `quad_split`) for interior edges, then a
        straight chord as the last resort. A stale or mismatched curve would
        no longer meet the edge's actual endpoints and the dense mesh would
        fail to weld there -- see CoarsePseudoQuadMesh.densification's
        edges_to_curves.
        """
        curves = {}
        for u, v in mesh.edges():
            a, b = mesh.vertex_coordinates(u), mesh.vertex_coordinates(v)
            curve = None
            if mesh.is_vertex_on_boundary(u) and mesh.is_vertex_on_boundary(v):
                curve = boundary_arc_between(a, b, poly)
                if curve is not None:
                    curve[0], curve[-1] = list(a), list(b)
            if curve is None:
                curve = (given or {}).get((u, v))
                if curve is None:
                    rev = (given or {}).get((v, u))
                    curve = list(reversed(rev)) if rev else None
            if (curve is None
                    or distance_point_point(curve[0], a) > 1e-6
                    or distance_point_point(curve[-1], b) > 1e-6):
                curve = [a, b]
            curves[u, v] = curve
        return curves

    work.set_strips_density_target(t=target_length)
    work.densification(edges_to_curves=curves_for_densification(work, work_curves))
    dense = work.get_quad_mesh()

    info = {
        'coarse': work,
        'guides': pts_list,
        'centrelines': pts_list,
        'rails': rails,
        'refine': refine,
        'note': 'snapped {} polyedge(s) {} ({} vertices); refined x{}{}'.format(
            n, [chosen[gi][0] for gi in range(n)], [len(r) for r in rails], refine, swung),
    }
    return dense, info


# =============================================================================
# GUIDE B -- inflate polyedges into bands and snap their rails
# =============================================================================

def guide_band_mesh(coarse, guides, target_length=0.6, half_width=None, refine=None,
                    boundary=None, min_rail_vertices=9, max_refine=3, relax=50,
                    swing=2, block_rows=1, edges_to_curves=None):
    """Densify ``coarse`` with an ``add_strip`` band snapped onto each guide.

    With the default ``block_rows=1`` the band densifies to a SINGLE row of
    quads, so each guide runs down the middle of its own row of faces rather
    than along an edge. That is what you want when the quads are blocks and the
    guide is a cable threaded through them: a cable on an edge is a cable in the
    joint between two blocks. ``block_rows=None`` restores the old behaviour, an
    edge on the guide with a proportionally-densified band either side.

    ``add_strip`` inflates one polyedge into a band of new faces, giving two
    fresh polyedges. Those are snapped onto true parallel offsets of the guide,
    so the band has constant width and its two rails stay parallel TO EACH OTHER
    however the guide is oriented or curved. That is why it beats
    ``guide_line_mesh`` off-axis: the patches inside the band keep a parallel
    partner side, which is what the discrete Coons patch interpolates between.

    Guides are processed SEQUENTIALLY in the order given, because each strip
    invalidates polyedge keys. Order matters -- put the guide you care most
    about first. Guides may not cross; use ``guide_line_mesh`` for that.

    Parameters
    ----------
    coarse : CoarseQuadMesh
        The coarse mesh to guide. Not modified -- it is copied first.
    guides : curve | list[curve]
        One guide curve or several. A curve is a ``Polyline``, a list of
        ``Point``, or a list of ``[x, y, z]``; used in XY. Each is the
        CENTRELINE of its band, extended along its end tangents to the mesh
        boundary or trimmed back to it.
    target_length : float, optional
        Target quad edge length for densification.
    half_width : float | list[float], optional
        Distance from each centreline to its rails; the band spans twice this,
        so with ``block_rows=1`` it is HALF THE BLOCK WIDTH along the guide. A
        list sets it per guide. ``None`` (default) uses ``target_length / 2``, so
        the guide's blocks match the surrounding quads. Note that the outer rail
        sits this much further out than a single snapped line would, so a wide
        band hits the fold limit sooner -- shrink it to buy range back,
        especially on a bow, where ``half_width`` and the sagitta spend the same
        budget. The narrow bands this default gives place far more readily than
        the old 1.2: two guides on a 12x20 plate are rejected at 1.2 and seat
        cleanly at 0.35.
    refine : int, optional
        Topological quad splits applied BEFORE any strip is added. These supply
        the slots the guides need and the rail vertices a curved guide bends at.
        ``None`` (default) picks the smallest number that gives every guide a
        slot and every curved guide a long enough rail, capped at ``max_refine``.
    boundary : list[[x, y, z]], optional
        The domain outline the rails' ends must land on. Defaults to the coarse
        mesh's own boundary polygon.
    min_rail_vertices, max_refine : int, optional
        Bounds for the automatic ``refine``.
    relax : int, optional
        Centroid-relaxation iterations applied to the coarse interior when a
        placement would otherwise fold, with the bands and the boundary held
        fixed. 0 disables it. Note that a band spends this reach twice over --
        both rails have to travel -- so ``half_width`` must stay below the local
        coarse row spacing, which refinement keeps shrinking.
    swing : int, optional
        Coarse rows OUTSIDE each band whose courses are turned onto the rail's
        local normal, decaying linearly with distance. Between the rails a band
        is already correct -- they are parallel offsets, so the courses joining
        them run along the normal by construction -- but beyond them, snapping
        leaves the crossing courses exactly where they were. 2 (default) is the
        measured sweet spot; 0 disables it and reproduces the pre-``swing``
        behaviour. Skipped, with a note, if it would fold the outline.
    block_rows : int | None, optional
        Dense quad rows across each band. 1 (default) puts the guide down the
        CENTRELINE of a single row of faces -- the cable-through-blocks reading.
        Use 2 to get an edge on the guide with one block either side, or
        ``None`` to let the band densify proportionally like every other strip.
    edges_to_curves : dict, optional
        ``{(u, v): [[x, y, z], ...]}`` keyed by ``coarse``'s OWN vertex indices
        (e.g. from a helper like ``boundary_edges_to_curves(coarse,
        decomposition)``), mapping a coarse edge to the real curve it was
        chorded from -- see ``CoarsePseudoQuadMesh.densification``. Without it,
        the final densification chords every coarse edge as a straight line, so
        a curved domain boundary comes out faceted even though the source mesh
        was smooth. Used only for edges of the FINAL coarse mesh whose two
        endpoints still sit exactly where ``coarse`` had them -- an edge inside
        a ``refine`` quad split has no original counterpart at all (the split
        replaces it with several shorter ones), a band's new rail vertices are
        brand new, and anything snap/relax/swing moved has no known curve any
        more; all of those fall back to a straight chord. On a small coarse
        mesh ``swing``'s reach (default 2) alone can cover most of the
        boundary, so most of the benefit shows up with ``refine=0`` and
        ``swing=0``, or on a larger/more refined mesh where untouched boundary
        edges remain.

    Returns
    -------
    (QuadMesh, dict)
        The dense quad mesh, and a dict with ``coarse`` (the edited coarse
        mesh), ``guides`` (the rail polylines plus the clipped centrelines, for
        drawing or measuring), ``rails`` (the (left, right) vertex chains per
        guide), ``refine`` and ``note``.

    Raises
    ------
    ValueError
        Not enough interior wall-to-wall polyedges to inflate; a guide or one of
        its offsets does not cross the domain; a guide doubles back on its own
        chord; two guides would cross; or snapping folded the coarse outline.
    """

    # -- input -----------------------------------------------------------------

    def as_curve_list(obj):
        """One curve or several -> list of curves, without guessing wrongly.

        A ``Polyline`` has ``.points``; a bare list whose first item is a number
        pair/triple is a single curve; anything else is a list of curves.
        """
        if hasattr(obj, 'points'):
            return [obj]
        seq = list(obj)
        if not seq:
            raise ValueError('no guide curves given')
        if hasattr(seq[0], 'points'):
            return seq
        try:
            if isinstance(seq[0][0], (int, float)):
                return [seq]                            # a single curve of points
        except (TypeError, IndexError):
            pass
        return seq

    def as_points(curve):
        """Polyline | list[Point] | list[[x, y, z]]  ->  list of [x, y, 0.0]."""
        raw = getattr(curve, 'points', curve)
        pts = [[float(p[0]), float(p[1]), 0.0] for p in raw]
        clean = [pts[0]]
        for p in pts[1:]:
            if distance_point_point(p, clean[-1]) > 1e-9:
                clean.append(p)
        if len(clean) < 2:
            raise ValueError('a guide curve degenerates to a single point')
        return clean

    # -- geometry helpers ------------------------------------------------------

    def closest(p, poly):
        """(distance, unit tangent, foot) of the closest point on a polyline."""
        best = (float('inf'), None, None)
        for a, b in pairwise(poly):
            ab = subtract_vectors(b, a)
            length2 = dot_vectors(ab, ab)
            if length2 == 0.0:
                continue
            t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / length2))
            q = [a[i] + ab[i] * t for i in range(3)]
            d = distance_point_point(p, q)
            if d < best[0]:
                best = (d, normalize_vector(ab), q)
        return best

    def boundary_polygon(mesh):
        """The mesh's outer boundary as an ordered closed list of points."""
        loops = mesh.boundaries()
        if not loops:
            raise ValueError('the coarse mesh has no boundary')
        loop = max(loops, key=len)                      # outer, if there are holes
        return [mesh.vertex_coordinates(v) for v in loop]

    def clip_to_boundary(pts, poly):
        """Extend the guide along its end tangents to the boundary, or trim it.

        The crossings taken are the ones BRACKETING the middle of the guide, not
        the outermost ones, so a ray that re-enters a non-convex domain further
        out does not stretch the guide across the gap.
        """
        far = 1.0e4
        head = add_vectors(pts[0], scale_vector(normalize_vector(
            subtract_vectors(pts[0], pts[1])), far))
        tail = add_vectors(pts[-1], scale_vector(normalize_vector(
            subtract_vectors(pts[-1], pts[-2])), far))
        ext = [head] + pts + [tail]

        hits = []
        for i, (a, b) in enumerate(pairwise(ext)):
            ab = subtract_vectors(b, a)
            length2 = dot_vectors(ab, ab)
            if length2 == 0.0:
                continue
            for p, q in pairwise(list(poly) + [poly[0]]):
                if distance_point_point(p, q) == 0.0:
                    continue                            # duplicated corner sample
                x = intersection_segment_segment_xy((a, b), (p, q))
                if x is None:
                    continue
                x = [x[0], x[1], 0.0]
                t = dot_vectors(subtract_vectors(x, a), ab) / length2
                hits.append((i + t, x))

        anchor = 1.0 + 0.5 * (len(pts) - 1)             # ext-index of the guide's middle
        before = [h for h in hits if h[0] <= anchor]
        after = [h for h in hits if h[0] >= anchor]
        if not before or not after:
            raise ValueError('a guide does not cross the domain -- check that its '
                             'middle lies inside the coarse mesh boundary')
        s0, first = max(before, key=lambda h: h[0])
        s1, last = min(after, key=lambda h: h[0])

        out = [first] + [ext[j] for j in range(len(ext)) if s0 < j < s1] + [last]
        clean = [out[0]]
        for p in out[1:]:
            if distance_point_point(p, clean[-1]) > 1e-9:
                clean.append(p)
        if len(clean) < 2:
            raise ValueError('a guide clips to a degenerate segment')
        return clean

    def offset_curve(pts, distance):
        """A true parallel offset of the guide, in XY.

        Each sample moves along the normal of its LOCAL tangent (the average of
        the incoming and outgoing segment directions), so the offset keeps a
        constant distance from the curve. This is what gives the band constant
        width -- deriving each rail independently from a base point and a
        direction, as a base+direction+bulge formulation has to, leaves the two
        rails with different chord lengths and a width that breathes.

        Tight inner curves with a large offset can self-intersect; the fold
        guard catches the consequence.
        """
        out = []
        n = len(pts)
        for i, p in enumerate(pts):
            if i == 0:
                tangent = subtract_vectors(pts[1], pts[0])
            elif i == n - 1:
                tangent = subtract_vectors(pts[-1], pts[-2])
            else:
                tangent = add_vectors(
                    normalize_vector(subtract_vectors(p, pts[i - 1])),
                    normalize_vector(subtract_vectors(pts[i + 1], p)))
            tangent = normalize_vector(tangent)
            normal = [-tangent[1], tangent[0], 0.0]
            out.append([p[0] + normal[0] * distance, p[1] + normal[1] * distance, 0.0])
        return out

    def chord_parametrise(pts):
        """(start, end, point_at) with t = the CHORD PROJECTION, not arc length.

        This is what makes snapping a minimal move: a vertex keeps its position
        along the span and only travels sideways onto the curve. It requires the
        guide to be monotone along its own chord, which is exactly the class of
        guides a single wall-to-wall polyedge can carry.
        """
        start, end = pts[0], pts[-1]
        chord = subtract_vectors(end, start)
        chord2 = dot_vectors(chord, chord)
        if chord2 == 0.0:
            raise ValueError('a guide starts and ends at the same point')
        us = [dot_vectors(subtract_vectors(p, start), chord) / chord2 for p in pts]
        for a, b in pairwise(us):
            if b <= a + 1e-12:
                raise ValueError('a guide doubles back on its own chord; a single '
                                 'polyedge cannot carry it')

        def point_at(t):
            t = min(1.0, max(0.0, t))
            for i in range(len(us) - 1):
                if us[i] <= t <= us[i + 1]:
                    span = us[i + 1] - us[i]
                    f = 0.0 if span == 0.0 else (t - us[i]) / span
                    return [pts[i][k] + (pts[i + 1][k] - pts[i][k]) * f for k in range(3)]
            return list(pts[-1])

        return start, end, point_at

    def is_curved(pts):
        """Whether the guide bows off its chord enough to need rail vertices."""
        start, end = pts[0], pts[-1]
        length = distance_point_point(start, end)
        if length == 0.0:
            return False
        d = scale_vector(subtract_vectors(end, start), 1.0 / length)
        normal = [-d[1], d[0], 0.0]
        bow = max(abs(dot_vectors(subtract_vectors(p, start), normal)) for p in pts)
        return bow > 1e-6 * length

    # -- coarse-mesh helpers ---------------------------------------------------

    def split_curve_at_half(curve):
        """Arc-length midpoint of a polyline, and the two halves either side of
        it (both include the split point, so rejoining reproduces the input)."""
        seglens = [distance_point_point(curve[i], curve[i + 1]) for i in range(len(curve) - 1)]
        total = sum(seglens)
        if total == 0.0:
            i = len(curve) // 2
            return curve[: i + 1], curve[i:]
        target = total / 2.0
        acc = 0.0
        for i, seglen in enumerate(seglens):
            if acc + seglen >= target or i == len(seglens) - 1:
                t = 0.0 if seglen == 0.0 else max(0.0, min(1.0, (target - acc) / seglen))
                a, b = curve[i], curve[i + 1]
                point = [a[k] + (b[k] - a[k]) * t for k in range(3)]
                return curve[: i + 1] + [point], [point] + curve[i + 1:]
            acc += seglen
        i = len(curve) // 2
        return curve[: i + 1], curve[i:]

    def quad_split(vertices, faces, edges_to_curves=None):
        """Catmull-Clark TOPOLOGICAL split: every n-gon becomes n quads.

        Original vertices keep their position -- no smoothing. On an all-quad
        mesh every new vertex (face centre and edge midpoint alike) comes out
        valence 4, so this buys slots and rail resolution without adding a
        singularity.

        edges_to_curves : dict, optional
            ``{(u, v): [[x, y, z], ...]}`` keyed by the SAME positional indices
            as ``vertices``/``faces``. When a split edge has a curve, its new
            midpoint is placed at the curve's own arc-length midpoint (not the
            plain average of the two endpoints) and the curve is chopped in two
            for the resulting sub-edges. Without this, refining a curved
            boundary before densification throws the curve away one split
            early -- the new vertices would sit on straight chords no matter
            what ``edges_to_curves`` densification is later given. Returns a
            third value: the same mapping, re-keyed to the output vertices,
            covering only edges that had a curve to begin with.
        """
        new_vertices = [list(p) for p in vertices]
        edges_to_curves = edges_to_curves or {}
        edge_mid = {}
        new_curves = {}

        def curve_for(u, v):
            if (u, v) in edges_to_curves:
                return edges_to_curves[u, v]
            if (v, u) in edges_to_curves:
                return list(reversed(edges_to_curves[v, u]))
            return None

        def mid(u, v):
            key = (min(u, v), max(u, v))
            if key not in edge_mid:
                idx = len(new_vertices)
                edge_mid[key] = idx
                curve = curve_for(u, v)
                if curve and len(curve) >= 2:
                    first, second = split_curve_at_half(curve)
                    new_vertices.append(list(first[-1]))
                    new_curves[u, idx] = first
                    new_curves[idx, v] = second
                else:
                    a, b = vertices[u], vertices[v]
                    new_vertices.append([(a[i] + b[i]) / 2.0 for i in range(3)])
            return edge_mid[key]

        new_faces = []
        for fv in faces:
            n = len(fv)
            centre = len(new_vertices)
            new_vertices.append([sum(vertices[v][i] for v in fv) / n for i in range(3)])
            mids = [mid(fv[i], fv[(i + 1) % n]) for i in range(n)]
            for i in range(n):
                new_faces.append([fv[i], mids[i], centre, mids[i - 1]])
        return new_vertices, new_faces, new_curves

    def prepared(times, edges_to_curves=None):
        """A quad-split copy of the input mesh, with strips collected, plus
        ``edges_to_curves`` carried through the splits and re-keyed to the
        output mesh's own vertex indices (see ``quad_split``).

        NOTE: rebuilding drops mesh attributes, including pseudo-quad pole data
        -- do not refine a pole-bearing mesh.
        """
        mesh = coarse.copy()
        curves = dict(edges_to_curves) if edges_to_curves else {}
        for _ in range(times):
            index = {vkey: i for i, vkey in enumerate(mesh.vertices())}
            vertices = [mesh.vertex_coordinates(v) for v in mesh.vertices()]
            faces = [[index[v] for v in mesh.face_vertices(f)] for f in mesh.faces()]
            curves = {(index[u], index[v]): c for (u, v), c in curves.items()
                      if u in index and v in index}
            vertices, faces, curves = quad_split(vertices, faces, curves)
            mesh = CoarseQuadMesh.from_vertices_and_faces(vertices, faces)
        mesh.collect_strips()
        return mesh, curves

    def all_slots(mesh, blocked):
        """Every polyedge that can carry a guide: open, interior, wall to wall,
        and not touching a rail an earlier guide already claimed.

        Interior AND wall-to-wall: both ends on the boundary but not the whole
        chain, or inflating it would drag the outline of the slab around. Open
        only -- a closed ring cannot carry an open guide. These are also exactly
        add_strip's preconditions.

        ``blocked`` is why bands cannot cross: add_strip deletes and duplicates
        the vertices of the polyedge it inflates, so a band crossing an earlier
        band would tear that band's rails off their curve.
        """
        out = []
        for pkey, polyedge in mesh.collect_polyedges():
            polyedge = list(polyedge)
            if len(polyedge) < 3 or polyedge[0] == polyedge[-1]:
                continue
            if not (mesh.is_vertex_on_boundary(polyedge[0])
                    and mesh.is_vertex_on_boundary(polyedge[-1])):
                continue
            if all(mesh.is_vertex_on_boundary(v) for v in polyedge):
                continue                                        # the boundary itself
            if any(v in blocked for v in polyedge):
                continue                                        # would cross a band
            out.append((pkey, polyedge))
        return out

    def score(mesh, polyedge, pts):
        """(alignment, score) of a slot against a guide.

        Alignment is |cos| between the two chords, 1 == parallel; it is returned
        separately so a slot pointing the wrong way can be rejected outright
        rather than merely ranked low. The score subtracts the slot's mean
        distance from the curve, normalised by guide length so it is scale-free.
        """
        chord = normalize_vector(subtract_vectors(pts[-1], pts[0]))
        length = distance_point_point(pts[0], pts[-1])
        cps = [mesh.vertex_coordinates(v) for v in polyedge]
        align = abs(dot_vectors(normalize_vector(subtract_vectors(cps[-1], cps[0])), chord))
        offset = sum(closest(p, pts)[0] for p in cps) / len(cps)
        return align, align - 1.5 * offset / length

    def rail_side(mesh, chain, pts, exclude):
        """Which side of the guide a rail's own faces sit on.

        add_strip's (left, right) convention is topological -- it says nothing
        about which one is on which side of YOUR guide. Look at the neighbours
        that are not part of the new strip and sum their signed offset from the
        curve.
        """
        total = 0.0
        for vkey in chain:
            for nbr in mesh.vertex_neighbors(vkey):
                if nbr in exclude:
                    continue
                q = mesh.vertex_coordinates(nbr)
                _, tangent, foot = closest(q, pts)
                if tangent is None:
                    continue
                normal = [-tangent[1], tangent[0], 0.0]
                total += dot_vectors(subtract_vectors(q, foot), normal)
        return total

    def snap_rail(mesh, rail, start, end, point_at):
        """Move a chain of coarse vertices onto a rail curve.

        Each vertex keeps its projection parameter on the straight chord and
        only moves across onto the curve, so a straight guide moves everything
        the minimum distance and a curved one just adds the bow. The two ends
        are boundary vertices and must stay on the boundary, so they go to the
        rail's endpoints -- which are its boundary crossings.
        """
        chord = subtract_vectors(end, start)
        chord2 = dot_vectors(chord, chord)
        last = len(rail) - 1
        for i, vkey in enumerate(rail):
            p = mesh.vertex_coordinates(vkey)
            if i == 0 or i == last:
                t = 0.0 if distance_point_point(p, start) <= distance_point_point(p, end) else 1.0
            else:
                t = min(1.0, max(0.0, dot_vectors(subtract_vectors(p, start), chord) / chord2))
            q = point_at(t)
            mesh.vertex[vkey]['x'], mesh.vertex[vkey]['y'], mesh.vertex[vkey]['z'] = q

    def walk_column(mesh, source, first, steps):
        """The quad polyedge leaving ``source`` through ``first``, ``steps`` long.

        Continues straight across each four-valent vertex -- of the neighbours of
        the vertex reached, the continuation is the one sharing no face with the
        vertex arrived from. Stops early at the boundary or at a singularity,
        which is what should happen: neither has a column to turn.
        """
        out, previous, current = [first], source, first
        for _ in range(steps - 1):
            if mesh.is_vertex_on_boundary(current):
                break
            neighbours = mesh.vertex_neighbors(current)
            if len(neighbours) != 4:
                break
            behind = set(mesh.vertex_faces(previous))
            ahead = [v for v in neighbours
                     if v != previous and not (set(mesh.vertex_faces(v)) & behind)]
            if len(ahead) != 1:
                break
            previous, current = current, ahead[0]
            out.append(current)
        return out

    def swing_columns(mesh, chains, steps, poly, protected):
        """Turn the coarse courses leaving each rail onto that rail's normal.

        This is the OUTSIDE of the band. Between the rails, ``place`` already
        re-seats the second rail on its partner's normal; without that the band
        is no better than a single line there either.

        Outside, a band has the same defect as the single line: snapping moves
        each rail vertex along the chord normal only, so it keeps the x of the
        coarse vertex beyond it and the patch side joining them is untouched, and
        the fan of quads on the far side of each rail still leaves the guide at
        90 deg minus the local tilt. Measured, rectangle + bowed guide, most
        tilted band: outward 73.6 -> 89.3 deg.

        Only the outward columns are turned here; the inward ones start on the
        opposite rail, which is in ``protected``, and are skipped. Boundary
        vertices are projected back onto the outline -- a wall has to stay a wall.
        """
        ring = list(poly) + [poly[0]]

        def project(p):
            """Back onto the domain outline, so a wall vertex stays on the wall."""
            best, best_d = p, float('inf')
            for a, b in pairwise(ring):
                ab = subtract_vectors(b, a)
                length2 = dot_vectors(ab, ab)
                if length2 == 0.0:
                    continue
                t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / length2))
                q = [a[i] + ab[i] * t for i in range(3)]
                d = distance_point_point(p, q)
                if d < best_d:
                    best, best_d = q, d
            return best

        proposed = {}
        for rail, pts in chains:
            on_rail = set(rail)
            for vkey in rail:
                anchor = mesh.vertex_coordinates(vkey)
                _, tangent, _ = closest(anchor, pts)
                if tangent is None:
                    continue
                normal = [-tangent[1], tangent[0], 0.0]
                for first in mesh.vertex_neighbors(vkey):
                    if first in on_rail or first in protected:
                        continue                # along the rail, or across the band
                    side = 1.0 if dot_vectors(subtract_vectors(
                        mesh.vertex_coordinates(first), anchor), normal) > 0.0 else -1.0
                    arc, previous = 0.0, anchor
                    for k, v in enumerate(walk_column(mesh, vkey, first, steps)):
                        p = mesh.vertex_coordinates(v)
                        arc += distance_point_point(previous, p)
                        previous = p
                        if v in protected:
                            continue            # never move a vertex that is on a rail
                        weight = 1.0 - k / float(steps)
                        target = [anchor[i] + side * arc * normal[i] for i in range(3)]
                        moved = [p[i] * (1.0 - weight) + target[i] * weight for i in range(3)]
                        if mesh.is_vertex_on_boundary(v):
                            moved = project(moved)
                        proposed.setdefault(v, []).append(moved)

        for vkey, candidates in proposed.items():
            p = [sum(c[i] for c in candidates) / len(candidates) for i in range(3)]
            mesh.vertex[vkey]['x'], mesh.vertex[vkey]['y'], mesh.vertex[vkey]['z'] = p

    def folded_faces(mesh):
        """Faces with non-positive signed area in XY -- a folded coarse outline.
        Densifying one produces a self-overlapping slab, silently."""
        bad = []
        for fkey in mesh.faces():
            pts = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
            area2 = sum(p[0] * q[1] - q[0] * p[1] for p, q in zip(pts, pts[1:] + pts[:1]))
            if area2 <= 1e-9:
                bad.append(fkey)
        return bad

    def relaxed(mesh, fixed, poly, iterations):
        """Centroid-relax the coarse interior, bands and boundary held fixed.

        Snapping moves a rail but leaves its neighbouring rows where they were,
        so a rail travelling much more than the local row spacing folds the
        faces in between -- and refining makes that WORSE, because it packs the
        rows closer. Letting the rows that are not on a band redistribute is
        what buys the reach back.

        This is not the dense-mesh smoothing a guide exists to avoid. The rails
        stay exactly on their curves; only coarse patch corners away from the
        bands move, and densification still builds its Coons patches between
        them.
        """
        ring = list(poly) + [poly[0]]

        def project(p):
            """Back onto the domain outline, so a sliding vertex stays on it."""
            best, best_d = p, float('inf')
            for a, b in pairwise(ring):
                ab = subtract_vectors(b, a)
                length2 = dot_vectors(ab, ab)
                if length2 == 0.0:
                    continue
                t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / length2))
                q = [a[i] + ab[i] * t for i in range(3)]
                d = distance_point_point(p, q)
                if d < best_d:
                    best, best_d = q, d
            return best

        # domain corners must not slide, or the outline gets rounded off
        corners = []
        for i in range(len(ring) - 1):
            u = normalize_vector(subtract_vectors(ring[i], ring[i - 2 if i == 0 else i - 1]))
            v = normalize_vector(subtract_vectors(ring[i + 1], ring[i]))
            if dot_vectors(u, v) < 0.985:                       # over ~10 deg of turn
                corners.append(ring[i])

        interior, sliding = [], set()
        for v in mesh.vertices():
            if v in fixed:
                continue
            if mesh.is_vertex_on_boundary(v):
                p = mesh.vertex_coordinates(v)
                if not any(distance_point_point(p, c) < 1e-6 for c in corners):
                    sliding.add(v)
            else:
                interior.append(v)

        free = interior + sorted(sliding)
        for _ in range(iterations):
            moved = {}
            for v in free:
                nbrs = mesh.vertex_neighbors(v)
                if not nbrs:
                    continue
                ps = [mesh.vertex_coordinates(u) for u in nbrs]
                c = [sum(p[i] for p in ps) / len(ps) for i in range(3)]
                moved[v] = project(c) if v in sliding else c
            for v, p in moved.items():
                mesh.vertex[v]['x'], mesh.vertex[v]['y'], mesh.vertex[v]['z'] = p

    def boundary_scrambled(mesh, poly):
        """True if snapping has slid a boundary vertex past its neighbour.

        The fold guard cannot catch this once the interior has been relaxed:
        relaxation pulls the inverted faces back to positive area while the
        OUTLINE itself stays crossed over, which densifies into a slab that
        overlaps along its own edge. The invariant that survives relaxation is
        that the boundary vertices keep their cyclic order along the domain
        outline, so check that directly.
        """
        ring = list(poly) + [poly[0]]
        cum = [0.0]
        for a, b in pairwise(ring):
            cum.append(cum[-1] + distance_point_point(a, b))

        def arc(p):
            best, best_d = 0.0, float('inf')
            for i, (a, b) in enumerate(pairwise(ring)):
                ab = subtract_vectors(b, a)
                length2 = dot_vectors(ab, ab)
                if length2 == 0.0:
                    continue
                t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / length2))
                q = [a[j] + ab[j] * t for j in range(3)]
                d = distance_point_point(p, q)
                if d < best_d:
                    best, best_d = cum[i] + t * distance_point_point(a, b), d
            return best

        try:
            loops = mesh.boundaries()
        except (KeyError, TypeError):
            return True                         # the loop no longer even walks
        if not loops:
            return True
        loop = max(loops, key=len)
        s = [arc(mesh.vertex_coordinates(v)) for v in loop]
        # one descent going round if the loop runs with the outline, all but one
        # if it runs against it; anything else means two vertices have swapped
        descents = sum(1 for a, b in pairwise(s + s[:1]) if b < a - 1e-9)
        return descents not in (1, len(s) - 1)

    # -- run -------------------------------------------------------------------

    curves = as_curve_list(guides)
    poly = boundary_polygon(coarse) if boundary is None else [
        [float(p[0]), float(p[1]), 0.0] for p in boundary]

    pts_list = [clip_to_boundary(as_points(c), poly) for c in curves]
    for p in pts_list:
        chord_parametrise(p)                    # validate every centreline up front
    curved = [is_curved(p) for p in pts_list]
    n = len(pts_list)

    if half_width is None:
        # a block as wide as the quads around it, so the guide's row does not
        # read as a different piece of mesh
        half_width = target_length / 2.0
    widths = [float(half_width)] * n if isinstance(half_width, (int, float)) else \
        [float(w) for w in half_width]
    if len(widths) != n:
        raise ValueError('half_width has {} values for {} guides'.format(len(widths), n))

    def place(k):
        """Seat a band on every guide at refinement level ``k``.

        Each guide is tried against its aligned slots in score order, on a COPY
        of the mesh, and the first band that does not fold the outline is kept.
        A copy is needed rather than an undo because add_strip changes topology.
        Backtracking matters here for the same reason as for the line: refining
        buys more slots but packs the rows closer, so each one travels less far.

        Returns ((mesh, mesh_curves, rails, drawn, notes), None) or (None, reason).
        """
        mesh, mesh_curves = prepared(k, edges_to_curves)
        blocked = set()
        rails, drawn, notes = [], [], []

        for gi in range(n):
            pts = pts_list[gi]
            cands = []
            for pkey, polyedge in all_slots(mesh, blocked):
                if curved[gi] and len(polyedge) < min_rail_vertices:
                    continue                    # too few vertices: would snap to a chevron
                align, s = score(mesh, polyedge, pts)
                if align < 0.5:
                    continue                    # over 60 deg off; not this guide's slot
                cands.append((s, pkey, polyedge))
            if not cands:
                return None, 'guide {} has no free slot aligned with it{}'.format(
                    gi, ' and long enough to bend' if curved[gi] else '')
            cands.sort(key=lambda c: -c[0])

            seated = None
            for _, pkey, slot in cands:
                trial = mesh.copy()
                trial.collect_strips()
                try:
                    # add_strip pops from the list it is given, so hand it a copy
                    skey, old_to_new = add_strip(trial, list(slot))
                    left = [old_to_new[v][0] for v in slot]
                    right = [old_to_new[v][1] for v in slot]

                    # the new vertices are created ON TOP of the vertex they
                    # replace, so the strip has zero width until it is
                    # positioned -- densifying it in that state raises
                    # ZeroDivisionError in Polyline.point_at
                    fresh = set(left) | set(right)
                    sign = 1.0 if rail_side(trial, left, pts, fresh) >= 0.0 else -1.0

                    made = []
                    for chain, distance in ((left, sign * widths[gi]),
                                            (right, -sign * widths[gi])):
                        rail_pts = clip_to_boundary(offset_curve(pts, distance), poly)
                        r_start, r_end, r_point_at = chord_parametrise(rail_pts)
                        snap_rail(trial, chain, r_start, r_end, r_point_at)
                        made.append(rail_pts)

                    # The rails are SHAPED as parallel offsets, but each is then
                    # clipped and parametrised by projection onto ITS OWN chord,
                    # and for a bowed guide those two chords do not agree. So
                    # partner vertices end up at different stations and the
                    # courses crossing the band are no more normal to the guide
                    # than anything else -- measured, 73.6 deg in the most tilted
                    # band, exactly what an unguided mesh gives. Re-seat the
                    # second rail on its partner's normal, which is the whole
                    # point of using a band. The two end vertices stay put: they
                    # are boundary vertices and have to remain on the boundary.
                    for i in range(1, len(left) - 1):
                        _, tangent, foot = closest(trial.vertex_coordinates(left[i]), pts)
                        if tangent is None:
                            continue
                        normal = [-tangent[1], tangent[0], 0.0]
                        q = [foot[j] - normal[j] * sign * widths[gi] for j in range(3)]
                        trial.vertex[right[i]]['x'] = q[0]
                        trial.vertex[right[i]]['y'] = q[1]
                        trial.vertex[right[i]]['z'] = q[2]
                except (ValueError, ZeroDivisionError, KeyError):
                    continue                    # this slot cannot carry the band
                bad = folded_faces(trial)
                if bad and relax:
                    relaxed(trial, blocked | fresh, poly, relax)
                    bad = folded_faces(trial)
                if bad or boundary_scrambled(trial, poly):
                    continue                    # relaxation can hide the boundary case
                trial.collect_strips()          # refresh strip data after the topology edit
                seated = (trial, left, right, made,
                          'strip {} on polyedge {} ({} v/rail)'.format(skey, pkey, len(left)))
                break

            if seated is None:
                return None, ('guide {} cannot be seated: every aligned slot either folds '
                              'the outline or crosses a band already placed (bands cannot '
                              'cross -- use guide_line_mesh)'.format(gi))

            mesh, left, right, made, note = seated
            blocked |= set(left) | set(right)
            rails.append((left, right))
            drawn.extend(made)
            notes.append(note)

        return (mesh, mesh_curves, rails, drawn, notes), None

    reason = None
    if refine is None:
        result = None
        for k in range(max_refine + 1):
            result, reason = place(k)
            if result:
                refine = k
                break
    else:
        result, reason = place(refine)

    if not result:
        raise ValueError(
            'could not place the guides ({}). Move them towards positions the mesh '
            'already has lines at, reduce half_width, raise max_refine, or lower '
            'min_rail_vertices -- note that refining supplies more slots but shortens '
            'how far each can travel.'.format(reason))

    work, work_curves, rails, drawn, notes = result

    swung = ''
    if swing:
        before = {v: work.vertex_coordinates(v) for v in work.vertices()}
        # drawn holds the two offset curves per guide, left then right
        chains = [(rails[gi][side], drawn[2 * gi + side])
                  for gi in range(n) for side in (0, 1)]
        swing_columns(work, chains, swing, poly,
                      {v for left, right in rails for v in left + right})
        if folded_faces(work) or boundary_scrambled(work, poly):
            for v, p in before.items():
                work.vertex[v]['x'], work.vertex[v]['y'], work.vertex[v]['z'] = p
            swung = '; swing SKIPPED (would fold the outline)'
        else:
            swung = '; swung x{}'.format(swing)

    work.set_strips_density_target(t=target_length)

    blocks = ''
    if block_rows:
        # One band is exactly one strip, and its density is the number of rows
        # ACROSS it -- so pinning it to 1 makes the guide the centreline of a
        # single row of faces instead of a line of edges. Look the strip up by a
        # rung rather than trusting the key add_strip returned: collect_strips
        # has re-run since, and renumbers.
        for left, right in rails:
            skey = None
            for i in range(1, len(left) - 1):
                for edge in ((left[i], right[i]), (right[i], left[i])):
                    if work.has_edge(edge):
                        skey = work.edge_strip(edge)
                        break
                if skey is not None:
                    break
            if skey is None:
                blocks = '; block_rows FAILED (no rung found)'
                break
            work.set_strip_density(skey, block_rows)
        else:
            blocks = '; {} row(s)/band'.format(block_rows)

    def boundary_arc_between(pa, pb, poly):
        """Arc of the closed boundary polygon `poly` between the points on it
        closest to `pa` and `pb`, oriented pa -> pb. `poly` need not already be
        closed (both `boundary_polygon`'s open loop and a caller-supplied
        explicitly-closed curve are accepted).

        This is what lets an edge whose endpoints add_strip/snap/relax/swing
        MOVED still densify along the true curve: both endpoints are still ON
        the boundary (that is what boundary reprojection guarantees), just
        not where they started, so re-deriving the arc from their CURRENT
        positions works where reusing the ORIGINAL curve (which no longer
        meets them) cannot.
        """
        pts = list(poly)
        if distance_point_point(pts[0], pts[-1]) > 1e-9:
            pts = pts + [pts[0]]
        n = len(pts) - 1                            # unique vertices
        if n < 2:
            return None
        doubled = pts[:-1] + pts[:-1]                # two laps, no wraparound needed
        segs = list(pairwise(doubled))
        cum = [0.0]
        for a, b in segs:
            cum.append(cum[-1] + distance_point_point(a, b))
        total = cum[n]                               # length of exactly one lap

        def project(p):
            best_d, best_s = float('inf'), None
            for (a, b), s0 in zip(segs[:n], cum[:n]):
                L = distance_point_point(a, b)
                if L == 0:
                    continue
                ab = subtract_vectors(b, a)
                t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / (L * L)))
                q = [a[k] + ab[k] * t for k in range(3)]
                d = distance_point_point(p, q)
                if d < best_d:
                    best_d, best_s = d, s0 + t * L
            return best_d, best_s

        _, sa = project(pa)
        _, sb = project(pb)
        if sa is None or sb is None:
            return None

        fwd = (sb - sa) % total
        bwd = (sa - sb) % total
        reverse = fwd > bwd
        s0, length = (sb, bwd) if reverse else (sa, fwd)
        s1 = s0 + length                             # <= s0 + total < 2 * total: safe, no wrap

        def point_at(s):
            for i, (c0, c1) in enumerate(zip(cum, cum[1:])):
                if c0 - 1e-9 <= s <= c1 + 1e-9:
                    a, b = segs[i]
                    L = c1 - c0
                    t = 0.0 if L == 0 else (s - c0) / L
                    return [a[k] + (b[k] - a[k]) * t for k in range(3)]
            return list(doubled[-1])

        arc = [point_at(s0)]
        for i, c in enumerate(cum):
            if s0 < c < s1:
                arc.append(list(doubled[i]))
        arc.append(point_at(s1))
        if reverse:
            arc.reverse()
        return arc

    def curves_for_densification(mesh, given):
        """For a boundary edge, prefer re-deriving the curve from `poly`
        between the endpoints' CURRENT positions (see `boundary_arc_between`)
        -- this is what makes the boundary densify along the true curve even
        where a guide moved it. Falls back to `given` (curves carried through
        refine by `prepared`, keyed by `mesh`'s OWN vertex indices where a
        split happened -- see `quad_split`) for interior edges, then a
        straight chord as the last resort. A stale or mismatched curve would
        no longer meet the edge's actual endpoints and the dense mesh would
        fail to weld there -- see CoarsePseudoQuadMesh.densification's
        edges_to_curves.
        """
        curves = {}
        for u, v in mesh.edges():
            a, b = mesh.vertex_coordinates(u), mesh.vertex_coordinates(v)
            curve = None
            if mesh.is_vertex_on_boundary(u) and mesh.is_vertex_on_boundary(v):
                curve = boundary_arc_between(a, b, poly)
                if curve is not None:
                    curve[0], curve[-1] = list(a), list(b)
            if curve is None:
                curve = (given or {}).get((u, v))
                if curve is None:
                    rev = (given or {}).get((v, u))
                    curve = list(reversed(rev)) if rev else None
            if (curve is None
                    or distance_point_point(curve[0], a) > 1e-6
                    or distance_point_point(curve[-1], b) > 1e-6):
                curve = [a, b]
            curves[u, v] = curve
        return curves

    work.densification(edges_to_curves=curves_for_densification(work, work_curves))
    dense = work.get_quad_mesh()

    info = {
        'coarse': work,
        'guides': drawn + pts_list,
        'centrelines': pts_list,
        'rails': rails,
        'refine': refine,
        'note': '; '.join(notes) + '; refined x{}{}{}'.format(refine, swung, blocks),
    }
    return dense, info


# =============================================================================
# GUIDE C -- let the guide shape the layout, as a curve feature
# =============================================================================

def guide_feature_mesh(outer_boundary, guides, inner_boundaries=(), point_features=(),
                       target_length=0.6, edges_to_curves=None):
    """Build the coarse layout AROUND the guides, by passing them as curve features.

    The other two functions in this module take a coarse mesh that already
    exists and move it onto the guides. This one is the third option, and the
    one Oval's thesis describes (SS4.3.2, Figs 4.17-4.19): hand the guides to
    :func:`~compas_singular.algorithms.boundary_triangulation` as
    ``polyline_features``, which cuts the Delaunay along them so the medial axis
    cannot cross them, and let the skeleton decomposition lay the patches out
    around the cut.

    Which of the three to reach for
    -------------------------------
    * :func:`guide_line_mesh` / :func:`guide_band_mesh` -- you HAVE a layout you
      want to keep and a guide you want it to follow. The patch layout is yours;
      the guide bends to it (line) or buys a strip (band).
    * :func:`guide_feature_mesh` -- the guide is a HARD constraint on the design
      and you are willing to let it decide the layout. The guide is embedded as
      a continuous course by construction rather than by snapping, and its
      extremities become singularities of the layout: three-valent where a guide
      ends on the boundary, a pole where it ends in the interior.

    The price is that you no longer choose the patch layout, and that a guide
    which ends in the interior costs a singularity that the other two routes do
    not spend.

    Parameters
    ----------
    outer_boundary : list
        The outer boundary as a list of XYZ coordinates.
    guides : Polyline | list
        One curve or a list of curves, in the forms this module accepts
        everywhere: a ``compas.geometry.Polyline``, a list of ``Point``, or a
        list of ``[x, y, z]``. Guides that meet are welded at the junction; see
        :func:`~compas_singular.algorithms.weld_polyline_features`.
    inner_boundaries : list, optional
        Holes, as lists of XYZ coordinates.
    point_features : list, optional
        Points to become poles.
    target_length : float, optional
        Target edge length of the pattern.
    edges_to_curves : dict, optional
        Passed through to ``densification`` so coarse edges densify along their
        curve rather than their chord.

    Returns
    -------
    (dense, info)
        ``info`` carries ``coarse``, ``guides``, ``centrelines``, ``rails``
        (empty -- this route snaps nothing), ``repair_notes`` from the
        decomposition, and a ``note``.

    """
    from compas.geometry import is_point_in_polygon_xy

    from compas_singular.algorithms import as_curves
    from compas_singular.algorithms import as_points
    from compas_singular.algorithms import boundary_triangulation
    from compas_singular.algorithms import SkeletonDecomposition

    def spacing_of(boundary):
        """The domain's own discretisation, as the median boundary segment."""
        lengths = sorted(distance_point_point(a, b) for a, b in pairwise(boundary + boundary[:1]))
        return lengths[len(lengths) // 2] if lengths else 1.0

    def resample(curve, spacing):
        """Walk the curve and emit a point every `spacing`, keeping both ends.

        A curve feature has to be sampled at the density of everything else, or
        its segments are not edges of the Delaunay and the topological cut has
        nothing to cut along. This is what ``discrete_mapping`` does to a Rhino
        curve before handing it over.
        """
        out = [list(curve[0])]
        for a, b in pairwise(curve):
            length = distance_point_point(a, b)
            if length <= 0.0:
                continue
            for i in range(1, max(1, int(round(length / spacing))) + 1):
                t = i / float(max(1, int(round(length / spacing))))
                out.append([a[k] + t * (b[k] - a[k]) for k in range(3)])
        return out

    def clip(curve, boundary):
        """Keep the longest run of the curve that lies inside the domain.

        A guide drawn past the wall -- the usual way to say "all the way across"
        -- has its outside points thrown away by ``boundary_triangulation``,
        which leaves the cut short of the wall and the guide with no effect at
        all. Trim it here instead, and snap the trimmed end onto the nearest
        boundary point so the extremity lands ON the wall: the thesis's
        three-valent boundary vertex rather than an interior pole.
        """
        inside = [is_point_in_polygon_xy(pt, boundary) for pt in curve]
        runs, run = [], []
        for pt, ok in zip(curve, inside):
            if ok:
                run.append(pt)
            elif run:
                runs.append(run)
                run = []
        if run:
            runs.append(run)
        if not runs:
            return [], True
        best = max(runs, key=lambda r: sum(distance_point_point(a, b) for a, b in pairwise(r)) if len(r) > 1 else 0.0)
        clipped = len(best) != len(curve)
        if clipped:
            # snap each trimmed end to the nearest boundary point, so the cut
            # reaches the wall and shares its vertex with the boundary
            first, last = curve.index(best[0]), curve.index(best[-1])
            if first > 0:
                best[0] = list(min(boundary, key=lambda q: distance_point_point(q, best[0])))
            if last < len(curve) - 1:
                best[-1] = list(min(boundary, key=lambda q: distance_point_point(q, best[-1])))
        return best, clipped

    curves = as_curves(guides)
    outer = as_points(outer_boundary, close=False)
    inners = as_curves(inner_boundaries, close=False)
    poles = [as_points([p])[0] for p in (point_features or [])]

    spacing = spacing_of(outer)
    prepared, skipped = [], []
    for i, curve in enumerate(curves):
        curve, was_clipped = clip(resample(curve, spacing), outer)
        if len(curve) < 2:
            skipped.append(i)
            continue
        prepared.append(curve)
    curves = prepared

    trimesh = boundary_triangulation(outer, inners, curves, poles)
    decomposition = SkeletonDecomposition.from_mesh(trimesh)
    coarse = decomposition.decomposition_mesh(poles)

    coarse.collect_strips()
    coarse.set_strips_density_target(target_length)
    if edges_to_curves:
        coarse.densification(edges_to_curves=edges_to_curves)
    else:
        coarse.densification()
    dense = coarse.get_quad_mesh()

    notes = list(getattr(decomposition, 'repair_notes', []))
    if skipped:
        notes.append('guide(s) {} lie entirely outside the domain and were skipped'.format(skipped))
    info = {
        'coarse': coarse,
        'guides': curves,
        'centrelines': curves,
        'rails': [],
        'refine': 0,
        'repair_notes': notes,
        'note': '{} guide(s) embedded as curve features; {} coarse patches{}'.format(
            len(curves), coarse.number_of_faces(),
            '; {} repair note(s)'.format(len(notes)) if notes else ''),
    }
    return dense, info


# =============================================================================
# Worked example -- the pentagon floor plate
# =============================================================================

if __name__ == '__main__':

    import sys
    from math import pi, sin, acos, degrees

    from compas.geometry import Line, Point, Polyline, Translation

    from compas_singular.algorithms import boundary_triangulation, SkeletonDecomposition

    # -- the design domain ----------------------------------------------------

    corners = [Point(6, 5, 0), Point(5, -5, 0), Point(-5, -5, 0),
               Point(-6, 5, 0), Point(0, 7, 0)]
    corners += [corners[0]]

    outer = []
    for a, b in pairwise(corners):
        outer.extend([[p.x, p.y, 0.0] for p in Line(a, b).to_polyline(n=10).points])

    class PlainDecomposition(SkeletonDecomposition):
        """No point features are used here, so there are no poles to record."""

        def store_pole_data(self, poles):
            self.mesh.attributes['face_pole'] = {}

    def baseline():
        trimesh = boundary_triangulation(outer_boundary=outer, inner_boundaries=[],
                                         polyline_features=[], point_features=[])
        mesh = PlainDecomposition.from_mesh(trimesh).decomposition_mesh([])
        mesh.collect_strips()
        return mesh

    # -- guide curves, as plain input curves ----------------------------------

    def bow(start, end, bulge, samples=48):
        """A half-sine bow off the chord -- flat in curvature where it meets the
        boundary rather than kinking into it. bulge=0 gives the straight chord."""
        chord = subtract_vectors(end, start)
        normal = normalize_vector([-chord[1], chord[0], 0.0])
        return Polyline([Point(*[start[i] + chord[i] * (k / samples)
                                 + normal[i] * bulge * sin(pi * k / samples)
                                 for i in range(3)]) for k in range(samples + 1)])

    def line(a, b):
        return Polyline([Point(*a), Point(*b)])

    # half_width is left at its default throughout -- target_length / 2, so each
    # band is one row of blocks the same size as the quads around it
    cases = [
        ('one straight', [line([-8, -1.5, 0], [8, -1.5, 0])], {}),
        ('one oblique', [line([-8, -2.7, 0], [8, -0.3, 0])], {}),            # ~8.5 deg
        ('one curved', [bow([-5.5, -1.5, 0], [5.5, -1.5, 0], 1.1)], {}),
        ('two parallel', [line([-8, -2.5, 0], [8, -2.5, 0]),
                          line([-8, 1.5, 0], [8, 1.5, 0])], {}),
        ('three parallel', [line([-8, -3.0, 0], [8, -3.0, 0]),
                            line([-8, 0.0, 0], [8, 0.0, 0]),
                            line([-8, 3.0, 0], [8, 3.0, 0])], {}),
        # the crossing guide has to run somewhere the plate already has a
        # near-vertical line: its decomposition has none near x = 0
        ('two crossing', [line([-8, -1.0, 0], [8, -1.0, 0]),
                          line([2.0, -8, 0], [2.0, 8, 0])], {}),
    ]

    # -- measurement ----------------------------------------------------------

    def off_centre(mesh, guide_polylines):
        """How far each guide runs from the middle of the faces it passes through.

        In units of half the face width across the guide, so 0 is dead centre
        and 1 is on the face's edge. That is the whole difference between the
        two readings: a guide on an edge scores ~1, a guide down the middle of a
        row of blocks scores ~0.
        """
        def median(values):
            s = sorted(values)
            k = len(s)
            return s[k // 2] if k % 2 else 0.5 * (s[k // 2 - 1] + s[k // 2])

        worst, all_ratios = 0.0, []
        for guide in guide_polylines:
            for fkey in mesh.faces():
                ps = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
                c = [sum(p[i] for p in ps) / len(ps) for i in range(3)]
                best = (float('inf'), None)
                for a, b in pairwise(guide):
                    ab = subtract_vectors(b, a)
                    length2 = dot_vectors(ab, ab)
                    if length2 == 0.0:
                        continue
                    t = max(0.0, min(1.0, dot_vectors(subtract_vectors(c, a), ab) / length2))
                    q = [a[i] + ab[i] * t for i in range(3)]
                    d = distance_point_point(c, q)
                    if d < best[0]:
                        best = (d, normalize_vector(ab))
                d, tangent = best
                if tangent is None:
                    continue
                # only faces the guide actually passes through
                if d >= max(distance_point_point(c, p) for p in ps):
                    continue
                normal = [-tangent[1], tangent[0], 0.0]
                proj = [dot_vectors(subtract_vectors(p, c), normal) for p in ps]
                half = (max(proj) - min(proj)) / 2.0
                if half > 1e-9:
                    all_ratios.append(d / half)
                    worst = max(worst, d / half)
        if not all_ratios:
            return ''
        return ' | off-centre {:4.2f} med, {:4.2f} worst'.format(median(all_ratios), worst)

    def measure(mesh, guide_polylines, on_tol=0.05, margin=0.0, steep=12.0):
        """Courses ALONG a guide want 0 deg to the tangent; interfaces ACROSS it
        want 90.

        ACROSS is also reported for the stations where the guide is most tilted
        (over ``steep`` degrees off horizontal), because the plain median hides
        the only failure worth seeing. Where a guide runs flat, the mesh scores
        ~90 whether or not it did anything -- the unguided mesh is already square
        there. It is where the guide tilts that a course has to be turned to
        follow it, and for a bowed guide that is at the ends.

        Which is also why ``margin`` now defaults to 0. Excluding a fan around
        the boundary looks reasonable and, on a bow, deletes exactly the tilted
        part: same mesh, margin 0 -> 76.4 deg / 62% within 15; 1.5 -> 78.7 / 79%;
        2.0 -> 80.6 / 88%. None of those three numbers is about the boundary.
        """

        def near(p, poly):
            best = (float('inf'), None)
            for a, b in pairwise(poly):
                ab = subtract_vectors(b, a)
                length2 = dot_vectors(ab, ab)
                if length2 == 0.0:
                    continue
                t = max(0.0, min(1.0, dot_vectors(subtract_vectors(p, a), ab) / length2))
                q = [a[i] + ab[i] * t for i in range(3)]
                d = distance_point_point(p, q)
                if d < best[0]:
                    best = (d, normalize_vector(ab))
            return best

        def acute(u, v):
            return degrees(acos(max(-1.0, min(1.0, abs(dot_vectors(u, v))))))

        def median(values):
            s = sorted(values)
            k = len(s)
            return s[k // 2] if k % 2 else 0.5 * (s[k // 2 - 1] + s[k // 2])

        on = {}
        for vkey in mesh.vertices():
            d, tangent = min((near(mesh.vertex_coordinates(vkey), g)
                              for g in guide_polylines), key=lambda dt: dt[0])
            if d <= on_tol and tangent is not None:
                on[vkey] = tangent

        bnd = [mesh.vertex_coordinates(v) for loop in mesh.boundaries() for v in loop]

        along, across = [], []
        for u, v in mesh.edges():
            pu, pv = mesh.vertex_coordinates(u), mesh.vertex_coordinates(v)
            if distance_point_point(pu, pv) == 0.0:
                continue
            edge_dir = normalize_vector(subtract_vectors(pv, pu))
            on_u, on_v = u in on, v in on
            if on_u and on_v:
                along.append(acute(edge_dir, on[u]))
            elif on_u != on_v:
                gv = u if on_u else v
                if margin and any(distance_point_point(mesh.vertex_coordinates(gv), q) < margin
                                  for q in bnd):
                    continue                                     # endpoint fan
                tilt = acute(on[gv], [1.0, 0.0, 0.0])
                across.append((tilt, acute(edge_dir, on[gv])))

        out = 'on-guide {:3d}'.format(len(on))
        if along:
            out += ' | along {:5.1f} deg'.format(median(along))
        if across:
            angles = [a for _, a in across]
            within = 100.0 * sum(1 for a in angles if a >= 75.0) / len(angles)
            out += ' | across {:5.1f} deg, {:3.0f}% within 15'.format(median(angles), within)
            tilted = [a for t, a in across if t >= steep]
            out += ' | tilted {}'.format(
                '{:5.1f} deg'.format(median(tilted)) if tilted else '    --   ')
        return out

    # -- run ------------------------------------------------------------------

    scenes = []
    print('=' * 78)
    print('PENTAGON FLOOR PLATE   guide_line_mesh (F) vs guide_band_mesh (B)')
    print('=' * 78)

    for name, curves, kwargs in cases:
        for label, fn in (('F line', guide_line_mesh), ('B band', guide_band_mesh)):
            opts = dict(kwargs) if label == 'B band' else {}
            title = '{:14s} {}'.format(name, label)
            try:
                dense, info = fn(baseline(), curves, **opts)
            except ValueError as e:
                print('\n{}\n  REJECTED: {}'.format(title, e))
                continue
            print('\n{}\n  {}'.format(title, info['note']))
            print('  dense: {} faces | {}'.format(
                dense.number_of_faces(), measure(dense, info['guides'])))
            print('  {:>11s}          | {}'.format(
                '', off_centre(dense, info['centrelines']).lstrip(' |')))
            # the same thing without swing, to show what snapping alone gives
            try:
                flat, flat_info = fn(baseline(), curves, swing=0, **opts)
                print('  {:>11s}   swing=0 | {}'.format(
                    '', measure(flat, flat_info['guides'])))
            except ValueError:
                pass
            scenes.append((title, info['coarse'], dense, info['guides']))

    if '--view' in sys.argv:
        from compas_viewer import Viewer

        viewer = Viewer()
        step = 16.0
        for i, (title, coarse_mesh, dense, guide_polylines) in enumerate(scenes):
            T = Translation.from_vector([i * step, 0.0, 0.0])
            group = viewer.scene.add_group(name=title)
            group.add(dense.transformed(T), name='dense quad mesh')

            patches = viewer.scene.add_group(name='coarse patches', parent=group)
            for u, v in coarse_mesh.edges():
                patches.add(Line(Point(*coarse_mesh.vertex_coordinates(u)),
                                 Point(*coarse_mesh.vertex_coordinates(v))).transformed(T),
                            linewidth=3)

            # lift the guides clear of the mesh plane, otherwise they z-fight
            # with the very course they produced and are invisible
            drawn = viewer.scene.add_group(name='guide', parent=group)
            for g in guide_polylines:
                drawn.add(Polyline([Point(p[0], p[1], 0.15) for p in g]).transformed(T),
                          linewidth=6)

        span = (len(scenes) - 1) * step + 14.0
        centre = [0.5 * (len(scenes) - 1) * step, 1.0, 0.0]
        viewer.renderer.camera.target = centre
        viewer.renderer.camera.position = [centre[0], centre[1] - 0.92 * span, 0.58 * span]
        viewer.show()
