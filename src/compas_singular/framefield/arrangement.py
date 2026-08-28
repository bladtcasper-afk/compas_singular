"""Step 6a.5 -- the planar arrangement ``from_polylines`` assumes it is given.

``CoarsePseudoQuadMesh.from_polylines`` hands every polyline SEGMENT to
``Mesh.from_lines`` and recovers faces from the resulting graph. ``from_lines``
keys its vertices off segment endpoints and nothing else: it never computes an
intersection. So two separatrices that cross in their interiors contribute four
segment ends and no node at the crossing, and the graph the face search walks is
not a planar embedding at all. The walk goes in one side of the crossing and out
of a branch it was never supposed to reach, and the face it recovers VISITS THE
SAME VERTEX TWICE.

That is the whole failure. Measured on the raw ``from_polylines`` output, before
any repair:

    disc @0.4     3 interior crossings   ->  8 of 14 faces repeat a vertex,
                                             one of them 14-sided
    hexagon @0.5  2 interior crossings   ->  1 face of 7, 19-sided
    ellipse @0.5  2 interior crossings   ->  1 face of 3, 14-sided
    stadium @0.5  8 interior crossings   ->  2 faces of 5, one 32-sided,
                                             and the mesh is non-manifold

``solve_non_quad_faces`` then quad-splits those n-gons, and a quad split of a
face that doubles back on itself is a fan of inverted quads: 58 of 140 faces on
the disc come out with non-positive area, and ``densifiable`` rejects the layout
at ``area2 <= 0.0``. The layout is then discarded for the polygon fallback --
which is why a curved domain silently stops following its own field.

The tell that this is a discrete bug and not a geometric limit: the hexagon
fails at target_length 0.5 and 0.4 and SUCCEEDS at 0.35. A genuine geometric
obstruction does not come and go with background resolution. The number of
crossings does.

So compute the arrangement properly first. Four things, in order:

1. **Land the loose ends.** A separatrix endpoint sitting within ``tol`` of
   another polyline's interior -- typically a trace that stopped just short of
   the wall -- is moved onto that polyline and the polyline is split there. That
   turns a degree-1 dangling end into a proper T-junction node.
2. **Split at every crossing.** Every pair of polylines, every pair of their
   segments, one shared node per crossing.
3. **Make the node set watertight.** Split points near an existing node reuse
   that node's exact coordinates rather than adding a second one beside it --
   ``TOL.geometric_key`` matches to three decimals, so "nearly the same point"
   is a different point and a sliver arc is worse than no split.
4. **Prune what is still dangling**, and drop doubled edges.

The output is the same ``(boundary, others)`` pair ``build_network`` returns, so
this module drops in between that and ``from_polylines`` and nothing downstream
changes shape.
"""
from math import atan2

from compas.geometry import distance_point_point
from compas.itertools import pairwise


__all__ = ['planar_arrangement', 'faces_with_repeated_vertices',
           'count_interior_crossings', 'count_dangling_ends',
           'faces_from_arrangement']


# ``TOL.geometric_key`` rounds to three decimals by default, and that is what
# ``from_polylines`` matches endpoints with. Any arc shorter than this is a node
# as far as the mesh is concerned, so producing one is producing a self-loop.
GKEY_RESOLUTION = 1e-3


def _seg_seg(p1, p2, p3, p4):
    """Parameters ``(t, u)`` where two segments meet, or ``None``.

    Endpoints count as meeting: a crossing that happens to fall on a polyline
    vertex still has to become a node. Collinear overlaps return ``None`` --
    there is no single crossing point to insert, and duplicate-edge removal is
    what handles those.
    """
    d1x, d1y = p2[0] - p1[0], p2[1] - p1[1]
    d2x, d2y = p4[0] - p3[0], p4[1] - p3[1]
    den = d1x * d2y - d1y * d2x
    if abs(den) < 1e-14:
        return None
    ox, oy = p3[0] - p1[0], p3[1] - p1[1]
    t = (ox * d2y - oy * d2x) / den
    u = (ox * d1y - oy * d1x) / den
    if -1e-9 <= t <= 1.0 + 1e-9 and -1e-9 <= u <= 1.0 + 1e-9:
        return max(0.0, min(1.0, t)), max(0.0, min(1.0, u))
    return None


def _closest_on(points, p):
    """``(segment index, t, point, distance)`` of the closest point on a polyline."""
    best = (0, 0.0, list(points[0]), float('inf'))
    for i, (a, b) in enumerate(pairwise(points)):
        abx, aby = b[0] - a[0], b[1] - a[1]
        length2 = abx * abx + aby * aby
        if length2 == 0.0:
            continue
        t = ((p[0] - a[0]) * abx + (p[1] - a[1]) * aby) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + abx * t, a[1] + aby * t, 0.0]
        d = distance_point_point(p, q)
        if d < best[3]:
            best = (i, t, q, d)
    return best


def _length_of(points):
    return sum(distance_point_point(a, b) for a, b in pairwise(points))


def _distance_to_loop(p, loop):
    """Shortest distance from a point to a closed loop, in XY."""
    best = float('inf')
    for a, b in pairwise(list(loop) + list(loop[:1])):
        abx, aby = b[0] - a[0], b[1] - a[1]
        length2 = abx * abx + aby * aby
        if length2 == 0.0:
            continue
        t = ((p[0] - a[0]) * abx + (p[1] - a[1]) * aby) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + abx * t, a[1] + aby * t, 0.0]
        best = min(best, distance_point_point(p, q))
    return best


def _cumulative(loop):
    """Arc-length parameter of each loop vertex, plus the loop's total length."""
    cum = [0.0]
    for a, b in pairwise(list(loop) + list(loop[:1])):
        cum.append(cum[-1] + distance_point_point(a, b))
    return cum


def _loop_param(loop, cum, p):
    """Arc-length position of the loop point closest to ``p``."""
    best_d, best_s = float('inf'), 0.0
    for i, (a, b) in enumerate(pairwise(list(loop) + list(loop[:1]))):
        abx, aby = b[0] - a[0], b[1] - a[1]
        length2 = abx * abx + aby * aby
        if length2 == 0.0:
            continue
        t = ((p[0] - a[0]) * abx + (p[1] - a[1]) * aby) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + abx * t, a[1] + aby * t, 0.0]
        d = distance_point_point(p, q)
        if d < best_d:
            best_d, best_s = d, cum[i] + t * (cum[i + 1] - cum[i])
    return best_d, best_s


def _partition_loop(loop, marks, tol):
    """Cut a closed loop into arcs at ``marks``: a true PARTITION, once round.

    ``build_network._split_loop`` walks the loop between consecutive marks by
    VERTEX INDEX, and when two marks sort into the wrong order relative to each
    other its walk runs the long way round and the arcs it emits overlap.
    Measured on the stadium: 17 arcs covering the wall about one and a half
    times, with ``(0,0) -> (12.93,0.15)`` and ``(0,0.98) -> (12.93,0.15)`` both
    present -- the second containing the first. Every downstream stage then sees
    a boundary that crosses itself, and no amount of splitting at crossings can
    repair a wall that is drawn twice.

    Sorting by ARC LENGTH instead makes the order total and unambiguous, so the
    arcs emitted here tile the loop exactly once and share only their endpoints.
    """
    cum = _cumulative(loop)
    total = cum[-1]
    if total <= 0.0:
        return []

    placed = []
    for p in marks:
        d, s = _loop_param(loop, cum, p)
        if d > tol:
            continue
        placed.append((s % total, [float(p[0]), float(p[1]), 0.0]))

    if len(placed) < 2:
        # fewer than two cuts leaves a closed arc, which bounds a degenerate
        # patch; the same case SkeletonDecomposition guards with
        # branches_splitting_collapsed_boundaries.
        placed = [(total * k / 3.0, None) for k in range(3)]

    placed.sort(key=lambda m: m[0])
    unique = [placed[0]]
    for s, p in placed[1:]:
        if s - unique[-1][0] > max(tol * 0.25, GKEY_RESOLUTION):
            unique.append((s, p))
    if len(unique) > 1 and (total - unique[-1][0]) + unique[0][0] <= \
            max(tol * 0.25, GKEY_RESOLUTION):
        unique.pop()
    if len(unique) < 2:
        unique = [(total * k / 3.0, None) for k in range(3)]

    def point_at(s):
        s = s % total
        for i, (c0, c1) in enumerate(zip(cum, cum[1:])):
            if c0 - 1e-12 <= s <= c1 + 1e-12:
                a = loop[i]
                b = loop[(i + 1) % len(loop)]
                span = c1 - c0
                t = 0.0 if span == 0 else (s - c0) / span
                return [a[k] + (b[k] - a[k]) * t for k in range(3)]
        return [float(c) for c in loop[0]]

    arcs = []
    n = len(unique)
    for k in range(n):
        s0, p0 = unique[k]
        s1, p1 = unique[(k + 1) % n]
        a0 = list(p0) if p0 is not None else point_at(s0)
        a1 = list(p1) if p1 is not None else point_at(s1)
        s1w = s1 if s1 > s0 else s1 + total

        # the loop's OWN vertices in between, so a curved wall keeps its shape.
        # Sorted by arc length, not by index: an arc that wraps past the loop's
        # seam takes its high-index vertices FIRST and its low-index ones after.
        between = []
        for i, sv in enumerate(cum[:-1]):
            for s in (sv, sv + total):
                if s0 + 1e-9 < s < s1w - 1e-9:
                    between.append((s, [float(c) for c in loop[i]]))
        between.sort(key=lambda m: m[0])

        arc = [a0]
        for _, p in between:
            if distance_point_point(p, arc[-1]) > GKEY_RESOLUTION:
                arc.append(p)
        if distance_point_point(a1, arc[-1]) <= GKEY_RESOLUTION and len(arc) > 1:
            arc.pop()
        arc.append(a1)
        if len(arc) >= 2 and _length_of(arc) > GKEY_RESOLUTION:
            arcs.append(arc)
    return arcs


class _Nodes(object):
    """Positions that several polylines must agree on EXACTLY.

    Every shared point in the arrangement is looked up here and written back
    from here, so two polylines meeting at a crossing get byte-identical
    coordinates rather than two values that merely round to the same
    ``geometric_key``.
    """

    def __init__(self, tol, loops=None):
        self.tol = tol
        self.loops = list(loops or [])
        self.points = []

    def _on_wall(self, p):
        return any(_distance_to_loop(p, loop) < self.tol * 0.25 for loop in self.loops)

    def add(self, p, tol=None, exclude=None):
        """Return the canonical position for ``p``, creating a node if new.

        Two points that both sit on a WALL are held to a quarter of ``tol``.
        See ``repair._cluster`` for why: between two landings this tolerance
        stops being numerical and starts deciding the layout, and merging two
        landings deletes the patch between them. ``exclude`` names a node this
        point must not merge with -- the other end of its own chain, which by
        construction is a different node however short the chain is.
        """
        tol = self.tol if tol is None else tol
        wall = self._on_wall(p) if self.loops else False
        for q in self.points:
            if exclude is not None and q is exclude:
                continue
            limit = self.tol * 0.25 if (wall and self._on_wall(q)) else tol
            if distance_point_point(p, q) < min(limit, tol):
                return q
        q = [float(p[0]), float(p[1]), 0.0]
        self.points.append(q)
        return q


def _normalise(cut, n_points):
    """Move a cut sitting on a segment END onto the next segment's start.

    Without this a crossing at ``t = 1`` of segment ``i`` and the same crossing
    reported at ``t = 0`` of segment ``i + 1`` sort as two different cuts and the
    split between them is empty.
    """
    i, t = cut
    if t > 1.0 - 1e-9 and i + 1 < n_points - 1:
        return i + 1, 0.0
    return i, t


def _split_chain(points, cuts):
    """Cut a polyline at ``(segment, t, position)`` marks. Ends are kept."""
    if not cuts:
        return [[list(p) for p in points]]
    cuts = sorted(cuts, key=lambda c: (c[0], c[1]))

    pieces = []
    current = [list(points[0])]
    k = 0
    for i in range(len(points) - 1):
        while k < len(cuts) and cuts[k][0] == i:
            q = cuts[k][2]
            k += 1
            if distance_point_point(current[-1], q) > 1e-12:
                current.append(list(q))
            if len(current) >= 2:
                pieces.append(current)
                current = [list(q)]
            else:
                current = [list(q)]
        if distance_point_point(current[-1], points[i + 1]) > 1e-12:
            current.append(list(points[i + 1]))
    if len(current) >= 2:
        pieces.append(current)
    return pieces


def planar_arrangement(boundary, others, tol, report=None, loops=None):
    """Split every polyline at every crossing so the network is a planar graph.

    Parameters
    ----------
    boundary : list[polyline]
        Boundary arcs, in loop order. Order is preserved: ``_single_patch``
        walks them as a ring.
    others : list[polyline]
        Separatrices and walls.
    tol : float
        How far apart two points may be and still be the same node. The
        background target length is the right scale -- a traced separatrix that
        stopped short of the wall missed it by a fraction of one step.
    report : dict, optional
        Updated in place with what the arrangement had to do.
    loops : list[list[[x, y, z]]], optional
        The domain's own boundary loops, outer first. Given these, the wall arcs
        are rebuilt as a true partition of each loop rather than trusted as
        handed in -- see :func:`_partition_loop` for why that is necessary.

    Returns
    -------
    (list[polyline], list[polyline])
        The same two lists, split.
    """
    report = report if report is not None else {}

    # -- 0. the wall, tiled exactly once -----------------------------------
    overlaps = 0
    if loops:
        rebuilt = []
        for loop in loops:
            cum = _cumulative(loop)
            marks = []
            # every arc end build_network chose -- its corner detection and its
            # landing snapping are kept; only the WALK between them is redone --
            # and every separatrix end that reached this wall, which is a cut in
            # it whether or not build_network recorded a landing for it
            for pl in list(boundary) + list(others):
                for p in (pl[0], pl[-1]):
                    d, _ = _loop_param(loop, cum, p)
                    if d <= tol * 0.5:
                        marks.append(p)
            rebuilt.extend(_partition_loop(loop, marks, tol * 0.5))
        if rebuilt:
            # NOT a count of overlaps removed -- an overlapping arc set can have
            # the same size as the partition that replaces it. It is how many
            # arcs the wall is now made of, which is the number that has to
            # match the patches around it.
            overlaps = len(rebuilt)
            boundary = rebuilt

    chains = ([('b', [list(p) for p in pl]) for pl in boundary]
              + [('o', [list(p) for p in pl]) for pl in others])
    if not chains:
        return list(boundary), list(others)

    nodes = _Nodes(tol, loops)

    # -- 1. every existing endpoint is a node ------------------------------
    # Seeded first and in one pass, so that a crossing landing near one of them
    # snaps to it rather than the other way round. The arms of a singularity
    # were already clustered by ``build_network``; this keeps that agreement.
    for _, pts in chains:
        head = nodes.add(pts[0])
        pts[0] = list(head)
        pts[-1] = list(nodes.add(pts[-1], exclude=head))

    # -- 2. loose ends onto the polyline they nearly touch -----------------
    # This is the second of the two root causes: a separatrix that lands on the
    # boundary without the boundary being split there leaves a degree-1 node,
    # and the face walk goes UP the dangling arm and back down it -- which is
    # another way to visit a vertex twice.
    #
    # Done per NODE, not per chain end. Four arms of one singularity share a
    # node, and projecting each arm's end onto its own nearest polyline pulls
    # that node into four -- measured on the stadium, whose round end is dense
    # enough that 43 of 58 chain ends found something within ``tol`` and the
    # layout came apart. Whatever happens here happens to all the chains at once.
    #
    # A node that already joins two or more chains is a real junction and is
    # left alone unless it is close enough to be an overshoot. A DEGREE-1 node
    # is not a junction at all -- it has to attach to something or the region
    # around it cannot close -- so it reaches much further.
    near_tol = tol * 0.35
    far_tol = tol * 2.5
    landed = 0
    for _ in range(4):
        moves = []
        for node_p in list(nodes.points):
            owners = [(i, end) for i, (_, pts) in enumerate(chains)
                      for end in (0, -1)
                      if distance_point_point(pts[end], node_p) < 1e-9]
            if not owners:
                continue
            reach = far_tol if len(owners) < 2 else near_tol
            own = set(i for i, _ in owners)
            best = None
            for j, (kind_j, pts_j) in enumerate(chains):
                if j in own:
                    continue
                _, _, q, d = _closest_on(pts_j, node_p)
                if d >= reach or d < 1e-9:
                    continue
                # already a shared node -- an arm meeting another arm end on
                # does not need the other one split
                if min(distance_point_point(q, pts_j[0]),
                       distance_point_point(q, pts_j[-1])) < near_tol:
                    continue
                if best is None or d < best[0]:
                    best = (d, q)
            if best is not None:
                moves.append((node_p, best[1], owners))
        if not moves:
            break
        for node_p, q, owners in moves:
            # The BOUNDARY does not move: snapping a wall vertex to wherever a
            # trace happened to stop would pull the domain out of shape. The
            # loose end moves onto it instead, and the wall is split there
            # (step 3 sees the moved end as a crossing and cuts it).
            node_p[0], node_p[1] = q[0], q[1]
            for i, end in owners:
                chains[i][1][end] = list(node_p)
            landed += 1

    # -- 3. every crossing becomes a shared node ---------------------------
    # Run AFTER step 2, so the endpoint positions the crossing parameters are
    # computed from are final. An end that was moved onto another polyline now
    # touches it, and this is the step that actually splits it there.
    cuts = [[] for _ in chains]
    crossings = 0
    for i in range(len(chains)):
        kind_i, pts_i = chains[i]
        for j in range(i + 1, len(chains)):
            kind_j, pts_j = chains[j]
            if kind_i == 'b' and kind_j == 'b':
                # arcs of a domain wall meet at their ends and nowhere else
                continue
            for si, (a1, a2) in enumerate(pairwise(pts_i)):
                for sj, (b1, b2) in enumerate(pairwise(pts_j)):
                    hit = _seg_seg(a1, a2, b1, b2)
                    if hit is None:
                        continue
                    t, u = hit
                    p = [a1[0] + (a2[0] - a1[0]) * t,
                         a1[1] + (a2[1] - a1[1]) * t, 0.0]
                    # A crossing at an END of a polyline is already a node
                    # there: cut only the polyline it is interior to. Cutting
                    # both would leave an arc shorter than ``geometric_key``
                    # can distinguish, which is a self-loop, not an edge.
                    near_i = min(distance_point_point(p, pts_i[0]),
                                 distance_point_point(p, pts_i[-1])) < near_tol
                    near_j = min(distance_point_point(p, pts_j[0]),
                                 distance_point_point(p, pts_j[-1])) < near_tol
                    if near_i and near_j:
                        continue
                    shared = nodes.add(p, near_tol)
                    if not near_i:
                        cuts[i].append(_normalise((si, t), len(pts_i)) + (shared,))
                    if not near_j:
                        cuts[j].append(_normalise((sj, u), len(pts_j)) + (shared,))
                    crossings += 1

    # -- 4. split, then rebuild the graph ----------------------------------
    pieces = []
    for (kind, pts), marks in zip(chains, cuts):
        # two crossings landing on one node give the same cut twice
        unique = []
        for m in marks:
            if not any(distance_point_point(m[2], n[2]) < 1e-9
                       and m[0] == n[0] for n in unique):
                unique.append(m)
        for piece in _split_chain(pts, unique):
            # An EXACT lookup, not a snap: every piece end is already a node
            # position, either an original endpoint or an inserted cut. Snapping
            # at ``tol`` here would merge two genuinely distinct nodes that
            # happen to sit close, which collapses a short but real arc.
            piece[0] = list(nodes.add(piece[0], 1e-9))
            piece[-1] = list(nodes.add(piece[-1], 1e-9))
            pieces.append((kind, piece))

    # An arc shorter than the key resolution IS its own endpoint as far as
    # ``from_polylines`` is concerned. Drop it and let the two ends be one node.
    kept = []
    collapsed = 0
    for kind, piece in pieces:
        if (_length_of(piece) <= GKEY_RESOLUTION
                or distance_point_point(piece[0], piece[-1]) <= GKEY_RESOLUTION
                and _length_of(piece) < tol * 0.05):
            collapsed += 1
            continue
        kept.append((kind, piece))
    pieces = kept

    # -- 5. doubled edges --------------------------------------------------
    # ``from_lines`` cannot tell which side of a doubled edge a face is on. Two
    # pieces spanning the same node pair along the same route are one edge.
    def key_of(piece):
        a = tuple(round(c, 9) for c in piece[0][:2])
        b = tuple(round(c, 9) for c in piece[-1][:2])
        return (a, b) if a <= b else (b, a)

    seen = {}
    unique_pieces = []
    doubled = 0
    for kind, piece in pieces:
        k = key_of(piece)
        mid = piece[len(piece) // 2]
        hit = False
        for other_mid in seen.get(k, ()):
            if distance_point_point(mid, other_mid) < tol:
                hit = True
                break
        if hit:
            doubled += 1
            continue
        seen.setdefault(k, []).append(mid)
        unique_pieces.append((kind, piece))
    pieces = unique_pieces

    # -- 6. prune what is still dangling -----------------------------------
    # A trace that stopped in open space, far enough from anything that step 2
    # found nothing to land on. Its arm bounds no region, and leaving it in
    # makes the face around it walk up and back.
    pruned = 0
    while True:
        degree = {}
        for _, piece in pieces:
            for p in (piece[0], piece[-1]):
                k = tuple(round(c, 9) for c in p[:2])
                degree[k] = degree.get(k, 0) + 1
        drop = set()
        for idx, (kind, piece) in enumerate(pieces):
            if kind == 'b':
                continue            # the wall is a closed loop; never cut it
            for p in (piece[0], piece[-1]):
                if degree[tuple(round(c, 9) for c in p[:2])] < 2:
                    drop.add(idx)
        if not drop:
            break
        pieces = [pc for idx, pc in enumerate(pieces) if idx not in drop]
        pruned += len(drop)

    report.update({
        'wall_arcs_partitioned': overlaps,
        'crossings_split': crossings,
        'ends_landed': landed,
        'arcs_collapsed': collapsed,
        'doubled_dropped': doubled,
        'dangling_pruned': pruned,
    })

    out_boundary = [p for kind, p in pieces if kind == 'b']
    out_others = [p for kind, p in pieces if kind == 'o']
    return out_boundary, out_others


# ----------------------------------------------------------------------------
# the guard
# ----------------------------------------------------------------------------

def faces_with_repeated_vertices(mesh):
    """Faces that visit the same vertex -- or the same POSITION -- twice.

    This is the direct symptom of a face recovered from a graph that was not a
    planar embedding, and it is what makes the difference between a layout that
    densifies and one that gets thrown away. Checking it costs nothing and it
    fails loudly, whereas the inverted quads it turns into show up several
    frames later as ``area2 <= 0.0`` with no indication of where they came from.

    Positions are compared as well as keys because ``from_polylines`` can emit
    two distinct vertices at one coordinate -- measured on the ellipse, a quad
    with the corner (11.88, 5.19) twice.

    Returns
    -------
    list[face key]
    """
    out = []
    for fkey in mesh.faces():
        fv = mesh.face_vertices(fkey)
        if len(set(fv)) != len(fv):
            out.append(fkey)
            continue
        keys = [tuple(round(c, 6) for c in mesh.vertex_coordinates(v)) for v in fv]
        if len(set(keys)) != len(keys):
            out.append(fkey)
    return out


# ----------------------------------------------------------------------------
# measurement, used by the checks
# ----------------------------------------------------------------------------

def count_interior_crossings(polylines, tol):
    """Crossings that are in the INTERIOR of both polylines. Should be zero."""
    found = []
    for i in range(len(polylines)):
        for j in range(i + 1, len(polylines)):
            a, b = polylines[i], polylines[j]
            for si, (a1, a2) in enumerate(pairwise(a)):
                for sj, (b1, b2) in enumerate(pairwise(b)):
                    hit = _seg_seg(a1, a2, b1, b2)
                    if hit is None:
                        continue
                    t, u = hit
                    p = [a1[0] + (a2[0] - a1[0]) * t,
                         a1[1] + (a2[1] - a1[1]) * t, 0.0]
                    if any(distance_point_point(p, e) < tol
                           for e in (a[0], a[-1], b[0], b[-1])):
                        continue
                    if not any(distance_point_point(p, q) < tol for q in found):
                        found.append(p)
    return found


def count_dangling_ends(boundary, others, tol):
    """Endpoints that no second polyline reaches. Should be zero."""
    nodes, degree = [], []

    def node(p):
        for k, q in enumerate(nodes):
            if distance_point_point(p, q) < tol:
                return k
        nodes.append(p)
        degree.append(0)
        return len(nodes) - 1

    for pts in list(boundary) + list(others):
        degree[node(pts[0])] += 1
        degree[node(pts[-1])] += 1
    return [nodes[k] for k, d in enumerate(degree) if d < 2]


# ----------------------------------------------------------------------------
# face recovery
# ----------------------------------------------------------------------------

def _point_in_ring(p, ring):
    """Ray casting. ``ring`` is a closed loop of [x, y, z]."""
    x, y = p[0], p[1]
    inside = False
    n = len(ring)
    for i in range(n):
        ax, ay = ring[i][0], ring[i][1]
        bx, by = ring[(i + 1) % n][0], ring[(i + 1) % n][1]
        if (ay > y) != (by > y):
            t = (y - ay) / (by - ay)
            if x < ax + t * (bx - ax):
                inside = not inside
    return inside


def _interior_sample(points):
    """A point inside the polygon ``points``, or ``None``.

    The centroid is tried first and is right for anything convex. A half-annulus
    is not convex and its centroid can land in the hole, so fall back to
    midpoints between pairs of vertices and take the first that is inside.
    """
    n = len(points)
    cx = sum(p[0] for p in points) / n
    cy = sum(p[1] for p in points) / n
    if _point_in_ring([cx, cy], points):
        return [cx, cy, 0.0]
    for i in range(n):
        for j in range(i + 2, n):
            m = [(points[i][0] + points[j][0]) / 2.0,
                 (points[i][1] + points[j][1]) / 2.0, 0.0]
            if _point_in_ring(m, points):
                return m
    return None


def faces_from_arrangement(boundary, others, loops):
    """Recover the domain's faces from the polyline network directly.

    ``CoarsePseudoQuadMesh.from_polylines`` cannot do this for every domain. It
    keeps a face only when at least one of its vertices is NOT on a boundary
    polyline::

        if len(notonboundary):
            faces.append(indices)

    That is a reasonable proxy for "this face is not the outside" as long as some
    interior line exists, and it is why a square -- whose single patch has all
    four corners on the wall -- comes back with no faces at all and is built by
    hand in ``decomposition._single_patch``. **An annulus opened by two cuts is
    the same case, twice over:** the cuts run boundary to boundary, so every
    vertex of both half-annuli lies on a boundary loop and both real faces are
    discarded, leaving artifacts covering 16% of the domain.

    So the faces are traversed here instead, with no heuristic about which is
    which. Standard planar-subdivision walk: sort the half-edges leaving each
    node by angle, and from an arriving half-edge take the next one clockwise
    from its reverse. Every face of the embedding comes out exactly once. The
    ones that are not domain -- the unbounded face, and the inside of each hole
    -- are then rejected by testing a point genuinely inside each against the
    domain, rather than by guessing from the vertices.

    Requires a clean arrangement: run :func:`planar_arrangement` first, or a
    crossing with no node at it will merge two faces into one.

    Parameters
    ----------
    boundary, others : list[list[[x, y, z]]]
        The network, already arranged.
    loops : list[list[[x, y, z]]]
        Domain boundary loops, outer first. Used only to test containment.

    Returns
    -------
    list[list[[x, y, z]]]
        One entry per face, its corner points in order.
    """
    def node(p):
        return (round(p[0] / GKEY_RESOLUTION), round(p[1] / GKEY_RESOLUTION))

    # Each polyline is ONE edge with two half-edges, ``(index, +1)`` running
    # first-to-last and ``(index, -1)`` the other way. Keying half-edges by their
    # node pair instead would silently drop PARALLEL EDGES -- and the annulus is
    # made of them: its outer circle, split at the two cut landings, is two
    # distinct arcs between the same two nodes. One overwrites the other in a
    # dict and the traversal then walks a graph that is missing half its wall.
    edges = []
    for pts in list(boundary) + list(others):
        a, b = node(pts[0]), node(pts[-1])
        if a != b:
            edges.append((a, b, list(pts)))

    def ends(h):
        i, d = h
        a, b, pts = edges[i]
        return (a, b, pts) if d > 0 else (b, a, pts[::-1])

    out = {}
    for i in range(len(edges)):
        for d in (1, -1):
            a, _b, pts = ends((i, d))
            ang = atan2(pts[1][1] - pts[0][1], pts[1][0] - pts[0][0])
            out.setdefault(a, []).append((ang, (i, d)))
    for a in out:
        out[a].sort()

    faces = []
    seen = set()
    for i in range(len(edges)):
        for d in (1, -1):
            start = (i, d)
            if start in seen:
                continue
            cycle = []
            h = start
            while h not in seen:
                seen.add(h)
                cycle.append(h)
                _a, b, _pts = ends(h)
                ring = out.get(b)
                rev = (h[0], -h[1])
                pos = next((k for k, (_ang, e) in enumerate(ring) if e == rev), None)
                if pos is None:
                    cycle = None
                    break
                h = ring[(pos - 1) % len(ring)][1]
            if not cycle:
                continue

            points = []
            for e in cycle:
                _a, _b, pts = ends(e)
                points.extend(pts[:-1])
            if len(points) < 3:
                continue
            area2 = sum(p[0] * q[1] - q[0] * p[1]
                        for p, q in zip(points, points[1:] + points[:1]))
            if area2 <= 0.0:
                continue                   # unbounded face, or wound backwards

            sample = _interior_sample(points)
            if sample is None:
                continue
            if not _point_in_ring(sample, loops[0]):
                continue                   # outside the domain
            if any(_point_in_ring(sample, hole) for hole in loops[1:]):
                continue                   # inside a hole, not in the domain

            faces.append([list(ends(e)[2][0]) for e in cycle])
    return faces
