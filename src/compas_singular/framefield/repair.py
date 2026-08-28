"""Step 6a -- turn traced separatrices into a network ``from_polylines`` accepts.

``CoarsePseudoQuadMesh.from_polylines(boundary_polylines, other_polylines)`` is
the whole integration seam with compas_singular (see the design doc), but it is
unforgiving about what it is given. It matches polyline endpoints by
``TOL.geometric_key``, so "nearly the same point" is a different point, and it
recovers faces from a planar embedding, so polylines that cross anywhere other
than at a shared endpoint produce garbage.

A raw traced network satisfies neither condition, and this module is what closes
the gap:

* separatrices start half a step off their singularity, not on it;
* several separatrices reaching the same place land microns apart;
* the domain boundary is one loop, not the arcs between landing points.

The last one is not a numerical detail. Boundary arcs must ALSO be split at the
domain's own corners, or every patch comes out a triangle -- a pentagon's five
separatrices cut it into five regions bounded by two separatrices and one wall
each. Splitting that wall at the pentagon corner turns each into a proper
four-sided patch. ``SkeletonDecomposition.branches_boundary`` does the same
thing for the same reason.
"""
from math import pi

from compas.geometry import distance_point_point
from compas.geometry import angle_vectors
from compas.geometry import length_vector
from compas.geometry import subtract_vectors
from compas.itertools import pairwise


__all__ = ['build_network', 'boundary_corners', 'solve_non_quad_faces',
           'topological_quad_split']


def _arc_midpoint(loops, a, b, tol=1e-6):
    """Midpoint of the boundary arc between two points, if both lie on one loop.

    A quad split places each new vertex at the CHORD midpoint of the edge it
    splits. On a curved wall that midpoint sits inside the domain, so every split
    cuts the corner and the mesh shrinks away from its own boundary -- measured
    on an ellipse, coverage fell to 65% of the domain area, and 74% on a disc.
    Putting the new vertex on the arc instead keeps the boundary where it is.
    """
    for loop in loops or []:
        ring = list(loop) + list(loop[:1])
        best = [None, None]
        for k, p in enumerate((a, b)):
            bd, bs = float('inf'), None
            acc = 0.0
            for u, v in pairwise(ring):
                abx, aby = v[0] - u[0], v[1] - u[1]
                length = (abx * abx + aby * aby) ** 0.5
                if length == 0.0:
                    continue
                t = ((p[0] - u[0]) * abx + (p[1] - u[1]) * aby) / (length * length)
                t = max(0.0, min(1.0, t))
                q = [u[0] + abx * t, u[1] + aby * t, 0.0]
                d = distance_point_point(p, q)
                if d < bd:
                    bd, bs = d, acc + t * length
                acc += length
            if bd > tol:
                best = None
                break
            best[k] = (bs, acc)
        if best is None or best[0] is None or best[1] is None:
            continue

        sa, sb = best[0][0], best[1][0]
        total = best[0][1]
        forward = (sb - sa) % total
        span = min(forward, total - forward)

        # Both endpoints being on the loop is NOT enough: an edge can be a CHORD
        # across a hole, with both ends on the hole's rim and its interior well
        # inside the domain. Snapping such an edge's midpoint onto the rim moves
        # a vertex somewhere it does not belong -- measured on a round-holed
        # plate, that produced 114% coverage and a non-manifold mesh.
        #
        # A genuine boundary edge runs ALONG the wall, so its arc and its chord
        # have nearly the same length. A chord across a hole does not.
        chord = distance_point_point(a, b)
        if span > max(chord * 1.5, 1e-9):
            continue

        s = (sa + forward * 0.5) % total if forward <= total - forward \
            else (sb + (total - forward) * 0.5) % total

        acc = 0.0
        for u, v in pairwise(ring):
            abx, aby = v[0] - u[0], v[1] - u[1]
            length = (abx * abx + aby * aby) ** 0.5
            if length == 0.0:
                continue
            if acc <= s <= acc + length:
                t = (s - acc) / length
                return [u[0] + abx * t, u[1] + aby * t, 0.0]
            acc += length
    return None


def topological_quad_split(vertices, faces, loops=None):
    """Catmull-Clark TOPOLOGICAL split: every n-gon becomes n quads.

    Original vertices keep their position -- no smoothing, no geometry change,
    only connectivity. On any polygonal mesh the result is ALL QUADS, and it is
    conforming: neighbouring faces share the same new edge midpoints, so no
    T-vertices appear.

    That guarantee is the point. A separatrix network does not always partition
    into four-sided regions, and without a fallback those domains produce no mesh
    at all -- which is useless however good the field was. This always produces
    one.

    The price is honest and visible: an n-gon's centre becomes a valence-n
    vertex, so a 14-sided patch leaves a 14-valent singularity. A bad layout
    stays bad; it just stops being *nothing*.

    Returns
    -------
    (list[[x, y, z]], list[list[int]])
    """
    new_vertices = [list(p) for p in vertices]
    edge_mid = {}

    # Two distinct edges can snap to the SAME point on the wall, and a mesh with
    # two vertices at one position is not manifold. Measured on a round-holed
    # plate: 3 duplicated coordinates and a non-manifold result. Keep the
    # occupied positions and fall back to the chord when one is taken.
    def key_of(p):
        return (round(p[0], 9), round(p[1], 9))

    taken = set(key_of(p) for p in vertices)

    def mid(u, v):
        key = (min(u, v), max(u, v))
        if key not in edge_mid:
            edge_mid[key] = len(new_vertices)
            a, b = vertices[u], vertices[v]
            point = _arc_midpoint(loops, a, b)
            if point is None or key_of(point) in taken:
                point = [(a[i] + b[i]) / 2.0 for i in range(3)]
            taken.add(key_of(point))
            new_vertices.append(point)
        return edge_mid[key]

    new_faces = []
    for fv in faces:
        n = len(fv)
        centre = len(new_vertices)
        new_vertices.append([sum(vertices[v][i] for v in fv) / n for i in range(3)])
        mids = [mid(fv[i], fv[(i + 1) % n]) for i in range(n)]
        for i in range(n):
            new_faces.append([fv[i], mids[i], centre, mids[i - 1]])
    return new_vertices, new_faces


def _signed_area(pts):
    return 0.5 * sum(p[0] * q[1] - q[0] * p[1] for p, q in zip(pts, pts[1:] + pts[:1]))


def _min_angle(pts):
    """Smallest interior angle of a polygon, in radians."""
    n = len(pts)
    worst = pi
    for i in range(n):
        u = subtract_vectors(pts[(i - 1) % n], pts[i])
        v = subtract_vectors(pts[(i + 1) % n], pts[i])
        if length_vector(u) < 1e-12 or length_vector(v) < 1e-12:
            return 0.0
        worst = min(worst, angle_vectors(u, v))
    return worst


def _fan(fv):
    """Cut an n-gon into quads plus at most one triangle, along its DIAGONALS.

    Takes the first four corners as a quad and starts again from the fourth,
    which is the idiom ``RobustSkeletonDecomposition`` uses for the same job on
    the medial-axis front end. Only existing vertices are used, so the cut is a
    diagonal of the face and no neighbour sees anything change: the layout stays
    conforming, which is the whole reason a per-face split is normally not an
    option.
    """
    out = []
    fv = list(fv)
    while len(fv) > 4:
        out.append(fv[0:4])
        fv = fv[3:] + fv[0:1]
    out.append(fv)
    return out


def _best_fan(fv, coords):
    """The best rotation of :func:`_fan`, or ``None`` if none is valid.

    A diagonal of a NON-CONVEX face can pass outside it, and the pieces either
    side of such a diagonal are not a partition of the face -- they overlap, or
    they cover ground the face never did. Both show up as the pieces' areas
    failing to add up to the face's, which is what is checked here; every
    rotation that fails it is discarded, and the survivor with the largest
    minimum interior angle is kept.
    """
    n = len(fv)
    whole = _signed_area([coords[v] for v in fv])
    best = None
    for start in range(n):
        rotated = fv[start:] + fv[:start]
        pieces = _fan(rotated)
        areas = [_signed_area([coords[v] for v in piece]) for piece in pieces]
        if any(a <= 1e-12 for a in areas):
            continue
        if abs(sum(areas) - whole) > 1e-6 * max(abs(whole), 1.0):
            continue
        score = min(_min_angle([coords[v] for v in piece]) for piece in pieces)
        if best is None or score > best[0]:
            best = (score, pieces)
    return best[1] if best else None


def _choose_pole(fv, coords):
    """Which corner of a triangle to collapse the pseudo-quad's fourth side at.

    A pseudo-quad ``(p, a, b)`` is the quad ``(p, a, b, p)``: ``(p, a)`` and
    ``(b, p)`` are one strip -- ``face_opposite_edge`` maps each to the other --
    and ``(a, b)`` is the side opposite the collapsed corner. Densification
    therefore fans the face from ``p`` onto ``(a, b)`` with the two rails at one
    density, so the fan is even when the two rails are the same length. Pick the
    corner whose two edges are closest to equal.
    """
    best = None
    for i in range(3):
        p, a, b = fv[i], fv[(i + 1) % 3], fv[(i + 2) % 3]
        ra = distance_point_point(coords[p], coords[a])
        rb = distance_point_point(coords[b], coords[p])
        if ra + rb <= 0.0:
            continue
        imbalance = abs(ra - rb) / (ra + rb)
        if best is None or imbalance < best[0]:
            best = (imbalance, p)
    return best[1] if best else fv[0]


def solve_non_quad_faces(mesh, cls=None, loops=None, poles=()):
    """Make ``mesh`` a quad mesh, whatever it started as.

    Four stages, cheapest and most local first. Only the last one changes any
    face that was already fine, and it only runs if the others left something
    that will not densify:

    1. **Triangles against the wall.** Ported from
       ``SkeletonDecomposition.solve_triangular_faces``: a triangle with two
       boundary vertices at the same point is a degenerate quad, so merging them
       removes the face.

    2. **n-gons along their own diagonals** (:func:`_best_fan`). A pentagon
       becomes a quad and a triangle, a hexagon two quads, and so on. No vertex
       is added, so no neighbour gains a T-junction and the layout stays
       conforming -- which is what rules out the obvious per-face split, where
       midpoints on the offending face's edges leave its neighbours
       non-conforming and ``collect_strips`` comes apart.

    3. **Whatever triangles are left become POLES.** A triangle is a
       pseudo-quad: the quad ``(p, a, b, p)`` with one side collapsed.
       ``CoarsePseudoQuadMesh`` is built for exactly this -- ``face_pole``,
       ``face_opposite_edge`` and ``densification`` all handle it -- and the
       same trick resolves the interior triangles the medial-axis front end
       leaves behind (see ``INTERIOR_TRIANGLE_POLE_FIX.md``). The pole is an
       honest valence-3 singularity and it costs no extra patches.

    4. **The global** :func:`topological_quad_split`, only if 1-3 did not leave
       something ``densifiable``. It cannot fail, and it doubles the coarse
       resolution -- one stray triangle used to take the hexagon from 11 patches
       to 43 and one stray pentagon the ellipse from 7 to 30 -- so it is the
       last resort, not the first.

    Parameters
    ----------
    mesh : CoarsePseudoQuadMesh
    cls : type, optional
    loops : list[list[[x, y, z]]], optional
        Boundary loops, so a split vertex lands on the wall and not on a chord.
    poles : list[[x, y, z]], optional
        Preferred pole positions. A triangle with a corner at one of these
        collapses there rather than wherever stage 3 would have chosen.

    Returns
    -------
    (mesh, note)
        The mesh, and a one-line description of what was needed -- empty when
        the input was already all quads.
    """
    cls = cls or type(mesh)
    sides = [len(mesh.face_vertices(f)) for f in mesh.faces()]
    if not sides:
        return mesh, 'empty mesh'

    # ``from_polylines`` recovers faces from a planar embedding and does not
    # guarantee they all wind the same way. A backwards face breaks strip
    # propagation, so the two edges of a quad that ought to share a strip end up
    # in different ones with different densities -- and densification then hands
    # ``discrete_coons_patch`` two opposite sides of unequal length and dies with
    # an IndexError far from the cause. Measured on a hexagon: 2 faces wound
    # backwards out of 44, 11 faces left with mismatched opposite densities.
    mesh, flipped = _unify(mesh, cls)

    sides = [len(mesh.face_vertices(f)) for f in mesh.faces()]
    if set(sides) == {4}:
        return mesh, ('unified {} backwards face(s)'.format(flipped) if flipped else '')

    before = {}
    for n in sides:
        before[n] = before.get(n, 0) + 1

    # -- stage 1: triangles touching the boundary ---------------------------
    merged = 0
    for fkey in list(mesh.faces()):
        if fkey not in list(mesh.faces()) or len(mesh.face_vertices(fkey)) != 3:
            continue
        on_wall = [v for v in mesh.face_vertices(fkey) if mesh.is_vertex_on_boundary(v)]
        if len(on_wall) != 2:
            continue
        a, b = on_wall
        if distance_point_point(mesh.vertex_coordinates(a),
                                mesh.vertex_coordinates(b)) > 1e-6:
            continue                    # genuinely distinct corners, not a doubling
        mesh.delete_face(fkey)
        for other in list(mesh.faces()):
            fv = mesh.face_vertices(other)
            if b in fv:
                mesh.face[other] = [a if v == b else v for v in fv]
        merged += 1
    if merged:
        mesh = _rebuild(mesh, cls)

    sides = [len(mesh.face_vertices(f)) for f in mesh.faces()]
    if set(sides) == {4}:
        return mesh, 'merged {} doubled boundary vertex/ices'.format(merged)

    stray = {n: c for n, c in before.items() if n != 4}

    # -- stages 2 and 3: local, conforming, and no extra patches -------------
    local, fanned, poled = _local_repair(mesh, cls, poles)
    if local is not None:
        ok, _ = densifiable(local)
        if ok:
            note = []
            if fanned:
                note.append('fanned {} n-gon(s) along their own diagonals'.format(fanned))
            if poled:
                note.append('{} triangle(s) kept as pseudo-quad pole(s)'.format(poled))
            return local, 'resolved {} non-quad patch(es) {}: {}'.format(
                sum(stray.values()), stray, '; '.join(note))

    # -- stage 4: the guaranteed step ---------------------------------------
    index = {vkey: i for i, vkey in enumerate(mesh.vertices())}
    vertices = [mesh.vertex_coordinates(v) for v in mesh.vertices()]
    faces = [[index[v] for v in mesh.face_vertices(f)] for f in mesh.faces()]
    vertices, faces = topological_quad_split(vertices, faces, loops)
    out = cls.from_vertices_and_faces(vertices, faces)

    note = ('quad-split {} non-quad patch(es) {} -- the local repair did not '
            'leave a densifiable layout, so every n-gon centre is now a '
            'valence-n vertex, worst {}'.format(
                sum(stray.values()), stray, max(before)))
    return out, note


def _local_repair(mesh, cls, poles=()):
    """Stages 2 and 3: fan the n-gons, pole the triangles. No new vertices.

    Returns
    -------
    (mesh or None, int, int)
        The repaired mesh -- ``None`` if some face could not be fanned validly
        -- and how many faces were fanned and how many poles registered.
    """
    index = {vkey: i for i, vkey in enumerate(mesh.vertices())}
    coords = [mesh.vertex_coordinates(v) for v in mesh.vertices()]

    faces = []
    fanned = 0
    for fkey in mesh.faces():
        fv = [index[v] for v in mesh.face_vertices(fkey)]
        if len(fv) <= 4:
            faces.append(fv)
            continue
        pieces = _best_fan(fv, coords)
        if pieces is None:
            return None, 0, 0             # a face no diagonal partitions
        faces.extend(pieces)
        fanned += 1

    out = cls.from_vertices_and_faces(coords, faces)

    pole_keys = set()
    for p in poles or ():
        pole_keys.add(tuple(round(c, 6) for c in p[:2]))

    face_poles = {}
    for fkey in out.faces():
        fv = out.face_vertices(fkey)
        if len(fv) != 3:
            continue
        chosen = None
        for vkey in fv:
            if tuple(round(c, 6) for c in out.vertex_coordinates(vkey)[:2]) in pole_keys:
                chosen = vkey
                break
        if chosen is None:
            chosen = _choose_pole(fv, {v: out.vertex_coordinates(v) for v in fv})
        face_poles[fkey] = chosen
    out.attributes['face_pole'] = face_poles
    return out, fanned, len(face_poles)


def densifiable(mesh):
    """Whether ``densification`` will actually succeed on this coarse mesh.

    Three conditions, and all three are needed. Checking only "is it all quads"
    is what let a hexagon through with 44 quad patches that still killed
    ``discrete_coons_patch`` with an ``IndexError`` several frames away:

    * every face is a quad, or a triangle REGISTERED as a pseudo-quad pole --
      Coons patches take exactly four sides, and a pole supplies the fourth as
      a collapsed one. An unregistered triangle is still fatal: ``face_pole`` is
      looked up unconditionally while walking a strip, so it raises ``KeyError``
      several frames away from the cause;
    * no face is inverted in XY -- a self-overlapping patch densifies into a
      self-overlapping slab, silently;
    * each face's OPPOSITE edges share a strip. Densification gives a face's four
      sides the point counts of their strips, so if opposite sides sit in
      different strips they get different counts and the Coons interpolation
      indexes off the end of the shorter one. A pseudo-quad's two rails are the
      pair that has to match; its third side is opposite the collapsed corner
      and answers to nothing.

    Returns
    -------
    (bool, str)
        Verdict and, when false, the first reason.
    """
    faces = list(mesh.faces())
    if not faces:
        return False, 'no faces'

    if not mesh.is_manifold():
        return False, 'coarse mesh is not manifold'

    face_pole = mesh.attributes.get('face_pole') or {}
    for fkey in faces:
        fv = mesh.face_vertices(fkey)
        if len(fv) == 3 and fkey in face_pole:
            if face_pole[fkey] not in fv:
                return False, 'face {} has a pole that is not one of its corners'.format(fkey)
        elif len(fv) != 4:
            return False, 'face {} has {} sides'.format(fkey, len(fv))
        pts = [mesh.vertex_coordinates(v) for v in fv]
        area2 = sum(p[0] * q[1] - q[0] * p[1] for p, q in zip(pts, pts[1:] + pts[:1]))
        if area2 <= 0.0:
            return False, 'face {} is inverted or degenerate'.format(fkey)

    try:
        mesh.collect_strips()
        edge_strip = {}
        for skey, edges in mesh.strips(data=True):
            for u, v in edges:
                edge_strip[u, v] = skey
                edge_strip[v, u] = skey
        for fkey in faces:
            he = list(mesh.face_halfedges(fkey))
            if len(he) == 3:
                pole = face_pole[fkey]
                rails = [e for e in he if pole in e]
                if edge_strip[rails[0]] != edge_strip[rails[1]]:
                    return False, 'pseudo-quad {} has its two rails in different ' \
                                  'strips'.format(fkey)
                continue
            s = [edge_strip[e] for e in he]
            if s[0] != s[2] or s[1] != s[3]:
                return False, 'face {} has opposite edges in different strips'.format(fkey)
    except Exception as exc:
        return False, 'strip collection failed: {}'.format(type(exc).__name__)

    return True, ''


def _unify(mesh, cls):
    """Make every face wind counter-clockwise in XY. Returns how many flipped.

    ``Mesh.unify_cycles`` makes the winding CONSISTENT but picks its seed
    arbitrarily, so it can settle on all-clockwise. For a planar layout the
    correct answer is known outright -- positive signed area -- so the faces are
    simply checked and reversed individually, which needs no seed and cannot get
    the global sense wrong.
    """
    index = {vkey: i for i, vkey in enumerate(mesh.vertices())}
    vertices = [mesh.vertex_coordinates(v) for v in mesh.vertices()]
    faces = [[index[v] for v in mesh.face_vertices(f)] for f in mesh.faces()]

    def signed(fv):
        pts = [vertices[v] for v in fv]
        return sum(p[0] * q[1] - q[0] * p[1] for p, q in zip(pts, pts[1:] + pts[:1]))

    if all(signed(fv) > 0 for fv in faces):
        return mesh, 0

    # Reversing individual negative-area faces is WRONG, even though it makes
    # every signed area positive. Orientation is a topological property shared
    # across edges: flipping one face and not its neighbour leaves them
    # traversing their common edge the same way, which is exactly what
    # non-manifold means. Doing that turned a manifold backstop mesh
    # non-manifold while every face reported positive area.
    #
    # So unify TOPOLOGICALLY first -- that makes the winding consistent but
    # picks its seed arbitrarily, so it can settle on all-clockwise -- then flip
    # the WHOLE mesh if the total area came out negative.
    # Written out rather than calling ``Mesh.unify_cycles`` -- that did not
    # return in ten minutes on a few-thousand-face backstop mesh. This is a
    # single breadth-first pass over face adjacency.
    edge_faces = {}
    for i, fv in enumerate(faces):
        for u, v in pairwise(fv + fv[:1]):
            edge_faces.setdefault((min(u, v), max(u, v)), []).append(i)

    seen = [False] * len(faces)
    for seed in range(len(faces)):
        if seen[seed]:
            continue
        seen[seed] = True
        stack = [seed]
        while stack:
            i = stack.pop()
            fv = faces[i]
            directed = set(pairwise(fv + fv[:1]))
            for u, v in directed:
                for j in edge_faces.get((min(u, v), max(u, v)), ()):
                    if j == i or seen[j]:
                        continue
                    seen[j] = True
                    # neighbours must traverse the shared edge the OTHER way
                    if (u, v) in set(pairwise(faces[j] + faces[j][:1])):
                        faces[j] = list(reversed(faces[j]))
                    stack.append(j)

    if sum(signed(fv) for fv in faces) < 0:
        faces = [list(reversed(fv)) for fv in faces]

    return cls.from_vertices_and_faces(vertices, faces), sum(
        1 for fv in faces if signed(fv) <= 0)


def _rebuild(mesh, cls):
    """Re-index a mesh after face surgery, dropping orphaned vertices."""
    keep = {v for f in mesh.faces() for v in mesh.face_vertices(f)}
    index = {vkey: i for i, vkey in enumerate(sorted(keep))}
    vertices = [mesh.vertex_coordinates(v) for v in sorted(keep)]
    faces = [[index[v] for v in mesh.face_vertices(f)] for f in mesh.faces()]
    return cls.from_vertices_and_faces(vertices, faces)


def _arc_lengths(loop):
    """Per-segment lengths and cumulative arc length of a closed loop."""
    ring = list(loop) + list(loop[:1])
    seg = [distance_point_point(a, b) for a, b in pairwise(ring)]
    cum = [0.0]
    for s in seg:
        cum.append(cum[-1] + s)
    return seg, cum


#: A turn this sharp is a CORNER, whatever its neighbours are doing. See
#: :func:`boundary_corners` for why the threshold is an angle and not a length,
#: and ``15_baseline.py --corners`` for the measurement that keeps it honest.
SHARP_TURN = pi / 4.0


def boundary_corners(loop, limit=pi / 12.0, spacing=None, sharp=SHARP_TURN):
    """Indices of the loop's corners: places where the TANGENT JUMPS.

    The per-vertex turning angle alone is not that, and the difference is worth
    a lot of coarse patches. A curved wall arrives here already discretised, and
    every vertex of a discretised arc turns -- the stadium's semicircle is
    sampled every 18 degrees, so nine of its vertices turn 18 degrees each, all
    of them past any useful ``limit``. Taking those as corners split the wall
    nine times where the domain has no corner at all: 11 corners for a shape
    with two, 17 wall arcs, three sliver triangles against the round end, and a
    "coarse" layout of 73 patches.

    TWO REGIMES, because the ambiguity is real but it is BOUNDED
    ------------------------------------------------------------

    Above ``sharp`` there is no ambiguity at all. A 90 degree tangent jump is not
    a sample of anything: to sample a curve that coarsely you would need four
    points per full turn, and a shape with four points per turn is a square, not
    a circle. So a turn above ``sharp`` is a corner unconditionally and is never
    grouped with anything.

    Below ``sharp`` the ambiguity is genuine and no purely local test resolves
    it, because after ``_densify_loop`` -- which inserts its points ON the chord
    -- a sampled arc really IS a polygon with 18-degree corners. The same
    vertices, the same angles. What differs is SPACING: a corner of the domain
    stands alone, while a sampled curve turns again a chord later, and again,
    all the way round. So those candidates are grouped into runs by arc-length
    gap, a run of several is read as one curve, and a run still yields a corner
    where the turn stands out from the rest of the run by more than ``limit`` --
    an arc meeting a straight edge at a genuine angle keeps its corner, since
    that vertex carries the arc's own turn PLUS the jump.

    WHY ``sharp`` HAD TO EXIST -- the comb, and the class of bug it belongs to
    -------------------------------------------------------------------------

    Run-grouping alone was the stadium fix, and applied alone it BROKE the comb
    plate. ``spacing`` defaults to 1/25 of the loop; the comb's perimeter is 74,
    so the gap is 2.96, and its teeth are built from 2-unit segments. Those fall
    under the gap, so the plate's right angles grouped into runs and were read as
    a sampled curve -- and because every turn in such a run is EXACTLY 90
    degrees, none of them exceeds ``typical + limit`` and the run yielded no
    corner at all. 2 corners for a 16-corner plate, and the domain then raised
    outright for want of anything to bound a face with. The rectangle-with-slot
    lost the two reflex corners at the bottom of its slot the same way -- its
    perimeter is 50, so ``spacing`` is exactly 2.0 and the slot's bottom segment
    is exactly 2.0 -- and produced 2 coarse patches with a 0.00/180.00 degree
    face.

    That is the fourth instance of one species of bug in this module: **a
    tolerance scaled to a GLOBAL quantity, applied to a domain whose LOCAL
    feature is smaller than that scale.** ``_cluster``'s landing tolerance, the
    duplicate test's ``0.35 * span`` and this function's ``spacing`` were all the
    same mistake. Tuning ``spacing`` down is not a fix, it is the next instance:
    the stadium's arc chords are 0.94 long and the square plate's hole has
    2-unit sides, so a gap small enough to spare the comb re-breaks the domain
    the constant was introduced for.

    ``sharp`` is not vulnerable to it, and that is the point of choosing an
    ANGLE. An angle has no length scale, so no small feature can be smaller than
    it, and ``_densify_loop`` only ever subdivides -- it never makes a turn
    coarser -- so a domain that passes at one background passes at all of them.

    Measured across every outline in the suite, the two populations do not come
    close to touching::

        sampled curves   disc 7.5, round hole 7.5, ellipse <=13.9, stadium 18.0
        real corners     hexagon 60.0, pentagon 66.8-81.2, rectilinear 90.0

    ``sharp`` at 45 degrees sits 2.5x above the coarsest sample and 1.33x below
    the shallowest true corner. ``15_baseline.py --corners`` asserts that gap on
    every run, so it stops being latent.

    WHAT BREAKS NEXT, and why it is a different and smaller risk
    ------------------------------------------------------------

    An outline sampled coarser than ``sharp`` -- fewer than 8 points per full
    turn -- has its samples read as corners. That is now a statement about
    SAMPLING RESOLUTION rather than about domain size, which is the whole
    improvement: it cannot be triggered by a small feature on a large domain,
    which is what the previous three instances all were. A domain with genuine
    corners shallower than ``sharp`` (interior angle above 135 degrees) falls
    through to the run-grouping and behaves exactly as it does today, so that
    band is unchanged rather than newly at risk.

    Parameters
    ----------
    loop : list[[x, y, z]]
    limit : float, optional
        Turn above which a vertex is a corner candidate, and the excess over
        its run above which a candidate inside a run survives.
    spacing : float, optional
        Gap below which two turning vertices belong to one curve. Defaults to
        1/25 of the loop. Applies ONLY to candidates below ``sharp``; see above
        for why that restriction is what makes the default safe.
    sharp : float, optional
        Turn above which a vertex is a corner regardless of its neighbours.
        Defaults to :data:`SHARP_TURN`, 45 degrees.

    Returns
    -------
    list[int]
    """
    n = len(loop)
    seg, cum = _arc_lengths(loop)
    total = cum[-1]
    if total <= 0.0 or n < 3:
        return []
    gap = spacing if spacing is not None else total / 25.0

    # ``turn`` holds the CANDIDATES; ``every_turn`` holds all of them. The lone-run
    # rule below needs a candidate's neighbours even when those neighbours turn
    # too little to be candidates themselves -- which is the whole point, since
    # a neighbour turning just under ``limit`` is evidence of a sampled curve,
    # not evidence of nothing.
    turn = {}
    every_turn = {}
    for i in range(n):
        if seg[(i - 1) % n] < 1e-12 or seg[i] < 1e-12:
            continue
        b = loop[i]
        angle = angle_vectors(subtract_vectors(b, loop[(i - 1) % n]),
                              subtract_vectors(loop[(i + 1) % n], b))
        every_turn[i] = angle
        if angle > limit:
            turn[i] = angle

    # -- regime 1: unambiguous corners, never grouped ------------------------
    # Held out of the runs entirely, not merely added afterwards. Leaving them
    # in would let a 90 degree corner set the median of a run of arc samples,
    # and then no member of that run -- corner included -- clears
    # ``typical + limit``.
    # Compared with a relative epsilon, not exactly. A REGULAR polygon turns by
    # the same angle at every vertex, so when that angle lands on ``sharp`` the
    # float noise of ``angle_vectors`` decides each vertex separately: measured
    # on a regular octagon, whose 45.000000 degree turn is exactly SHARP_TURN, 6
    # of its 8 vertices came back as corners. Either reading of an octagon is
    # defensible -- see this function's note on what a 45 degree turn cannot be
    # told apart from -- but SIX of eight is not a reading at all.
    threshold = sharp * (1.0 - 1e-12)
    out = [i for i, a in turn.items() if a >= threshold]

    candidates = sorted(i for i, a in turn.items() if a < threshold)
    if not candidates:
        return sorted(out)
    if len(candidates) == 1:
        return sorted(out + candidates)

    # -- regime 2: ambiguous turns, grouped into runs, cyclically ------------
    # Two candidates join one run when they are within ``gap`` along the loop OR
    # when they are NEIGHBOURING VERTICES. The second test is not a tolerance and
    # has no length in it, which is the point: ``gap`` defaults to a fraction of
    # the loop, so on a loop carrying fewer than 25 points every chord is longer
    # than the gap, every turning vertex becomes a run of one, and a run of one
    # is returned as a corner below. A uniformly sampled circle then comes back
    # with a corner at EVERY vertex -- measured as a hard cliff on regular
    # n-gons: n <= 23 gives n corners, n >= 25 gives none.
    #
    # That is the fifth instance of this module's recurring species -- a
    # tolerance scaled to a GLOBAL quantity applied to a smaller LOCAL feature --
    # and here the feature is the loop's own SAMPLING. It reached the field front
    # end through ``rhino/helpers.curve_points``, which divides a curved Rhino
    # boundary by ``triangulation_spacing``: any curved loop shorter than ~24.5x
    # that spacing arrives under 25 points, so at the default 0.5 a round hole
    # needed radius 1.95 to be read as round at all. Below it, ``hole_launches``
    # sees a hole that is "already cornered", declines to cut it, and the layout
    # is discarded for a triangulation fallback.
    #
    # Two vertices that are NEIGHBOURS cannot be isolated from each other,
    # whatever the loop's length is, so this only ever merges what the gap test
    # would have merged on a better-sampled loop. It cannot join two distant
    # corners: that still needs ``gap``.
    def _same_run(prev, i):
        return cum[i] - cum[prev] <= gap or i == prev + 1

    runs = [[candidates[0]]]
    for prev, i in pairwise(candidates):
        if _same_run(prev, i):
            runs[-1].append(i)
        else:
            runs.append([i])
    wraps = ((total - cum[candidates[-1]]) + cum[candidates[0]] <= gap
             or (candidates[-1] == n - 1 and candidates[0] == 0))
    if len(runs) > 1 and wraps:
        runs[0] = runs.pop() + runs[0]

    for run in runs:
        if len(run) == 1:
            # A run of one is the only place this function emits a corner
            # without comparing the turn to ANYTHING, and that is how a
            # uniformly sampled loop leaks corners even after the adjacency
            # rule above: the candidate test is ``angle > limit``, so on a loop
            # whose samples turn by almost exactly ``limit`` only the handful of
            # vertices that clear it by a float epsilon become candidates, and
            # each of those is then isolated. Measured on a regular 24-gon,
            # whose 15.000000 degree turn sits exactly on the default limit: 3
            # vertices of 24 came back as corners while 23-gons and 25-gons came
            # back flat.
            #
            # So ask the same question the run does, against the two loop
            # NEIGHBOURS rather than against other candidates -- close to the
            # test ``SkeletonDecomposition`` uses for a boundary kink
            # (``relative_kink_angle_limit``), which is why the skeleton front
            # end never had this failure.
            #
            # BOTH neighbours must match, not their mean. A mean is meaningless
            # exactly where it matters most: at a corner split across two
            # vertices by the discretisation, one neighbour carries most of the
            # corner and the other is flat, and averaging 61.4 with 4.0 reads
            # like a 32.7 degree curve sample that the 42.6 degree candidate
            # barely stands out from. Measured on the ``decomposition`` and
            # ``RV`` domains, the mean dropped 10 baseline rows' worth of real
            # corners that way. Requiring both neighbours to turn within
            # ``limit`` of the candidate is the actual evidence for a sampled
            # curve -- my neighbours turn like me -- and a real corner, whose
            # neighbours are either flat or another corner, never provides it.
            i = run[0]
            here = turn[i]
            before = every_turn.get((i - 1) % n, 0.0)
            after = every_turn.get((i + 1) % n, 0.0)
            sampled = abs(here - before) <= limit and abs(here - after) <= limit
            if not sampled:
                out.append(i)
            continue
        angles = sorted(turn[i] for i in run)
        typical = angles[len(angles) // 2]
        out.extend(i for i in run if turn[i] - typical > limit)
    return sorted(out)


def _cluster(points, tol, loops=None, wall_tol=None, groups=None, symmetry=None):
    """Group points within ``tol`` and return ``index -> canonical point``.

    The canonical point is one MEMBER of the cluster, not the centroid: members
    are snapped to a value that already exists, so a landing point stays exactly
    on the boundary segment it was clipped to instead of drifting off it.

    Points sharing a ``groups`` entry are never merged with each other. The two
    ends of ONE trace are such a pair: this function is here to make ends that
    mean to meet actually meet, and a trace's own two ends never do. Without
    that, a short arm is deleted outright -- on the ellipse at target_length
    0.6 the left singularity sits 0.44 from the wall, its arm to the wall is
    0.44 long against a 0.48 tolerance, and collapsing it left a valence-3
    singularity with two arms and a pair of self-overlapping pentagons wrapped
    round the tip.

    Two points that BOTH sit on a wall are held to the tighter ``wall_tol``.
    ``tol`` is a numerical tolerance -- it exists so that arms which mean to
    meet actually do -- but between two landings it decides the layout instead,
    and merging two landings deletes a patch. Measured on the hexagon: at
    target_length 0.4 the two singularities land on the top wall 0.75 apart and
    the patch between them is a clean quad; at 0.5 they land 0.35 apart, inside
    a 0.40 tolerance, and that quad collapses to a triangle -- one stray
    triangle that then took the whole layout from 11 patches to 43. A landing is
    pinned to the wall and cannot drift; nothing numerical is being repaired by
    merging two of them.

    ``symmetry``, when given, projects the RESULT onto the group -- the loop
    below is untouched and still decides everything. It has to be projected
    because first-fit reads the list in trace order, which no symmetry respects:
    two points that are exact mirror images can join different clusters, and a
    cluster's representative is whichever member happened to come first. See
    ``symmetry.project_clusters``, which explains why the projected
    representative is still exactly on the wall when the original was.
    """
    wall_tol = tol * 0.25 if wall_tol is None else wall_tol

    def on_wall(p):
        return any(_distance_to_loop(p, loop) < tol * 0.25 for loop in (loops or []))

    flags = [on_wall(p) for p in points]

    canonical = {}
    reps = []
    members = {}
    for i, p in enumerate(points):
        for j, rep in reps:
            # No CLUSTER may hold both ends of one trace, not merely no pair.
            # Checking the pair is not enough: the arm's landing does not merge
            # with the arm's own start, it merges with a SIBLING arm's copy of
            # the same singularity centroid, and lands in the same cluster by
            # the back door.
            if groups is not None and groups[i] in members[j]:
                continue
            limit = wall_tol if (flags[i] and flags[j]) else tol
            if distance_point_point(p, rep) < limit:
                canonical[i] = j
                if groups is not None:
                    members[j].add(groups[i])
                break
        else:
            reps.append((i, p))
            canonical[i] = i
            members[i] = {groups[i]} if groups is not None else set()

    out = {i: points[j] for i, j in canonical.items()}
    if symmetry is not None and symmetry.enabled('network'):
        from .symmetry import project_clusters

        out = project_clusters(points, out, symmetry)
    return out


def _decimate(points, spacing):
    """Thin a traced polyline, keeping both ends.

    Every intermediate point becomes a vertex of the planar network inside
    ``from_polylines``, so an undecimated trace makes the face search needlessly
    slow without changing the result.
    """
    if len(points) < 3:
        return [list(p) for p in points]
    out = [list(points[0])]
    for p in points[1:-1]:
        if distance_point_point(p, out[-1]) >= spacing:
            out.append(list(p))
    if distance_point_point(points[-1], out[-1]) < spacing * 0.5 and len(out) > 1:
        out.pop()
    out.append(list(points[-1]))
    return out


def _resample(points, n):
    """``n`` points spread evenly by ARC LENGTH along a polyline, ends kept."""
    cum = [0.0]
    for a, b in pairwise(points):
        cum.append(cum[-1] + distance_point_point(a, b))
    total = cum[-1]
    if total == 0.0 or n < 2:
        return [list(points[0]), list(points[-1])]

    out = []
    for i in range(n):
        s = total * i / float(n - 1)
        for j, (c0, c1) in enumerate(zip(cum, cum[1:])):
            if c0 - 1e-12 <= s <= c1 + 1e-12:
                span = c1 - c0
                t = 0.0 if span == 0 else (s - c0) / span
                a, b = points[j], points[j + 1]
                out.append([a[k] + (b[k] - a[k]) * t for k in range(3)])
                break
        else:
            out.append(list(points[-1]))
    out[0], out[-1] = list(points[0]), list(points[-1])
    return out


def _length_of(points):
    """Total arc length of a polyline."""
    return sum(distance_point_point(a, b) for a, b in pairwise(points))


def _launch_direction(points, distance):
    """The direction a polyline leaves its first point in.

    Measured to the point ``distance`` along, not to the next vertex: a
    separatrix is noisiest exactly where it starts, because that is where the
    field's magnitude goes to zero, and the first decimated segment is mostly
    that noise.
    """
    acc = 0.0
    for a, b in pairwise(points):
        acc += distance_point_point(a, b)
        if acc >= distance:
            return subtract_vectors(b, points[0])
    return subtract_vectors(points[-1], points[0])


def _arc_midpoint_of(points):
    """The point half way along a polyline BY LENGTH."""
    return _resample(points, 3)[1]


def _average_polylines(a, b):
    """Mean of two polylines that share their endpoints.

    Both are resampled onto a common arc-length parameter first, since they are
    decimated independently and generally have different point counts.
    """
    n = max(len(a), len(b), 3)
    ra, rb = _resample(a, n), _resample(b, n)
    out = [[(p[k] + q[k]) / 2.0 for k in range(3)] for p, q in zip(ra, rb)]
    out[0], out[-1] = list(a[0]), list(a[-1])
    return out


def _param(loop, point):
    """(segment index, t) of the loop point closest to ``point``."""
    best, best_key = float('inf'), (0, 0.0)
    for i, (a, b) in enumerate(pairwise(loop + loop[:1])):
        ab = subtract_vectors(b, a)
        length2 = ab[0] * ab[0] + ab[1] * ab[1]
        if length2 == 0.0:
            continue
        t = ((point[0] - a[0]) * ab[0] + (point[1] - a[1]) * ab[1]) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + ab[0] * t, a[1] + ab[1] * t, 0.0]
        d = distance_point_point(point, q)
        if d < best:
            best, best_key = d, (i, t)
    return best_key


def _split_loop(loop, cuts, corner_limit):
    """Split a closed loop into arcs at ``cuts`` and at its own corners.

    Returns
    -------
    (list[polyline], list[[x, y, z]])
        The arcs, and the split points actually used.
    """
    marks = []
    for i in boundary_corners(loop, corner_limit):
        marks.append(((i, 0.0), list(loop[i])))
    for p in cuts:
        marks.append((_param(loop, p), list(p)))

    if len(marks) < 2:
        # a loop with fewer than two splits collapses to a degenerate patch;
        # SkeletonDecomposition.branches_splitting_collapsed_boundaries hits the
        # same case. Add evenly spaced marks so the loop still bounds something.
        n = len(loop)
        for k in range(3):
            i = (k * n) // 3
            marks.append(((i, 0.0), list(loop[i])))

    marks.sort(key=lambda m: m[0])

    # drop duplicates -- a separatrix landing exactly on a corner
    unique = [marks[0]]
    for m in marks[1:]:
        if distance_point_point(m[1], unique[-1][1]) > 1e-9:
            unique.append(m)
    if len(unique) > 1 and distance_point_point(unique[0][1], unique[-1][1]) < 1e-9:
        unique.pop()

    arcs = []
    used = [m[1] for m in unique]
    count = len(unique)
    for k in range(count):
        (i0, t0), p0 = unique[k]
        (i1, t1), p1 = unique[(k + 1) % count]
        arc = [list(p0)]
        j = i0
        # walk the loop's own vertices between the two marks
        guard = 0
        while guard <= len(loop):
            guard += 1
            j = (j + 1) % len(loop)
            if (j == i1 and t1 > 0.0) or (j == (i1 + 1) % len(loop) and t1 == 0.0):
                break
            arc.append(list(loop[j]))
            if j == i1:
                break
        if distance_point_point(arc[-1], p1) > 1e-9:
            arc.append(list(p1))
        else:
            arc[-1] = list(p1)
        if len(arc) >= 2:
            arcs.append(arc)
    return arcs, used


def build_network(field, separatrices, tol=None, corner_limit=pi / 12.0, spacing=None,
                  singularity_points=None, symmetry=None):
    """Snap, split and partition traced separatrices for ``from_polylines``.

    Parameters
    ----------
    field : CrossField
    separatrices : list[Separatrix]
    tol : float, optional
        Snapping distance. Defaults to 0.8 x the background target length --
        comfortably above ``TOL.geometric_key``'s precision, which is what
        ``from_polylines`` matches endpoints with.
    corner_limit : float, optional
        Turn angle above which a boundary vertex counts as a corner.
    spacing : float, optional
        Decimation spacing for traced polylines.
    singularity_points : dict[int, [x, y, z]], optional
        Where each singularity is, from ``Tracer._singularity_points``. Face
        centroids when omitted, which is what this function used to compute for
        itself -- and recomputing them is not harmless: it DISCARDS a position
        the tracer chose. On a symmetric domain that one line was the largest
        single source of asymmetry in the finished layout, taking traces that go
        in exactly symmetric (all 472 points of all 32 separatrices, deviation
        0.0000) and putting their front points back at 69%, deviation 0.3486.
    symmetry : :class:`framefield.symmetry.Symmetry`, optional
        When given and non-trivial, the end-clustering below is projected onto
        the group. See :func:`_cluster`.

    Returns
    -------
    (list[polyline], list[polyline], dict)
        Boundary polylines, other polylines, and a report.
    """
    background = field.background
    tol = tol if tol is not None else background.target_length * 0.8
    spacing = spacing if spacing is not None else background.target_length

    mesh = background.mesh
    if singularity_points:
        centroids = dict(singularity_points)
    else:
        centroids = {fkey: mesh.face_centroid(fkey) for fkey, _ in field.singularities()}

    # -- endpoints ----------------------------------------------------------
    # A separatrix starts half a step off its singularity, so put the centroid
    # back on the front: every arm of one singularity must share one exact point
    # or the patches around it will not close.
    corner_points = []
    for loop in [background.outer] + list(background.inners):
        for i in boundary_corners(loop, corner_limit):
            corner_points.append(list(loop[i]))

    traces = []
    for s in separatrices:
        pts = [list(p) for p in s.points]
        if s.source in centroids:
            pts.insert(0, list(centroids[s.source]))
        elif s.source is None and corner_points:
            # a corner launch starts one inset inside the wall; put it back on
            # the corner, which is already a split point of the boundary
            nearest = min(corner_points, key=lambda c: distance_point_point(c, pts[0]))
            pts.insert(0, list(nearest))
        if s.sink in centroids:
            pts.append(list(centroids[s.sink]))
        traces.append((s, _decimate(pts, spacing)))

    ends = []
    for _, pts in traces:
        ends.append(pts[0])
        ends.append(pts[-1])
    snapped = _cluster(ends, tol, loops=[background.outer] + list(background.inners),
                       groups=[k // 2 for k in range(len(ends))],
                       symmetry=symmetry)

    others = []
    landings = []
    for k, (s, pts) in enumerate(traces):
        pts[0] = list(snapped[2 * k])
        pts[-1] = list(snapped[2 * k + 1])
        if distance_point_point(pts[0], pts[-1]) < 1e-9:
            continue                        # collapsed to a point
        if s.reason == 'boundary':
            # snap a landing that arrived near a corner onto it exactly. The
            # corner is already a split mark of the boundary, so a landing a few
            # microns away would add a second mark beside it and cut off a sliver
            # arc that no patch can use.
            if corner_points:
                nearest = min(corner_points, key=lambda c: distance_point_point(c, pts[-1]))
                if distance_point_point(nearest, pts[-1]) < tol:
                    pts[-1] = list(nearest)
            landings.append(pts[-1])
            if s.source == 'hole':
                # A hole cut lands on one wall but STARTS on another -- the
                # smooth inner loop it was launched from. Only the landing was
                # ever recorded here, so that loop reached ``_split_loop`` with
                # NO cuts of its own and took the "fewer than two marks" branch,
                # which spreads three evenly spaced marks around it. Those do not
                # line up with the cuts, so each one sliced a sector in two: the
                # annulus came out as 4 sectors + 3 stray marks = 7 patches, 3 of
                # them triangles, instead of the 4 quads its 4 radial cuts define.
                landings.append(pts[0])
        others.append(pts)

    # -- drop duplicates ----------------------------------------------------
    # Two reflex corners facing each other launch the SAME line from both ends,
    # so it arrives twice, reversed. from_polylines would then see a doubled
    # edge and the planar face search cannot resolve which side is which.
    unique = []
    for pts in others:
        a, b = pts[0], pts[-1]
        match = None
        for k, kept in enumerate(unique):
            ka, kb = kept[0], kept[-1]
            forward = (distance_point_point(a, ka) < tol
                       and distance_point_point(b, kb) < tol)
            backward = (distance_point_point(a, kb) < tol
                        and distance_point_point(b, ka) < tol)
            if not (forward or backward):
                continue
            other = pts if forward else list(reversed(pts))
            # Compare ARC-LENGTH midpoints, not ``pts[len(pts)//2]``. The two
            # traces are decimated independently and run in opposite directions,
            # so the same index is a different position along the curve: on a
            # disc the index-midpoints of four duplicate pairs sat 0.55-1.10
            # apart against a 0.48 tolerance, and not one pair was recognised.
            # Every ring edge was then drawn and meshed twice.
            # Tolerance scaled to the connection's own LENGTH, not to the
            # background spacing. Two traces of one connection bow apart by some
            # fraction of how far they run -- both are noisiest where they leave
            # their singularity -- so a fixed tolerance recognises the short
            # pairs and misses the long ones. On a disc that left one of the
            # four ring edges duplicated while the other three merged.
            span = 0.5 * (_length_of(other) + _length_of(kept))
            if distance_point_point(_arc_midpoint_of(other),
                                    _arc_midpoint_of(kept)) >= max(tol, 0.35 * span):
                continue
            # Sharing both endpoints and passing near the same midpoint is NOT
            # enough. Two DIFFERENT arms of one singularity can run to the same
            # place and, over a long enough span, ``0.35 * span`` is wide enough
            # to swallow the gap between them. Measured on the ellipse: the two
            # arms leaving the left tip both reach the right tip, 14.6 apart
            # end to end, bowing 4.7 apart in the middle against a threshold of
            # 5.1 -- so all four traces of the pair collapsed to one line, the
            # tip singularity was left with two arms instead of three, and the
            # patches either side of it came out as pentagons.
            #
            # A singularity's arms leave it in DIFFERENT directions -- 120
            # degrees apart at index +1 -- while two traces of one connection
            # leave in the same one. That is the discriminating measurement, and
            # it is not close: across the disc at four resolutions, the hexagon,
            # the ellipse and the stadium, every genuine duplicate agrees to
            # within 22 degrees and every false one differs by more than 100.
            if angle_vectors(_launch_direction(other, tol * 2.5),
                             _launch_direction(kept, tol * 2.5)) > pi / 4.0:
                continue
            match = (k, other)
            break
        if match is None:
            unique.append(pts)
        else:
            # Neither trace is more correct than the other -- both wobble where
            # |u| -> 0 near the singularity they left. Averaging the two along a
            # common parameter is both smoother and symmetric, which is what a
            # connection between two singularities should be.
            k, other = match
            unique[k] = _average_polylines(unique[k], other)
    dropped = len(others) - len(unique)
    others = unique

    # -- boundary arcs ------------------------------------------------------
    boundary = []
    for loop in [background.outer] + list(background.inners):
        cuts = [p for p in landings if _distance_to_loop(p, loop) < tol]
        arcs, _ = _split_loop(loop, cuts, corner_limit)
        boundary.extend(arcs)

    report = {
        'separatrices': len(others),
        'boundary_arcs': len(boundary),
        'landings': len(landings),
        'duplicates_dropped': dropped,
        'snapped': len(ends) - len({id(v) for v in snapped.values()}),
    }
    return boundary, others, report


def _distance_to_loop(p, loop):
    best = float('inf')
    for a, b in pairwise(loop + loop[:1]):
        ab = subtract_vectors(b, a)
        length2 = ab[0] * ab[0] + ab[1] * ab[1]
        if length2 == 0.0:
            continue
        t = ((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1]) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + ab[0] * t, a[1] + ab[1] * t, 0.0]
        best = min(best, distance_point_point(p, q))
    return best
