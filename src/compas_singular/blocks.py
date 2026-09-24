"""**The singularity as a BLOCK, placed by a POINT.**

A point feature in this library becomes a POLE: a vertex where the mesh lines
gather, index 1, its incident faces fanned into pseudo-quads. A block system
cannot build that -- blocks meet along faces, never in a point.

This module is the other half of the pair. A **block point** is a picked point
that becomes the CENTRE of an n-gon FACE with quads all round it: the same
singularity, carried by an element instead of a joint.

THE ONE DIFFERENCE BETWEEN A POLE POINT AND A BLOCK POINT
----------------------------------------------------------

Whether the coarse vertex gets fanned. Pass the point to
:func:`~compas_singular.algorithms.boundary_triangulation` so the layout wraps
patches around it and it lands as a coarse vertex -- but WITHHOLD it from
``decomposition_mesh(poles)``. With no fan it survives densification as a
genuine valence-n joint with an all-quad 1-ring, at the picked location, and
:func:`block_points` truncates it into the block::

    point_features = [...]                       # poles, as before
    block_features = [P]                         # the new kind

    trimesh = boundary_triangulation(outer, inners, lines,
                                     point_features + block_features)
    coarse  = SkeletonDecomposition.from_mesh(trimesh).decomposition_mesh(
        point_features)                          # P withheld -- no fan
    dense   = densify(coarse)
    dense, report = block_points(dense, block_features)

Nothing in ``decomposition.py`` changes; the whole of it is which list the point
goes in.

WHAT A POINT CAN AND CANNOT ASK FOR
------------------------------------

The block's SIZE is free -- ``ratio``, how far along the incident edges its
corners sit. Its DEGREE is not: :func:`index_sum` pins it to the valence of the
joint it replaces, so a valence-3 joint gives a triangular block against 3
quads and a valence-5 joint a pentagon against 5. Asking for a hexagon leaves
+1 of index to park somewhere, i.e. a new valence-3 joint -- a point junction
again, just moved.

For the same reason a POLE can never become one block: a full pole carries raw
index 4, and no single face carries 4 (that would be a face of degree 0). This
is why the routing above withholds the point from the pole list rather than
converting a pole after the fact.

HOW THE TRUNCATION WORKS, AND WHAT IT COSTS
--------------------------------------------

Truncate the joint: put a new corner part-way along each of its n edges and
join them into the n-gon. That leaves the n incident quads as PENTAGONS, which
is not a block layout either, so each is repaired to two quads by splitting one
of its outer edges::

        \\ | /                  \\  |  /
         \\|/                    +--+--+
      ----o----      ->       ---|  N  |---        N = the n-gon block
         /|\\                    +--+--+
        / | \\                  /   |   \\

A split edge in a quad mesh is a quad-strip split: the face on the far side is
now a pentagon too, and so is the one after it. The split runs until its strip
ends, which on a plate means it runs to the wall. So a block costs n SEAMS
radiating from it -- ``report['strips']`` counts them and :func:`seam_edges`
draws them.

That cost is not an artefact of this construction, it is a conservation law.
For any sub-disc R of a quad mesh, ``4F = 2E_int + E_bnd``, so ``E_bnd`` is
EVEN. Re-mesh R to hold one n-gon and otherwise quads and
``4(F-1) + n = 2E_int + E_bnd``, so ``E_bnd`` has the parity of n. For odd n --
and 3 and 5 are the only valences that occur -- those disagree, so NO region,
however large, can contain the block. At least one split must escape. Two
blocks together have even parity, so a PAIR can absorb each other's seam, which
is why :func:`block_points` builds every point in ONE pass through
:func:`blocks_at` rather than one at a time.

TWO CHIRALITIES, AND THE CHOICE MATTERS
----------------------------------------

Each incident quad repairs itself by chording to ONE of the two block corners
it touches, and every block corner needs exactly one such chord. That is a
perfect matching between the n faces and n corners of a cycle, and a cycle has
exactly TWO -- the pinwheel spins one way or the other, and mixing them cannot
close. ``spin=0`` and ``spin=1`` are those two. They send their seams off in
different directions, so one can succeed where the other collides;
:func:`blockable` tries both and :func:`blocks_any_spin` searches assignments.

WHEN IT CANNOT BE DONE
----------------------

If a repair strip is CLOSED, the split has no end to run to: it comes back
round and lands on an edge of the joint itself, which the truncation cannot
survive (the incident face would need a 7-gon repair, and 7 is odd). If it runs
into a POLE, it stops at a face it cannot split. Call :func:`head_strips` to
see either before building -- it reports ``open`` per incident face, or which
of the three ways it failed. The 4 valence-3 joints of a disc each have one
returning strip in both spins at every density, so a disc takes no local blocks
at all; dualise it instead.

Note that a joint's own 1-ring must be all quads, so a POLE cannot be blocked
-- but a pole elsewhere in the mesh is fine, as long as no seam reaches it.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any

from compas.datastructures import Mesh
from compas.geometry import angle_points
from compas.geometry import centroid_points
from compas.geometry import distance_point_point


__all__ = [
    'Propagation', 'joints', 'joint_at', 'head_strips', 'blockable',
    'block_at', 'blocks_at', 'blocks_any_spin', 'block_points',
    'pole_at', 'pole_blocks',
    'relax', 'seam_edges', 'report_lines',
    'index_sum', 'face_degrees', 'interior_valences', 'quality',
]

DEG = 57.29577951308232


class Propagation(Exception):
    """A repair split had nowhere to go. Always a statement about the MESH."""


# ----------------------------------------------------------------------------
# reading a mesh
# ----------------------------------------------------------------------------

def joints(mesh: Mesh) -> list[int]:
    """Interior vertices where a number of mesh lines other than 4 gather."""
    return [v for v in mesh.vertices()
            if not mesh.is_vertex_on_boundary(v) and mesh.vertex_degree(v) != 4]


def face_degrees(mesh: Mesh) -> dict[int, int]:
    d = defaultdict(int)
    for f in mesh.faces():
        d[len(mesh.face_vertices(f))] += 1
    return dict(d)


def interior_valences(mesh: Mesh) -> dict[int, int]:
    d = defaultdict(int)
    for v in mesh.vertices():
        if not mesh.is_vertex_on_boundary(v):
            d[mesh.vertex_degree(v)] += 1
    return dict(d)


def index_sum(mesh: Mesh) -> int:
    """``Sum(4 - valence)`` over interior vertices ``+ Sum(4 - degree)`` over faces.

    The conserved quantity. A valence-5 joint and a pentagonal face carry the
    same -1, which is exactly why one can be traded for the other -- and why
    neither can be made to disappear. If a block insertion changes this number,
    it did not move the singularity, it invented one.
    """
    s = sum(4 - mesh.vertex_degree(v)
            for v in mesh.vertices() if not mesh.is_vertex_on_boundary(v))
    return s + sum(4 - len(mesh.face_vertices(f)) for f in mesh.faces())


def quality(mesh: Mesh) -> tuple[float, float, float]:
    """Worst corner angle (min, max, degrees) and worst edge aspect ratio."""
    lo, hi, asp = 180.0, 0.0, 1.0
    for f in mesh.faces():
        pts = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(f)]
        n = len(pts)
        for i in range(n):
            a = angle_points(pts[i], pts[(i - 1) % n], pts[(i + 1) % n]) * DEG
            lo, hi = min(lo, a), max(hi, a)
        lengths = [distance_point_point(pts[i], pts[(i + 1) % n])
                   for i in range(n)]
        if min(lengths) > 1e-9:
            asp = max(asp, max(lengths) / min(lengths))
    return lo, hi, asp


def _opposite_edge(mesh: Mesh, fkey: int, u: int, v: int) -> tuple[int, int] | None:
    """The edge of quad ``fkey`` opposite ``(u, v)``. None if not a quad."""
    vertices = mesh.face_vertices(fkey)
    if len(vertices) != 4:
        return None
    i = vertices.index(u)
    if vertices[(i + 1) % 4] != v:
        i = vertices.index(v)
    return (vertices[(i + 2) % 4], vertices[(i + 3) % 4])


def _across(mesh: Mesh, fkey: int, u: int, v: int) -> int:
    """The face on the other side of edge ``(u, v)`` from ``fkey``."""
    a, b = mesh.halfedge[u][v], mesh.halfedge[v][u]
    return b if a == fkey else a


def _head(mesh: Mesh, fkey: int, vkey: int, spin: int) -> tuple[int, int]:
    """Which outer edge of an incident quad the repair splits."""
    d = mesh.face_vertex_descendant(fkey, vkey)
    o = mesh.face_vertex_descendant(fkey, d)
    a = mesh.face_vertex_descendant(fkey, o)
    return (d, o) if spin == 0 else (o, a)


def _median_edge_length(mesh: Mesh) -> float:
    lengths = sorted(mesh.edge_length(edge) for edge in mesh.edges())
    if not lengths:
        return 0.0
    return lengths[len(lengths) // 2]


# ----------------------------------------------------------------------------
# will it work?
# ----------------------------------------------------------------------------

def head_strips(mesh: Mesh, vkey: int, spin: int = 0) -> list[tuple[str, int]]:
    """Classify the strip each repair split must run in, one per incident quad.

    Returns ``(verdict, faces)`` per incident quad. ``open`` means the split
    reaches a wall, which is what it needs. Anything else means it cannot get
    there, and then the block cannot be built at this joint in this spin --
    the three ways being that it comes back to ``vkey``, that its strip is
    closed, or that it runs into a face it cannot split (a POLE).
    """
    out = []
    for f in mesh.vertex_faces(vkey, ordered=True):
        edge = _head(mesh, f, vkey, spin)
        face, n, seen = f, 0, set()
        while True:
            g = _across(mesh, face, *edge)
            if g is None:
                out.append(('open', n))
                break
            nxt = _opposite_edge(mesh, g, *edge)
            if nxt is None:
                # ``g`` is not a quad, so the strip has no continuation there.
                # In a pseudo-quad mesh that is a pole, and the split would
                # land on a face the rebuild cannot repair.
                out.append(('runs into a non-quad face (a pole)', n))
                break
            if frozenset(nxt) in seen or n > 5000:
                out.append(('closed', n))
                break
            if vkey in nxt:
                out.append(('returns to the block', n))
                break
            seen.add(frozenset(nxt))
            face, edge, n = g, nxt, n + 1
    return out


def blockable(mesh: Mesh, vkey: int) -> tuple[int | None, str]:
    """``(spin, reason)`` -- which chirality works at this joint, if either.

    ``spin`` is None when neither does, and ``reason`` says why in one line.
    Checks the cheap disqualifications first: a joint whose own 1-ring is not
    all quads (a pole sits in it) cannot be truncated at all.
    """
    if mesh.is_vertex_on_boundary(vkey):
        return None, 'on the boundary'
    if any(len(mesh.face_vertices(f)) != 4
           for f in mesh.vertex_faces(vkey, ordered=True)):
        return None, 'a non-quad face (a pole) sits in the 1-ring'
    for spin in (0, 1):
        verdicts = head_strips(mesh, vkey, spin)
        if all(v == 'open' for v, _ in verdicts):
            return spin, 'spin={} -- all {} repair strips reach a wall'.format(
                spin, len(verdicts))
    bad = [v for v, _ in head_strips(mesh, vkey, 0) if v != 'open']
    return None, ('{} repair strip(s) cannot reach a wall in either spin -- '
                  '{}'.format(len(bad), '; '.join(sorted(set(bad)))))


def joint_at(
    mesh: Mesh, point: list[float], tol: float | None = None
) -> tuple[int | None, float | None, str]:
    """The blockable joint at ``point``: ``(vkey, distance, reason)``.

    ``vkey`` is None when there is no joint there to block, and ``reason`` says
    which of the ways it is missing -- the layout put a regular vertex at the
    point, or it put a POLE there (a 3-vertex face on it), or the nearest joint
    is further away than ``tol``.

    ``tol`` defaults to 1.5x the median edge length, i.e. "within about one
    element of where you clicked". Pass ``float('inf')`` to bind to the nearest
    joint whatever the distance.
    """
    if tol is None:
        tol = 1.5 * _median_edge_length(mesh)

    nearest, best = None, None
    nearest_any, best_any = None, None
    for v in mesh.vertices():
        d = distance_point_point(point, mesh.vertex_coordinates(v))
        if best_any is None or d < best_any:
            nearest_any, best_any = v, d
        if mesh.is_vertex_on_boundary(v) or mesh.vertex_degree(v) == 4:
            continue
        if any(len(mesh.face_vertices(f)) != 4
               for f in mesh.vertex_faces(v, ordered=True)):
            continue
        if best is None or d < best:
            nearest, best = v, d

    if nearest is not None and best <= tol:
        return nearest, best, 'valence {}'.format(mesh.vertex_degree(nearest))

    # Nothing blockable within reach -- say what IS at the point instead.
    if nearest_any is not None and best_any <= tol:
        v = nearest_any
        if any(len(mesh.face_vertices(f)) != 4
               for f in mesh.vertex_faces(v, ordered=True)):
            why = ('the layout made this point a POLE (a non-quad face sits on '
                   'it), which cannot be truncated')
        elif mesh.is_vertex_on_boundary(v):
            why = 'the nearest vertex is on the boundary'
        else:
            why = ('the layout put a regular valence-4 vertex here, so there '
                   'is no singularity to carry a block')
        return None, best_any, why

    if nearest is None:
        return None, None, 'the mesh has no blockable joint at all'
    return None, best, ('the nearest blockable joint is {:.4f} away, beyond '
                        'the {:.4f} tolerance'.format(best, tol))


# ----------------------------------------------------------------------------
# the operator
# ----------------------------------------------------------------------------

def block_at(
    mesh: Mesh,
    vkey: int,
    ratio: float = 0.4,
    spin: int = 0,
    _splits: dict[frozenset, bool] | None = None,
    _report: dict[str, Any] | None = None,
) -> tuple[Mesh | None, dict[str, Any]]:
    """A NEW mesh with the valence-n joint ``vkey`` replaced by an n-gon block.

    ``ratio`` is how far along each incident edge the block's corner sits, i.e.
    the block's size. It is the one free parameter -- the block's DEGREE is
    pinned to the valence by :func:`index_sum`, so a valence-5 joint gives a
    pentagon and nothing else.

    Returns ``(mesh, report)``; ``report['strips']`` is how many seams the block
    cost and ``report['strip_faces']`` how many faces those seams split.

    Note that ``vkey`` does not survive the call: the rebuild renumbers every
    vertex. Use :func:`block_points`, which is keyed on coordinates, to place
    several blocks or to drive this from a picked point.
    """
    if mesh.is_vertex_on_boundary(vkey):
        raise ValueError('vertex {} is on the boundary'.format(vkey))
    ring = mesh.vertex_faces(vkey, ordered=True)
    if any(len(mesh.face_vertices(f)) != 4 for f in ring):
        raise ValueError('the 1-ring of vertex {} is not all quads'.format(vkey))

    splits = {} if _splits is None else _splits
    trunc = {}
    report = ({'strips': 0, 'strip_faces': 0, 'ends': []}
              if _report is None else _report)

    for nbr in mesh.vertex_neighbors(vkey):
        trunc[(vkey, nbr)] = True

    for f in ring:
        edge = _head(mesh, f, vkey, spin)
        if frozenset(edge) not in splits:
            splits[frozenset(edge)] = True
            report['strips'] += 1
        # the split runs outward until its strip ends
        face, guard = f, 0
        while True:
            guard += 1
            if guard > 100000:
                raise Propagation('runaway strip at vertex {}'.format(vkey))
            g = _across(mesh, face, *edge)
            if g is None:
                report['ends'].append('wall')
                break
            nxt = _opposite_edge(mesh, g, *edge)
            if nxt is None:
                report['ends'].append('non-quad face')
                break
            if frozenset(nxt) in splits:
                report['ends'].append('met another seam')
                break
            splits[frozenset(nxt)] = True
            report['strip_faces'] += 1
            face, edge = g, nxt

    if _splits is not None:
        return None, report                      # caller is batching
    return _rebuild(mesh, {vkey: (ring, ratio, spin)}, splits, trunc), report


def blocks_at(
    mesh: Mesh, vkeys: list[int], ratio: float = 0.4, spins: dict[int, int] | None = None
) -> tuple[Mesh, dict[str, Any]]:
    """Every block in one pass, so their seams can end on each other."""
    splits, trunc, sing = {}, {}, {}
    spins = spins or {}
    report = {'strips': 0, 'strip_faces': 0, 'ends': []}
    for vkey in vkeys:
        ring = mesh.vertex_faces(vkey, ordered=True)
        if any(len(mesh.face_vertices(f)) != 4 for f in ring):
            raise ValueError(
                'the 1-ring of vertex {} is not all quads'.format(vkey))
        sing[vkey] = (ring, ratio, spins.get(vkey, 0))
        for nbr in mesh.vertex_neighbors(vkey):
            trunc[(vkey, nbr)] = True
    for vkey in vkeys:
        block_at(mesh, vkey, ratio, spins.get(vkey, 0),
                 _splits=splits, _report=report)
    return _rebuild(mesh, sing, splits, trunc), report


def blocks_any_spin(mesh: Mesh, vkeys: list[int], ratio: float = 0.4) -> tuple[Mesh, dict[str, Any]]:
    """Search the chirality assignments until one closes. Reports which."""
    import itertools
    n = len(vkeys)
    combinations = (itertools.product((0, 1), repeat=n) if n <= 8
                    else [(0,) * n, (1,) * n])
    last = None
    for bits in combinations:
        try:
            out, report = blocks_at(mesh, vkeys, ratio, dict(zip(vkeys, bits)))
            report['spins'] = bits
            return out, report
        except Propagation as exc:
            last = exc
    raise Propagation('no chirality assignment closes; last: {}'.format(last))


def _rebuild(
    mesh: Mesh,
    sing: dict[int, tuple[list[int], float, int]],
    splits: dict[frozenset, bool],
    trunc: dict[tuple[int, int], bool],
) -> Mesh:
    """Assemble the new mesh from the original faces and the recorded splits."""
    points, index = [], {}

    def add(tag: tuple, xyz: list[float]) -> int:
        if tag not in index:
            index[tag] = len(points)
            points.append(list(xyz))
        return index[tag]

    for v in mesh.vertices():
        if v not in sing:
            add(('v', v), mesh.vertex_coordinates(v))
    for key in splits:
        u, w = tuple(key)
        if u in sing or w in sing:
            raise Propagation(
                'a repair strip came back and split an edge AT block {} -- '
                'the block cannot be built here'.format(
                    u if u in sing else w))
        add(('s', key), mesh.edge_midpoint((u, w)))
    for (vkey, nbr) in trunc:
        ratio = sing[vkey][1]
        a = mesh.vertex_coordinates(vkey)
        b = mesh.vertex_coordinates(nbr)
        add(('t', vkey, nbr), [a[i] + ratio * (b[i] - a[i]) for i in range(3)])

    def V(v: int) -> int:
        return index[('v', v)]

    def S(u: int, w: int) -> int:
        return index[('s', frozenset((u, w)))]

    def T(vkey: int, nbr: int) -> int:
        return index[('t', vkey, nbr)]

    faces = []
    ring_faces = {f: v for v, (ring, _, _) in sing.items() for f in ring}

    for f in mesh.faces():
        vertices = mesh.face_vertices(f)
        cut = [e for e in _edges_of(vertices) if frozenset(e) in splits]

        if f in ring_faces:
            vkey = ring_faces[f]
            d = mesh.face_vertex_descendant(f, vkey)
            o = mesh.face_vertex_descendant(f, d)
            a = mesh.face_vertex_descendant(f, o)
            if len(cut) != 1:
                raise Propagation(
                    'a repair strip runs THROUGH the 1-ring of block {} '
                    '({} split edges on one incident quad) -- the block cannot '
                    'be built here'.format(vkey, len(cut)))
            # the incident quad is the pentagon (xa, xd, d, o, a); which of its
            # two outer edges is split is the chirality
            xd, xa = T(vkey, d), T(vkey, a)
            if sing[vkey][2] == 0:
                m = S(d, o)
                faces.append([xa, xd, V(d), m])
                faces.append([xa, m, V(o), V(a)])
            else:
                m = S(o, a)
                faces.append([xd, V(d), V(o), m])
                faces.append([xa, xd, m, V(a)])
            continue

        if not cut:
            faces.append([V(v) for v in vertices])
            continue
        if len(vertices) != 4:
            raise Propagation(
                'a repair strip ran into non-quad face {} (a pole) -- the '
                'block cannot be built here'.format(f))

        a, b, c, d = vertices
        marks = [frozenset(e) in splits for e in _edges_of(vertices)]
        if marks == [True, False, True, False]:
            m1, m2 = S(a, b), S(c, d)
            faces.append([V(a), m1, m2, V(d)])
            faces.append([m1, V(b), V(c), m2])
        elif marks == [False, True, False, True]:
            m1, m2 = S(b, c), S(d, a)
            faces.append([V(a), V(b), m1, m2])
            faces.append([m2, m1, V(c), V(d)])
        elif all(marks):
            # two seams crossing: four quads around a new centre. The centre is
            # a regular valence-4 joint, so a crossing costs nothing.
            mab, mbc, mcd, mda = S(a, b), S(b, c), S(c, d), S(d, a)
            ctr = add(('c', f), mesh.face_centroid(f))
            faces.append([V(a), mab, ctr, mda])
            faces.append([mab, V(b), mbc, ctr])
            faces.append([ctr, mbc, V(c), mcd])
            faces.append([mda, ctr, mcd, V(d)])
        else:
            raise Propagation(
                'face {} is split on {} edge(s) that are not opposite pairs -- '
                'a repair strip terminated in mid-mesh'.format(f, len(cut)))

    out = Mesh.from_vertices_and_faces(points, faces)
    for vkey in sing:
        out.add_face([T(vkey, p)
                      for p in mesh.vertex_neighbors(vkey, ordered=True)])
    out.unify_cycles()
    return out


def _edges_of(vertices: list[int]) -> list[tuple[int, int]]:
    n = len(vertices)
    return [(vertices[i], vertices[(i + 1) % n]) for i in range(n)]


# ----------------------------------------------------------------------------
# the point-driven entry point
# ----------------------------------------------------------------------------

def block_points(
    mesh: Mesh,
    points: list[list[float]],
    ratio: float = 0.4,
    tol: float | None = None,
    relax_iters: int = 0,
) -> tuple[Mesh, dict[str, Any]]:
    """Turn each picked point into an n-gon block face in a dense quad mesh.

    This is the entry point a picked point drives. It is keyed on COORDINATES
    rather than on vertex keys because :func:`block_at` renumbers every vertex
    when it rebuilds -- the picked point is the only handle that survives.

    Every accepted point is built in ONE pass through :func:`blocks_at`, so
    their seams can end on each other rather than each running to a wall; the
    parity law permits a pair of odd blocks to absorb one another's seam.
    Chirality is searched with :func:`blocks_any_spin` if the first assignment
    collides, and only then does it fall back to one block at a time.

    Parameters
    ----------
    mesh : Mesh
        A DENSE quad mesh. Poles elsewhere in it are fine; a pole AT a picked
        point is not, and is refused with that reason.
    points : list of point
        The block centres, in model coordinates.
    ratio : float, optional
        How far along its incident edges a block's corners sit -- the block's
        size. Its degree is not a parameter; the index pins it to the valence.
    tol : float, optional
        How far a point may be from the joint it binds to. Defaults to 1.5x the
        median edge length.
    relax_iters : int, optional
        Centroid smoothing passes afterwards, boundaries pinned. The block's
        corners land at ``ratio``, which is a guess until this runs.

    Returns
    -------
    (mesh, report)
        ``report['blocks']`` is one record per point, in the order given, each
        with ``ok``, ``reason``, ``valence``, ``degree``, ``distance`` (picked
        point to joint) and ``offset`` (picked point to the finished block's
        centroid). ``report['index_before']`` and ``report['index_after']``
        must be equal -- if they are not, a block was invented, not moved.
    """
    points = [list(p) for p in points]
    records = [{'point': p, 'ok': False, 'reason': '', 'valence': None,
                'degree': None, 'spin': None, 'distance': None, 'offset': None}
               for p in points]
    report = {'blocks': records, 'built': 0, 'refused': 0,
              'strips': 0, 'strip_faces': 0, 'ends': {},
              'index_before': index_sum(mesh), 'index_after': None,
              'spins': None, 'mode': 'none'}

    # -- bind each point to a joint ------------------------------------------
    taken = {}
    accepted = []
    for i, point in enumerate(points):
        vkey, distance, reason = joint_at(mesh, point, tol)
        records[i]['distance'] = distance
        if vkey is None:
            records[i]['reason'] = reason
            continue
        if vkey in taken:
            records[i]['reason'] = (
                'binds to the same joint as point {} -- one joint carries one '
                'block'.format(taken[vkey]))
            continue
        spin, why = blockable(mesh, vkey)
        if spin is None:
            records[i]['reason'] = why
            continue
        taken[vkey] = i
        records[i].update(valence=mesh.vertex_degree(vkey), spin=spin,
                          reason=why)
        accepted.append((i, vkey, spin, list(mesh.vertex_coordinates(vkey))))

    if not accepted:
        report['refused'] = len(records)
        report['index_after'] = report['index_before']
        return mesh, report

    # -- build them ----------------------------------------------------------
    vkeys = [vkey for _, vkey, _, _ in accepted]
    spins = {vkey: spin for _, vkey, spin, _ in accepted}
    built = list(accepted)
    try:
        out, built_report = blocks_at(mesh, vkeys, ratio, spins)
        report['mode'] = 'one pass'
    except Propagation:
        try:
            out, built_report = blocks_any_spin(mesh, vkeys, ratio)
            report['mode'] = 'one pass, chirality searched'
            report['spins'] = built_report.get('spins')
            for k, (i, vkey, _, xyz) in enumerate(accepted):
                records[i]['spin'] = built_report['spins'][k]
        except Propagation:
            out, built_report, built = _one_at_a_time(
                mesh, accepted, records, ratio, tol)
            report['mode'] = 'one at a time'

    report['strips'] = built_report['strips']
    report['strip_faces'] = built_report['strip_faces']
    ends = defaultdict(int)
    for e in built_report['ends']:
        ends[e] += 1
    report['ends'] = dict(ends)

    # -- report where each block actually landed -----------------------------
    for i, _, _, xyz in built:
        fkey = _block_face(out, xyz)
        records[i]['ok'] = True
        if fkey is not None:
            records[i]['degree'] = len(out.face_vertices(fkey))
            records[i]['offset'] = distance_point_point(
                records[i]['point'], out.face_centroid(fkey))

    if relax_iters:
        relax(out, relax_iters)

    report['built'] = sum(1 for r in records if r['ok'])
    report['refused'] = len(records) - report['built']
    report['index_after'] = index_sum(out)
    return out, report


def _one_at_a_time(
    mesh: Mesh,
    accepted: list[tuple[int, int, int, list[float]]],
    records: list[dict[str, Any]],
    ratio: float,
    tol: float | None,
) -> tuple[Mesh, dict[str, Any], list[tuple[int, int, int, list[float]]]]:
    """Fall back to one block per rebuild, re-locating each point as we go.

    The rebuild renumbers everything, so the remaining points have to be bound
    again against the new mesh -- which is exactly why the public API takes
    coordinates.
    """
    out = mesh
    total = {'strips': 0, 'strip_faces': 0, 'ends': []}
    built = []
    for i, _, _, xyz in accepted:
        vkey, distance, reason = joint_at(out, xyz, tol)
        if vkey is None:
            records[i]['reason'] = ('after an earlier block was built: '
                                    '{}'.format(reason))
            continue
        spin, why = blockable(out, vkey)
        if spin is None:
            records[i]['reason'] = ('after an earlier block was built: '
                                    '{}'.format(why))
            continue
        valence = out.vertex_degree(vkey)
        try:
            out, one = block_at(out, vkey, ratio, spin)
        except Propagation as exc:
            records[i]['reason'] = str(exc)
            continue
        records[i].update(valence=valence, spin=spin, reason=why)
        total['strips'] += one['strips']
        total['strip_faces'] += one['strip_faces']
        total['ends'] += one['ends']
        built.append((i, vkey, spin, xyz))
    return out, total, built


# ----------------------------------------------------------------------------
# the other route: collapse the pole's fan, the dual's trade done locally
# ----------------------------------------------------------------------------

def pole_at(
    mesh: Mesh, point: list[float], tol: float | None = None
) -> tuple[int | None, float | None, str]:
    """The pole at ``point``: ``(vkey, distance, reason)``.

    A pole here is any INTERIOR vertex carrying at least one non-quad face in
    its 1-ring -- which is what a point feature becomes, and what
    :func:`blockable` refuses. ``vkey`` is None when there is no pole to
    collapse within ``tol``.
    """
    if tol is None:
        tol = 1.5 * _median_edge_length(mesh)

    nearest, best = None, None
    nearest_any, best_any = None, None
    for v in mesh.vertices():
        d = distance_point_point(point, mesh.vertex_coordinates(v))
        if best_any is None or d < best_any:
            nearest_any, best_any = v, d
        if mesh.is_vertex_on_boundary(v):
            continue
        if all(len(mesh.face_vertices(f)) == 4
               for f in mesh.vertex_faces(v, ordered=True)):
            continue
        if best is None or d < best:
            nearest, best = v, d

    if nearest is not None and best <= tol:
        return nearest, best, 'fan of {} faces'.format(
            len(mesh.vertex_faces(nearest)))
    if nearest_any is not None and best_any <= tol and \
            mesh.is_vertex_on_boundary(nearest_any):
        return None, best_any, 'the vertex at this point is on the boundary'
    if nearest is None:
        return None, best_any, ('there is no pole here -- every face round the '
                                'nearest vertex is a quad')
    return None, best, ('the nearest pole is {:.4f} away, beyond the {:.4f} '
                        'tolerance'.format(best, tol))


def pole_blocks(
    mesh: Mesh, points: list[list[float]], tol: float | None = None, relax_iters: int = 0
) -> tuple[Mesh, dict[str, Any]]:
    """Collapse the FAN at each pole into ONE n-gon face, centred on the point.

    This is what ``mesh_dual_conway`` does to a pole, done at one vertex
    instead of to the whole plate: drop the pole and the fan of faces round it,
    and close the hole with a single face on the fan's ring. The block lands
    EXACTLY on the picked point and it costs **no seams at all** -- nothing
    propagates, every other face in the mesh is untouched.

    The price is on the ring instead. Each ring vertex loses its edge to the
    pole, so a regular one drops from valence 4 to valence 3: an n-gon block
    ringed by n valence-3 joints. That is the same trade the global dual makes
    at a pole, and :func:`index_sum` is conserved through it exactly -- the
    pole's index 4 does not vanish, it is spread around the ring.

    Contrast :func:`block_points`, which truncates a genuine valence-n JOINT:
    that one leaves the ring regular and pays in seams instead. Neither is free;
    they charge in different currency.

        collapse (here)   exact placement, 0 seams, n valence-3 joints on the
                          ring, degree = the fan's size (so it follows the
                          densities, and is 4 -- a plain quad -- on a coarse
                          4-patch pole)

        truncate          ring stays regular, n seams run to the wall, degree
                          pinned to the joint's valence, and it needs a joint
                          rather than a pole to start from

    Returns ``(mesh, report)`` with the same shape as :func:`block_points`'s,
    plus ``report['ring_joints']`` -- how many valence-3 joints the collapse
    left behind, which is the number to look at before choosing this route.

    Notes
    -----
    ``index_sum`` is conserved exactly when the fan lies clear of the wall, and
    NOT when it touches it -- and that is a property of the bookkeeping, not a
    fault in the collapse. ``index_sum`` counts interior vertices only, so a
    ring vertex ON the boundary loses a valence without its term being counted,
    and the total drops by one for each such vertex. The result is still
    manifold and still correct; ``report['index_before']`` and
    ``report['index_after']`` are both returned so the drop is visible rather
    than silent. Measured on the British Museum poles, where every pole's fan
    reaches the wall: 3 ring vertices on the wall costs -2, 1 costs -1.

    The fan is closed on its own boundary CYCLE, which for a partial pole (one
    with quads in the fan as well as triangles) is longer than the neighbour
    ring -- a quad's far corner is on the cycle but is not a neighbour of the
    pole. The block's degree is that cycle's length, so a 4-face partial fan
    can give a hexagon.
    """
    points = [list(p) for p in points]
    records = [{'point': p, 'ok': False, 'reason': '', 'valence': None,
                'degree': None, 'spin': None, 'distance': None, 'offset': None}
               for p in points]
    report = {'blocks': records, 'built': 0, 'refused': 0,
              'strips': 0, 'strip_faces': 0, 'ends': {},
              'index_before': index_sum(mesh), 'index_after': None,
              'spins': None, 'mode': 'collapse', 'ring_joints': 0}

    targets = {}
    for i, point in enumerate(points):
        vkey, distance, reason = pole_at(mesh, point, tol)
        records[i]['distance'] = distance
        if vkey is None:
            records[i]['reason'] = reason
            continue
        if vkey in targets:
            records[i]['reason'] = (
                'binds to the same pole as point {}'.format(targets[vkey]))
            continue
        targets[vkey] = i
        records[i].update(valence=mesh.vertex_degree(vkey), reason=reason)

    # The fan has to be a disc, and it is closed on its own boundary cycle --
    # NOT on the neighbour ring, which is only the same thing for a full pole.
    cycles = {}
    for vkey, i in list(targets.items()):
        cycle = _fan_boundary(mesh, mesh.vertex_faces(vkey))
        if cycle is None:
            records[i]['reason'] = ('the fan round this pole is not a simple '
                                    'disc, so it has no single ring to close on')
            records[i]['valence'] = None
            del targets[vkey]
            continue
        cycles[vkey] = cycle

    # Two poles whose fans touch would each eat the other's ring.
    for vkey, i in list(targets.items()):
        clash = [w for w in targets
                 if w != vkey and (w in cycles[vkey]
                                   or set(mesh.vertex_faces(w))
                                   & set(mesh.vertex_faces(vkey)))]
        if clash:
            records[i]['reason'] = (
                'this pole\'s fan touches the fan of point {} -- collapse '
                'them one at a time'.format(targets[clash[0]]))
            records[i]['valence'] = None
            del targets[vkey]

    if not targets:
        report['refused'] = len(records)
        report['index_after'] = report['index_before']
        return mesh, report

    dropped = set()
    for vkey in targets:
        dropped.update(mesh.vertex_faces(vkey))

    points_out, index = [], {}

    def add(v: int) -> int:
        if v not in index:
            index[v] = len(points_out)
            points_out.append(mesh.vertex_coordinates(v))
        return index[v]

    faces = []
    for f in mesh.faces():
        if f in dropped:
            continue
        faces.append([add(v) for v in mesh.face_vertices(f)])
    rings = {}
    for vkey, i in targets.items():
        ring = cycles[vkey]
        rings[i] = [mesh.vertex_coordinates(v) for v in ring]
        faces.append([add(v) for v in ring])

    out = Mesh.from_vertices_and_faces(points_out, faces)
    out.unify_cycles()

    for i, ring_xyz in rings.items():
        centre = centroid_points(ring_xyz)
        # match on the degree too -- a 4-fan collapses to a QUAD, which
        # ``_block_face`` would skip while looking for a non-quad.
        fkey = _face_of_degree_at(out, centre, len(ring_xyz))
        records[i]['ok'] = True
        if fkey is not None:
            records[i]['degree'] = len(out.face_vertices(fkey))
            records[i]['offset'] = distance_point_point(
                records[i]['point'], out.face_centroid(fkey))

    report['ring_joints'] = len(joints(out)) - len(joints(mesh))
    if relax_iters:
        relax(out, relax_iters)
    report['built'] = sum(1 for r in records if r['ok'])
    report['refused'] = len(records) - report['built']
    report['index_after'] = index_sum(out)
    return out, report


def _block_face(mesh: Mesh, xyz: list[float]) -> int | None:
    """The non-quad face nearest ``xyz`` -- the block that was just built."""
    best, fkey = None, None
    for f in mesh.faces():
        if len(mesh.face_vertices(f)) == 4:
            continue
        d = distance_point_point(xyz, mesh.face_centroid(f))
        if best is None or d < best:
            best, fkey = d, f
    return fkey


def _fan_boundary(mesh: Mesh, fkeys: list[int]) -> list[int] | None:
    """The ordered boundary cycle of a patch of faces, or None if it is not a disc.

    Not the same as ``vertex_neighbors(vkey, ordered=True)``, and the
    difference is the whole correctness of :func:`pole_blocks`. That shortcut
    is only right for a FULL pole, where every fan face is a triangle
    ``[P, r_i, r_i+1]`` and consecutive neighbours therefore share an edge. A
    PARTIAL pole has quads in its fan, and a quad's far corner is not a
    neighbour of the pole at all -- fill on the neighbour ring there and those
    corners are left dangling on edges that do not exist.
    """
    inner = set(fkeys)
    nxt = {}
    for f in inner:
        fv = mesh.face_vertices(f)
        n = len(fv)
        for i in range(n):
            u, v = fv[i], fv[(i + 1) % n]
            if mesh.halfedge[v].get(u) not in inner:
                if u in nxt:            # pinched: the patch is not a disc
                    return None
                nxt[u] = v
    if not nxt:
        return None
    start = next(iter(nxt))
    cycle, v = [start], nxt[start]
    while v != start:
        if v not in nxt:
            return None
        cycle.append(v)
        v = nxt[v]
    if len(cycle) != len(nxt):          # more than one cycle -- an annulus
        return None
    return cycle


def _face_of_degree_at(mesh: Mesh, xyz: list[float], degree: int) -> int | None:
    """The face of exactly ``degree`` sides nearest ``xyz``."""
    best, fkey = None, None
    for f in mesh.faces():
        if len(mesh.face_vertices(f)) != degree:
            continue
        d = distance_point_point(xyz, mesh.face_centroid(f))
        if best is None or d < best:
            best, fkey = d, f
    return fkey


def report_lines(report: dict[str, Any]) -> list[str]:
    """``block_points``'s report as printable lines, one per point plus totals."""
    lines = []
    for i, r in enumerate(report['blocks']):
        x, y, _ = r['point']
        head = 'point {} ({:7.3f}, {:7.3f})'.format(i, x, y)
        if r['ok'] and report['mode'] == 'collapse':
            what = 'BLOCK' if (r['degree'] or 4) != 4 else 'quad (NOT a block)'
            lines.append(
                '{}  {} degree {}  (fan of {}, centroid {:.4f} away)'.format(
                    head, what, r['degree'], r['valence'], r['offset'] or 0.0))
        elif r['ok']:
            lines.append(
                '{}  BLOCK degree {}  (valence {}, spin {}, joint {:.4f} away, '
                'centroid {:.4f} away)'.format(
                    head, r['degree'], r['valence'], r['spin'],
                    r['distance'] or 0.0, r['offset'] or 0.0))
        else:
            lines.append('{}  REFUSED: {}'.format(head, r['reason']))
    if report['mode'] == 'collapse':
        lines.append(
            'built {} of {}  (collapse); 0 seams -- nothing propagated, but '
            '{:+d} new joints on the ring'.format(
                report['built'], len(report['blocks']),
                report.get('ring_joints', 0)))
    else:
        lines.append(
            'built {} of {}  ({}); {} seams, {} faces split, ends {}'.format(
                report['built'], len(report['blocks']), report['mode'],
                report['strips'], report['strip_faces'],
                report['ends'] or '{}'))
    before, after = report['index_before'], report['index_after']
    if before == after:
        lines.append('index {} conserved -- the singularities were MOVED into '
                     'faces, not invented'.format(before))
    else:
        lines.append('! INDEX CHANGED {} -> {}: singularities were invented, '
                     'not moved'.format(before, after))
    return lines


# ----------------------------------------------------------------------------
# afterwards
# ----------------------------------------------------------------------------

def relax(mesh: Mesh, kmax: int = 30, slide: bool = False, kink: float = 30.0) -> Mesh:
    """Centroid smoothing with the boundary held, in place.

    The block's corners land on the incident edges at ``ratio``, which is a
    guess; this is what turns that guess into an element. Boundaries are read
    from ``vertices_on_boundaries()`` -- the PLURAL, because the singular form
    returns only the longest ring and would let holes round off.

    Parameters
    ----------
    slide : bool, optional
        How the boundary is held. ``False`` pins every boundary vertex where it
        stands, which keeps the outline exactly but also freezes the spacing of
        the boundary row, so a block near a wall cannot pull that row into line
        with itself. ``True`` lets a boundary vertex SLIDE along the outline:
        it moves to the centroid of its neighbours like an interior one, and is
        then projected back onto the polyline the boundary started as. The
        normal component of the move is discarded and only the tangential one
        survives, so the outline is preserved to the same tolerance as before
        while the spacing is free to even out.
    kink : float, optional
        Corner threshold in degrees. A boundary vertex whose outline turns by
        more than this is a CORNER and stays pinned whatever ``slide`` says --
        corners are the shape of the domain, and letting them drift rounds it
        off. The default keeps anything sharper than 150 degrees.

    Notes
    -----
    Sliding only ever moves a vertex along the ORIGINAL outline, captured once
    before the first pass, so repeated passes cannot creep inwards the way
    projecting onto the current boundary would.
    """
    from compas.geometry import closest_point_on_polyline

    rings = []
    for ring in mesh.vertices_on_boundaries():
        # the ring comes back closed -- the first vertex is repeated at the end
        rings.append(ring[:-1] if len(ring) > 1 and ring[0] == ring[-1]
                     else ring)

    fixed = set()
    for ring in rings:
        fixed.update(ring)

    if not slide:
        sliding, outlines = {}, {}
    else:
        sliding, outlines = {}, {}
        for r, ring in enumerate(rings):
            pts = [mesh.vertex_coordinates(v) for v in ring]
            outlines[r] = pts + pts[:1]              # closed polyline
            n = len(ring)
            for i, v in enumerate(ring):
                turn = angle_points(pts[i], pts[i - 1], pts[(i + 1) % n]) * DEG
                if turn >= 180.0 - kink:
                    sliding[v] = r                   # straight enough to slide
        fixed -= set(sliding)

    for _ in range(kmax):
        moved = {}
        for v in mesh.vertices():
            if v in fixed:
                continue
            target = centroid_points(
                [mesh.vertex_coordinates(u) for u in mesh.vertex_neighbors(v)])
            if v in sliding:
                # keep only the component along the outline
                target = closest_point_on_polyline(target, outlines[sliding[v]])
            moved[v] = target
        for v, xyz in moved.items():
            mesh.vertex_attributes(v, 'xyz', xyz)
    return mesh


def seam_edges(
    primal: Mesh, result: Mesh, tol: int = 5
) -> list[tuple[list[float], list[float]]]:
    """Edges of ``result`` that are not edges of ``primal`` -- what a block cost.

    Matched on COORDINATES, so call it before :func:`relax`: the whole claim is
    that the mesh away from the seams is untouched vertex for vertex, and
    smoothing is what stops that being literally true.
    """
    old = set()
    for u, v in primal.edges():
        a = tuple(round(c, tol) for c in primal.vertex_coordinates(u))
        b = tuple(round(c, tol) for c in primal.vertex_coordinates(v))
        old.add(frozenset((a, b)))
    out = []
    for u, v in result.edges():
        pa, pb = result.vertex_coordinates(u), result.vertex_coordinates(v)
        a = tuple(round(c, tol) for c in pa)
        b = tuple(round(c, tol) for c in pb)
        if frozenset((a, b)) not in old:
            out.append((pa, pb))
    return out
