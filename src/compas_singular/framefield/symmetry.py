"""Keep a symmetric domain's symmetry all the way through to the layout.

WIRED IN, and ON BY DEFAULT. ``FieldDecomposition.from_boundary`` detects the
symmetry of its own input and applies all of :data:`STEPS`; pass
``symmetry=None`` to get the pre-2026-08-14 behaviour back on any domain.
Adopted 2026-08-14 after the evaluation in ``24_symmetry.py``; the files as they
were before are in ``compas_singular/backup_2026-08-14_pre_symmetry``.

This module holds the machinery. The five places the pipeline consults it are
``background.BackgroundMesh.from_boundary``, ``constraints.from_curves``,
``trace.Tracer``, ``repair.build_network`` and ``repair._cluster`` -- each a
handful of lines, each pointing back here.

THE PROBLEM
-----------

A square with four symmetric guide arcs produces a visibly asymmetric coarse
layout. The solve is not at fault: on a symmetric mesh it is exactly
symmetry-preserving (measured below). The symmetry is lost before the field is
solved, and then the layout stage amplifies it.

Measured on a 10x10 square with four arcs that are exact D4 images of each
other, at spacing 0.5. ``exact`` is the share of points whose image under the
operation is also a point; ``dev`` is the worst distance to the nearest one::

    stage                       rot90            mirror
    boundary vertices           100%             100%
    interior vertices             0% / 0.65        0% / 0.65
    singular faces                0% / 0.49        0% / 0.36
    coarse vertices              15% / 0.49       15% / 0.49

0.65 at a 0.5 spacing is more than a whole triangle. FIVE independent causes,
in the order they bite. None of them is a tolerance, and none of them is the
solve:

1. **The interior grid is asymmetric by construction** -- ``background.py``.
   ``_jitter`` hashes the grid indices, and a hash respects no symmetry; the
   grid is also anchored at ``min(xs)`` rather than centred. The two are
   separable::

       spacing 0.5  (10/0.5 = 20.00)  jitter ON -> rot90   0%   OFF -> 100%
       spacing 0.3  (10/0.3 = 33.33)  jitter ON -> rot90   0%   OFF ->   0%

   The jitter breaks it always; the anchoring breaks it additionally whenever
   the spacing does not divide the domain. The jitter is NOT removable -- a
   regular grid makes every cell cocircular and Qhull then picks diagonals
   arbitrarily -- so :func:`interior_points` keeps it and makes it EQUIVARIANT
   instead.

2. **``from_curves`` is order-dependent.** Its ``claimed`` set gives each vertex
   to the FIRST curve in list order, so a vertex near two arcs and its rotated
   twin take their tangent from different arcs. With a perfectly symmetric mesh
   and perfectly symmetric arcs the constraint set still mismatches by 1.145
   under every non-identity element. ``constraints.from_curves`` no longer has
   it -- that one is a plain bug fix and is applied whether or not any symmetry
   was detected.

3. **A singularity is reported as a FACE, and the face is a cruder object than
   the field on it.** After 1 and 2 the field is symmetric to 2e-16 and the
   layout still is not, because the triangle that wins the winding tie around a
   zero is not equivariant -- and a symmetric field puts its zeros exactly ON
   the symmetry axes, where several triangles are equally entitled.
   :func:`snap_singularities` takes the position from the field's own ``|u|``
   minimum instead. The face-index tie in ``wrap_to_period`` is left exactly as
   it is; it stops mattering once nothing reads a position off it.

4. **``build_network`` recomputed those centroids** and put them back on the
   front of every trace, undoing 3. It now takes them from the tracer. This was
   the single largest residual: the traces go in exactly symmetric -- all 472
   points of all 32 separatrices, deviation 0.0000 -- and came out of that loop
   at 69%, deviation 0.3486.

5. **``repair._cluster`` is greedy first-fit in trace order**, so two points
   that are exact mirror images can join different clusters. It still is; the
   result is now projected onto the group by :func:`project_clusters`, which
   leaves every decision the original makes -- the ``groups`` rule, the tighter
   ``wall_tol`` -- to the original code and only makes it consistent across an
   orbit.

All five are now edits in place. They were developed as scoped monkeypatch
shims so the workflow could be left alone while the fix was measured; the shims
are gone.

WHAT THIS BUYS, measured on the four-arc square at spacing 0.5
--------------------------------------------------------------

::

                              patches  poles  launch pts   coarse rot   mirror
    baseline                      18      2      0%        15% / 0.59    15%
    + symmetric background        16      0     50%        40% / 0.42    32%
    + equivariant guides          16      0     50%        100%/ 0.44    36%
    + field symmetrisation        16      0     50%        84% / 0.44    68%
    + snapped singularities       16      0    100%        92% / 0.35    84%
    + network fixes               16      0    100%       100% / 0.0000  100%

Exact under all eight elements, deviation 0.0000, with Poincare-Hopf still
passing and no repair warnings. Element quality improves with it rather than
against it: minimum angle 29.75 -> 57.83 degrees, worst aspect 2.01 -> 1.60.

THE COST, because there is one
------------------------------

The layout stops splitting where the asymmetry used to split it, so a guide is
followed slightly less closely by a simpler layout. Measured as mean angle
between dense mesh edges and the guides::

    square+4 arcs   10.2 deg over 18 patches  ->  10.7 deg over 16
    rect+2 arcs     13.0 deg over 17 patches  ->  14.8 deg over  1

The rectangle row is the one that was argued about before adopting: the
original found 6
singularities there and this route finds none. Both are defensible -- two arcs
on a rectangle force no winding the walls cannot absorb, so the six were
artefacts of the asymmetric background -- but it is a real change in what the
layout does, not only in how symmetric it looks. Its minimum angle goes 12.54
-> 78.54 degrees.
"""
from __future__ import annotations

from math import ceil
from typing import Any

from compas.geometry import is_point_in_polygon_xy

from compas_singular.geometry.polyline import distance_to_loop
from compas_singular.framefield.background import _jitter
from compas_singular.framefield.field import CrossField


__all__ = ['Symmetry', 'STEPS', 'interior_points', 'symmetrise',
           'snap_singularities', 'project_clusters', 'invariance',
           'field_invariance', 'singularity_orbits']


#: ``(a, b, c, d, reflects, name)`` -- the linear part as a 2x2 matrix, whether
#: it reverses orientation, and a name to report it by. A cross field cares
#: about the reflection flag and nothing else about the matrix: a rotation by 90
#: degrees leaves ``exp(i*4*theta)`` UNCHANGED, and every reflection whose axis
#: passes through the centre at a multiple of 45 degrees conjugates it. See
#: :func:`symmetrise`.
ELEMENTS = [
    (1, 0, 0, 1, False, 'identity'),
    (0, -1, 1, 0, False, 'rot90'),
    (-1, 0, 0, -1, False, 'rot180'),
    (0, 1, -1, 0, False, 'rot270'),
    (-1, 0, 0, 1, True, 'mirror-vertical'),
    (1, 0, 0, -1, True, 'mirror-horizontal'),
    (0, 1, 1, 0, True, 'mirror-diagonal'),
    (0, -1, -1, 0, True, 'mirror-antidiagonal'),
]


#: The five places the pipeline has to be told about the symmetry, named so
#: they can be switched off one at a time. ``FieldDecomposition.from_boundary``
#: enables all of them; ``24_symmetry.py`` turns them off in sequence to show
#: what each is worth, which is the only reason this is a set and not five
#: hard-coded ``if`` statements.
STEPS = ('background', 'guides', 'field', 'singularities', 'network')


class Symmetry(object):
    """A subgroup of the square's symmetries, about a centre.

    Attributes
    ----------
    centre : [x, y, z]
    elements : list[tuple]
        A subset of :data:`ELEMENTS`, always containing the identity.
    steps : tuple[str]
        Which of :data:`STEPS` to apply. All of them unless a caller says
        otherwise.

    Notes
    -----
    Only the 8 symmetries of the square are considered. That is not a
    limitation of the method -- everything below works for any finite group of
    isometries -- but it covers every symmetric domain this front end has been
    put in front of, and detecting a 6-fold symmetry nobody drew would be a
    liability rather than a feature.
    """

    def __init__(
        self,
        centre: list[float],
        elements: list[tuple[int, int, int, int, bool, str]],
        steps: tuple[str, ...] = STEPS,
    ) -> None:
        self.centre = [float(centre[0]), float(centre[1]), 0.0]
        self.elements = list(elements)
        self.steps = tuple(steps)

    def enabled(self, step: str) -> bool:
        """Whether ``step`` should be applied. False for a trivial group.

        Every call site asks this rather than testing ``symmetry is not None``,
        so a domain with no symmetry takes exactly the path it took before this
        module existed -- which is the property the L-shape control in
        ``24_symmetry.py`` checks by comparing layouts byte for byte.
        """
        return not self.trivial and step in self.steps

    def __len__(self) -> int:
        return len(self.elements)

    def __repr__(self) -> str:
        return '<Symmetry order {} about ({:.3f}, {:.3f}): {}>'.format(
            len(self.elements), self.centre[0], self.centre[1],
            ', '.join(self.names()))

    def names(self) -> list[str]:
        return [g[5] for g in self.elements]

    @property
    def trivial(self) -> bool:
        return len(self.elements) < 2

    def apply(self, g: tuple[int, int, int, int, bool, str], point: list[float]) -> list[float]:
        """``g`` applied to a point, about :attr:`centre`."""
        a, b, c, d = g[0], g[1], g[2], g[3]
        x, y = point[0] - self.centre[0], point[1] - self.centre[1]
        return [a * x + b * y + self.centre[0],
                c * x + d * y + self.centre[1], 0.0]

    @staticmethod
    def linear(g: tuple[int, int, int, int, bool, str], vector: list[float]) -> list[float]:
        """``g``'s linear part applied to a free vector -- no centre involved."""
        a, b, c, d = g[0], g[1], g[2], g[3]
        return [a * vector[0] + b * vector[1], c * vector[0] + d * vector[1], 0.0]

    # ------------------------------------------------------------------
    # detection
    # ------------------------------------------------------------------

    @classmethod
    def detect(
        cls,
        loops: list[list[list[float]]],
        centre: list[float] | None = None,
        tol: float = 1e-6,
    ) -> Symmetry:
        """The subgroup that maps every one of ``loops`` onto the union of them.

        Parameters
        ----------
        loops : list[list[[x, y, z]]]
            Every curve the field will see -- outer boundary, holes, guides.
            All of them, because a symmetric plate with an off-centre cable is
            not a symmetric problem and must not be reported as one.
        centre : [x, y, z], optional
            Defaults to the centre of the outer boundary's bounding box. Passed
            explicitly when the bounding box is not centred on the symmetry --
            which cannot happen for any subgroup of the square's symmetries, so
            the default is safe and the parameter is for callers who know
            better.
        tol : float, optional
            Point-matching tolerance. Loose by geometric standards (1e-6 on a
            10-unit plate) because Rhino curve sampling is not bit-exact even
            for a mirrored curve, and a symmetry missed here is silently not
            enforced anywhere downstream.

        Returns
        -------
        Symmetry

        Notes
        -----
        Matching is by POINT SET, not curve by curve. A guide and its mirror
        image are two different curves and the map has to be allowed to swap
        them; requiring each curve to map to itself would find only the
        symmetries an individual arc happens to have.
        """
        points = [p for loop in loops for p in loop]
        if not points:
            raise ValueError('no geometry to detect a symmetry from')
        if centre is None:
            xs = [p[0] for p in loops[0]]
            ys = [p[1] for p in loops[0]]
            centre = [(max(xs) + min(xs)) / 2.0, (max(ys) + min(ys)) / 2.0, 0.0]

        self = cls(centre, [ELEMENTS[0]])
        have = {(round(p[0] / tol), round(p[1] / tol)) for p in points}
        found = []
        for g in ELEMENTS:
            ok = True
            for p in points:
                q = self.apply(g, p)
                if (round(q[0] / tol), round(q[1] / tol)) not in have:
                    ok = False
                    break
            if ok:
                found.append(g)
        return cls(centre, found)


# ----------------------------------------------------------------------
# 1 -- interior points invariant under the group
# ----------------------------------------------------------------------

def interior_points(
    symmetry: Symmetry,
    target_length: float,
    outer: list[list[float]],
    inners: Any = (),
    margin: float = 0.45,
    on_axis: bool = True,
) -> list[list[float]]:
    """Interior points, as a set exactly invariant under ``symmetry``."""
    centre = symmetry.centre
    half = target_length / 2.0
    parity = 0 if on_axis else 1          # doubled-lattice parity
    limit = margin * target_length

    xs = [p[0] for p in outer]
    ys = [p[1] for p in outer]
    reach = int(ceil(max(max(xs) - centre[0], centre[0] - min(xs),
                         max(ys) - centre[1], centre[1] - min(ys))
                     / target_length)) + 1

    def point_of(u: int, v: int, jitter: list[float]) -> list[float]:
        return [centre[0] + u * half + jitter[0],
                centre[1] + v * half + jitter[1], 0.0]

    def act(g: tuple[int, int, int, int, bool, str], u: int, v: int) -> tuple[int, int]:
        """The group acting on DOUBLED lattice indices. Exact -- integers only."""
        return (g[0] * u + g[1] * v, g[2] * u + g[3] * v)

    out = []
    seen_orbits = set()
    emitted = set()
    for i in range(-reach, reach + 1):
        for j in range(-reach, reach + 1):
            u, v = 2 * i + parity, 2 * j + parity
            orbit = [act(g, u, v) for g in symmetry.elements]
            rep = min(orbit)
            if rep in seen_orbits:
                continue
            seen_orbits.add(rep)

            # one jitter per orbit, from the representative, projected onto
            # what the representative's stabiliser fixes
            raw = [_jitter(rep[0], rep[1], 0.2 * target_length),
                   _jitter(rep[1], rep[0], 0.2 * target_length), 0.0]
            stabiliser = [g for g in symmetry.elements if act(g, *rep) == rep]
            jitter = [0.0, 0.0, 0.0]
            for g in stabiliser:
                moved = Symmetry.linear(g, raw)
                jitter = [jitter[k] + moved[k] / len(stabiliser) for k in range(3)]

            # keep or drop the whole orbit on the representative's verdict
            if not _inside(point_of(rep[0], rep[1], jitter), outer, inners, limit):
                continue

            for g in symmetry.elements:
                index = act(g, *rep)
                if index in emitted:
                    continue
                emitted.add(index)
                out.append(point_of(index[0], index[1], Symmetry.linear(g, jitter)))
    return out


def _inside(p: list[float], outer: list[list[float]], inners: Any, limit: float) -> bool:
    if not is_point_in_polygon_xy(p, outer):
        return False
    if distance_to_loop(p, outer) < limit:
        return False
    for loop in inners:
        if is_point_in_polygon_xy(p, loop):
            return False
        if distance_to_loop(p, loop) < limit:
            return False
    return True


# ----------------------------------------------------------------------
# 3 -- project the solved field onto the symmetric subspace
# ----------------------------------------------------------------------

def symmetrise(field: CrossField, symmetry: Symmetry, tol: float = 1e-9) -> CrossField:
    """Group-average ``u``, and return a new :class:`CrossField`.

    Cheap, and exact to machine precision: measured ``max|u - rho(g)u|`` of
    3.5e-16 afterwards, against 4.4e-01 before.

    The averaging is two lines because of the 4th-power representation. A
    rotation by 90 degrees maps ``theta -> theta + 90``, so ``exp(i*4*theta)``
    is UNCHANGED -- rotations act trivially. A reflection in a line at angle
    ``alpha`` maps ``theta -> 2*alpha - theta``, so ``u -> exp(i*8*alpha) *
    conj(u)``, and for every axis of the square (0, 45, 90, 135 degrees)
    ``exp(i*8*alpha) = 1``. So::

        u_sym(p) = mean over g of  ( conj(u(g p)) if g reflects else u(g p) )

    That the reflection part is antilinear does not stop the average being
    invariant: conjugation is additive, so ``rho(h)`` still commutes with the
    sum and ``rho(h) u_sym = u_sym``.

    WHY THIS IS NOT ENOUGH ON ITS OWN
    ---------------------------------

    An exactly symmetric field still gave an 84%-symmetric layout. The reason
    is worth understanding before touching anything else, because it looks
    exactly like a tolerance problem and is not one.

    An equivariant field is exactly REAL on a mirror axis (``u = conj(u)``
    there), so ``theta = arg(u)/4`` lands exactly on the +/-45 degree branch --
    exactly half a period. ``wrap_to_period`` resolves a half-period difference
    oddly, ``wrap(+45) = +45`` and ``wrap(-45) = -45``, so the winding around a
    triangle depends on the direction its vertex cycle is traversed -- and a
    reflection reverses that cycle. Two triangles that are exact mirror images,
    with exact mirror-image angles, measured::

        face    theta = 45.000, -26.273,  0.000  ->  index +1
        mirror  theta = 45.000,   0.000, 26.273  ->  index  0

    That is systematic rather than a knife-edge accident: a symmetric field
    puts its zeros ON the symmetry axes, and that is precisely where a
    face-based index is degenerate. The offending vertex above has
    ``|u| = 0.0028`` and sits at ``x = 5.0000`` -- the zero is AT a background
    vertex, on the axis, and no single face owns it.

    Do NOT fix this by biasing the tie-break. It makes the symptom vanish and
    the answer wrong: the singular faces then measure 100% invariant under all
    eight elements, with an index sum of 8 against a boundary winding of 0, and
    ``poincare_hopf`` reports ``ok: False``. That check is earning its keep.

    :func:`snap_singularities` resolves it instead by not reading a POSITION
    off the face at all. The index still comes from the winding, ties and all;
    only the point comes from the field. See its docstring for the numbers.
    """
    mesh = field.background.mesh
    grid = {}
    for vkey in mesh.vertices():
        p = mesh.vertex_coordinates(vkey)
        grid[(_snap(p[0], tol), _snap(p[1], tol))] = vkey

    values = {}
    orphans = 0
    for vkey in mesh.vertices():
        p = mesh.vertex_coordinates(vkey)
        acc, n = 0j, 0
        for g in symmetry.elements:
            q = symmetry.apply(g, p)
            partner = grid.get((_snap(q[0], tol), _snap(q[1], tol)))
            if partner is None:
                continue
            value = field.u[partner]
            acc += value.conjugate() if g[4] else value
            n += 1
        if n < len(symmetry.elements):
            orphans += 1
        values[vkey] = acc / n if n else field.u[vkey]

    out = CrossField(field.background, values, relaxed=field.relaxed,
                     iterations=field.iterations, residual=field.residual)
    out.symmetry_orphans = orphans
    return out


def _snap(value: float, tol: float) -> int:
    return int(round(value / tol))


# ----------------------------------------------------------------------
# 4 -- take the singularity off the arbitrary face it was reported on
# ----------------------------------------------------------------------

def snap_singularities(
    field: CrossField,
    relative_tol: float = 1e-9,
) -> tuple[CrossField, dict[int, list[float]], dict[str, Any]]:
    """Move each singularity onto the field's own ``|u|`` minimum near it.

    Returns a new :class:`CrossField` -- the caller's is not touched -- plus a
    report. The new field carries a rewritten singularity list, and
    :class:`SymmetricTracer` launches from the snapped POSITIONS rather than
    from face centroids.

    WHY THIS IS NEEDED EVEN WHEN THE FIELD IS ALREADY EXACT
    -------------------------------------------------------

    :func:`symmetrise` gets ``max|u - rho(g)u|`` to 2e-16, and the layout still
    came out 84% invariant rather than 100%. The reason is that a singularity
    is reported as a FACE, and a face is a much cruder object than the field
    on it. Measured on the four-arc square with everything else exact::

        the four +1 singularities, as faces
            (7.668, 5.174)  (2.332, 5.174)  (5.174, 7.668)  (4.826, 2.332)
        the |u| minimum under each of them
            (7.4522, 5.0000) (2.5478, 5.0000) (5.0000, 7.4522) (5.0000, 2.5478)

    The face centroids are not a rotational orbit -- rot90 of the first is
    ``(4.826, 7.668)``, which is not in the list. The minima are an exact
    orbit, and they sit exactly ON the mirror axes, which is where a symmetric
    field puts its singularities. The face reported is simply whichever of the
    several triangles around that vertex won the winding tie, and that choice
    is not equivariant.

    So the fix is not to compute the index differently -- it is to stop using
    the triangle as the position. The index still comes from the winding; only
    the point does not.

    HOW
    ---

    * the candidates are the singular face's own vertices plus its 1-ring, so a
      minimum one triangle away is still found;
    * ties are averaged rather than broken. An argmin is not equivariant when
      two vertices tie exactly, which on a symmetry axis is not rare -- it is
      what the axis IS. The centroid of the tied set is equivariant;
    * faces snapping to the same point are ONE singularity and are merged, with
      their indices summed. A pair summing to zero is dropped: two adjacent
      faces reporting ``+1`` and ``-1`` about the same zero are an artefact of
      the tie, not two singularities.

    Poincare-Hopf is re-checked afterwards and reported. If merging ever
    changes the index sum, that check is what says so -- do not skip it.
    """
    mesh = field.background.mesh
    singular = field.singularities()

    positions = {}
    for fkey, _ in singular:
        candidates = set(mesh.face_vertices(fkey))
        for vkey in list(candidates):
            candidates.update(mesh.vertex_neighbors(vkey))
        magnitudes = {v: abs(field.u[v]) for v in candidates}
        low = min(magnitudes.values())
        # equivariant argmin: average every vertex that ties for the minimum
        tied = [v for v, m in magnitudes.items()
                if m <= low + relative_tol * max(low, 1.0)]
        point = [sum(mesh.vertex_coordinates(v)[k] for v in tied) / len(tied)
                 for k in range(3)]
        positions[fkey] = point

    merged = {}
    for fkey, index in singular:
        key = tuple(_snap(c, 1e-9) for c in positions[fkey][:2])
        entry = merged.setdefault(key, {'faces': [], 'index': 0,
                                        'point': positions[fkey]})
        entry['faces'].append(fkey)
        entry['index'] += index

    points, indices = {}, []
    dropped = 0
    for entry in merged.values():
        if entry['index'] == 0:
            dropped += len(entry['faces'])
            continue
        fkey = entry['faces'][0]
        points[fkey] = entry['point']
        indices.append((fkey, entry['index']))

    out = CrossField(field.background, dict(field.u), relaxed=field.relaxed,
                     iterations=field.iterations, residual=field.residual)
    out._singularities = indices
    report = {
        'faces_in': len(singular),
        'singularities_out': len(indices),
        'merged': len(singular) - len(indices) - dropped,
        'dropped': dropped,
        'poincare_hopf': out.poincare_hopf(),
    }
    return out, points, report


# ----------------------------------------------------------------------
# 5 -- the end-snapping in repair.build_network
# ----------------------------------------------------------------------

def project_clusters(
    points: list[list[float]],
    base: dict[int, list[float]],
    symmetry: Symmetry,
    tol: float = 1e-6,
) -> dict[int, list[float]]:
    """Force ``base`` -- an ``{index: point}`` map -- to commute with the group.

    Keyed by POSITION, not by index. The list handed to ``_cluster`` is the two
    ends of every trace, and ends coincide constantly -- 44 of 64 entries share
    a position with another on the four-arc square, because every arm of a
    singularity ends at the singularity its neighbour started from. Keying by
    index made the orbit lookup return the first index holding a position and
    silently skip the rest, which is a projection that moves nothing.

    A position whose several indices were NOT all given the same representative
    is left exactly as the original left it. That happens only through the
    ``groups`` rule, which exists to stop one trace's two ends merging, and
    overriding it here would undo the thing it protects.
    """
    def key(p: list[float]) -> tuple[int, int]:
        return (_snap(p[0], tol), _snap(p[1], tol))

    by_position = {}
    for i, p in enumerate(points):
        by_position.setdefault(key(p), []).append(i)

    place = {}
    representative = {}
    for k, indices in by_position.items():
        place[k] = points[indices[0]]
        chosen = {key(base[i]) for i in indices}
        if len(chosen) == 1:
            representative[k] = base[indices[0]]

    out = dict(base)
    done = set()
    for k in sorted(representative):
        if k in done:
            continue
        orbit = []
        for g in symmetry.elements:
            q = symmetry.apply(g, place[k])
            if key(q) in representative:
                orbit.append((g, key(q)))
        if not orbit:
            continue
        # one member of the orbit decides for all of them
        source_g, source_k = min(orbit, key=lambda pair: pair[1])
        back = _inverse(source_g)
        for g, other in orbit:
            if other in done:
                continue
            moved = symmetry.apply(g, symmetry.apply(back, representative[source_k]))
            for i in by_position[other]:
                out[i] = moved
            done.add(other)
    return out


def _inverse(g: tuple[int, int, int, int, bool, str]) -> tuple[int, int, int, int, bool, str]:
    """The inverse of a signed permutation -- its transpose, both being orthogonal."""
    for other in ELEMENTS:
        if (other[0] * g[0] + other[1] * g[2] == 1
                and other[0] * g[1] + other[1] * g[3] == 0
                and other[2] * g[0] + other[3] * g[2] == 0
                and other[2] * g[1] + other[3] * g[3] == 1):
            return other
    raise ValueError('not an element of the group: {}'.format(g))


# ----------------------------------------------------------------------
# measurement
# ----------------------------------------------------------------------

def invariance(
    points: list[list[float]],
    symmetry: Symmetry,
    tol: float = 1e-6,
) -> dict[str, tuple[float, float]]:
    """How invariant a point set is, per group element.

    Returns
    -------
    dict[str, (float, float)]
        Name -> ``(share of points whose image is also a point, worst distance
        from an image to the nearest point)``.

    The two numbers say different things and both are needed. The share is what
    the layout cares about -- a patch corner either has a twin or it does not.
    The distance is what says how badly: ``build_network`` snaps at
    ``0.8 * target_length``, so a deviation near that size is the difference
    between a symmetric layout and one carrying an extra patch on one side.
    """
    pts = [(p[0], p[1]) for p in points]
    out = {}
    if not pts:
        return out
    have = {(_snap(x, tol), _snap(y, tol)) for x, y in pts}
    for g in symmetry.elements:
        hit = 0
        worst = 0.0
        for p in pts:
            q = symmetry.apply(g, p)
            if (_snap(q[0], tol), _snap(q[1], tol)) in have:
                hit += 1
            worst = max(worst, min(((q[0] - r[0]) ** 2 + (q[1] - r[1]) ** 2) ** 0.5
                                   for r in pts))
        out[g[5]] = (hit / float(len(pts)), worst)
    return out


def field_invariance(field: CrossField, symmetry: Symmetry, tol: float = 1e-9) -> float:
    """``max|u - rho(g)u|`` over vertices, worst over the group.

    Zero to machine precision means the field itself is symmetric, whatever the
    layout downstream then does with it.
    """
    mesh = field.background.mesh
    grid = {}
    for vkey in mesh.vertices():
        p = mesh.vertex_coordinates(vkey)
        grid[(_snap(p[0], tol), _snap(p[1], tol))] = vkey

    worst = 0.0
    for vkey in mesh.vertices():
        p = mesh.vertex_coordinates(vkey)
        for g in symmetry.elements:
            q = symmetry.apply(g, p)
            partner = grid.get((_snap(q[0], tol), _snap(q[1], tol)))
            if partner is None:
                continue
            value = field.u[partner]
            worst = max(worst, abs(field.u[vkey]
                                   - (value.conjugate() if g[4] else value)))
    return worst


def singularity_orbits(
    field: CrossField,
    symmetry: Symmetry,
    points: dict[int, list[float]] | None = None,
    tol: float = 1e-6,
) -> list[list[tuple[int, int]]]:
    """Singularities grouped into orbits, to see whether they come in sets.

    A layout cannot be symmetric if the singularities are not, and this says so
    directly -- an orbit of size 1 under a group of order 8 is a singularity
    the field put somewhere no symmetry maps it to.

    ``points`` are the positions to group by; the face centroids by default.
    Pass ``tracer._singularity_points()`` to measure what was actually TRACED
    from, which after :func:`snap_singularities` is not the centroid.
    """
    mesh = field.background.mesh
    sings = field.singularities()
    centroids = dict(points) if points else {f: mesh.face_centroid(f) for f, _ in sings}
    remaining = dict(sings)
    orbits = []
    while remaining:
        fkey = next(iter(remaining))
        index = remaining.pop(fkey)
        orbit = [(fkey, index)]
        for g in symmetry.elements:
            q = symmetry.apply(g, centroids[fkey])
            for other in list(remaining):
                c = centroids[other]
                if abs(c[0] - q[0]) < tol and abs(c[1] - q[1]) < tol:
                    orbit.append((other, remaining.pop(other)))
        orbits.append(orbit)
    return orbits
