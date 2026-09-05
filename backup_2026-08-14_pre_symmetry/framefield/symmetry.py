"""Keep a symmetric domain's symmetry all the way through to the layout.

NOT WIRED IN. Nothing imports this module; ``__init__.py`` does not export it
and no ``ff_0*`` command calls it. It is reached only through
:func:`decomposition`, which is a parallel route to
``FieldDecomposition.from_boundary``. That is deliberate -- everything here
changes the background triangulation, so adopting it moves every row of
``baseline.json``, and that decision is not this module's to take.

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
   arbitrarily -- so :func:`background` keeps it and makes it EQUIVARIANT
   instead.

2. **``from_curves`` is order-dependent.** Its ``claimed`` set gives each vertex
   to the FIRST curve in list order, so a vertex near two arcs and its rotated
   twin take their tangent from different arcs. With a perfectly symmetric mesh
   and perfectly symmetric arcs the constraint set still mismatches by 1.145
   under every non-identity element. :func:`guide_constraints` drops the set.

3. **A singularity is reported as a FACE, and the face is a cruder object than
   the field on it.** After 1 and 2 the field is symmetric to 2e-16 and the
   layout still is not, because the triangle that wins the winding tie around a
   zero is not equivariant -- and a symmetric field puts its zeros exactly ON
   the symmetry axes, where several triangles are equally entitled.
   :func:`snap_singularities` takes the position from the field's own ``|u|``
   minimum instead. The face-index tie in ``wrap_to_period`` is left exactly as
   it is; it stops mattering once nothing reads a position off it.

4. **``build_network`` recomputes those centroids** and puts them back on the
   front of every trace, undoing 3 -- see :class:`snapped_centroids`.

5. **``repair._cluster`` is greedy first-fit in trace order** -- see
   :class:`equivariant_clustering`.

4 and 5 are applied as SCOPED SHIMS around the untouched functions rather than
as edits. At adoption each is a small change in place; until then nothing on
disk changes for any other caller.

WHAT THIS BUYS, measured on the four-arc square at spacing 0.5
--------------------------------------------------------------

::

                              patches  poles  launch pts   coarse rot   mirror
    baseline                      18      2      0%        15% / 0.59    15%
    + symmetric background        16      0     50%        40% / 0.42    32%
    + equivariant guides          16      0     50%        100%/ 0.44    36%
    + field symmetrisation        16      0     50%        84% / 0.44    68%
    + snapped singularities       16      0    100%        92% / 0.35    84%
    + network shims               16      0    100%       100% / 0.0000  100%

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

The rectangle row is the one to look at before adopting: the original found 6
singularities there and this route finds none. Both are defensible -- two arcs
on a rectangle force no winding the walls cannot absorb, so the six were
artefacts of the asymmetric background -- but it is a real change in what the
layout does, not only in how symmetric it looks. Its minimum angle goes 12.54
-> 78.54 degrees.
"""
from math import ceil

from compas.geometry import cross_vectors
from compas.geometry import delaunay_triangulation
from compas.geometry import distance_point_point
from compas.geometry import is_point_in_polygon_xy
from compas.geometry import length_vector
from compas.geometry import normalize_vector
from compas.geometry import subtract_vectors
from compas.itertools import pairwise

from compas_singular.datastructures import Mesh

from .background import BackgroundMesh
from .background import _as_open_loop
from .background import _densify_loop
from .background import _distance_to_loop
from .background import _jitter
from .constraints import Constraint
from .constraints import from_boundary
from .constraints import representation
from .field import CrossField
from .trace import Tracer


__all__ = ['Symmetry', 'background', 'guide_constraints', 'symmetrise',
           'snap_singularities', 'SymmetricTracer', 'decomposition',
           'invariance', 'field_invariance']


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


class Symmetry(object):
    """A subgroup of the square's symmetries, about a centre.

    Attributes
    ----------
    centre : [x, y, z]
    elements : list[tuple]
        A subset of :data:`ELEMENTS`, always containing the identity.

    Notes
    -----
    Only the 8 symmetries of the square are considered. That is not a
    limitation of the method -- everything below works for any finite group of
    isometries -- but it covers every symmetric domain this front end has been
    put in front of, and detecting a 6-fold symmetry nobody drew would be a
    liability rather than a feature.
    """

    def __init__(self, centre, elements):
        self.centre = [float(centre[0]), float(centre[1]), 0.0]
        self.elements = list(elements)

    def __len__(self):
        return len(self.elements)

    def __repr__(self):
        return '<Symmetry order {} about ({:.3f}, {:.3f}): {}>'.format(
            len(self.elements), self.centre[0], self.centre[1],
            ', '.join(self.names()))

    def names(self):
        return [g[5] for g in self.elements]

    @property
    def trivial(self):
        return len(self.elements) < 2

    def apply(self, g, point):
        """``g`` applied to a point, about :attr:`centre`."""
        a, b, c, d = g[0], g[1], g[2], g[3]
        x, y = point[0] - self.centre[0], point[1] - self.centre[1]
        return [a * x + b * y + self.centre[0],
                c * x + d * y + self.centre[1], 0.0]

    @staticmethod
    def linear(g, vector):
        """``g``'s linear part applied to a free vector -- no centre involved."""
        a, b, c, d = g[0], g[1], g[2], g[3]
        return [a * vector[0] + b * vector[1], c * vector[0] + d * vector[1], 0.0]

    # ------------------------------------------------------------------
    # detection
    # ------------------------------------------------------------------

    @classmethod
    def detect(cls, loops, centre=None, tol=1e-6):
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
# 1 -- a background whose triangulation is invariant under the group
# ----------------------------------------------------------------------

def background(outer_boundary, inner_boundaries=None, target_length=None,
               symmetry=None, margin=0.45, on_axis=True):
    """``BackgroundMesh.from_boundary``, with an interior grid the group fixes.

    Same signature and same return type as the original, plus ``symmetry``.
    With ``symmetry`` trivial this still differs from the original -- the grid
    is centred rather than anchored at ``min(xs)`` -- so it is not a drop-in
    for reproducing an existing baseline.

    Parameters
    ----------
    symmetry : Symmetry, optional
        Detected from the boundaries when omitted. Guides are NOT detected from
        here because this function never sees them; pass a ``Symmetry`` built
        from everything, or use :func:`decomposition`.
    on_axis : bool, optional
        Whether grid lines pass THROUGH the centre (default) or cells straddle
        it. Measured: on-axis gives 16 clean quads on the four-arc square,
        straddling gives 23 patches and 2 poles. Kept as a parameter only
        because that measurement is one domain deep.

    Notes
    -----
    **The jitter is kept and made equivariant, not removed.** ``background.py``
    perturbs the grid because a regular one makes every cell's four corners
    cocircular, at which point Qhull's choice of diagonal is arbitrary and not
    reproducible across platforms. That reasoning is sound and survives here;
    what changes is only that the perturbation is now a function of the point's
    ORBIT rather than of its grid indices.

    The construction, which is why it is exact rather than nearly exact:

    * the lattice is indexed by integers, and every element of the group is a
      signed permutation, so the group acts on the INDICES exactly -- there is
      no floating-point tie to lose;
    * each orbit gets one jitter vector, generated from its canonical
      representative's indices by the original ``_jitter``;
    * that vector is averaged over the representative's STABILISER, which
      projects it onto the subspace the stabiliser fixes. A point at the centre
      gets no jitter, a point on a mirror axis may only slide ALONG the axis,
      and a point on a diagonal only along the diagonal. Without this the
      construction is ill-defined, because two group elements can map the
      representative to the same point and would otherwise disagree about where
      to put it;
    * the inside/near-the-wall filters are applied to the representative and
      the whole orbit is kept or dropped together, so a point sitting exactly
      ``margin * target_length`` from the wall cannot be admitted on one side
      and rejected on its mirror.
    """
    outer = _as_open_loop(outer_boundary)
    inners = [_as_open_loop(loop) for loop in (inner_boundaries or [])]

    xs = [p[0] for p in outer]
    ys = [p[1] for p in outer]
    if target_length is None:
        diagonal = ((max(xs) - min(xs)) ** 2 + (max(ys) - min(ys)) ** 2) ** 0.5
        target_length = diagonal / 25.0
    if symmetry is None:
        symmetry = Symmetry.detect([outer] + inners)

    outer = _densify_loop(outer, target_length)
    inners = [_densify_loop(loop, target_length) for loop in inners]

    points = list(outer)
    for loop in inners:
        points.extend(loop)
    points.extend(_equivariant_grid(symmetry, target_length, outer, inners,
                                    margin, on_axis))

    return BackgroundMesh(_triangulate(points, outer, inners), outer, inners,
                          target_length)


def _equivariant_grid(symmetry, target_length, outer, inners, margin, on_axis):
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

    def point_of(u, v, jitter):
        return [centre[0] + u * half + jitter[0],
                centre[1] + v * half + jitter[1], 0.0]

    def act(g, u, v):
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


def _inside(p, outer, inners, limit):
    if not is_point_in_polygon_xy(p, outer):
        return False
    if _distance_to_loop(p, outer) < limit:
        return False
    for loop in inners:
        if is_point_in_polygon_xy(p, loop):
            return False
        if _distance_to_loop(p, loop) < limit:
            return False
    return True


def _triangulate(points, outer, inners):
    """Delaunay, then the same cleanup ``BackgroundMesh.from_boundary`` does.

    Duplicated rather than shared because the original interleaves the cleanup
    with its own point generation, and this module's whole point is to be
    removable without having touched anything.
    """
    mesh = Mesh.from_vertices_and_faces(points, delaunay_triangulation(points))
    for fkey in list(mesh.faces()):
        a, b, c = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
        if length_vector(cross_vectors(subtract_vectors(b, a),
                                       subtract_vectors(c, a))) < 1e-12:
            mesh.delete_face(fkey)
            continue
        centroid = mesh.face_centroid(fkey)
        if not is_point_in_polygon_xy(centroid, outer):
            mesh.delete_face(fkey)
        elif any(is_point_in_polygon_xy(centroid, loop) for loop in inners):
            mesh.delete_face(fkey)
    for vkey in list(mesh.vertices()):
        if not mesh.vertex_faces(vkey):
            mesh.delete_vertex(vkey)
    return mesh


# ----------------------------------------------------------------------
# 2 -- guide constraints that do not depend on curve order
# ----------------------------------------------------------------------

def guide_constraints(bg, curves, mode='tangent', band=None, weight=1.0):
    """``from_curves`` without the first-come ``claimed`` set.

    ``from_curves`` claims each vertex for the first curve that reaches it, so
    a vertex within ``band`` of two guides takes its direction from whichever
    guide the caller listed first -- and its rotated twin can take a different
    one. Measured on a perfectly symmetric mesh with four perfectly symmetric
    arcs: constraint mismatch 1.145 under EVERY non-identity element, which is
    the full magnitude of a representation, i.e. a different constraint
    entirely.

    Here every guide within ``band`` contributes its own constraint and the
    solve accumulates them per vertex -- which ``field.solve`` already does for
    soft constraints, and which its own comment calls "the right behaviour
    where two guides overlap".

    HARD constraints are merged here instead, because ``field.solve`` puts them
    in a dict and the last one written wins -- order-dependent again, and not
    fixable from the caller's side. Two hard guides at one vertex are averaged
    in the 4th-power representation, which is the same averaging
    ``from_boundary`` uses at a corner.
    """
    if curves and not isinstance(curves[0][0], (list, tuple)):
        curves = [curves]

    radius = bg.target_length if band is None else band
    boundary = bg.boundary_vertices()

    out = []
    for curve in curves:
        pts = [[float(p[0]), float(p[1]), 0.0] for p in curve]
        for vkey in bg.mesh.vertices():
            if vkey in boundary:
                continue
            p = bg.mesh.vertex_coordinates(vkey)
            best, tangent = float('inf'), None
            for a, b in pairwise(pts):
                ab = subtract_vectors(b, a)
                length2 = ab[0] * ab[0] + ab[1] * ab[1]
                if length2 == 0.0:
                    continue
                t = ((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1]) / length2
                t = max(0.0, min(1.0, t))
                q = [a[0] + ab[0] * t, a[1] + ab[1] * t, 0.0]
                d = distance_point_point(p, q)
                if d < best:
                    best, tangent = d, normalize_vector(ab)
            if tangent is None or best > radius:
                continue
            direction = tangent
            if mode == 'perpendicular':
                direction = [-tangent[1], tangent[0], 0.0]
            elif mode != 'tangent':
                raise ValueError("mode must be 'perpendicular' or 'tangent'")
            out.append(Constraint(vkey, direction, weight))

    return _merge_hard(out) if weight is None else out


def _merge_hard(constraints):
    from cmath import phase as _phase
    from math import cos, sin

    from .constraints import PERIOD                       # noqa: F401

    grouped = {}
    for c in constraints:
        grouped.setdefault(c.vkey, []).append(c)
    out = []
    for vkey, group in grouped.items():
        if len(group) == 1:
            out.append(group[0])
            continue
        acc = sum((representation(c.direction) for c in group), 0j) / len(group)
        if abs(acc) < 1e-9:
            continue                    # they cancel: leave the vertex free
        theta = _phase(acc) / 4.0
        out.append(Constraint(vkey, [cos(theta), sin(theta), 0.0], None))
    return out


# ----------------------------------------------------------------------
# 3 -- project the solved field onto the symmetric subspace
# ----------------------------------------------------------------------

def symmetrise(field, symmetry, tol=1e-9):
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


def _snap(value, tol):
    return int(round(value / tol))


# ----------------------------------------------------------------------
# 4 -- take the singularity off the arbitrary face it was reported on
# ----------------------------------------------------------------------

def snap_singularities(field, relative_tol=1e-9):
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


class SymmetricTracer(Tracer):
    """:class:`Tracer` that launches from given points, not from face centroids.

    The base class computes the launch centre as ``face_centroid(fkey)`` in one
    place, :meth:`_singularity_points`, and everything else -- arm finding,
    integration, the self-hit guard, the proximity snapping in ``repair`` --
    consumes that dict. So this is the whole subclass.
    """

    def __init__(self, field, points, **kwargs):
        super(SymmetricTracer, self).__init__(field, **kwargs)
        self._points = dict(points)

    def _singularity_points(self):
        return {fkey: self._points.get(fkey, self.mesh.face_centroid(fkey))
                for fkey, _ in self.field.singularities()}


# ----------------------------------------------------------------------
# 5 -- the end-snapping in repair.build_network
# ----------------------------------------------------------------------

class equivariant_clustering(object):
    """Context manager: make ``repair._cluster`` order-independent, in scope.

    A SHIM, not an edit. ``repair._cluster`` is replaced for the duration of
    one call and put back afterwards, so the module on disk is untouched and
    every other caller sees the original. At adoption time this becomes six
    lines inside ``_cluster`` itself and the context manager goes away.

    WHY IT IS NEEDED, AND WHY IT IS LAST
    ------------------------------------

    With the first four fixes in place, everything the tracer produces is
    exactly symmetric -- measured on the four-arc square, all 472 points of all
    32 separatrices, deviation 0.0000, and the 32 arc lengths falling into
    exact groups of four and eight. The layout built from them is not: the
    polyline endpoints drop to 85%, deviation 0.3486.

    The whole of that gap is here. ``_cluster`` is greedy first-fit -- each
    point joins the FIRST representative within tolerance, scanning in input
    order -- and the input order is trace order, which no symmetry respects. So
    two points that are exact mirror images can join different clusters, and a
    cluster's representative (deliberately a member, not a centroid, so that a
    landing stays exactly on the wall it was clipped to) is whichever member
    happened to come first.

    THE WRAPPER DOES NOT REIMPLEMENT THE CLUSTERING. It runs the original,
    unmodified, and then projects the result onto the group: if input point
    ``m`` is ``g`` applied to input point ``c``, then ``m``'s representative is
    set to ``g`` applied to ``c``'s representative. Every subtlety the original
    docstring defends -- the ``groups`` rule that stops one trace's two ends
    merging, the tighter ``wall_tol`` between two landings -- is decided by the
    original code and only then made consistent across the orbit.

    The projected representative is still exactly on the wall when the original
    was: it is the image of a boundary point under a symmetry of that same
    boundary.

    With a trivial group, or on points with no orbit structure, this is the
    identity -- every orbit is a singleton and every representative is its own.
    """

    def __init__(self, symmetry, tol=1e-6):
        self.symmetry = symmetry
        self.tol = tol
        self._original = None

    def __enter__(self):
        from . import repair

        self._original = repair._cluster
        symmetry, tol = self.symmetry, self.tol
        original = self._original

        def clustered(points, cluster_tol, **kwargs):
            base = original(points, cluster_tol, **kwargs)
            return _project_clusters(points, base, symmetry, tol)

        repair._cluster = clustered
        return self

    def __exit__(self, *exc):
        from . import repair

        repair._cluster = self._original
        return False


def _project_clusters(points, base, symmetry, tol):
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
    def key(p):
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


class snapped_centroids(object):
    """Context manager: let ``build_network`` see the SNAPPED singularity points.

    A SHIM for the same reason :class:`equivariant_clustering` is one, and the
    more important of the two.

    :func:`snap_singularities` moves the launch point off the face centroid, and
    ``repair.build_network`` then puts the centroid straight back::

        centroids = {fkey: mesh.face_centroid(fkey) for fkey, _ in
                     field.singularities()}
        ...
        if s.source in centroids:
            pts.insert(0, list(centroids[s.source]))

    That insert is not optional -- every arm of one singularity has to share one
    exact point or the patches around it do not close -- but the point it
    inserts is recomputed rather than taken from the tracer. So the traces go in
    exactly symmetric (measured: all 472 points, deviation 0.0000) and come out
    of that loop at 69%, deviation 0.3486, which is the whole of the residual
    asymmetry in the final layout.

    At adoption this is one line in ``repair.build_network`` -- take the points
    from the tracer instead of recomputing them -- and this class goes away.
    Until then the override is scoped to the ``build_network`` call alone, not
    to the whole layout build, because ``face_centroid`` on the background mesh
    is a general-purpose query and ``densify`` uses it for something else.
    """

    def __init__(self, points):
        self.points = dict(points)
        self._original = None

    def __enter__(self):
        from . import decomposition as module

        self._original = module.build_network
        original = self._original
        points = self.points

        def patched(field, separatrices, **kwargs):
            mesh = field.background.mesh
            saved = mesh.__dict__.get('face_centroid')

            def face_centroid(fkey):
                if fkey in points:
                    return list(points[fkey])
                return type(mesh).face_centroid(mesh, fkey)

            mesh.face_centroid = face_centroid
            try:
                return original(field, separatrices, **kwargs)
            finally:
                if saved is None:
                    mesh.__dict__.pop('face_centroid', None)
                else:
                    mesh.face_centroid = saved

        module.build_network = patched
        return self

    def __exit__(self, *exc):
        from . import decomposition as module

        module.build_network = self._original
        return False


def _inverse(g):
    """The inverse of a signed permutation -- its transpose, both being orthogonal."""
    for other in ELEMENTS:
        if (other[0] * g[0] + other[1] * g[2] == 1
                and other[0] * g[1] + other[1] * g[3] == 0
                and other[2] * g[0] + other[3] * g[2] == 0
                and other[2] * g[1] + other[3] * g[3] == 1):
            return other
    raise ValueError('not an element of the group: {}'.format(g))


# ----------------------------------------------------------------------
# the whole route, as a parallel to FieldDecomposition.from_boundary
# ----------------------------------------------------------------------

def decomposition(outer_boundary, inner_boundaries=None, guides=None,
                  mode='perpendicular', target_length=None, guide_weight=1.0,
                  guide_band=None, relax=False, field_tau=None,
                  symmetry='auto',
                  steps=('background', 'guides', 'field', 'singularities',
                         'network'),
                  on_axis=True):
    """``FieldDecomposition.from_boundary`` with the symmetric pieces swapped in.

    Returns the same :class:`FieldDecomposition`, so everything downstream --
    ``decomposition_mesh``, ``quad_mesh``, the strip grammar, the Rhino steps --
    is unchanged and unaware.

    Parameters
    ----------
    symmetry : Symmetry or 'auto' or None, optional
        ``'auto'`` detects it from the boundaries AND the guides together.
        ``None`` disables every step, which makes this the original route with
        a centred grid.
    steps : tuple[str], optional
        Which of the four fixes to apply, for measuring them separately:
        ``'background'``, ``'guides'``, ``'field'``, ``'singularities'``,
        ``'network'``. Their effects are cumulative: ``'singularities'`` snaps
        to the minimum of a field it assumes is already symmetric, and
        ``'network'`` projects a clustering whose input it assumes is already
        symmetric.

    Notes
    -----
    ``relax=True`` and symmetrisation are independent and compose, but note the
    order: the relaxation is a NONLINEAR flow, so unlike the default solve it
    has no uniqueness guarantee and may genuinely prefer an asymmetric steady
    state. Symmetrising afterwards projects that back, which is a claim about
    what you want rather than about what the flow found. The default solve
    needs no such caveat: it is a convex problem with a unique minimiser, so a
    symmetric discretisation gives a symmetric field -- measured at 1e-15 on
    walls alone.
    """
    from .decomposition import FieldDecomposition

    loops = [outer_boundary] + list(inner_boundaries or []) + list(guides or [])
    if symmetry == 'auto':
        symmetry = Symmetry.detect(loops)
    if symmetry is None or symmetry.trivial:
        steps = ()

    if 'background' in steps:
        bg = background(outer_boundary, inner_boundaries, target_length,
                        symmetry=symmetry, on_axis=on_axis)
    else:
        bg = BackgroundMesh.from_boundary(outer_boundary, inner_boundaries,
                                          target_length=target_length)

    constraints = from_boundary(bg)
    if guides:
        if 'guides' in steps:
            constraints += guide_constraints(bg, guides, mode=mode,
                                             weight=guide_weight, band=guide_band)
        else:
            from .constraints import from_curves
            constraints += from_curves(bg, guides, mode=mode,
                                       weight=guide_weight, band=guide_band)

    field = CrossField.solve(bg, constraints, relax=relax, tau=field_tau)
    if 'field' in steps:
        field = symmetrise(field, symmetry)

    snap_report = {}
    points = {}
    if 'singularities' in steps:
        field, points, snap_report = snap_singularities(field)
        tracer = SymmetricTracer(field, points)
    else:
        tracer = Tracer(field)

    separatrices, report = tracer.separatrices()
    out = FieldDecomposition(bg, field, separatrices, tracer)
    out.guides = list(guides or [])
    out.trace_report = report
    out.symmetry = symmetry
    out.snap_report = snap_report
    if 'network' in steps:
        _wrap_network(out, symmetry, points)
    return out


def _wrap_network(decomp, symmetry, points):
    """Make this decomposition build its network under :class:`equivariant_clustering`.

    The clustering happens inside ``decomposition_polylines`` and
    ``decomposition_mesh``, which run LATER than this function -- the
    decomposition is returned first and the layout is asked for afterwards. So
    the shim is attached to the three public entry points of this ONE object,
    rather than held open around a call that has not happened yet. Instance
    attributes shadow the class, so internal ``self.decomposition_mesh(...)``
    calls go through the wrapper too.
    """
    def wrap(name):
        bound = getattr(decomp, name)

        def wrapped(*args, **kwargs):
            with equivariant_clustering(symmetry), snapped_centroids(points):
                return bound(*args, **kwargs)
        return wrapped

    for name in ('decomposition_polylines', 'decomposition_mesh', 'quad_mesh'):
        setattr(decomp, name, wrap(name))


# ----------------------------------------------------------------------
# measurement
# ----------------------------------------------------------------------

def invariance(points, symmetry, tol=1e-6):
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


def field_invariance(field, symmetry, tol=1e-9):
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


def singularity_orbits(field, symmetry, points=None, tol=1e-6):
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
