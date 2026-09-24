"""Patch patterns for coarse quad mesh densification.

A pattern is a *template* on the unit square plus a *morph* onto a coarse face.
That split is what ``discrete_coons_patch`` already does internally for the
plain grid; naming it lets every other pattern reuse the same machinery.

    template            morph                 result
    (u, w) in [0,1]^2   4 boundary polylines  one dense patch

Templates
---------
A template is a pure list of ``(u, w)`` parameter pairs plus faces indexing
into it. It knows nothing about geometry -- no coordinates, no curves, no
poles. ``u`` runs along the a->b side of the patch, ``w`` along the a->d side,
following the orientation ``compas.geometry.discrete_coons_patch`` documents::

    b -----> c          w
    ^        ^          ^
    |        |          |
    |        |          |
    a -----> d          +---> u

Morph
-----
``pattern_morph`` evaluates the transfinite (Coons) map at each template
parameter::

    P(u, w) =   (1-w) AB(u) + w DC(u)          # blend of the two u-curves
              + (1-u) AD(w) + u BC(w)          # blend of the two w-curves
              - bilinear(A, B, C, D)           # minus the double-counted hypar

where ``AB(u)`` evaluates the *already densified* side polyline at ``u``. The
side polylines carry whatever ``edges_to_curves`` supplied, so a curved coarse
edge stays curved for every pattern, not just the grid.

Two properties make this safe to share:

* At a grid node the map reduces to exactly the arithmetic
  ``discrete_coons_patch`` performs, so ``ortho`` reproduces the existing
  densification bit for bit.
* A patch corner that collapsed to a pole is just a side polyline that is
  ``None``; the fill below turns it into a repeated corner point and the map
  degenerates on its own. For the GRID templates -- ``ortho`` and
  ``diagonal`` -- poles cost no pattern-specific code: the collapse is one row
  of the grid.

The fan is the exception. Its polar centres are the patch corners, and a
triangular patch has three of them, not four with two on top of each other.
It gets its own template and its own morph -- ``_pattern_fan_triangle`` and
``pattern_morph_triangle`` -- and ``densify`` picks them for any ``fan`` patch
with a collapsed side.

Why divisions are not free
--------------------------
Adjacent patches weld only if they put the *same* points on their shared edge,
so the division counts come from the strip densities, not from the pattern.
``ortho`` accepts any ``(nu, nw)``. ``diagonal`` and ``fan`` need ``nu == nw``,
and even -- see ``_pattern_diagonal`` and ``_pattern_fan`` -- and
``reconcile_strip_densities`` is what makes a layout legal for them before any
patch is built.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import annotations

from typing import Any
from typing import Sequence
from typing import TYPE_CHECKING

from compas.itertools import pairwise

if TYPE_CHECKING:
    from compas_singular.datastructures import CoarseQuadMesh


__all__ = [
    'PATTERNS',
    'create_pattern',
    'pattern_morph',
    'pattern_morph_triangle',
    'patch_divisions',
    'reconcile_strip_densities',
]


# ==============================================================================
# create patterns
# ==============================================================================


PATTERNS = [
    'ortho',
    'diagonal',
    'fan',
]


def create_pattern(type: str, nu: int, nw: int) -> tuple[list[tuple[float, float]], list[list[int]]]:
    """Build the template of a pattern at a given number of divisions.

    Parameters
    ----------
    kind : str
        A key of :attr:`PATTERNS`.
    nu, nw : int
        Divisions along u (the a->b side) and along w (the a->d side).

    Returns
    -------
    list[(float, float)]
        Template vertices as parameter pairs on the unit square.
    list[list[int]]
        Template faces, indexing into those.

    Raises
    ------
    ValueError
        If ``kind`` is not a known pattern, or the divisions do not suit it.
    """
    PATTERNS = {
        'ortho': _pattern_ortho,
        'diagonal': _pattern_diagonal,
        'fan': _pattern_fan,
        # internal: what 'fan' becomes on a triangular patch. Its parameters
        # are barycentric, so it pairs with pattern_morph_triangle, not
        # pattern_morph -- which is why it is not in the public PATTERNS.
        'fan_triangle': _pattern_fan_triangle,
    }

    if type not in PATTERNS:
        raise ValueError('unknown pattern {!r} -- pick one of {}'.format(type, sorted(PATTERNS)))
    return PATTERNS[type](max(1, nu), max(1, nw))

# ==============================================================================
# patterns
# ==============================================================================

def _pattern_ortho(nu: int, nw: int) -> tuple[list[tuple[float, float]], list[list[int]]]:
    """The plain grid: ``nu`` by ``nw`` quads on the unit square.

    Vertex order (``u`` major, ``w`` minor) and face winding are copied from
    ``discrete_coons_patch`` so that morphing this template is not merely
    equivalent to the old densification but identical to it.

    Parameters
    ----------
    nu, nw : int
        Divisions along u and along w.

    Returns
    -------
    list[(float, float)], list[list[int]]
    """
    uw = [(i / nu, j / nw) for i in range(nu + 1) for j in range(nw + 1)]

    m = nw + 1
    faces = [
        [i * m + j, i * m + j + 1, (i + 1) * m + j + 1, (i + 1) * m + j]
        for i in range(nu)
        for j in range(nw)
    ]
    return uw, faces


def _pattern_diagonal(nu: int, nw: int) -> tuple[list[tuple[float, float]], list[list[int]]]:
    """The grid with both diagonals of the patch drawn through it.

    Each diagonal is a chain of cell diagonals -- cells ``(i, i)`` for a->c,
    ``(i, n-1-i)`` for b->d -- and every cell on it is cut into two triangles
    along it.

    Cutting a cell does not move any vertex, so the boundary of this template
    is the plain grid's boundary: a diagonal patch welds against an ortho one.

    Why ``nu == nw``, and even
    --------------------------
    Both are about the diagonal staying one unbroken line:

    * **Square.** A cell is the patch scaled by ``1/nu`` along u and ``1/nw``
      along w, so it has the patch's own height/width ratio only when
      ``nu == nw``. Only then is a cell diagonal a piece of the patch
      diagonal, and do consecutive ones share a vertex. Otherwise the cells
      the diagonal crosses form a staircase whose cuts do not touch, and
      NEITHER diagonal is connected -- measured on 4x2, 6x3 and 8x5.
    * **Even.** The two diagonals cross at the patch centre. For ``nu`` even
      that is a grid vertex; for ``nu`` odd it is the middle of a cell, which
      can be cut along only one of them, so b->d breaks there -- measured on
      3x3, 5x5 and 7x7.

    ``nu`` and ``nw`` are the densities of the two strips crossing the patch,
    so both are constraints on the layout, not on the patch.
    ``reconcile_strip_densities`` enforces them.

    Parameters
    ----------
    nu, nw : int
        Divisions along u and along w. Must be equal and even.

    Returns
    -------
    list[(float, float)], list[list[int]]

    Raises
    ------
    ValueError
        If ``nu`` and ``nw`` are not equal, or not even.
    """
    if nu != nw or nu % 2:
        raise ValueError(
            'the diagonal pattern needs the same even number of divisions on '
            'all four sides of a patch, got {} by {} -- run '
            'reconcile_strip_densities() first'.format(nu, nw)
        )

    uw, quads = _pattern_ortho(nu, nw)

    n = nu
    faces = []
    for idx, (i, j) in enumerate((i, j) for i in range(n) for j in range(n)):
        # v0 = (i, j), v1 = (i, j+1), v2 = (i+1, j+1), v3 = (i+1, j)
        v0, v1, v2, v3 = quads[idx]
        if i == j:
            # a -> c, cut v0-v2
            faces.append([v0, v1, v2])
            faces.append([v0, v2, v3])
        elif i + j == n - 1:
            # b -> d, cut v1-v3; never the same cell as a -> c, since n is even
            faces.append([v1, v2, v3])
            faces.append([v1, v3, v0])
        else:
            faces.append(quads[idx])
    return uw, faces


# the four quadrants of a fan, as (pole, first mid-edge, second mid-edge) in
# (u, w). Ordered so the angular sweep runs the same way around every pole,
# which keeps the face winding consistent with the grid templates.
_FAN_QUADRANTS = [
    ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5)),   # at a, between sides ab and ad
    ((1.0, 0.0), (1.0, 0.5), (0.5, 0.0)),   # at b, between sides bc and ab
    ((1.0, 1.0), (0.5, 1.0), (1.0, 0.5)),   # at c, between sides dc and bc
    ((0.0, 1.0), (0.0, 0.5), (0.5, 1.0)),   # at d, between sides ad and dc
]

_FAN_CENTRE = (0.5, 0.5)


def _pattern_fan(nu: int, nw: int) -> tuple[list[tuple[float, float]], list[list[int]]]:
    """Four polar fans, one per corner, meeting at the patch centre.

    Each quadrant is a polar grid running from its corner (the pole) out to the
    L-shaped path ``mid-edge -> centre -> mid-edge``. ``rings`` counts the
    radial steps, ``rays`` the angular ones. The innermost ring sits on the
    pole, so its quads collapse to the triangles that give the pattern its
    look; ``pattern_morph`` does that collapse.

    Why ``nu == nw``, and even
    --------------------------
    Two constraints, both about points having to coincide:

    * **Rays.** The leg from a mid-edge to the centre belongs to two quadrants
      at once, walked from opposite ends. Each puts half its rays on it, so
      neighbouring quadrants need equal ray counts -- which chains around all
      four.
    * **Rings.** A quadrant's first ray lies *on* a patch side, so it lays down
      that side's near half. The far half comes from the next corner's
      quadrant. A side therefore carries ``rings(corner) + rings(corner)``
      divisions, which makes the ring count a property of the coarse *vertex*,
      not the face -- and the neighbouring face across that side has to pick
      the same number for it.

    Solving that in general is a constraint system over the whole layout
    (``density(strip) == rings(u) + rings(v)`` for every edge of every strip),
    with no solution for some layouts. The uniform solution always exists and
    is the one taken here: one ring count per patch, so its two strips carry
    ``2 * rings`` each and the patch is ``d`` by ``d`` with ``d`` even.
    ``reconcile_strip_densities`` enforces it, and since a strip crosses many
    patches, on a layout that is all fan it comes out as one density
    everywhere.

    Parameters
    ----------
    nu, nw : int
        Divisions along u and along w. Must be equal and even.

    Returns
    -------
    list[(float, float)], list[list[int]]

    Raises
    ------
    ValueError
        If ``nu`` and ``nw`` are not equal, or not even.
    """
    if nu != nw or nu % 2:
        raise ValueError(
            'the fan pattern needs the same even number of divisions on all '
            'four sides of a patch, got {} by {} -- run '
            'reconcile_strip_densities() first'.format(nu, nw)
        )

    rings = nu // 2      # radial steps, pole -> L
    rays = nu            # angular steps across the whole L
    half = rays // 2     # angular steps per leg of the L

    uw = []
    index = {}

    def add(p: tuple[float, float]) -> int:
        # one vertex per parameter pair, so the four quadrants share their
        # seams, their poles and the centre instead of duplicating them
        key = (round(p[0], 12), round(p[1], 12))
        if key not in index:
            index[key] = len(uw)
            uw.append((p[0], p[1]))
        return index[key]

    faces = []
    for pole, m1, m2 in _FAN_QUADRANTS:
        grid = {}
        for i in range(rays + 1):
            # angular: walk the far boundary m1 -> centre -> m2, half the rays
            # on each leg, so the sweep is at constant speed
            if i <= half:
                far = _lerp(m1, _FAN_CENTRE, i, half)
            else:
                far = _lerp(_FAN_CENTRE, m2, i - half, half)
            for j in range(rings + 1):
                # radial: slide from the pole out to the far boundary
                grid[i, j] = add(_lerp(pole, far, j, rings))

        for i in range(rays):
            for j in range(rings):
                faces.append([grid[i, j], grid[i + 1, j], grid[i + 1, j + 1], grid[i, j + 1]])

    return uw, faces


# the three corner fans of a triangular patch, as (corner, first mid-edge,
# second mid-edge) in barycentric coordinates on the corners (V0, V1, V2),
# counter-clockwise, V0 the pole. A cyclic shift of the corners maps the table
# onto itself -- that is what makes every corner, the pole included, the same.
_FAN_TRIANGLE_REGIONS = [
    ((1.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5)),   # at V0, between sides V0V1 and V2V0
    ((0.0, 1.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.5, 0.0)),   # at V1, between sides V1V2 and V0V1
    ((0.0, 0.0, 1.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)),   # at V2, between sides V2V0 and V1V2
]

_FAN_TRIANGLE_CENTRE = (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)


def _pattern_fan_triangle(nu: int, nw: int) -> tuple[list[tuple[float, float, float]], list[list[int]]]:
    """Three polar fans, one per corner, meeting at the centre of a triangle.

    The fan of a triangular patch (a pseudo-quad). Each corner owns the
    quadrilateral ``corner -> mid-edge -> centre -> mid-edge`` and fills it with
    the polar grid of :func:`_pattern_fan`, with the same ``rings`` and
    ``rays``, so a side still carries ``2 * rings`` divisions and the patch
    welds against any neighbour.

    The parameters are barycentric ``(l0, l1, l2)`` on the corners, not
    ``(u, w)``, and :func:`pattern_morph_triangle` maps them. Two things were
    measured that rule out reusing the quad machinery:

    * **The quad fan on a zero-length side.** The Coons morph collapses the
      WHOLE side ``w = 0`` onto the pole, so near it ``w`` is the radius and
      ``u`` the angle. The quad fan centres its polar grids on the parameter
      POINTS (0, 0) and (1, 0) instead, so its rings land on lines through the
      pole: d - 2 zero-area faces survive the weld, and the pole gets valency
      3d - 3 with edges lying on top of each other, against d + 1 for ortho.
    * **Three fans laid out in (u, w).** The map from (u, w) to the triangle is
      bilinear (``l_pole = 1 - w, l_c = w u, l_d = w (1 - u)``), so the
      pole's spokes come out straight and the other two corners' curved: 0.08
      off its own 120 degree rotation on an equilateral triangle of side 2.

    In barycentric coordinates the construction is the same at every corner.
    The one asymmetry is bookkeeping: V0's innermost ring gets a vertex per ray
    rather than one shared vertex, so each innermost cell keeps two coincident
    corners. The weld collapses them and ``densify`` registers the pole faces,
    exactly as it does for a grid template. V1 and V2 share one vertex, as the
    corners of :func:`_pattern_fan` do. The geometry is identical either way.

    Parameters
    ----------
    nu, nw : int
        Divisions along the sides. Must be equal and even, as for
        :func:`_pattern_fan`.

    Returns
    -------
    list[(float, float, float)], list[list[int]]
        Template vertices as barycentric coordinates on (pole, V1, V2), and
        template faces.

    Raises
    ------
    ValueError
        If ``nu`` and ``nw`` are not equal, or not even.
    """
    if nu != nw or nu % 2:
        raise ValueError(
            'the fan pattern needs the same even number of divisions on all '
            'three sides of a triangular patch, got {} by {} -- run '
            'reconcile_strip_densities() first'.format(nu, nw)
        )

    rings = nu // 2      # radial steps, corner -> L
    rays = nu            # angular steps across the whole L
    half = rays // 2     # angular steps per leg of the L

    lam = []
    index = {}

    def add(key: Any, p: tuple[float, float, float]) -> int:
        # one vertex per key, so the three corners share their seams and the
        # centre -- and V1 and V2 their innermost ring -- instead of
        # duplicating them
        if key not in index:
            index[key] = len(lam)
            lam.append(tuple(p))
        return index[key]

    faces = []
    for n, (corner, m1, m2) in enumerate(_FAN_TRIANGLE_REGIONS):
        grid = {}
        for i in range(rays + 1):
            if i <= half:
                far = _lerp(m1, _FAN_TRIANGLE_CENTRE, i, half)
            else:
                far = _lerp(_FAN_TRIANGLE_CENTRE, m2, i - half, half)
            for j in range(rings + 1):
                p = _lerp(corner, far, j, rings)
                if n == 0 and j == 0:
                    key = ('pole', i)
                else:
                    key = tuple(round(x, 12) for x in p)
                grid[i, j] = add(key, p)

        for i in range(rays):
            for j in range(rings):
                faces.append([grid[i, j], grid[i + 1, j], grid[i + 1, j + 1], grid[i, j + 1]])

    return lam, faces


def _lerp(p: Sequence[float], q: Sequence[float], k: int, n: int) -> list[float]:
    """Point ``k/n`` of the way from ``p`` to ``q``.

    Written as a weighted average of the two ends -- ``(p(n-k) + qk)/n`` --
    rather than the usual ``p + t(q - p)``, for two reasons:

    * the ends come back exactly (``k == 0`` and ``k == n`` are special-cased),
    * the expression is symmetric, so walking a segment forwards from ``p`` and
      backwards from ``q`` produces bit-identical points.

    The second one is not pedantry: the fan visits every seam of a patch twice,
    once from each of the two quadrants that share it, and visits every patch
    side twice more from the neighbouring face. Symmetric arithmetic means
    those visits agree exactly instead of within a welding tolerance.

    Parameters
    ----------
    p, q : sequence[float]
        The two ends.
    k, n : int
        Step index and step count.

    Returns
    -------
    list[float]
    """
    if k == 0:
        return [float(a) for a in p]
    if k == n:
        return [float(b) for b in q]
    return [(a * (n - k) + b * k) / n for a, b in zip(p, q)]


# ==============================================================================
# morph
# ==============================================================================

def pattern_morph(uw: list[tuple[float, float]], faces: list[list[int]], sides: list[list[list[float]] | None]) -> tuple[list[list[float]], list[list[int]]]:
    """Map a template onto a coarse face through its four side polylines.

    Parameters
    ----------
    uw : list[(float, float)]
        Template vertices, as parameters on the unit square.
    faces : list[list[int]]
        Template faces.
    sides : list
        The four side polylines ``[ab, bc, dc, ad]``, each a list of points, in
        the orientation ``discrete_coons_patch`` documents. A side that is
        ``None`` is a collapsed one -- the patch has a pole at that corner.

    Returns
    -------
    list[[float, float, float]]
        The morphed vertices.
    list[list[int]]
        The faces, with vertices repeated inside a face dropped -- that is what
        turns the fan's innermost ring of quads into triangles.
    """
    ab, bc, dc, ad = sides

    # a missing side is the corner it collapsed to, repeated (as in compas)
    if not ab:
        ab = [ad[0]] * len(dc)
    if not bc:
        bc = [ab[-1]] * len(ad)
    if not dc:
        dc = [bc[-1]] * len(ab)
    if not ad:
        ad = [dc[0]] * len(bc)

    a, b, c, d = ab[0], bc[0], dc[-1], ad[-1]

    vertices = []
    for u, w in uw:
        p_ab = _polyline_point_at(ab, u)
        p_dc = _polyline_point_at(dc, u)
        p_ad = _polyline_point_at(ad, w)
        p_bc = _polyline_point_at(bc, w)

        point = []
        for k in range(3):
            # blend the two curves running along u, then the two along w ...
            along_u = p_ab[k] * (1 - w) + p_dc[k] * w
            along_w = p_ad[k] * (1 - u) + p_bc[k] * u
            # ... and take out the bilinear corner patch, which both blends
            # already contain.
            #
            # Summed with the builtin, not with a running total: compas adds
            # the four corner terms through ``sum_vectors``, and since 3.12
            # CPython's ``sum`` compensates the rounding of a float series. A
            # hand-rolled accumulation lands a unit in the last place away
            # from it, which is enough to shift a welded vertex.
            corner = sum((
                a[k] * ((1 - u) * (1 - w)),
                b[k] * (u * (1 - w)),
                c[k] * (u * w),
                d[k] * ((1 - u) * w),
            ))
            point.append(along_u + along_w - corner)
        vertices.append(point)

    faces = [[p for p, q in pairwise(face + face[:1]) if p != q] for face in faces]
    return vertices, faces


def pattern_morph_triangle(lam: list[tuple[float, float, float]], faces: list[list[int]], sides: list[list[list[float]] | None]) -> tuple[list[list[float]], list[list[int]]]:
    """Map a barycentric template onto a triangular patch through its three sides.

    The side-vertex interpolant (Nielson, 1979). From each corner ``Vi`` a
    straight spoke runs through the point to the side opposite ``Vi``::

        S_i = l_i Vi + (1 - l_i) E_i(s_i),   s_i = l_(i+2) / (l_(i+1) + l_(i+2))

    and the three spokes are blended with ``w_i = l_(i+1) l_(i+2)``, normalised.
    On side ``i`` only ``w_i`` survives, so every side polyline is reproduced;
    for straight sides the whole map is affine, so a straight line in the
    template is a straight line in the patch.

    Why not :func:`pattern_morph`: a pseudo-quad's Coons map singles out the
    collapsed side. With curved sides it put the same template 0.07 to 0.14 off
    its own 120 degree rotation, and moved it by as much when a different corner
    was the pole. This map is symmetric in the three corners by construction, so
    which corner is the pole changes nothing but the bookkeeping.

    Parameters
    ----------
    lam : list[(float, float, float)]
        Template vertices, as barycentric coordinates on (pole, V1, V2).
    faces : list[list[int]]
        Template faces.
    sides : list
        The four side polylines ``[ab, bc, dc, ad]`` of a pseudo-quad, as for
        :func:`pattern_morph`, exactly one of them ``None``.

    Returns
    -------
    list[[float, float, float]]
        The morphed vertices.
    list[list[int]]
        The faces, with vertices repeated inside a face dropped.
    """
    ab, bc, dc, ad = sides

    # back to one loop a -> b -> c -> d -> a; the collapsed side names the pole,
    # and the three sides after it run V0 -> V1 -> V2 -> V0
    loop = [ab, bc, dc[::-1] if dc else None, ad[::-1] if ad else None]
    k = [side is None for side in loop].index(True)
    t0, t1, t2 = loop[(k + 1) % 4], loop[(k + 2) % 4], loop[(k + 3) % 4]
    corners = [t0[0], t1[0], t2[0]]
    opposite = [t1, t2, t0]   # E_i, running V(i+1) -> V(i+2)

    vertices = []
    for l in lam:
        # on a corner or a side, return the input point untouched rather than
        # the blend, which would divide by zero on a corner and could come back
        # a unit in the last place away from the side polyline on a side --
        # and side points are what the neighbouring patch welds to
        zero = [i for i in range(3) if l[i] == 0.0]
        if len(zero) == 2:
            vertices.append(list(corners[3 - sum(zero)]))
            continue
        if zero:
            i = zero[0]
            j, h = (i + 1) % 3, (i + 2) % 3
            vertices.append(_polyline_point_at(opposite[i], l[h] / (l[j] + l[h])))
            continue

        spokes, weights = [], []
        for i in range(3):
            j, h = (i + 1) % 3, (i + 2) % 3
            e = _polyline_point_at(opposite[i], l[h] / (l[j] + l[h]))
            spokes.append([l[i] * v + (1.0 - l[i]) * x for v, x in zip(corners[i], e)])
            weights.append(l[j] * l[h])
        total = sum(weights)
        vertices.append([sum(weights[i] * spokes[i][c] for i in range(3)) / total for c in range(3)])

    faces = [[p for p, q in pairwise(face + face[:1]) if p != q] for face in faces]
    return vertices, faces


def _polyline_point_at(points: list[list[float]], t: float) -> list[float]:
    """Evaluate a polyline at ``t`` in [0, 1], parametrised by point *index*.

    Index parametrisation -- not arc length -- is what ``discrete_coons_patch``
    uses (``normalize_values(range(n))``), and what the side polylines were
    sampled with. Using anything else here would silently reparametrise every
    curved boundary.

    ``t`` landing on a node returns that node untouched, so a template whose
    parameters are grid nodes reproduces the discrete Coons patch exactly.

    Parameters
    ----------
    points : list[[float, float, float]]
        The polyline points.
    t : float
        Parameter in [0, 1].

    Returns
    -------
    [float, float, float]
    """
    n = len(points) - 1
    if n <= 0:
        return list(points[0])

    s = t * n
    node = int(round(s))
    if abs(s - node) < 1e-9:
        return list(points[min(max(node, 0), n)])

    i = int(s)
    i = min(max(i, 0), n - 1)
    f = s - i
    p, q = points[i], points[i + 1]
    return [a + f * (b - a) for a, b in zip(p, q)]


# ==============================================================================
# divisions
# ==============================================================================

def patch_divisions(sides: list[list[list[float]] | None]) -> tuple[int, int]:
    """The number of divisions a set of side polylines asks for.

    Parameters
    ----------
    sides : list
        ``[ab, bc, dc, ad]``, any one of which may be ``None``.

    Returns
    -------
    (int, int)
        Divisions along u and along w.
    """
    ab, bc, dc, ad = sides
    return len(ab if ab else dc) - 1, len(bc if bc else ad) - 1


# the patterns whose template needs a patch to be d by d with d even -- see
# ``_pattern_diagonal`` and ``_pattern_fan`` for why
_SQUARE_EVEN = ('diagonal', 'fan')


def reconcile_strip_densities(coarse: "CoarseQuadMesh", patterns: "str | dict[int, str]") -> dict[int, tuple[int, int]]:
    """Force the strip densities of a layout to something its patterns can tile.

    Densities live on strips, not on faces, because two patches sharing a
    coarse edge have to put the same points on it or the dense mesh will not
    weld. A pattern that constrains its divisions therefore constrains the
    whole layout, and has to do so *before* any patch is built.

    ``ortho`` takes any densities. ``diagonal`` and ``fan`` need their patch
    to be ``d`` by ``d`` with ``d`` even. A patch's ``nu`` and ``nw`` are the
    densities of the two strips crossing it, so such a patch ties those two
    strips to one density -- and each of them runs on through other patches,
    which may tie it to a third. The strips fall into groups that have to
    share a density, and each group takes the largest one in it, rounded up
    to even: refining is recoverable, coarsening is not.

    An ``ortho`` patch ties nothing, so in a layout mixing patterns only the
    strips the constrained patches chain together are touched.

    Parameters
    ----------
    coarse : :class:`CoarseQuadMesh`
        The layout. Its strip densities are modified in place.
    patterns : str or dict
        A key of :attr:`PATTERNS` for every face, or one per face as
        ``{fkey: pattern}``.

    Returns
    -------
    dict
        ``{skey: (old, new)}`` for every strip whose density was changed.
    """
    if not isinstance(patterns, dict):
        patterns = {fkey: patterns for fkey in coarse.faces()}

    edge_strip = {}
    for skey, edges in coarse.strips(data=True):
        for u, v in edges:
            edge_strip[u, v] = skey
            edge_strip[v, u] = skey

    # union-find over the strips; a square patch merges the two crossing it
    parent = {skey: skey for skey in coarse.strips()}

    def find(skey: int) -> int:
        while parent[skey] != skey:
            parent[skey] = parent[parent[skey]]
            skey = parent[skey]
        return skey

    tied = []
    for fkey in coarse.faces():
        if patterns.get(fkey) not in _SQUARE_EVEN:
            continue
        # the strips ``_patch_sides`` reads nu and nw from: two, or one if a
        # strip crosses itself here. Read off the halfedges rather than
        # ``face_strips`` so a pseudo-quad needs no special case.
        roots = {find(edge_strip[u, v]) for u, v in coarse.face_halfedges(fkey)}
        root = roots.pop()
        for other in roots:
            parent[other] = root
        tied.append(root)

    groups = {}
    for skey in parent:
        groups.setdefault(find(skey), []).append(skey)

    changed = {}
    for root in {find(skey) for skey in tied}:
        d = max(coarse.get_strip_density(skey) for skey in groups[root])
        d = max(2, d + d % 2)
        for skey in groups[root]:
            old = coarse.get_strip_density(skey)
            if old != d:
                coarse.set_strip_density(skey, d)
                changed[skey] = (old, d)

    if changed:
        print('reconcile_strip_densities: {} strip(s) raised to fit the '
              'diagonal/fan pattern(s): {}'.format(
                  len(changed),
                  ', '.join('{} {}->{}'.format(skey, old, new)
                            for skey, (old, new) in sorted(changed.items()))))

    return changed
