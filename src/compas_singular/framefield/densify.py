"""Step 7 -- **field-aware densification**. The field decides patch INTERIORS.

The missing half of the objective. ``cable -> field`` has worked and been
measured since ``12_cables`` part 1: a hard-constrained cable puts the field
0.0 degrees off itself. ``field -> mesh`` never worked, because
``CoarsePseudoQuadMesh.densification`` fills every patch with
``discrete_coons_patch`` -- a bilinear blend of the four coarse edges that
never consults the field. A square with a diagonal cable came out as a plain
axis-aligned grid with the cable running through it, identically at every
background resolution and however hard the cable pulled.

WHAT IS AND IS NOT FREE TO MOVE
-------------------------------

The four boundary polylines of a patch are **fixed**, exactly as
``discrete_coons_patch`` receives them. Only the interior nodes move. That is
not a simplification, it is the load-bearing decision:

* the shared edge between two patches comes from one ``edges_to_curves``
  polyline sampled at one strip density, so leaving boundaries alone keeps
  the patches welding into a manifold mesh, and keeps opposite edges of a patch
  at matching densities -- ``collect_strips``, ``set_strips_density`` and
  ``add_strip`` all keep working on the result (``11_downstream``);
* the coarse edges ALREADY follow the field. They are separatrices, handed
  through ``edges_to_curves``. It is only the interiors that ignored it.

THE ENERGY
----------

Write the patch as a grid ``x[i][j]``, ``i`` along ``ab``, ``j`` along ``ad``,
the same indexing ``discrete_coons_patch`` uses. Let ``e1`` be the cross arm
the ``i`` family follows and ``e2 = perp(e1)`` the one the ``j`` family
follows. A perfectly field-aligned grid has every ``i`` step parallel to
``e1`` and every ``j`` step parallel to ``e2``, so minimise what is left over::

    E = sum over i-edges  ( n1 . (x[i+1][j] - x[i][j]) )^2
      + sum over j-edges  ( n2 . (x[i][j+1] - x[i][j]) )^2

with ``n1 = perp(e1)``, ``n2 = perp(e2) = e1``, both sampled at the edge
midpoint, each row scaled by ``1 / |edge|`` so the residual is an angle rather
than a distance and a long edge cannot outvote a short one.

This is linear in the free nodes, so each round is one sparse least-squares
solve. The only nonlinearity is which arm of the cross each family follows, and
that is re-decided from the current geometry between solves.

Alignment ALONE is not enough, and the failure is loud rather than subtle. A
square patch under a field constant at 36.87 degrees -- ``15_baseline``'s
``square+cable`` -- wants a rotated grid inside a boundary pinned to the walls,
which is impossible. The unregularised solve tears it apart: 25 of 100 quads
folded, minimum angle 2.68 degrees, and it never converges. So the energy
carries a second term, ``stiffness`` times the Laplacian of the CORRECTION
``x - x_coons``, and the stiffness is chosen per patch by
:func:`relax_patch` -- weakest that :func:`_accepts` will take. The conflict
between a fixed boundary and a field that disagrees with it is genuine, and the
answer is to price it rather than to pretend it is not there.

Note that the two families constrain *different* components: an ``i`` edge says
nothing about spacing ALONG ``e1``, only about drift across it. The ``j`` edges
supply exactly that missing component, and vice versa. Together they pin both
coordinates of every interior node; neither family alone would.

WHY A CONSTANT FIELD GIVES THE GRID BACK, EXACTLY
-------------------------------------------------

The rectilinear plates -- square, L, T, U, plus, comb, slotted rectangle,
square with a square hole -- have an exactly constant field (measured: maximum
angular deviation 0.000 degrees). On a rectangular patch with uniformly sampled
straight sides, the Coons grid has every ``i`` step parallel to ``e1`` and every
``j`` step parallel to ``e2``, so **every residual above is identically zero**.
A sum of squares cannot go below zero, so the Coons grid is already the global
minimum and the solver has nothing to do. The Tikhonov term is anchored to the
Coons grid rather than to the running iterate for the same reason: it is zero
there too. Those eight domains come out bit-for-bit unchanged, which
``--constant-field`` in ``16_field_densify.py`` checks directly rather than
inferring from the baseline.

POLES ARE NOT DONE HERE, DELIBERATELY
-------------------------------------

A pseudo-quad patch keeps its Coons interior. Three reasons, and they are
reasons rather than an omission:

* a pole is a collapsed side, so one whole boundary row of the grid is a single
  repeated point and the least-squares system loses rank along it;
* the field at a pole is at or beside a singularity, which is precisely where
  ``theta`` has no continuous branch -- the arm-matching that makes the ``i``
  family well defined across a patch relies on the patch interior being
  singularity-free, which is true of every ordinary patch by construction of
  the separatrix layout and false here;
* the pole faces are where the suite's remaining degeneracies already live
  (``12_cables``' ring cable), so a scheme that guessed at them would be
  measured against noise.

:func:`field_densification` reports how many patches it skipped for this
reason, so a caller is never left inferring it.

QUALITY IS SPENT ONLY WHERE A GUIDE ASKED FOR IT
------------------------------------------------

Every relaxed interior is guarded, per patch rather than per mesh, so one bad
patch costs one patch rather than the whole domain's alignment. But the guard
is not one rule: a patch in a domain carrying a GUIDE CURVE may be spent down
to ``PATCH_MIN_ANGLE`` / ``PATCH_MAX_ANGLE`` / ``PATCH_MAX_ASPECT``, and a
patch in a domain without one may not be made worse in any respect at all.

:func:`_accepts` sets out the reasoning and the measurements. The short version
is that a cable is an instruction and distortion is its price, whereas on an
unguided domain the leftover misalignment is discretisation noise and paying
twelve degrees of minimum angle for half a degree of it is a trade nobody
wanted.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from math import atan2
from math import cos
from math import degrees
from math import pi
from math import sin

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from compas.geometry import discrete_coons_patch
from compas.itertools import pairwise
from compas.tolerance import TOL
from compas_singular.datastructures import meshes_join_and_weld
from compas_singular.datastructures import PseudoQuadMesh


__all__ = ['field_densification', 'FieldSampler', 'relax_patch',
           'PATCH_MIN_ANGLE', 'PATCH_MAX_ANGLE', 'PATCH_MAX_ASPECT']


#: The cross's period. Four arms, 90 degrees apart.
PERIOD = pi / 2.0

#: How far a patch may be spent to satisfy a GUIDE curve -- and only a guide;
#: see :func:`_accepts`. Angles in degrees, aspect as longest edge over
#: shortest. Floors, not targets: a patch already worse than these keeps its
#: own values as the limit rather than being allowed down to them.
PATCH_MIN_ANGLE = 30.0
PATCH_MAX_ANGLE = 150.0
PATCH_MAX_ASPECT = 5.0

#: Solve / re-arm rounds within one stiffness. The element quality is final
#: long before the last digit of the position is: on ``square+cable`` at
#: stiffness 1.0 the minimum angle is 42.7 degrees after 8 rounds and 42.65
#: after 40, while the per-round movement is still 0.15 at round 8. Twelve is
#: where the geometry has stopped changing; chasing the fixed point costs
#: another 30 rounds per stiffness and buys two hundredths of a degree.
ITERATIONS = 12

#: Stiffnesses tried, weakest first: how hard the smoothing term resists the
#: alignment term. The first one that settles AND yields an acceptable patch
#: wins, so a patch buys the most alignment its own boundary can carry and no
#: more. Measured on ``square+cable`` at bg 0.50, which is the suite's hardest
#: case for this -- a field constant at 36.87 degrees inside a square boundary
#: pinned axis-aligned::
#:
#:      stiffness    min     max  folds  ARmax  cable  verdict
#:           0.25  21.99  171.83      5  52.29   21.2  rejected: 5 folds
#:           0.50  26.25  163.75      0  22.29   23.4  rejected: 163.75 > 150
#:           1.00  42.69  144.23      0   4.00   27.0  ACCEPTED
#:           2.00  57.52  126.38      0   1.88   30.2
#:           4.00  71.97  109.74      0   1.38   32.8
#:          16.00  88.05   92.06      0   1.03   36.5
#:          64.00  89.87   90.13      0   1.00   36.8  ~ the Coons grid
#:
#: ``16_field_densify.py --tradeoff`` regenerates that table. Note what the top
#: of it costs: 0.50 buys 3.6 degrees of alignment over 1.00 and pays with an
#: aspect ratio of 22.3 and a 163.75 degree corner. The schedule does not skip
#: it -- :func:`_accepts` rejects it, which is the rule doing its job rather
#: than a constant chosen to suit one domain.
STIFFNESS = (0.25, 0.5, 1.0, 2.0, 4.0, 16.0, 64.0)

#: Weight of the pull back toward the Coons position, relative to the mean
#: alignment row weight. Conditions the system and bounds how far a patch can
#: drift; it is exactly zero on a constant field, where Coons IS the minimiser.
TIKHONOV = 1e-3

#: Stop early once no node moved further than this fraction of the patch's
#: mean edge length.
TOLERANCE = 1e-4

#: How far the per-round movement must have contracted, as a fraction of the
#: first round's, before the result is treated as a solution rather than as
#: wherever an oscillation happened to stop.
#:
#: An absolute threshold is the wrong tool here and was tried first: at
#: stiffness 1.0 the movement is still 0.15 of an edge length after 8 rounds
#: while the angles have not changed in six, so any threshold tight enough to
#: exclude a genuine oscillation also excludes every stiffness worth having.
#: Contraction separates them cleanly -- measured on ``square+cable``, the
#: alignment-only solve oscillates at 2.2, 1.5, 1.5, 1.1, 1.96, 1.31 and never
#: contracts at all, while stiffness 0.5 falls 1.646 -> 0.231 over eight.
CONTRACTION = 0.25


# ------------------------------------------------------------------
# sampling the field at an arbitrary point
# ------------------------------------------------------------------

class FieldSampler(object):
    """``theta`` at any point of the domain, with the tracer's point location.

    Wraps :meth:`Tracer.locate` and :meth:`CrossField.angle_in_face` so a patch
    interior -- which has no vertices of its own on the background mesh -- can
    ask the field what it is doing there. ``hint`` is threaded through because
    grid neighbours land in the same or an adjacent triangle almost always.

    The angle comes back folded into the cross's period; choosing a branch is
    :func:`_unwrap`'s job, not this one's.
    """

    def __init__(self, field, tracer):
        self.field = field
        self.tracer = tracer
        self.mesh = field.background.mesh

    def theta(self, point, hint=None):
        """``(theta, fkey)``, or ``(None, None)`` if the point is unreachable."""
        found = self.tracer.locate(point, hint)
        if found is None:
            found = self._nearest(point)
        if found is None:
            return None, None
        fkey, bary = found
        return self.field.angle_in_face(fkey, bary), fkey

    def _nearest(self, point):
        """Closest background face, for a node that fell outside the domain.

        A Coons interior can bulge a hair past a concave boundary arc. Rather
        than drop the node's constraint -- which would leave a hole in the
        least-squares system -- use the nearest triangle's centre. The error is
        bounded by the background spacing and only ever applies to nodes that
        are already outside.
        """
        tracer = self.tracer
        if tracer._grid is None:
            tracer._build_grid()
        i = int((point[0] - tracer._origin[0]) // tracer._cell)
        j = int((point[1] - tracer._origin[1]) // tracer._cell)
        best, best_d = None, float('inf')
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                for fkey in tracer._grid.get((i + di, j + dj), ()):
                    c = self.mesh.face_centroid(fkey)
                    d = (c[0] - point[0]) ** 2 + (c[1] - point[1]) ** 2
                    if d < best_d:
                        best, best_d = fkey, d
        if best is None:
            return None
        return best, (1 / 3.0, 1 / 3.0, 1 / 3.0)


# ------------------------------------------------------------------
# the patch solver
# ------------------------------------------------------------------

def _unwrap(grid, sampler):
    """A CONTINUOUS ``i``-family angle over the patch, or ``None``.

    ``angle_in_face`` returns the cross angle folded into one period, which
    makes ``e1`` and ``e2`` interchangeable at every node independently. Left
    like that, neighbouring rows of one patch can pick perpendicular arms and
    the grid kinks. So pick a branch once, from the patch's own ``i`` direction
    at one corner, and propagate it: along the ``j = 0`` column first, then
    along each row from it. Every step chooses the branch within a quarter turn
    of its predecessor, which is the same arm-matching that keeps a traced
    separatrix going straight instead of turning 90 degrees at an arbitrary
    triangle.

    Path-independent because a patch interior contains no singularity -- the
    separatrix layout puts every singularity on a patch CORNER. Pseudo-quads
    are the exception and they never reach here.
    """
    n, m, _ = grid.shape
    raw = np.zeros((n, m))
    ok = np.zeros((n, m), dtype=bool)
    hint = None
    for i in range(n):
        for j in range(m):
            t, hint = sampler.theta(grid[i, j], hint)
            if t is None:
                hint = None
                continue
            raw[i, j] = t
            ok[i, j] = True

    if not ok.any():
        return None

    theta = np.zeros((n, m))

    def branch(value, reference):
        return value + PERIOD * round((reference - value) / PERIOD)

    seed = atan2(grid[1, 0][1] - grid[0, 0][1], grid[1, 0][0] - grid[0, 0][0])
    theta[0, 0] = branch(raw[0, 0], seed) if ok[0, 0] else seed
    for i in range(1, n):
        theta[i, 0] = branch(raw[i, 0], theta[i - 1, 0]) if ok[i, 0] else theta[i - 1, 0]
    for i in range(n):
        for j in range(1, m):
            theta[i, j] = (branch(raw[i, j], theta[i, j - 1]) if ok[i, j]
                           else theta[i, j - 1])
    return theta


def _solve(grid, coons, theta, stiffness, tikhonov):
    """One least-squares solve for the interior, boundary held fixed.

    Three kinds of row:

    * one **alignment** residual per interior-touching grid edge, scaled by
      ``1 / |edge|``;
    * one **smoothing** residual per free node per axis, ``stiffness`` times the
      discrete Laplacian OF THE CORRECTION ``x - x_coons``, not of ``x``. That
      distinction is the whole reason a curved patch can be left alone: the
      Laplacian of a Coons grid over a curved patch is not zero, so smoothing
      ``x`` would pull every such patch about for reasons that have nothing to
      do with the field. Smoothing the correction is zero at ``x = x_coons`` on
      any patch whatever, so with no field influence nothing moves at all;
    * one **Tikhonov** row per free coordinate, likewise anchored to Coons.

    Alignment alone is not enough and the failure is not subtle. A square patch
    under a field constant at 36.87 degrees -- ``15_baseline``'s
    ``square+cable`` -- wants a rotated grid inside a boundary that is pinned to
    the walls, which is impossible, and the unregularised solve tears it: 25 of
    100 quads folded, minimum angle 2.68 degrees, and it never converges. The
    conflict is genuine, so the answer is to trade it off rather than to pretend
    it is not there.

    Returns a new grid, or ``None`` if the system is singular.
    """
    n, m, _ = grid.shape
    free = {}
    for i in range(1, n - 1):
        for j in range(1, m - 1):
            free[i, j] = len(free) * 2
    if not free:
        return grid

    rows, cols, vals, rhs = [], [], [], []
    weights = []

    def add(row, i, j, sign, nx, ny, w):
        """+/- (nx, ny) . x[i][j], to the free columns or to the right side."""
        if (i, j) in free:
            k = free[i, j]
            rows.extend([row, row])
            cols.extend([k, k + 1])
            vals.extend([sign * nx * w, sign * ny * w])
        else:
            rhs[row] -= sign * w * (nx * grid[i, j][0] + ny * grid[i, j][1])

    def edge(ia, ja, ib, jb):
        """One alignment row for the step (ia, ja) -> (ib, jb)."""
        # The family this step belongs to decides which arm it follows: an i
        # step follows e1 = theta, a j step follows e2 = theta + 90 degrees.
        # Only the perpendicular of that arm enters the residual, and perp is
        # insensitive to sign, so no orientation bookkeeping is needed.
        ang = 0.5 * (theta[ia, ja] + theta[ib, jb])
        if ia != ib:                      # i step -> perp(e1)
            nx, ny = -sin(ang), cos(ang)
        else:                             # j step -> perp(e2) = e1
            nx, ny = cos(ang), sin(ang)
        d = grid[ib, jb] - grid[ia, ja]
        length = (d[0] * d[0] + d[1] * d[1]) ** 0.5
        w = 1.0 / max(length, 1e-9)
        if (ia, ja) not in free and (ib, jb) not in free:
            return
        row = len(rhs)
        rhs.append(0.0)
        weights.append(w)
        add(row, ib, jb, +1.0, nx, ny, w)
        add(row, ia, ja, -1.0, nx, ny, w)

    for i in range(n - 1):
        for j in range(m):
            edge(i, j, i + 1, j)
    for i in range(n):
        for j in range(m - 1):
            edge(i, j, i, j + 1)

    if not rhs:
        return grid

    mean_w = sum(weights) / len(weights)

    # smoothing: stiffness * Laplacian(x - x_coons) = 0, per free node per axis.
    # A FIXED neighbour is at its Coons position by construction, so its
    # correction is exactly zero and it drops out of the row entirely -- neither
    # a column nor a right-hand-side term. That is what makes the correction
    # decay smoothly to nothing at the patch boundary.
    mu = stiffness * mean_w
    for (i, j), k in free.items():
        movable = [nb for nb in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1))
                   if nb in free]
        share = mu / 4.0
        for axis in (0, 1):
            row = len(rhs)
            target = mu * coons[i, j][axis]
            rows.append(row)
            cols.append(k + axis)
            vals.append(mu)
            for ni, nj in movable:
                rows.append(row)
                cols.append(free[ni, nj] + axis)
                vals.append(-share)
                target -= share * coons[ni, nj][axis]
            rhs.append(target)

    lam = tikhonov * mean_w
    for (i, j), k in free.items():
        for axis in (0, 1):
            row = len(rhs)
            rhs.append(lam * coons[i, j][axis])
            rows.append(row)
            cols.append(k + axis)
            vals.append(lam)

    A = coo_matrix((vals, (rows, cols)), shape=(len(rhs), 2 * len(free))).tocsr()
    b = np.array(rhs)
    normal = (A.T @ A).tocsc()
    try:
        solution = spsolve(normal, A.T @ b)
    except Exception:
        return None
    if solution is None or not np.all(np.isfinite(solution)):
        return None

    out = grid.copy()
    for (i, j), k in free.items():
        out[i, j][0] = solution[k]
        out[i, j][1] = solution[k + 1]
    return out


def _measure(grid, coons):
    """``(min_angle, max_angle, aspect_max, folds)`` of a relaxed patch.

    ``folds`` counts quads whose orientation flipped relative to their Coons
    counterpart -- compared against Coons rather than against an absolute
    winding, so a patch whose Coons parameterisation is itself reversed is
    judged on whether the relaxation CHANGED it.
    """
    n, m, _ = grid.shape
    lo, hi, aspect, folds = 180.0, 0.0, 0.0, 0
    for i in range(n - 1):
        for j in range(m - 1):
            quad = [grid[i, j], grid[i, j + 1], grid[i + 1, j + 1], grid[i + 1, j]]
            ref = [coons[i, j], coons[i, j + 1], coons[i + 1, j + 1], coons[i + 1, j]]
            if _signed_area(quad) * _signed_area(ref) <= 0.0:
                folds += 1
            angles = _angles(quad)
            lo = min(lo, min(angles))
            hi = max(hi, max(angles))
            lengths = [((quad[k][0] - quad[(k + 1) % 4][0]) ** 2
                        + (quad[k][1] - quad[(k + 1) % 4][1]) ** 2) ** 0.5
                       for k in range(4)]
            shortest = min(lengths)
            aspect = (float('inf') if shortest <= 0.0
                      else max(aspect, max(lengths) / shortest))
    return lo, hi, aspect, folds


def _accepts(grid, coons, spend):
    """**The rule that lets tier 2 and tier 3 coexist.**

    Nothing may fold, ever. Past that there are two regimes, and which one a
    patch is in is decided by whether the domain carries a GUIDE CURVE.

    * ``spend=False`` -- no guide. **Improvement only.** The relaxation may
      make the patch better aligned but may not make a single element worse:
      minimum angle may not fall, maximum angle may not rise, aspect ratio may
      not grow. If no stiffness manages that, the patch keeps its Coons
      interior and the domain does not move at all.
    * ``spend=True`` -- a guide curve. The patch may be spent down to
      ``PATCH_MIN_ANGLE`` / ``PATCH_MAX_ANGLE`` / ``PATCH_MAX_ASPECT``, or to
      the Coons patch's own values where those are already worse.

    That distinction is the design decision in this module, so it is worth
    being explicit about why it is not merely a way of passing two sets of
    tests. **Element quality is spent only to satisfy something the user
    actually asked for.** A cable is an explicit instruction that a curve
    matters and that the mesh should follow it; distortion is the price and the
    user has asked to pay it. A domain with no guide has issued no such
    instruction -- its field is just the smoothest cross field the boundary
    admits, and there is nothing to buy.

    That "nothing to buy" is measured, not assumed. Per-patch mean
    misalignment of the Coons interior, over the suite::

        disc         1.3 - 11.0 deg      hexagon      3.1 - 16.3 deg
        ellipse      2.3 - 15.0 deg      stadium      1.8 -  7.3 deg
        square+cable       23.3 deg      (36.9 deg along the cable itself)

    and the alignment the weakest stiffness buys on the unguided ones is 0.1 to
    5.5 degrees of that. Paying 12 degrees of minimum angle -- which is what
    the disc did at bg 0.50 before this rule existed -- to recover half a
    degree of alignment that was discretisation noise to begin with is not a
    trade anybody wanted.
    """
    lo, hi, aspect, folds = _measure(grid, coons)
    if folds:
        return False
    ref_lo, ref_hi, ref_aspect, _ = _measure(coons, coons)
    if not spend:
        return lo >= ref_lo and hi <= ref_hi and aspect <= ref_aspect
    return (lo >= min(PATCH_MIN_ANGLE, ref_lo)
            and hi <= max(PATCH_MAX_ANGLE, ref_hi)
            and aspect <= max(PATCH_MAX_ASPECT, ref_aspect))


def _signed_area(points):
    total = 0.0
    for a, b in pairwise(list(points) + list(points[:1])):
        total += a[0] * b[1] - b[0] * a[1]
    return 0.5 * total


def _angles(points):
    out = []
    k = len(points)
    for i in range(k):
        a, b, c = points[i - 1], points[i], points[(i + 1) % k]
        ux, uy = a[0] - b[0], a[1] - b[1]
        vx, vy = c[0] - b[0], c[1] - b[1]
        lu = (ux * ux + uy * uy) ** 0.5
        lv = (vx * vx + vy * vy) ** 0.5
        if lu < 1e-12 or lv < 1e-12:
            out.append(0.0)
            continue
        dot = max(-1.0, min(1.0, (ux * vx + uy * vy) / (lu * lv)))
        out.append(degrees(np.arccos(dot)))
    return out


def relax_patch(ab, bc, dc, ad, sampler, iterations=ITERATIONS,
                tikhonov=TIKHONOV, stiffness=STIFFNESS, guard=True,
                spend=False):
    """**The field-aware replacement for** ``discrete_coons_patch``.

    Same signature, same vertex order -- ``i * m + j`` -- so it drops into the
    densification loop unchanged and a patch that does not move comes out
    bit-identical.

    Works down the ``stiffness`` schedule, weakest first, and takes the first
    result :func:`_accepts` passes. Weakest-first is the point: a patch buys as
    much alignment as its own fixed boundary can carry, and a patch that can
    carry none is left exactly as Coons made it. There is no global tuning
    constant deciding that -- each patch answers for itself.

    Parameters
    ----------
    ab, bc, dc, ad : list[[x, y, z]]
        The four boundary polylines, exactly as ``discrete_coons_patch`` wants
        them. Held FIXED; only the interior moves.
    sampler : FieldSampler
    iterations, tikhonov, stiffness : optional
    guard : bool, optional
        ``False`` takes the first stiffness that SETTLES without asking whether
        the patch it produced is usable. Only for measuring the trade-off --
        ``16_field_densify.py --tradeoff`` -- never for producing a mesh.
    spend : bool, optional
        Whether element quality may be spent to buy alignment. True only where
        a guide curve asked for it. See :func:`_accepts`.

    Returns
    -------
    (list, list, str, dict)
        Vertices, faces, a one-word verdict -- ``'relaxed'``, ``'flat'`` (no
        interior node, or none that moved) or ``'guarded'`` (nothing on the
        schedule was acceptable, so the Coons interior was kept) -- and what
        was measured: ``stiffness``, ``min``, ``max``, ``folds``.
    """
    vertices, faces = discrete_coons_patch(ab, bc, dc, ad)
    n, m = len(ab), len(bc)
    coons = np.array([[v[0], v[1], v[2]] for v in vertices], dtype=float)
    coons = coons.reshape(n, m, 3)
    info = {'stiffness': None, 'min': None, 'max': None, 'folds': None,
            'aspect': None}

    if n < 3 or m < 3:
        return vertices, faces, 'flat', info

    scale = float(np.mean(np.abs(np.diff(coons[:, :, :2], axis=0)))) or 1.0
    still = TOLERANCE * scale

    for mu in stiffness:
        grid = coons.copy()
        first = shift = float('inf')
        for _ in range(iterations):
            theta = _unwrap(grid, sampler)
            if theta is None:
                return vertices, faces, 'guarded', info
            nxt = _solve(grid, coons, theta, mu, tikhonov)
            if nxt is None:
                return vertices, faces, 'guarded', info
            shift = float(np.max(np.abs(nxt - grid)))
            first = shift if first == float('inf') else first
            grid = nxt
            if shift < still:
                break

        if shift > CONTRACTION * first and shift > still:
            continue                      # never contracted -- not a solution

        if float(np.max(np.abs(grid - coons))) < still:
            # The Coons grid was ALREADY the minimiser and the solve only added
            # its own round-off -- a constant field over a rectangular patch,
            # which is every rectilinear plate in the suite. Hand back the input
            # untouched so those domains come out bit-identical rather than
            # identical to fourteen decimal places. TRAP 3, module docstring.
            return vertices, faces, 'flat', info

        lo, hi, aspect, folds = _measure(grid, coons)
        info = {'stiffness': mu, 'min': lo, 'max': hi, 'folds': folds,
                'aspect': aspect}
        if not guard or _accepts(grid, coons, spend):
            out = [[float(p[0]), float(p[1]), float(p[2])]
                   for p in grid.reshape(n * m, 3)]
            return out, faces, 'relaxed', info

    return vertices, faces, 'guarded', info


# ------------------------------------------------------------------
# the densification itself
# ------------------------------------------------------------------

def field_densification(coarse, field, tracer, edges_to_curves=None,
                        iterations=ITERATIONS, tikhonov=TIKHONOV,
                        stiffness=STIFFNESS, guard=True, spend=False,
                        field_aware=True):
    """**Densify a coarse layout with the field steering patch interiors.**

    A drop-in for ``CoarsePseudoQuadMesh.densification`` -- same strip
    densities, same ``edges_to_curves``, same welding, same ``face_pole``
    bookkeeping, and it sets the result on ``coarse`` so ``get_quad_mesh``
    keeps working. The single difference is that an ordinary patch's interior
    comes from :func:`relax_patch` instead of ``discrete_coons_patch``.

    ``compas_singular`` itself is not touched; this is a parallel
    implementation of the same loop.

    Parameters
    ----------
    coarse : CoarsePseudoQuadMesh
        Strips already collected and densities already set.
    field : CrossField
    tracer : Tracer
    edges_to_curves : dict, optional
    iterations, tikhonov, stiffness, guard, spend : optional
        Passed to :func:`relax_patch`.
    field_aware : bool, optional
        ``False`` reproduces ``discrete_coons_patch`` exactly, for A/B tests.

    Returns
    -------
    (mesh, dict)
        The dense mesh, and counts: ``patches``, ``relaxed``, ``poles``
        (skipped because they are pseudo-quads), ``guarded`` (relaxed and
        rejected) and ``flat`` (no interior node to move), plus ``stiffness``
        -- the stiffnesses actually used, worst first -- and ``folds``.
    """
    edge_strip = {}
    for strip, edges in coarse.strips(data=True):
        for u, v in edges:
            edge_strip[u, v] = strip
            edge_strip[v, u] = strip

    pole_map = [TOL.geometric_key(coarse.vertex_coordinates(pole))
                for pole in coarse.poles()]

    sampler = FieldSampler(field, tracer) if field_aware else None
    stats = {'patches': 0, 'relaxed': 0, 'poles': 0, 'guarded': 0, 'flat': 0,
             'stiffness': [], 'folds': 0}
    # ``geometric_key -> (u, v, index)``: which coarse edge each patch-boundary node
    # came from, and its position ALONG that edge. Geometric key is the one identity
    # that survives ``meshes_join_and_weld``, which is keyed the same way.
    #
    # A post-pass that lets a seam node SLIDE along its own separatrix needs both
    # halves. The curve, because recovering it afterwards by proximity picks the
    # wrong polyline wherever two seams meet. The index, because a node free to slide
    # anywhere along a curve can slide PAST its neighbour and collapse the edge
    # between them -- measured, that turns the ring cable into a 0.00/180.00 degree
    # mesh. Order along the curve is the guard, and here is the only place it is
    # still known exactly.
    stats['seam_edge'] = {}

    meshes = []
    for fkey in coarse.faces():
        stats['patches'] += 1
        polylines = []
        for u, v in coarse.face_halfedges(fkey):
            d = coarse.get_strip_density(edge_strip[u, v])
            # Falls back to a straight chord for any edge missing from
            # ``edges_to_curves`` instead of raising -- the mapping no longer
            # has to be complete for the whole layout, only for the edges the
            # caller (``CoarseQuadMesh.densification``'s ``boundary_curvature``
            # / ``skeleton_curvature`` split) actually wants curved.
            polyline = coarse._create_patch_edge(u, v, d, edges_to_curves)
            for index, point in enumerate(polyline):
                stats['seam_edge'].setdefault(TOL.geometric_key(point), (u, v, index))
            polylines.append(polyline)

        pseudo = coarse.is_face_pseudo_quad(fkey)
        if pseudo:
            pole = coarse.attributes['face_pole'][fkey]
            idx = coarse.face_vertices(fkey).index(pole)
            polylines.insert(idx, None)

        ab, bc, cd, da = polylines
        dc = cd[::-1] if cd else None
        ad = da[::-1] if da else None

        if sampler is None or pseudo:
            # A pole keeps its Coons interior -- see the module docstring. The
            # collapsed side costs the system its rank, and the field has no
            # continuous branch at a singularity to integrate anyway.
            vertices, faces = discrete_coons_patch(ab, bc, dc, ad)
            stats['poles' if pseudo else 'coons'] = stats.get(
                'poles' if pseudo else 'coons', 0) + 1
        else:
            vertices, faces, verdict, info = relax_patch(
                ab, bc, dc, ad, sampler, iterations=iterations,
                tikhonov=tikhonov, stiffness=stiffness, guard=guard,
                spend=spend)
            stats[verdict] += 1
            if verdict == 'relaxed':
                stats['stiffness'].append(info['stiffness'])
                stats['folds'] += info['folds'] or 0

        faces = [[u for u, v in pairwise(face + face[:1]) if u != v]
                 for face in faces]
        meshes.append(PseudoQuadMesh.from_vertices_and_faces_with_face_poles(
            vertices, faces))

    face_pole_map = {}
    for mesh in meshes:
        for fkey in mesh.faces():
            fv = mesh.face_vertices(fkey)
            for u, v in pairwise(fv + fv[:1]):
                gk = TOL.geometric_key(mesh.vertex_coordinates(u))
                if gk in pole_map and gk == TOL.geometric_key(mesh.vertex_coordinates(v)):
                    face_pole_map[TOL.geometric_key(mesh.face_center(fkey))] = gk
                    break

    coarse.set_quad_mesh(meshes_join_and_weld(meshes))
    dense = coarse.get_quad_mesh()

    face_pole = {}
    for fkey in dense.faces():
        gk = TOL.geometric_key(dense.face_center(fkey))
        if gk in face_pole_map:
            for vkey in dense.face_vertices(fkey):
                if TOL.geometric_key(dense.vertex_coordinates(vkey)) == face_pole_map[gk]:
                    face_pole[fkey] = vkey
                    break
    dense.attributes['face_pole'] = face_pole

    return dense, stats
