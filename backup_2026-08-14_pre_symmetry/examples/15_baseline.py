"""EXAMPLE 6 -- the LOCKED QUALITY BASELINE, and the diff that makes a change reviewable.

Why this exists
---------------

There was no automated quality gate. ``d.route()`` was the de-facto one, and it
stopped discriminating the moment ``arrangement.py`` put all 16 domains on route
``'field'``. ``_acceptable`` checked all-quad, manifold and 97-103% coverage --
the same three properties that certified a mesh carrying a literal 180 degree
angle. So quality was being checked BY HAND, and hand-checking failed in the way
hand-checking always fails: a report claimed the ellipse improved to 65.3 degrees
minimum angle when the number belonged to a different row of the same table.

This script is the replacement. It measures, it commits the numbers, and it
diffs. **The diff is the deliverable** -- an absolute threshold cannot catch
"this change made the ellipse worse", because minimum angle across the field
route legitimately spans 24 to 90 degrees over this suite. Only a committed
per-row baseline can.

Two tiers, and both matter
--------------------------

* the **hard floor** in ``framefield/quality.py`` -- absolute and permanent,
  degeneracies only, checked on every run regardless of the baseline;
* the **regression check** here -- every row, every metric, against
  ``baseline.json``.

Run:
    python 15_baseline.py                 # compare against the committed baseline
    python 15_baseline.py --regenerate    # rewrite baseline.json
    python 15_baseline.py --domain disc   # one domain, all resolutions
    python 15_baseline.py --poles         # pseudo-quad metrics self-test
    python 15_baseline.py --corners       # repair.SHARP_TURN margin self-test
    python 15_baseline.py --relax         # same 64 rows under CrossField(relax=True)

``--relax`` diffs the OPT-IN field solver against the committed baseline, which
was measured with the default one -- so every row it moves is a difference
between the two solvers, not a regression. It refuses ``--regenerate``: the
baseline records the default solver, and a baseline written from the other one
would silently redefine what every future diff means. See
``19_field_accuracy.py`` for why the opt-in solver exists.

Exit code 0 means every row matches the committed baseline and no NEW hard-floor
breach appeared. Non-zero means something moved, and the report names the row,
the metric, the old value and the new one.

Standing defects recorded in the baseline
-----------------------------------------

One of the sixteen domains is still broken, and the baseline records it as
broken rather than hiding it:

* **square+ring cable** -- a 180.00 degree face at backgrounds 0.5, 0.4 and 0.3,
  and a route that degrades from ``'field'`` to ``'polygon'`` at 0.4 and below.
  Also the domain carrying the suite's only POLES (14 in the dense mesh at bg
  0.5), which is why it is in the matrix at all: a pole is what a naive metric
  divides by zero on. ``--poles`` checks that.

Two more were caught by this harness's FIRST run and have since been fixed. They
stay written down because the SHAPE of the bug matters more than the bug:

* **comb-plate** used to raise at every resolution, and **rect with slot** used
  to give 2 coarse patches for an 8-corner plate, 100.5% coverage and a
  0.00/180.00 degree face while reporting route ``'field'``, all-quad and
  manifold. Both are now 90.0/90.0/AR 1.00 at 100% coverage, 11 and 5 patches.

  One cause. ``boundary_corners``' ``spacing`` default of perimeter/25 = 2.96
  swallowed the comb's 2-unit teeth, so its right angles grouped into a "run",
  were read as a sampled curve, and -- every turn in that run being exactly 90
  degrees, so none exceeding ``typical + limit`` -- yielded no corner at all.
  The rect-with-slot lost the two reflex corners at its slot bottom the same
  way: its perimeter is exactly 50, so ``spacing`` is exactly 2.0, and the
  slot's bottom segment is exactly 2.0.

  That was the **fourth instance of one species**: a tolerance scaled to a
  GLOBAL quantity, applied to a domain whose LOCAL feature is smaller than that
  scale. The fix is ``repair.SHARP_TURN`` -- an ANGLE rather than a length,
  because an angle has no length scale and so no small feature can be smaller
  than it. See ``boundary_corners``' docstring for what breaks next, and
  ``--corners``, which fails if the threshold ever loses its margin here.

What the resolution sweep says
------------------------------

Coarse patch count is resolution-independent for almost every domain -- ellipse
10/10/10/10, stadium 4/4/4/4, hexagon 11/11/11/11, pentagon 5/5/5/5 -- and the
question was whether quality follows it. It does, with exactly one exception,
and the exception has a mechanism.

Minimum angle across backgrounds 0.6/0.5/0.4/0.3, worst spread first::

    ellipse           24/65/61/66   spread 41.5   poles 4/0/0/0
    disc+round hole   18/21/27/22   spread  9.0   poles 0/0/0/0
    disc              59/57/56/52   spread  6.8   poles 0/0/0/0
    stadium           58/60/55/56   spread  5.8   poles 0/0/0/0
    pentagon          66/68/65/71   spread  5.7   poles 0/0/0/0
    hexagon           68/65/62/66   spread  5.6   poles 0/0/0/0
    the 7 rectilinear plates        spread  0.0   poles 0/0/0/0

The ellipse's 41.5 degree spread is 4.5x the next largest, and it is **one row**:
background 0.6, the only ellipse resolution whose layout carries POLES. Drop that
row and its spread is 4.5 degrees -- smaller than the disc, stadium, pentagon or
hexagon. So this is not gradual resolution sensitivity at all. It is a single
discrete state change, and the state that changes is pole count, not patch count:
the ellipse holds 10 patches at every resolution while its poles go 4/0/0/0.

**No other domain does this.** The only other domain whose pole count varies is
``square+ring cable`` (5/14/0/0), and its quality is poor at every resolution, so
there is nothing for the pole to make discontinuous.

That is as far as this goes deliberately -- it may simply be what an ellipse
does. It is written down because field-aware densification is next, it will move
every one of these numbers, and "the ellipse is jumpy" is a much worse thing to
discover then than "the ellipse has a pole at 0.6 and does not below it".
"""
import json
import os
import sys
from math import cos, pi, sin

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from framefield.decomposition import FieldDecomposition             # noqa: E402
from framefield.quality import hard_floor                           # noqa: E402


BASELINE = os.path.join(HERE, 'baseline.json')
DATA_DIR = os.path.join(HERE, '..', 'data')


# --- the domains -------------------------------------------------------------
# Every domain the suite exercises, including the two that only ever appeared in
# prose (comb, rectangle with a slot) and the two cable cases from 12_cables.
# The ring cable earns its place: it is the only domain in the suite that
# produces POLES, and a pole is precisely what a naive quality metric divides by
# zero on.

SQUARE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]
L_SHAPE = [[0, 0, 0], [12, 0, 0], [12, 4, 0], [5, 4, 0], [5, 10, 0], [0, 10, 0]]
PENTAGON = [[0, 0, 0], [10, 0, 0], [13, 7, 0], [5, 12, 0], [-2, 7, 0]]
SQUARE_HOLE = [[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]]
T_SHAPE = [[0, 0, 0], [12, 0, 0], [12, 4, 0], [8, 4, 0],
           [8, 10, 0], [4, 10, 0], [4, 4, 0], [0, 4, 0]]
U_SHAPE = [[0, 0, 0], [12, 0, 0], [12, 10, 0], [9, 10, 0],
           [9, 4, 0], [3, 4, 0], [3, 10, 0], [0, 10, 0]]
PLUS = [[4, 0, 0], [8, 0, 0], [8, 4, 0], [12, 4, 0], [12, 8, 0], [8, 8, 0],
        [8, 12, 0], [4, 12, 0], [4, 8, 0], [0, 8, 0], [0, 4, 0], [4, 4, 0]]

#: three 2-unit slots cut into a plate -- the shape ``decomposition.py`` cites
#: for why the backstop may never triangulate coarser than the background.
#: Slots AND teeth are both 2 units: at 1 unit the teeth are thinner than two
#: background elements at tl=0.6, no separatrix survives in them, and the domain
#: raises rather than meshing. That is a resolution limit of the tooth, not a
#: property of combs, and it is not what this row is here to measure.
COMB = [[0, 0, 0], [14, 0, 0], [14, 8, 0],
        [12, 8, 0], [12, 3, 0], [10, 3, 0], [10, 8, 0],
        [8, 8, 0], [8, 3, 0], [6, 3, 0], [6, 8, 0],
        [4, 8, 0], [4, 3, 0], [2, 3, 0], [2, 8, 0],
        [0, 8, 0]]
RECT_SLOT = [[0, 0, 0], [12, 0, 0], [12, 8, 0], [7, 8, 0],
             [7, 3, 0], [5, 3, 0], [5, 8, 0], [0, 8, 0]]

DISC = [[5 + 5 * cos(2 * pi * i / 48), 5 + 5 * sin(2 * pi * i / 48), 0.0]
        for i in range(48)]
ROUND_HOLE = [[5 + 2 * cos(-2 * pi * i / 48), 5 + 2 * sin(-2 * pi * i / 48), 0.0]
              for i in range(48)]
HEXAGON = [[5 + 5 * cos(2 * pi * i / 6), 5 + 5 * sin(2 * pi * i / 6), 0.0]
           for i in range(6)]
ELLIPSE = [[5 + 7 * cos(2 * pi * i / 60), 5 + 3 * sin(2 * pi * i / 60), 0.0]
           for i in range(60)]
STADIUM = ([[0, 0, 0], [12, 0, 0]]
           + [[12 + 3 * cos(a), 3 + 3 * sin(a), 0.0]
              for a in [-pi / 2 + pi * i / 10 for i in range(1, 10)]]
           + [[12, 6, 0], [0, 6, 0]])

CABLE = [[1.0, 2.0, 0.0], [5.0, 5.0, 0.0], [9.0, 8.0, 0.0]]
RING = [[5 + 3 * cos(2 * pi * i / 24), 5 + 3 * sin(2 * pi * i / 24), 0.0]
        for i in range(25)]


def _decomposition_boundary():
    """Outer + hole loops of ``examples/data/01_decomposition.json``.

    That file is ``[outer_boundary, inner_boundaries, polyline_features,
    point_features]`` -- the input to the OLD ``SkeletonDecomposition`` route in
    ``02_decomposition_discrete_planar.py``. ``polyline_features`` is empty here
    and ``point_features`` (5 points) has no equivalent in ``FieldDecomposition``
    -- ``from_boundary`` takes curve ``guides``, not point constraints -- so only
    the two boundary loops carry over. This is a real freeform outline (111
    points) with two freeform holes, not a parametric shape, which is the point
    of adding it: every other curved domain in this suite is a circle, ellipse
    or stadium.
    """
    path = os.path.join(DATA_DIR, '01_decomposition.json')
    with open(path) as f:
        outer, holes, _polyline_features, _point_features = json.load(f)
    return outer, holes


DECOMPOSITION_OUTER, DECOMPOSITION_HOLES = _decomposition_boundary()


# --- Rhino curve-OBJ boundaries (RV, RV+hole) ---------------------------------
# Rhino exports a closed planar curve as a ``cstype rat bspline`` entity --
# degree-2 rational B-spline control points as ``v x y z w`` (weight in the 4th
# column), knot vector as ``parm u`` -- NOT as a mesh. There is no compas
# backend that evaluates this here (no ``compas_occ`` in this env; ``NurbsCurve``
# without a plugin raises ``PluginNotInstalledError``), so it is parsed and
# evaluated by hand with De Boor's algorithm (Piegl & Tiller, "The NURBS Book").

def _find_span(n, p, u, knots):
    """Knot span index containing ``u``, for a clamped B-spline of degree ``p``
    with ``n + 1`` control points."""
    if u >= knots[n + 1]:
        return n
    if u <= knots[p]:
        return p
    lo, hi = p, n + 1
    mid = (lo + hi) // 2
    while u < knots[mid] or u >= knots[mid + 1]:
        if u < knots[mid]:
            hi = mid
        else:
            lo = mid
        mid = (lo + hi) // 2
    return mid


def _curve_point(curve, u):
    """Point on a rational B-spline at parameter ``u``, De Boor's algorithm."""
    degree, knots, points = curve['degree'], curve['knots'], curve['points']
    span = _find_span(len(points) - 1, degree, u, knots)
    d = [[x * w, y * w, z * w, w] for x, y, z, w in points[span - degree: span + 1]]
    for r in range(1, degree + 1):
        for j in range(degree, r - 1, -1):
            i = span - degree + j
            denom = knots[i + degree - r + 1] - knots[i]
            alpha = 0.0 if denom == 0.0 else (u - knots[i]) / denom
            d[j] = [(1.0 - alpha) * d[j - 1][c] + alpha * d[j][c] for c in range(4)]
    x, y, z, w = d[degree]
    return x / w, y / w, z / w


def _read_curve_obj(path):
    """Every ``cstype rat bspline`` entity in a Rhino curve-OBJ, as
    ``{'degree', 'u0', 'u1', 'knots', 'points'}`` dicts. Control point indices
    in ``curv`` are global across the file, 1-based, matching plain ``v`` order."""
    vertices, curves = [], []
    degree = indices = u0 = u1 = None
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == 'v':
                vertices.append(tuple(float(x) for x in parts[1:5]))
            elif parts[0] == 'deg':
                degree = int(parts[1])
            elif parts[0] == 'curv':
                u0, u1 = float(parts[1]), float(parts[2])
                indices = [int(i) for i in parts[3:]]
            elif parts[0] == 'parm' and parts[1] == 'u':
                curves.append({'degree': degree, 'u0': u0, 'u1': u1,
                               'knots': [float(x) for x in parts[2:]],
                               'points': [vertices[i - 1] for i in indices]})
    return curves


def _sample_curve(curve, spacing):
    """Closed loop of ``[x, y, 0]`` points, roughly ``spacing`` apart.

    Rhino's export parameterises each curve close to arc length -- these
    shapes came from a PolyCurve of lines and arcs converted to one NURBS, and
    the u-range (147.07 for RV's outer wall) is of the same order as its own
    perimeter -- so sampling uniformly in u is close enough to uniform arc
    length; measured segment lengths come out 0.32-0.52 at ``spacing=0.5``.
    The last sample (``u1``) is dropped: that is where the curve returns to its
    start control point (how the OBJ closes it), so keeping it would duplicate
    the first point. The curve is flat in Rhino's XZ plane here (every vertex's
    middle coordinate is 0), so ``(x, z)`` becomes this suite's 2D ``(x, y)``.
    """
    u0, u1 = curve['u0'], curve['u1']
    n = max(8, round((u1 - u0) / spacing))
    pts = [_curve_point(curve, u0 + (u1 - u0) * i / n) for i in range(n)]
    return [[x, z, 0.0] for x, _y, z in pts]


def _signed_area(loop):
    return 0.5 * sum(a[0] * b[1] - b[0] * a[1]
                     for a, b in zip(loop, loop[1:] + loop[:1]))


def _curve_obj_boundary(filename, spacing=0.5):
    """Outer + hole loops of a Rhino curve-OBJ, largest loop first.

    ``RV.obj`` has one curve (the outer wall); ``RV_with_hole.obj`` has two
    (the same outer wall plus a small closed hole) -- whichever loop has the
    larger bounding area is the outer boundary, the rest are holes. Winding
    does not need fixing here: ``FieldDecomposition`` derives outer-CCW /
    hole-CW itself from the background mesh, not from input order (see
    ``field.py``'s ``_oriented_boundary_loops``).
    """
    path = os.path.join(DATA_DIR, filename)
    loops = [_sample_curve(curve, spacing) for curve in _read_curve_obj(path)]
    loops.sort(key=lambda loop: abs(_signed_area(loop)), reverse=True)
    return loops[0], loops[1:]


RV_OUTER, RV_HOLES = _curve_obj_boundary('RV.obj')
RV_WITH_HOLE_OUTER, RV_WITH_HOLE_HOLES = _curve_obj_boundary('RV_with_hole.obj')


#: ``(label, outer, holes, guides, extra kwargs)``. 21 domains x 4 background
#: resolutions = 84 rows.
DOMAINS = [
    ('square',            SQUARE,  None,           None,     {}),
    ('L-shape',           L_SHAPE, None,           None,     {}),
    ('T-plate',           T_SHAPE, None,           None,     {}),
    ('U-plate',           U_SHAPE, None,           None,     {}),
    ('plus-plate',        PLUS,    None,           None,     {}),
    ('comb-plate',        COMB,    None,           None,     {}),
    ('rect with slot',    RECT_SLOT, None,         None,     {}),
    ('square+sq hole',    SQUARE,  [SQUARE_HOLE],  None,     {}),
    ('pentagon',          PENTAGON, None,          None,     {}),
    ('hexagon',           HEXAGON, None,           None,     {}),
    ('disc',              DISC,    None,           None,     {}),
    ('ellipse',           ELLIPSE, None,           None,     {}),
    ('stadium',           STADIUM, None,           None,     {}),
    ('disc+round hole',   DISC,    [ROUND_HOLE],   None,     {}),
    ('square+cable',      SQUARE,  None,           [CABLE],
     dict(guide_band=3.0, guide_weight=None)),
    ('square+ring cable', SQUARE,  None,           [RING],
     dict(guide_band=1.5, guide_weight=5.0)),
    ('square+circle',     SQUARE,  [ROUND_HOLE],    None,     {}),
    ('circle+square',     DISC,    [SQUARE_HOLE],   None,     {}),
    ('decomposition',     DECOMPOSITION_OUTER, DECOMPOSITION_HOLES, None, {}),
    ('RV',                RV_OUTER, RV_HOLES,  None,     {}),
    ('RV+hole',           RV_WITH_HOLE_OUTER, RV_WITH_HOLE_HOLES, None, {}),
]

#: background triangulation spacing -- NOT the quad size
RESOLUTIONS = (0.6, 0.5, 0.4, 0.3)

#: quad size, held fixed across the sweep so the only thing varying is the
#: background. 1.0 is what the nine-row reference table was measured at.
QUAD_TARGET = 1.0

#: the metrics compared, and how close counts as unchanged. Integers must match
#: exactly; the floats get a tolerance that is far tighter than any real
#: regression and far looser than last-bit drift.
FIELDS = [
    ('route', None),
    ('coarse_faces', 0),
    ('poles', 0),
    ('coverage', 5e-4),
    ('faces', 0),
    ('min_angle', 5e-3),
    ('max_angle', 5e-3),
    ('aspect_max', 5e-3),
    ('share_below_20', 5e-4),
    ('irregular_interior', 0),
    # derived from the metrics above, but compared anyway: a row crossing the
    # hard floor in either direction is the single most important thing a diff
    # can say, and it should never be something the reader has to infer.
    ('hard_floor', None),
]


def measure_one(label, outer, holes, guides, kwargs, target_length, relax=False):
    """One row. Never raises -- a domain that blows up is recorded as blown up.

    A harness that dies on the first bad domain reports nothing about the other
    63, which is the opposite of what a baseline is for.
    """
    row = {'domain': label, 'target_length': target_length}
    try:
        d = FieldDecomposition.from_boundary(
            outer, holes, guides=guides, target_length=target_length,
            relax=relax, **kwargs)
        dense = d.quad_mesh(target_length=QUAD_TARGET)
        q = d.quality(dense)
    except Exception as exc:
        row['route'] = 'EXCEPTION'
        row['error'] = '{}: {}'.format(type(exc).__name__, str(exc)[:120])
        return row

    row['route'] = q['route']
    row['coarse_faces'] = q['coarse_faces']
    row['poles'] = q['poles']
    row['coverage'] = round(q['coverage'], 6)
    row['faces'] = q['faces']
    row['min_angle'] = round(q['min_angle'], 4)
    row['max_angle'] = round(q['max_angle'], 4)
    row['aspect_max'] = round(q['aspect_max'], 4)
    row['share_below_20'] = round(q['share_below'], 6)
    row['irregular_interior'] = q['irregular_interior']

    ok, why = hard_floor(q)
    row['hard_floor'] = 'pass' if ok else why
    return row


def measure(only=None, relax=False):
    rows = []
    for label, outer, holes, guides, kwargs in DOMAINS:
        if only and label != only:
            continue
        for tl in RESOLUTIONS:
            rows.append(measure_one(label, outer, holes, guides, kwargs, tl, relax=relax))
    return rows


# --- reporting ---------------------------------------------------------------

HEAD = ('{:18s} {:>4s} {:>13s} {:>7s} {:>5s} {:>6s} {:>6s} {:>7s} {:>7s} '
        '{:>7s} {:>7s} {:>4s}')
LINE = ('{:18s} {:4.2f} {:>13s} {:7d} {:5d} {:6.3f} {:6d} {:7.2f} {:7.2f} '
        '{:7.2f} {:7.1%} {:4d}')


def show(rows):
    print(HEAD.format('domain', 'bg', 'route', 'patches', 'poles', 'cover',
                      'faces', 'min', 'max', 'ARmax', '<20deg', 'irr'))
    print('-' * 108)
    for r in rows:
        if r['route'] == 'EXCEPTION':
            print('{:18s} {:4.2f}  EXCEPTION -- {}'.format(
                r['domain'], r['target_length'], r['error']))
            continue
        print(LINE.format(
            r['domain'], r['target_length'], r['route'], r['coarse_faces'],
            r['poles'], r['coverage'], r['faces'], r['min_angle'],
            r['max_angle'], r['aspect_max'], r['share_below_20'],
            r['irregular_interior']))


def report_hard_floor(rows, old=None):
    """TIER 1, on every run, independent of the baseline and not negotiable.

    Seven rows breach the floor TODAY -- see "Standing defects" in the module
    docstring. They are real degeneracies, not fixable from here, so they are
    recorded in ``baseline.json`` and reported on every run as a standing defect
    list rather than being quietly tolerated.

    What must never happen is a NEW one, and that is what the return value
    counts. A known breach that changes character is caught too: ``hard_floor``
    is a compared field, so its text moving shows up in the diff.
    """
    bad = [r for r in rows if r.get('hard_floor', 'pass') != 'pass']
    print()
    if not bad:
        print('HARD FLOOR: all {} rows clear -- no angle at or near 0 or 180 '
              'degrees, no non-finite aspect ratio.'.format(len(rows)))
        return 0

    known = {}
    if old:
        known = {key_of(r): r.get('hard_floor', 'pass') for r in old}

    fresh = [r for r in bad if known.get(key_of(r), 'pass') == 'pass']
    standing = [r for r in bad if r not in fresh]

    if standing:
        print('HARD FLOOR -- {} KNOWN breach(es), recorded in baseline.json. '
              'Degeneracies, not tastes; each is a standing defect:'.format(
                  len(standing)))
        for r in standing:
            print('  {} @ bg {:.2f}: {}'.format(r['domain'], r['target_length'],
                                                r['hard_floor']))
    if fresh:
        print('HARD FLOOR -- {} NEW breach(es). This change introduced a '
              'degenerate element:'.format(len(fresh)))
        for r in fresh:
            print('  {} @ bg {:.2f}: {}'.format(r['domain'], r['target_length'],
                                                r['hard_floor']))
    return len(fresh)


def key_of(row):
    return (row['domain'], row['target_length'])


def compare(new, old):
    """TIER 2. Every row that moved, with the metric named and both values."""
    old_by_key = {key_of(r): r for r in old}
    new_by_key = {key_of(r): r for r in new}

    moved = []
    for k in sorted(set(old_by_key) | set(new_by_key)):
        a, b = old_by_key.get(k), new_by_key.get(k)
        if a is None:
            moved.append((k, [('row', 'absent from baseline', 'present')]))
            continue
        if b is None:
            moved.append((k, [('row', 'present in baseline', 'absent')]))
            continue
        deltas = []
        for name, tol in FIELDS:
            if name not in a and name not in b:
                continue
            va, vb = a.get(name), b.get(name)
            if tol is None or isinstance(va, str) or isinstance(vb, str):
                if va != vb:
                    deltas.append((name, va, vb))
            elif va is None or vb is None:
                if va != vb:
                    deltas.append((name, va, vb))
            elif abs(float(va) - float(vb)) > tol:
                deltas.append((name, va, vb))
        if deltas:
            moved.append((k, deltas))
    return moved


def report_diff(moved, total):
    if not moved:
        print()
        print('BASELINE: all {} rows match baseline.json exactly.'.format(total))
        return 0
    print()
    print('BASELINE: {} of {} rows MOVED.'.format(len(moved), total))
    print()
    for (domain, tl), deltas in moved:
        print('  {} @ bg {:.2f}'.format(domain, tl))
        for name, was, now in deltas:
            def fmt(v):
                return '{:.4f}'.format(v) if isinstance(v, float) else str(v)
            print('      {:20s} {:>14s}  ->  {:>14s}'.format(
                name, fmt(was), fmt(now)))
    return len(moved)


def check_poles():
    """VERIFY: every metric survives a pseudo-quad, and is doing real work.

    ``12_cables``' ring cable is the live case -- the only domain in the suite
    that produces poles. Both its meshes are checked, because the pole means
    different things in each: in the COARSE layout a pole face is a genuine
    pseudo-quad, the quad ``(p, a, b, p)`` stored with its fourth side collapsed;
    in the DENSE mesh it is the triangle fan that pseudo-quad densifies into.

    The naive reading is printed alongside, because a metric that passes here by
    accident is worth nothing. Reconstituting the phantom fourth corner is what
    a metric written without poles in mind does, and it turns every one of these
    faces into a 0 degree angle and a division by zero.
    """
    from framefield.quality import HARD_MIN_ANGLE, face_angles, mesh_quality
    from framefield.quality import _face_edges

    print('--- pseudo-quads: 12_cables ring cable ---------------------------')
    d = FieldDecomposition.from_boundary(SQUARE, guides=[RING],
                                         target_length=0.5, guide_band=1.5,
                                         guide_weight=5.0)
    coarse = d.decomposition_mesh()
    dense = d.quad_mesh(target_length=QUAD_TARGET)

    failures = []
    for name, mesh in (('coarse', coarse), ('dense', dense)):
        poles = mesh.attributes.get('face_pole') or {}
        q = mesh_quality(mesh)
        print('  {:6s}: {} faces, {} pole face(s), min {:.2f} deg, max {:.2f} deg, '
              'ARmax {:.2f}'.format(name, mesh.number_of_faces(), len(poles),
                                    q['min_angle'], q['max_angle'], q['aspect_max']))
        if not poles:
            failures.append('{}: expected pole faces, found none'.format(name))
            continue

        worst_angle, worst_aspect = 180.0, 0.0
        naive_zero = 0
        for fkey, pole in poles.items():
            pts = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
            angles = face_angles(pts)
            edges = _face_edges(pts)
            worst_angle = min(worst_angle, min(angles))
            worst_aspect = max(worst_aspect, max(edges) / min(edges))

            # what a metric written without poles in mind would have done:
            # put the collapsed corner back and measure the "quad".
            fv = list(mesh.face_vertices(fkey))
            if len(fv) == 3 and pole in fv:
                i = fv.index(pole)
                naive = [mesh.vertex_coordinates(v)
                         for v in fv[:i + 1] + [pole] + fv[i + 1:]]
                if min(face_angles(naive)) <= 0.0:
                    naive_zero += 1

        print('         pole faces read as they ARE: min angle {:.2f} deg, '
              'worst aspect {:.2f} -- all finite'.format(worst_angle, worst_aspect))
        print('         the same faces read as (p, a, b, p): {} of {} report '
              '0.00 deg and divide by zero'.format(naive_zero, len(poles)))

        if worst_angle <= 0.0:
            failures.append('{}: a pole face reports a 0 degree angle'.format(name))
        if worst_aspect != worst_aspect or worst_aspect == float('inf'):
            failures.append('{}: a pole face reports a non-finite aspect'.format(name))
        if worst_angle < HARD_MIN_ANGLE:
            failures.append('{}: pole min angle {:.3f} is under the hard '
                            'floor'.format(name, worst_angle))

    # the other domain that carries poles, per the baseline sweep
    print()
    print('  ellipse @ bg 0.60 -- the suite\'s other pole case:')
    de = FieldDecomposition.from_boundary(ELLIPSE, target_length=0.6)
    ed = de.quad_mesh(target_length=QUAD_TARGET)
    eq = de.quality(ed)
    print('         {} poles in the coarse layout, dense min {:.2f} deg, '
          'ARmax {:.2f}'.format(
              (de.mesh.attributes.get('face_pole') or {}) and
              len(de.mesh.attributes.get('face_pole') or {}),
              eq['min_angle'], eq['aspect_max']))
    if eq['aspect_max'] == float('inf') or eq['min_angle'] <= 0.0:
        failures.append('ellipse @0.60: pole handling failed')

    print()
    if failures:
        for f in failures:
            print('FAIL: {}'.format(f))
        return 1
    print('all pole checks passed -- no spurious 0 degree angle, no division '
          'by zero, on either mesh of either domain.')
    return 0


#: Required separation between `SHARP_TURN` and the nearest real turn angle in
#: the suite, as a ratio. Not a tuning knob -- it is the alarm that says the
#: constant has stopped having margin and somebody must look at it again.
CORNER_MARGIN = 1.25


def check_corners():
    """VERIFY: `repair.SHARP_TURN` still separates sampled curves from corners.

    ``boundary_corners`` splits turns into two regimes at ``SHARP_TURN``: above
    it a vertex is a corner unconditionally, below it the arc-length run-grouping
    decides. That split is only sound while the suite's two populations stay well
    clear of the threshold on both sides, and NOTHING made that check visible
    before -- which is exactly how ``spacing`` came to swallow the comb's teeth
    while nobody noticed.

    So the separation is measured on every run and asserted. If a future domain
    puts a turn near 45 degrees, this fails and says so, instead of silently
    handing that domain to whichever regime it happens to land in.
    """
    from framefield.background import _densify_loop
    from framefield.repair import SHARP_TURN, boundary_corners
    from compas.geometry import angle_vectors, subtract_vectors
    from framefield.repair import _arc_lengths
    from math import degrees, pi

    limit = pi / 12.0
    sharp_deg = degrees(SHARP_TURN)
    print('--- SHARP_TURN separation: {:.1f} deg ---------------------------'.format(
        sharp_deg))

    loops = []
    for label, outer, holes, guides, kwargs in DOMAINS:
        loops.append((label, outer))
        for j, h in enumerate(holes or []):
            loops.append(('{} hole{}'.format(label, j), h))

    seen = set()
    samples, corners_seen = [], []
    print('  {:20s} {:>9s} {:>9s}  {}'.format(
        'domain', 'max<45', 'min>=45', 'reading'))
    for label, loop in loops:
        key = tuple(tuple(p) for p in loop)
        if key in seen:
            continue
        seen.add(key)
        dense = _densify_loop([list(p) for p in loop], min(RESOLUTIONS))
        seg, _cum = _arc_lengths(dense)
        turns = []
        for i in range(len(dense)):
            if seg[(i - 1) % len(dense)] < 1e-12 or seg[i] < 1e-12:
                continue
            b = dense[i]
            a = angle_vectors(subtract_vectors(b, dense[(i - 1) % len(dense)]),
                              subtract_vectors(dense[(i + 1) % len(dense)], b))
            if a > limit:
                turns.append(degrees(a))
        below = [t for t in turns if t < sharp_deg]
        above = [t for t in turns if t >= sharp_deg]
        samples.extend(below)
        corners_seen.extend(above)
        reading = []
        if below:
            reading.append('{} sampled'.format(len(below)))
        if above:
            reading.append('{} corner(s)'.format(len(above)))
        print('  {:20s} {:>9s} {:>9s}  {}'.format(
            label,
            '{:.1f}'.format(max(below)) if below else '-',
            '{:.1f}'.format(min(above)) if above else '-',
            ', '.join(reading) or 'flat'))

    worst_sample = max(samples) if samples else 0.0
    best_corner = min(corners_seen) if corners_seen else 180.0
    margin_lo = sharp_deg / worst_sample if worst_sample else float('inf')
    margin_hi = best_corner / sharp_deg

    print()
    print('  coarsest turn read as a CURVE SAMPLE : {:6.1f} deg  '
          '({:.2f}x below the threshold)'.format(worst_sample, margin_lo))
    print('  shallowest turn read as a CORNER     : {:6.1f} deg  '
          '({:.2f}x above the threshold)'.format(best_corner, margin_hi))
    print('  required margin on both sides        : {:6.2f}x'.format(CORNER_MARGIN))

    # the comb and the slot: the two the threshold exists for
    print()
    for name in ('comb-plate', 'rect with slot', 'stadium'):
        loop = [o for lbl, o, _h, _g, _k in DOMAINS if lbl == name][0]
        dense = _densify_loop([list(p) for p in loop], 0.5)
        print('  {:16s} boundary_corners -> {} corner(s)'.format(
            name, len(boundary_corners(dense))))

    print()
    failures = []
    if margin_lo < CORNER_MARGIN:
        failures.append('a curve sample at {:.1f} deg is only {:.2f}x below '
                        'SHARP_TURN'.format(worst_sample, margin_lo))
    if margin_hi < CORNER_MARGIN:
        failures.append('a real corner at {:.1f} deg is only {:.2f}x above '
                        'SHARP_TURN'.format(best_corner, margin_hi))
    if failures:
        for f in failures:
            print('FAIL: {} -- SHARP_TURN has stopped having margin on this '
                  'suite. Do not tune it silently; read boundary_corners\' '
                  'docstring first.'.format(f))
        return 1
    print('SHARP_TURN separates the two populations with margin to spare on '
          'both sides.')
    return 0


def main():
    if '--poles' in sys.argv:
        return check_poles()
    if '--corners' in sys.argv:
        return check_corners()

    only = None
    if '--domain' in sys.argv:
        only = sys.argv[sys.argv.index('--domain') + 1]

    relax = '--relax' in sys.argv
    rows = measure(only, relax=relax)
    show(rows)

    if relax and '--regenerate' in sys.argv:
        print()
        print('refusing to regenerate under --relax. baseline.json records the '
              'DEFAULT solver; a baseline written from the opt-in one would '
              'silently redefine what every future diff is measured against.')
        return 2

    if '--regenerate' in sys.argv:
        if only:
            report_hard_floor(rows)
            print()
            print('refusing to regenerate from a single domain -- baseline.json '
                  'would lose every other row. Drop --domain.')
            return 2
        report_hard_floor(rows)
        with open(BASELINE, 'w') as f:
            json.dump({'quad_target': QUAD_TARGET,
                       'resolutions': list(RESOLUTIONS),
                       'rows': rows}, f, indent=1, sort_keys=True)
            f.write('\n')
        print()
        print('wrote {} rows to {}'.format(len(rows), os.path.basename(BASELINE)))
        return 0

    if not os.path.exists(BASELINE):
        report_hard_floor(rows)
        print()
        print('no baseline.json yet -- run with --regenerate to create it.')
        return 2

    with open(BASELINE) as f:
        old = json.load(f)['rows']
    if only:
        old = [r for r in old if r['domain'] == only]

    fresh_breaches = report_hard_floor(rows, old)
    moved = compare(rows, old)
    changed = report_diff(moved, len(rows))
    return 1 if (changed or fresh_breaches) else 0


if __name__ == '__main__':
    sys.exit(main())
