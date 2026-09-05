"""**Pre-flight: is this input sampled finely enough for a cross field?**

There is one way to hand a frame field a problem it will solve perfectly and
wrongly, and no solver parameter can rescue it.

**A cross is 90-degree symmetric.** ``field.wrap_to_period`` folds every angle
difference into plus or minus 45 degrees, because past 45 a direction is better
matched by the next arm round and the two readings are the *same object*. Not an
approximation -- it is what a cross is. So **a boundary turn above 45 degrees
cannot be read as an arc continuing; only as a corner.**

For a real corner that is correct, and it is why a square has no interior
singularity: its 90-degree turns fold to 0, the boundary-tangent field is
constant. For a SAMPLED ARC it is a lie the discretisation told, and nothing
downstream catches it:

* each aliased tip loses exactly one full period of boundary winding;
* :meth:`CrossField.poincare_hopf` reports ``ok=True``, residual 0, at every
  background resolution and under both solvers -- it compares against the
  MEASURED boundary winding, so the field is not violating anything. It is
  faithfully solving the wrong problem;
* ``repair.SHARP_TURN`` is ``pi/4`` too, so ``boundary_corners`` agrees with the
  field -- consistently and wrongly -- and ``hole_launches`` then skips the loop
  as already-cornered.

The measured case: an ellipse hole with semi-axes 4.2 by 0.5 in a disc of radius
5, tip curvature radius ``b^2/a = 0.0595``. At 48 points each tip turns 57.7
degrees, and all three of the above flip at once.

**So the remedy is always to resample the curve finer, and never a solver
knob.** An automated caller that does not know this will spend its whole budget
sweeping ``relax`` and ``guide_band`` against a discretisation artefact. This
module exists to say so first, in one call, before any of that.

**A sharp turn is not automatically a fault, and the test for which is not the
obvious one.** "A corner has straight neighbours" is wrong: a square has four
points, every one of them turning 90 degrees, so every corner's neighbours are
as sharp as it is. The honest position is the one ``boundary_corners`` already
takes -- above ``sharp`` the sample data *cannot* distinguish a polygon from an
arc, because to sample a curve that coarsely you would need four points per
turn, and a shape with four points per turn is a square.

What CAN be distinguished is UNIFORM turning from a PEAK:

* every vertex turning about the same amount is a polygon, and a coarse circle
  at six points really is a hexagon -- reported as ``corners``;
* one vertex turning far more than its neighbours, which are themselves turning,
  is a curvature peak on a curve that is otherwise smooth. That is the aliasing
  fault, and it is what the ellipse tips are: 57.7 degrees between neighbours of
  30.3.

So a square reports four corners and no faults, a hexagon six corners and no
faults, and the ellipse two faults.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from math import acos
from math import ceil
from math import degrees
from math import sqrt


__all__ = ['LIMIT_DEGREES', 'SAFE_DEGREES', 'STRAIGHT_DEGREES', 'UNIFORM_RATIO',
           'turn_angles', 'check_curve', 'check_inputs']


#: Degrees. ``degrees(repair.SHARP_TURN)`` and the cross's half-period, which are
#: the same number for the same reason. Above this a turn IS a corner, to the
#: field and to the corner detector both.
LIMIT_DEGREES = 45.0

#: Recommended headroom. A turn just under the limit is a coin flip against
#: rounding and against a slightly different sampling, so a refinement is sized
#: to land here rather than at 45.
SAFE_DEGREES = 30.0

#: Below this a segment counts as continuing straight, so a sharp turn beside it
#: is an isolated corner on an otherwise flat run.
STRAIGHT_DEGREES = 5.0

#: A sharp turn whose neighbours turn at least this share of it is UNIFORM
#: turning -- a polygon vertex, not a peak. A regular polygon scores 1.0; the
#: measured ellipse tip scores 30.3 / 57.7 = 0.53.
UNIFORM_RATIO = 0.85


def _unit(a, b):
    v = [b[i] - a[i] for i in range(3)]
    n = sqrt(sum(c * c for c in v))
    if n == 0.0:
        return None
    return [c / n for c in v]


def _turn(a, b, c):
    """Degrees the direction turns at ``b``. ``None`` on a zero-length segment."""
    u = _unit(a, b)
    v = _unit(b, c)
    if u is None or v is None:
        return None
    dot = sum(u[i] * v[i] for i in range(3))
    dot = max(-1.0, min(1.0, dot))
    return degrees(acos(dot))


def _normalise(points, closed=None):
    """``(points, closed)`` -- drop a repeated closing point, detect closure."""
    pts = [list(p) for p in points]
    while len(pts) > 1 and all(abs(pts[0][i] - pts[-1][i]) < 1e-9
                               for i in range(min(3, len(pts[0])))):
        pts.pop()
        if closed is None:
            closed = True
    for p in pts:
        while len(p) < 3:
            p.append(0.0)
    return pts, bool(closed)


def turn_angles(points, closed=None):
    """``[(index, degrees), ...]`` -- how far the direction turns at each point.

    A closed curve turns at every point; an open one does not turn at its two
    ends, and they are absent from the result rather than reported as zero.
    """
    pts, closed = _normalise(points, closed)
    n = len(pts)
    if n < 3:
        return []

    out = []
    span = range(n) if closed else range(1, n - 1)
    for i in span:
        a = pts[(i - 1) % n]
        b = pts[i]
        c = pts[(i + 1) % n]
        angle = _turn(a, b, c)
        if angle is not None:
            out.append((i, angle))
    return out


def check_curve(points, name='curve', kind='curve', closed=None,
                limit=LIMIT_DEGREES, safe=SAFE_DEGREES,
                straight=STRAIGHT_DEGREES, uniform=UNIFORM_RATIO):
    """**Is this one curve sampled finely enough?**

    Parameters
    ----------
    points : list[[x, y, z]]
    name : str
        What to call it in the report.
    kind : {'outer', 'hole', 'guide', 'curve'}
    closed : bool, optional
        Detected from a repeated closing point when not given. Boundaries are
        closed; guides usually are not.
    limit, safe, straight : float, optional
        See the module constants.

    Returns
    -------
    dict
        ``ok`` -- no aliased points. ``faults`` -- the aliased ones, each with
        ``index``, ``point``, ``turn`` and ``refine`` (the factor this curve's
        sampling must be multiplied by to bring that turn under ``safe``).
        ``corners`` -- sharp turns that stand between two straight runs, which
        are corners and are fine. ``near_limit`` -- turns between ``safe`` and
        ``limit``, which will alias under a small change of sampling.
        ``max_turn``, ``points``, ``closed``, ``refine`` (the worst factor).
    """
    pts, closed = _normalise(points, closed)
    turns = turn_angles(pts, closed=closed)
    by_index = dict(turns)
    n = len(pts)

    report = {
        'name': name, 'kind': kind, 'points': n, 'closed': closed,
        'ok': True, 'faults': [], 'corners': [], 'near_limit': [],
        'max_turn': 0.0, 'refine': 1,
    }
    if not turns:
        report['reason'] = 'fewer than 3 distinct points -- nothing to measure'
        return report

    report['max_turn'] = max(angle for _i, angle in turns)

    for index, angle in turns:
        if angle <= safe:
            continue
        if angle <= limit:
            report['near_limit'].append({'index': index, 'turn': angle,
                                         'point': pts[index]})
            continue

        # Sharp. Corner or alias? Not decidable from the turn alone -- see the
        # module docstring. Two things ARE corners:
        #   isolated -- neighbours flat, so this is a kink on a straight run;
        #   uniform  -- neighbours turning as much, so this is a polygon vertex
        #               and a six-point circle really is a hexagon.
        # What is left is a PEAK: a vertex turning far more than neighbours that
        # are themselves turning, i.e. high curvature on a smooth curve, sampled
        # too coarsely to survive the cross's half-period.
        before = by_index.get((index - 1) % n, 0.0)
        after = by_index.get((index + 1) % n, 0.0)
        peak = max(before, after)
        entry = {'index': index, 'turn': angle, 'point': pts[index],
                 'neighbour_turns': [before, after]}

        if peak <= straight:
            entry['why'] = 'isolated -- neighbours are straight'
            report['corners'].append(entry)
        elif peak >= angle * uniform:
            entry['why'] = 'uniform -- neighbours turn as much, so this is a polygon vertex'
            report['corners'].append(entry)
        else:
            entry['refine'] = int(ceil(angle / safe))
            entry['why'] = (
                'peak -- turns {:.1f} between neighbours turning {:.1f} and '
                '{:.1f}, so the curve is smooth here and undersampled'.format(
                    angle, before, after))
            report['faults'].append(entry)

    if report['faults']:
        report['ok'] = False
        report['refine'] = max(f['refine'] for f in report['faults'])
    return report


def check_inputs(outer=None, inners=None, guides=None,
                 limit=LIMIT_DEGREES, safe=SAFE_DEGREES,
                 straight=STRAIGHT_DEGREES, uniform=UNIFORM_RATIO):
    """**Check every curve the field will be solved from.** Run this first.

    Parameters
    ----------
    outer : list[[x, y, z]], optional
    inners : list, optional
    guides : list, optional

    Returns
    -------
    dict
        ``ok`` -- nothing aliased anywhere. ``curves`` -- one
        :func:`check_curve` report each. ``faults``, ``near_limit``,
        ``corners`` -- totals. ``worst`` -- the single worst fault, or ``None``.
        ``advice`` -- what to do, in words, because the remedy is never a solver
        parameter and a caller that does not know that will go looking for one.
    """
    reports = []
    if outer is not None:
        reports.append(check_curve(outer, name='outer', kind='outer',
                                   closed=True, limit=limit, safe=safe,
                                   straight=straight, uniform=uniform))
    for i, loop in enumerate(inners or []):
        reports.append(check_curve(loop, name='hole {}'.format(i), kind='hole',
                                   closed=True, limit=limit, safe=safe,
                                   straight=straight, uniform=uniform))
    for i, guide in enumerate(guides or []):
        reports.append(check_curve(guide, name='guide {}'.format(i),
                                   kind='guide', limit=limit, safe=safe,
                                   straight=straight, uniform=uniform))

    faults = [dict(f, curve=r['name']) for r in reports for f in r['faults']]
    near = sum(len(r['near_limit']) for r in reports)
    corners = sum(len(r['corners']) for r in reports)
    worst = max(faults, key=lambda f: f['turn']) if faults else None

    out = {
        'ok': not faults,
        'curves': reports,
        'faults': len(faults),
        'near_limit': near,
        'corners': corners,
        'worst': worst,
        'limit': limit,
        'safe': safe,
    }

    if faults:
        out['advice'] = (
            '{} point(s) on {} curve(s) turn more than {:.0f} degrees inside a '
            'run that is curving. A cross field reads each of those as a CORNER '
            'and cannot do otherwise -- it will solve the resulting problem '
            'exactly, and the problem will not be the one you drew. RESAMPLE '
            'those curves at least {}x finer. No solver parameter -- not relax, '
            'not guide_band, not target_length -- changes this.'.format(
                len(faults), len(set(f['curve'] for f in faults)), limit,
                max(r['refine'] for r in reports)))
    elif near:
        out['advice'] = (
            '{} point(s) turn between {:.0f} and {:.0f} degrees. Nothing is '
            'aliased yet, but a small change of sampling would do it. Worth '
            'refining if the layout comes out with singularities nobody '
            'drew.'.format(near, safe, limit))
    else:
        out['advice'] = 'Sampling is fine. Every turn is either a corner or below {:.0f} degrees.'.format(safe)

    return out
