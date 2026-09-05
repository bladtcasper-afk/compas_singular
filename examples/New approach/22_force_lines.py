"""EXAMPLE 22 -- the CONSTRAINT INPUT: cables and force-line patterns.

``12_cables.py`` establishes that the chain works, using ONE cable on a square.
This asks the next question -- what can actually be fed in, and which patterns
survive the trip from constraint to mesh -- by putting ten force-line patterns
through the same plate at the same settings and measuring every one of them.

The short version, and it is not the flattering one:

    THE FIELD TAKES EVERYTHING.  Ten patterns, 0.0 to 3.8 degrees off the curves
    that were fed in, against 36.9 for a mesh that was given no cable at all.
    There is no pattern in this file the cross field fails to represent,
    including the ones that force singularities.

    THE LAYOUT TAKES ALL TEN, AND PAYS FOR HALF.  Patterns COMPATIBLE with the
    walls come out perfect -- minimum angle untouched at 90, one patch. Patterns
    that force singularities keep the field route but pay in element quality:
    minimum angle 90 -> 12.2, 10.8, 7.2, 7.1 and, worst, 1.1 on the fan.

    Until 2026-08-11 three of the ten lost the field route outright and came
    back as a quad-split triangulation. That was a single bad coarse EDGE, not
    a limit of the method -- see WHAT GOES WRONG 1. The quality figures above
    are the honest remainder.

The dividing line is not cable count, curvature or weight. It is whether the
pattern forces a SINGULARITY. A cable the walls can absorb is free; a cable that
forces winding the walls cannot absorb has to be paid for in the layout, and the
layout is where this front end is weakest.

WHAT GOES RIGHT
---------------

* an orthogonal cable grid is exactly free -- 0.0 degrees off in the field AND
  0.0 in the mesh, minimum angle still 90.0, still one patch. A cross field is
  already the right object for two orthogonal families, so nothing is spent.
* a cable with a free end inside the domain is absorbed smoothly: 0.0 degrees in
  the field, 15.3 in the mesh, no singularity, one patch, minimum angle 59.5.
  A force line does not have to reach a wall to be usable.
* a parallel bundle behaves like one cable, not like three: 1.4 degrees in the
  field, 16.1 in the mesh, minimum angle 64.4, still one patch. Repeating a
  direction costs nothing, because there is no new winding in it.
* the narrow default band beats the wide one everywhere, on both metrics.
  ``12_cables`` quotes 26.9 degrees for its headline cable using
  ``guide_band=3.0``; the same cable at the DEFAULT band measures **18.8**, and
  the sweep in part 3 is monotone the whole way. The wide band is the worst
  setting tested.
* ``relax=True`` is worth more than any other knob on a guided domain, and it
  clears most of the suite's one known-broken row -- see below.

WHAT GOES WRONG
---------------

1.  **FIXED 2026-08-11 -- three patterns used to lose the field route.** The
    arch, the tied arch and the fan came back on route ``'triangulation'`` --
    2520 faces of quad-split triangles that said nothing about the cables. All
    ten now hold the field route.

    The cause was ONE coarse edge, and it was never the layout. On the arch, the
    bottom-left patch is the triangle (0.00, 2.11) (0.00, 0.00) (1.66, 0.00),
    whose third side has both ends on the outer wall and no traced separatrix
    between them. ``edges_to_curves`` therefore handed it the boundary ARC,
    which ran along the bottom wall to the corner and up the left wall -- 46 of
    its 82 points on ``x = 0``. The patch then retraced its own other two sides,
    enclosed nothing, and densified into two ZERO-AREA quads: minimum angle
    0.00, maximum 180.00. ``_acceptable`` rejected the whole 195-face mesh for
    one of them (reported three steps from the cause, as
    ``not all quads (face 12 has 3 sides)``, because a zero-area quad collapses
    to a triangle that is not a registered pole) and fell back.

    The fix is a guard in :meth:`FieldDecomposition._arc_is_one_edge`: a
    boundary arc is only a piece of wall if no OTHER coarse corner lies on it.
    An arc through one has gone round a corner of the domain and come back, so
    the edge takes its chord instead. All 84 rows of ``baseline.json`` are
    unchanged -- no domain in the baseline suite has the pathology.

2.  **The bad rows are still unstable, and that is NOT fixed.** The same fan
    across four backgrounds (part 2) now holds route ``field`` at all four, but
    gives 4 or 6 singularities, 13 to 19 patches, and minimum angles of 0.16,
    1.13, 10.54 and 9.34 degrees -- the hard floor still breached at the
    COARSEST background and clear at the two finest. Nothing is monotone in
    resolution, so nothing here is fixable by meshing finer. Holding the field
    route made these rows honest, not good.

3.  **Root cause of the lost separatrices, measured.** A guide pushes
    singularities against the walls, and ``Tracer.launch_directions`` probes a
    circle of radius ``target_length * 0.45`` at its smallest (``trace.py``
    line 204). A singularity closer to a wall than that has every probe circle
    leave the domain, so all six radii are skipped and it launches NO
    separatrices at all. Confirmed twice in part 4, and in both runs the
    singularities inside ``r_min`` are exactly the ones that failed: fan at
    background 0.60, singularity 0.227 from the wall against a smallest radius
    of 0.270; crossing cables at 0.50, singularity 0.159 against 0.225.

4.  **Two cable families 45 degrees apart cannot both be had, and the loser is
    decided by LIST ORDER.** A cross is 4-fold symmetric, so two directions 45
    degrees apart are the maximally incompatible pair -- that much is geometry,
    not an implementation limit. What IS an implementation limit is the silence:
    ``from_curves`` claims each vertex for the first curve that reaches it, so
    ``guides=[H, D]`` puts H at 0.0 degrees and D at 7.5, and ``guides=[D, H]``
    puts D at 0.0 and H at 7.5. Same domain, same settings, different mesh, no
    warning. 6 of 91 constraints are dropped where they cross.

5.  **Three ways to hand in a guide that is silently ignored** (part 1): a curve
    wholly outside the domain, a degenerate single-point curve, and the boundary
    vertices of any curve that reaches a wall -- ``skip_boundary`` is hard-wired
    ``True`` in ``FieldDecomposition.from_boundary``. The first two return 0
    constraints, no warning, and a mesh that looks fine and ignored the input.

6.  **A band wider than the domain raises.** ``guide_band=8.0`` on a 10-unit
    plate hard-constrains every vertex and ``CrossField.solve`` raises
    ``ValueError: every vertex is hard-constrained; nothing to solve``. Correct,
    and it is the only guide input in this file that fails loudly.

7.  **On the one REALISTIC pattern, the guides made the mesh worse.** Part 6
    traces principal stress trajectories of the Kirsch solution -- a plate with
    a hole in uniaxial tension -- and feeds them in as force lines. Against the
    walls-only mesh, mean mesh-to-load-path angle near the hole went 6.1 -> 8.2
    degrees (soft) and -> 14.4 (hard), and the soft run breached the hard floor
    at 0.288 degrees. The boundary already implied the load paths; the guides
    only added constraint the layout had to pay for.

    ``relax=True`` reverses it: 6.1 -> **5.1** near the hole, 3.9 -> **3.1** in
    the far field, hard floor clear. That is the only setting in this file where
    real force lines beat the walls alone.

WHY ``relax=True`` MATTERS SO MUCH HERE
---------------------------------------

The default solver leaves ``|u|`` free, and the singularities are where it goes
to zero. Part 5 measures what a hard guide does to it: on the unguided square
``|u|`` is 1.000 everywhere, and with a hard diagonal cable it is 1.000 ON the
cable and 0.604 two bands away -- a factor of 1.66. The guide builds a magnitude
RIDGE, the singularities slide off it into the walls, and that is mechanism 3
above. Under ``relax=True`` the magnitude is pinned to 1 everywhere and the
ridge does not exist, so the guide moves only the ANGLE -- which is all it was
ever meant to move.

So the field-model inaccuracy ``19_field_accuracy.py`` measures on the unguided
disc is not a curiosity that only shows up against a published number: adding
guide constraints is precisely the thing that makes it bite.

``relax`` still defaults to ``False`` and this file does not change that -- it
moves rows in ``baseline.json`` that have nothing to do with cables, which is a
separate decision. But on a GUIDED domain it is the first thing to turn on.
Part 2b runs the suite's one known-broken row, ``square+ring cable``, at
``15_baseline``'s own settings under both solvers:

    background      relax=False              relax=True
        0.60    field    3.55 / 144.31    field    3.02 / 177.06
        0.50    field    7.38 / 180.00    field   18.07 / 149.72   floor cleared
        0.40    polygon  5.38 / 180.00    field   12.84 / 137.29   floor cleared
        0.30    polygon  5.01 / 180.00    polygon  3.73 / 180.00   unchanged

Two of that domain's three hard-floor breaches clear and the field route comes
back at 0.40. It is a large partial repair, not a cure.

Run:
    python 22_force_lines.py
    python 22_force_lines.py --no-view          # numbers only
    python 22_force_lines.py fan ring           # a subset, same layout
"""
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas_singular.framefield import from_curves                       # noqa: E402
from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402
from compas_singular.framefield.quality import curve_alignment           # noqa: E402
from compas_singular.framefield.quality import curve_alignment_profile   # noqa: E402
from compas_singular.framefield.quality import hard_floor                # noqa: E402


#: The plate every pattern in part 2 sits on. One domain for all ten, so the
#: PATTERN is the only variable -- a different outline per pattern would make
#: the table a comparison of shapes rather than of constraints.
PLATE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]

#: Background spacing (the resolution the FIELD is solved at) and target quad
#: edge length. Same values as ``12_cables`` so the numbers are comparable.
SPACING = 0.5
QUAD = 1.0

#: The smallest probe circle ``Tracer.launch_directions`` tries, as a multiple
#: of ``target_length`` -- the last entry of its ``attempts`` list,
#: ``trace.py`` line 204. Part 4 compares singularity-to-wall distances against
#: it; if that list ever changes, this must change with it.
PROBE_MIN_SCALE = 0.45


# ---------------------------------------------------------------------------
# force-line patterns
# ---------------------------------------------------------------------------

def line(a, b, n=12):
    """A straight guide, sampled. Sampled rather than given as two points
    because ``from_curves`` measures distance to SEGMENTS and a long segment is
    fine, but ``curve_alignment`` samples the polyline and wants resolution."""
    return [[a[0] + (b[0] - a[0]) * i / n, a[1] + (b[1] - a[1]) * i / n, 0.0]
            for i in range(n + 1)]


def arc(cx, cy, r, a0, a1, n=32):
    return [[cx + r * math.cos(a0 + (a1 - a0) * i / n),
             cy + r * math.sin(a0 + (a1 - a0) * i / n), 0.0] for i in range(n + 1)]


def parabola(x0, x1, y0, rise, n=20):
    """A parabolic cable profile -- the shape a uniformly loaded cable takes."""
    return [[x0 + (x1 - x0) * (i / n), y0 + rise * 4 * (i / n) * (1 - i / n), 0.0]
            for i in range(n + 1)]


def ray(origin, degrees, inset=0.2):
    """A ray from ``origin``, stopped just inside the plate.

    Guides are clipped to the domain rather than run past it. The part outside
    constrains nothing -- there are no vertices there -- so an unclipped ray
    changes no number in this file, but it draws a force line hanging in space
    outside the plate, which is a picture of something that is not happening.
    """
    angle = math.radians(degrees)
    dx, dy = math.cos(angle), math.sin(angle)
    limit = float('inf')
    for position, direction, low, high in ((origin[0], dx, 0.0, 10.0),
                                           (origin[1], dy, 0.0, 10.0)):
        if abs(direction) < 1e-12:
            continue
        edge = (high - inset) if direction > 0 else (low + inset)
        limit = min(limit, (edge - position) / direction)
    return line(origin, [origin[0] + dx * limit, origin[1] + dy * limit, 0.0])


#: ``(label, guides, what it tests)``. Ordered easiest to hardest, which turns
#: out to be the same order as "how much winding does it force".
PATTERNS = [
    ('none', None,
     'the reference -- walls only, which is what every other example uses'),
    ('orthogonal grid',
     [line([0.5, 3, 0], [9.5, 3, 0]), line([0.5, 7, 0], [9.5, 7, 0]),
      line([3, 0.5, 0], [3, 9.5, 0]), line([7, 0.5, 0], [7, 9.5, 0])],
     'two families at 90 deg -- exactly what a cross field already is'),
    ('free end', [line([2, 5, 0], [7, 8, 0])],
     'a cable that stops inside the domain, touching no wall'),
    ('single cable', [line([1, 2, 0], [9, 8, 0])],
     "12_cables' diagonal, for continuity with its numbers"),
    ('parallel x3',
     [line([0.5, 1.5 + k * 3.0, 0], [9.5, 4.5 + k * 3.0, 0]) for k in range(3)],
     'a prestress bundle -- one direction, repeated, no winding'),
    ('crossing 45',
     [line([0.5, 5, 0], [9.5, 5, 0]), line([1, 1, 0], [9, 9, 0])],
     'the pair a cross CANNOT hold: 45 deg apart is maximally incompatible'),
    ('arch', [parabola(0.5, 9.5, 2.0, 5.0)],
     'a curved thrust line -- turning the field through a full arc'),
    ('arch+tie', [parabola(0.5, 9.5, 2.0, 5.0), line([0.5, 2.0, 0], [9.5, 2.0, 0])],
     'a tied arch: a curved member and a straight one meeting at both ends'),
    ('fan x5', [ray([5, 0.2, 0], d) for d in (20, 50, 90, 130, 160)],
     'load paths converging on a support -- a hub IS a singularity'),
    ('ring', [arc(5, 5, 3, 0, 2 * math.pi)],
     'a closed hoop -- forces winding no wall can absorb'),
]


def pattern(label):
    """The guides of one named pattern. ``PATTERNS`` holds triples, so it is
    not a mapping and ``dict(PATTERNS)`` does not build one."""
    for name, guides, _why in PATTERNS:
        if name == label:
            return guides
    raise KeyError(label)


def selected(label, wanted):
    """Whether a pattern was asked for on the command line.

    SUBSTRING match, not equality. Six of the ten names contain a space --
    ``fan x5``, ``orthogonal grid`` -- so an exact match would need shell
    quoting to select any of them, and ``22_force_lines.py fan ring`` would
    silently select nothing at all rather than failing.
    """
    return not wanted or any(w.lower() in label.lower() for w in wanted)


# ---------------------------------------------------------------------------
# measurement
# ---------------------------------------------------------------------------

def field_offset(decomposition, curves):
    """Mean angle between the FIELD and the curves it was given, along them.

    The same convention ``12_cables`` part 1 uses and the same one
    ``curve_alignment`` uses on a mesh, so the field number and the mesh number
    below it are directly comparable. 0 means the field runs on the curve; 45
    is the worst a cross field can do.
    """
    bg = decomposition.background
    offs = []
    for curve in curves or []:
        for a, b in zip(curve, curve[1:]):
            mid = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, 0.0]
            fkey = min(bg.mesh.faces(),
                       key=lambda f: (bg.mesh.face_centroid(f)[0] - mid[0]) ** 2
                       + (bg.mesh.face_centroid(f)[1] - mid[1]) ** 2)
            theta = decomposition.field.angle_in_face(fkey, (1 / 3., 1 / 3., 1 / 3.))
            tangent = math.atan2(b[1] - a[1], b[0] - a[0])
            d = (tangent - theta) % (math.pi / 2)
            offs.append(math.degrees(min(d, math.pi / 2 - d)))
    return statistics.mean(offs) if offs else float('nan')


def mesh_offset(mesh, curves):
    """Mean of :func:`quality.curve_alignment` over several curves."""
    vals = [v for v in (curve_alignment(mesh, c) for c in curves or []) if v == v]
    return statistics.mean(vals) if vals else float('nan')


def build(guides, relax=False, **kw):
    """One pass through the pipeline at the settings part 3 recommends.

    HARD constraints at the DEFAULT band -- ``guide_band=None`` means
    ``background.target_length``, i.e. the curve steers its own neighbourhood
    and nothing else. Part 3 is the measurement behind that choice.
    """
    options = dict(guide_weight=None, guide_band=None)
    options.update(kw)
    d = FieldDecomposition.from_boundary(PLATE, guides=guides,
                                         target_length=SPACING, relax=relax,
                                         **options)
    return d, d.quad_mesh(target_length=QUAD)


# ---------------------------------------------------------------------------
# 1. what the constraint layer accepts
# ---------------------------------------------------------------------------

def part1():
    print('--- 1. what can be fed in, and what is silently dropped ----------')
    print()
    print('    from_boundary(background, weight=None, corner_tolerance=0.1)')
    print('        walls, tangent, HARD. Always applied; not exposed on')
    print('        FieldDecomposition.from_boundary at all.')
    print('    from_curves(background, curves, mode, band, weight, skip_boundary)')
    print('        cables and force lines. mode/band/weight reach the caller as')
    print('        mode=, guide_band=, guide_weight=. skip_boundary does NOT --')
    print('        it is hard-wired True.')
    print()

    diag = line([1, 2, 0], [9, 8, 0])

    print('    mode= is a NO-OP today, and correctly so:')
    for mode in ('tangent', 'perpendicular'):
        d, dense = build([diag], mode=mode)
        print('        {:14s} field {:.4f} deg   mesh {:.4f} deg   {} faces'.format(
            mode, field_offset(d, [diag]), curve_alignment(dense, diag),
            dense.number_of_faces()))
    print('        Identical to four decimals. exp(i*4*theta) is unchanged by a')
    print('        90 deg rotation, so "along the cable" and "across it" are one')
    print('        constraint. They separate only once warp.py makes the frame')
    print('        non-orthogonal -- the parameter is reserved, not dead.')
    print()

    print('    TWO ways to hand in a guide that is silently ignored ENTIRELY,')
    print('    and one that is silently thinned:')
    d0 = FieldDecomposition.from_boundary(PLATE, target_length=SPACING)
    for label, curve in [
            ('wholly outside the domain', line([-8, -8, 0], [-2, -2, 0])),
            ('degenerate single point', [[5, 5, 0], [5, 5, 0]]),
            ('lying along a wall', line([0.0, 0.0, 0], [10.0, 0.0, 0])),
            ('the same curve, mid-plate', line([0.0, 5.0, 0], [10.0, 5.0, 0])),
    ]:
        n = len(from_curves(d0.background, [curve], weight=None))
        print('        {:28s} {:3d} constraints'.format(label, n))
    print('        The first two are worth a warning and do not get one: a mesh')
    print('        comes back that looks fine and used none of the input.')
    print('        The last pair is skip_boundary -- the same curve keeps {} of'
          .format(len(from_curves(d0.background,
                                  [line([0.0, 5.0, 0], [10.0, 5.0, 0])],
                                  weight=None))))
    print('        its constraints mid-plate and {} on the wall. Leaving the'
          .format(len(from_curves(d0.background,
                                  [line([0.0, 0.0, 0], [10.0, 0.0, 0])],
                                  weight=None))))
    print('        wall to from_boundary is right -- the wall must win there --')
    print('        but a cable RUNNING to a wall loses its last vertices the')
    print('        same way, which is exactly where an anchorage would be.')
    print()

    print('    ONE input that fails loudly -- a band wider than the domain:')
    for band in (5.0, 8.0):
        try:
            d, _ = build([diag], guide_band=band)
            n = len(from_curves(d.background, [diag], weight=None, band=band))
            print('        guide_band={:<5} {:3d} of {} vertices constrained'.format(
                band, n, d.background.mesh.number_of_vertices()))
        except Exception as exc:                                  # noqa: BLE001
            print('        guide_band={:<5} {}: {}'.format(
                band, type(exc).__name__, exc))
    print()

    print('    AND ONE that depends on the order of the list. from_curves')
    print('    "claims" each vertex for the first curve that reaches it, so')
    print('    where two guides overlap the later one is simply not applied:')
    horizontal, diagonal = pattern('crossing 45')
    for label, guides in [('guides=[H, D]', [horizontal, diagonal]),
                          ('guides=[D, H]', [diagonal, horizontal])]:
        d, _dense = build(guides)
        print('        {}   H is {:4.1f} deg off,  D is {:4.1f} deg off'.format(
            label, field_offset(d, [horizontal]), field_offset(d, [diagonal])))
    alone = sum(len(from_curves(d0.background, [c], weight=None))
                for c in (horizontal, diagonal))
    together = len(from_curves(d0.background, [horizontal, diagonal], weight=None))
    print('        {} constraints separately, {} together -- {} dropped where'
          .format(alone, together, alone - together))
    print('        they cross. Whichever curve is listed FIRST is satisfied')
    print('        exactly and the other absorbs the whole disagreement.')
    print()
    print('        That two guides 45 deg apart cannot both be had is geometry:')
    print('        a cross is 4-fold symmetric, so 45 deg is the maximally')
    print('        incompatible pair and no cross satisfies both. Deciding it')
    print('        by list position, silently, is not.')


# ---------------------------------------------------------------------------
# 2. ten patterns
# ---------------------------------------------------------------------------

def run_patterns(wanted=None, relax=False):
    """Every pattern through the pipeline once. Returns the rows for part 2."""
    rows = []
    for label, guides, _why in PATTERNS:
        if not selected(label, wanted):
            continue
        try:
            d, dense = build(guides, relax=relax)
        except Exception as exc:                                  # noqa: BLE001
            rows.append((label, guides, None, None, exc))
            continue
        rows.append((label, guides, d, dense, None))
    return rows


def part2(rows, relax=False):
    print('--- 2. ten force-line patterns, one plate, one setting -----------')
    print()
    print('    Hard constraints at the default band, background {}, quads {}.'
          .format(SPACING, QUAD))
    print('    fld / msh: mean degrees off the guides, field and mesh, same')
    print('    scale. The "none" row measures against the single cable, as a')
    print('    reference for what an UNGUIDED mesh scores on a curve.')
    print()
    print('    {:17s} {:>5s} {:>5s} {:>4s} {:>4s} {:>6s} {:>14s} {:>7s} {:>5s}'.format(
        'pattern', 'fld', 'msh', 'sng', 'sep', 'patch', 'route', 'min ang', 'floor'))

    reference = [line([1, 2, 0], [9, 8, 0])]
    for label, guides, d, dense, exc in rows:
        if exc is not None:
            print('    {:17s} RAISED {}: {}'.format(
                label, type(exc).__name__, str(exc)[:50]))
            continue
        q = d.quality()
        ok, _why = hard_floor(q)
        curves = guides or reference
        print('    {:17s} {:5.1f} {:5.1f} {:4d} {:4d} {:6d} {:>14s} {:7.1f} {:>5s}'
              .format(label, field_offset(d, curves), mesh_offset(dense, curves),
                      len(d.field.singularities()), len(d.separatrices),
                      d.mesh.number_of_faces(), q['route'], q['min_angle'],
                      'ok' if ok else 'BAD'))
        for warning in d.warnings():
            if 'launch directions' in warning or 'QUALITY FAILURE' in warning \
                    or 'DISCARDED' in warning:
                print('      ! {}'.format(warning[:96]))

    print()
    print('    Read the min-angle column against the "none" row. Every pattern')
    print('    that forces a singularity costs most of it; every pattern the')
    print('    walls can absorb costs none of it. That is the whole dividing')
    print('    line, and it is not about how many cables there are.')
    print()
    print('    And the bad rows are not merely bad, they are UNSTABLE. The same')
    print('    fan across four background resolutions, nothing else changed:')
    print()
    print('    {:>6s} {:>5s} {:>7s} {:>14s} {:>8s} {:>6s}'.format(
        'bg', 'sing', 'patches', 'route', 'min ang', 'floor'))
    for spacing in (0.6, 0.5, 0.4, 0.35):
        d = FieldDecomposition.from_boundary(
            PLATE, guides=pattern('fan x5'), target_length=spacing,
            guide_weight=None, guide_band=None)
        d.quad_mesh(target_length=QUAD)
        q = d.quality()
        ok, _ = hard_floor(q)
        print('    {:6.2f} {:5d} {:7d} {:>14s} {:8.2f} {:>6s}'.format(
            spacing, len(d.field.singularities()), d.mesh.number_of_faces(),
            q['route'], q['min_angle'], 'ok' if ok else 'BAD'))
    print()
    print('    Route, patch count and minimum angle all jump around, and the')
    print('    COARSEST background is the one that breaches the hard floor.')
    print('    None of it is monotone in resolution, so none of it can be tuned')
    print('    away by meshing finer. Part 4 is why.')


def part2b(wanted=None):
    """The same ten, resolved by the other solver."""
    print()
    print('--- 2b. the same ten under relax=True ----------------------------')
    print()
    print('    {:17s} {:>14s} {:>7s} {:>5s}   {:>14s} {:>7s} {:>5s}'.format(
        '', 'relax=False', 'min ang', 'floor', 'relax=True', 'min ang', 'floor'))
    for label, guides, _why in PATTERNS:
        if not selected(label, wanted):
            continue
        cells = []
        for relax in (False, True):
            try:
                d, _dense = build(guides, relax=relax)
                q = d.quality()
                ok, _ = hard_floor(q)
                cells.append((q['route'], q['min_angle'], 'ok' if ok else 'BAD'))
            except Exception as exc:                              # noqa: BLE001
                cells.append((type(exc).__name__, float('nan'), '-'))
        print('    {:17s} {:>14s} {:7.1f} {:>5s}   {:>14s} {:7.1f} {:>5s}'.format(
            label, cells[0][0], cells[0][1], cells[0][2],
            cells[1][0], cells[1][1], cells[1][2]))
    print()
    print('    relax=True is the single most valuable knob on a guided domain.')
    print('    Part 5 says why: it removes the |u| ridge a hard guide builds.')
    print()
    print('    And it is not only true at THIS file\'s settings. The suite\'s one')
    print('    known-broken row is square+ring cable, which 15_baseline runs at')
    print('    guide_band=1.5, guide_weight=5.0 across four backgrounds. Its')
    print('    committed defects are a 180.00 deg face at 0.5/0.4/0.3 and a')
    print('    route that falls to polygon at 0.4/0.3. Both solvers, same rows:')
    print()
    print('    {:>6s}   {:>14s} {:>7s} {:>7s} {:>5s}   {:>14s} {:>7s} {:>7s} {:>5s}'
          .format('bg', 'relax=False', 'min', 'max', 'floor',
                  'relax=True', 'min', 'max', 'floor'))
    ring = [arc(5, 5, 3, 0, 2 * math.pi, n=24)]
    for spacing in (0.6, 0.5, 0.4, 0.3):
        cells = []
        for relax in (False, True):
            d = FieldDecomposition.from_boundary(
                PLATE, guides=ring, target_length=spacing, relax=relax,
                guide_band=1.5, guide_weight=5.0)
            d.quad_mesh(target_length=QUAD)
            q = d.quality()
            ok, _ = hard_floor(q)
            cells.append((q['route'], q['min_angle'], q['max_angle'],
                          'ok' if ok else 'BAD'))
        print('    {:6.2f}   {:>14s} {:7.2f} {:7.2f} {:>5s}   {:>14s} {:7.2f} '
              '{:7.2f} {:>5s}'.format(spacing, *(cells[0] + cells[1])))
    print()
    print('    Two of the three hard-floor breaches clear, and the field route')
    print('    comes back at 0.4. The 0.3 row does not move. So relax=True is a')
    print('    large partial repair of the suite\'s worst domain, not a cure --')
    print('    and it is still not the default, because turning it on moves')
    print('    rows in baseline.json that have nothing to do with cables.')


# ---------------------------------------------------------------------------
# 3. the knobs
# ---------------------------------------------------------------------------

def part3():
    print('--- 3. band and weight, measured rather than guessed --------------')
    diag = line([1, 2, 0], [9, 8, 0])

    print()
    print('    BAND -- how far from the curve vertices are constrained.')
    print('    {:>7s} {:>6s} {:>8s} {:>8s} {:>8s}'.format(
        'band', 'cons', 'fld off', 'msh off', 'min ang'))
    for band in (None, 1.0, 2.0, 3.0, 5.0, 8.0):
        try:
            d, dense = build([diag], guide_band=band)
            n = len(from_curves(d.background, [diag], weight=None, band=band))
            print('    {:>7s} {:6d} {:8.1f} {:8.1f} {:8.1f}'.format(
                'default' if band is None else '{:.1f}'.format(band), n,
                field_offset(d, [diag]), curve_alignment(dense, diag),
                d.quality()['min_angle']))
        except Exception as exc:                                  # noqa: BLE001
            print('    {:>7s} {}: {}'.format(
                '{:.1f}'.format(band), type(exc).__name__, str(exc)[:52]))
    print()
    print('    The field reaches 0.0 at every band that solves, so a wider band')
    print('    buys NOTHING in the field -- and it costs in the mesh, steadily.')
    print('    The reason is the ceiling 12_cables names: patch boundaries are')
    print('    held fixed, and here they ARE the four walls. A narrow band bends')
    print('    the field only where the cable is; a wide one bends it across the')
    print('    whole plate, where the fixed boundary then has to unbend it.')

    print()
    print('    WEIGHT -- hard vs least-squares, at the default band.')
    print('    {:>7s} {:>8s} {:>8s} {:>8s}'.format(
        'weight', 'fld off', 'msh off', 'min ang'))
    for w in (0.1, 0.5, 1.0, 5.0, None):
        d, dense = build([diag], guide_weight=w)
        print('    {:>7s} {:8.1f} {:8.1f} {:8.1f}'.format(
            'hard' if w is None else '{:.1f}'.format(w),
            field_offset(d, [diag]), curve_alignment(dense, diag),
            d.quality()['min_angle']))
    print()
    print('    Weight SATURATES: by 5.0 it is within half a degree of hard, and')
    print('    nothing above that is worth setting. Below it the guide is a')
    print('    suggestion -- at 0.1 the field is still 33.8 degrees off the')
    print('    cable, barely better than no cable at all, and the plate keeps a')
    print('    minimum angle of 88.2.')
    print()
    print('    This is the ONE knob that trades rather than merely costs. Band')
    print('    spends quality and buys nothing; weight spends quality and buys')
    print('    alignment, monotonically, all the way along. If a cable is')
    print('    advisory rather than structural, weight is where to say so.')


# ---------------------------------------------------------------------------
# 4. why the separatrices go missing
# ---------------------------------------------------------------------------

def part4():
    print('--- 4. the lost separatrices, root cause --------------------------')
    print()
    print('    "singularity N: found 0 launch directions, expected 5" is the')
    print('    warning decomposition.warnings() calls the one that matters. On a')
    print('    guided domain it has one cause, and it is measurable.')
    print()
    print('    Tracer.launch_directions walks a circle around the singularity')
    print('    and counts where the field points radially. It tries six radii,')
    print('    the smallest target_length * {} -- trace.py line 204. A circle'
          .format(PROBE_MIN_SCALE))
    print('    that leaves the domain is skipped. A singularity nearer a wall')
    print('    than the SMALLEST radius has every circle skipped, so it returns')
    print('    the empty best-effort set and launches nothing.')
    print()

    #: Band 1.0 rather than the default here, and the two backgrounds are the
    #: ones that reproduce the failure. The mechanism is not settings-specific,
    #: but WHICH run shows it is -- a singularity has to land inside r_min, and
    #: at other settings the same pattern puts it just outside.
    cases = [
        ('fan x5', pattern('fan x5'), 0.6),
        ('crossing 45', pattern('crossing 45'), 0.5),
    ]
    for label, guides, spacing in cases:
        d = FieldDecomposition.from_boundary(
            PLATE, guides=guides, target_length=spacing,
            guide_weight=None, guide_band=1.0)
        d.quad_mesh(target_length=QUAD)
        mismatch = {f: (got, want)
                    for f, got, want in (d.trace_report.get('arm_mismatch') or [])}
        r_min = spacing * PROBE_MIN_SCALE
        print('    {} at background {:.2f}, band 1.0 -- smallest probe radius '
              '{:.3f}'.format(label, spacing, r_min))
        for fkey, k in d.field.singularities():
            x, y, _ = d.background.mesh.face_centroid(fkey)
            wall = min(x, 10.0 - x, y, 10.0 - y)
            tag = ('lost all {} arms'.format(mismatch[fkey][1])
                   if fkey in mismatch else '')
            print('        {:+d} at ({:5.2f},{:5.2f})  {:.3f} from the wall  {}{}'
                  .format(k, x, y, wall,
                          '<-- INSIDE r_min  ' if wall < r_min else '',
                          tag))
    print()
    print('    Both mismatches, and only those, sit inside the smallest probe')
    print('    radius. Guides are what put them there: a cable pushes the field')
    print("    against the wall, and a cable's hub -- the fan -- lands ON one.")


# ---------------------------------------------------------------------------
# 5. the magnitude ridge
# ---------------------------------------------------------------------------

def part5():
    print('--- 5. what a hard guide does to |u| ------------------------------')
    print()
    print('    field.py leaves |u| completely free when relax=False, and the')
    print('    singularities are where it collapses to zero. So anything that')
    print('    changes the magnitude MOVES the singularities. A hard constraint')
    print('    pins |u| = 1 at the vertices it touches -- and only there.')
    print()
    diag = line([1, 2, 0], [9, 8, 0])
    # The buckets are the CONSTRAINED set and everything more than four bands
    # clear of it -- not an arbitrary radius. Measuring "within 1.0" of a cable
    # constrained to 0.5 mixes pinned and free vertices into one average and
    # understates the ridge it is trying to show.
    print('    on cable: within the band ({}); away: more than {} clear.'.format(
        SPACING, 4 * SPACING))
    print('    {:24s} {:>12s} {:>10s} {:>7s}'.format(
        '', '|u| on cable', '|u| away', 'ratio'))
    for label, guides, kw in [
            ('no guide', None, {}),
            ('hard cable', [diag], {}),
            ('hard cable, relax=True', [diag], dict(relax=True)),
    ]:
        d, _dense = build(guides, **kw)
        bg = d.background
        near, far = [], []
        for vkey in bg.mesh.vertices():
            if bg.mesh.is_vertex_on_boundary(vkey):
                continue
            p = bg.mesh.vertex_coordinates(vkey)
            closest = min(((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) ** 0.5
                          for q in diag)
            if closest <= SPACING:
                near.append(abs(d.field.u[vkey]))
            elif closest > 4 * SPACING:
                far.append(abs(d.field.u[vkey]))
        on, off = statistics.mean(near), statistics.mean(far)
        print('    {:24s} {:12.3f} {:10.3f} {:7.2f}'.format(label, on, off, on / off))
    print()
    print('    The cable builds a RIDGE in the magnitude -- pinned to exactly 1')
    print('    on itself and collapsed to about half of that away from it. The')
    print('    singularities slide off it, and where they land is part 4.')
    print('    relax=True normalises |u| to 1 at every step, so the ridge does')
    print('    not exist and the guide only moves the ANGLE, which is all it was')
    print('    ever meant to move.')
    print()
    print('    19_field_accuracy.py measures this same defect on the UNGUIDED')
    print('    disc against a published singularity radius. This is why it is')
    print('    not an academic point: constraints are what make it bite.')


# ---------------------------------------------------------------------------
# 6. real load paths
# ---------------------------------------------------------------------------

#: Kirsch: an infinite plate with a circular hole, uniaxial tension along x.
#: The hole radius and the plate the trajectories are traced across.
HOLE_RADIUS = 2.0
KIRSCH_PLATE = [[-10, -6, 0], [10, -6, 0], [10, 6, 0], [-10, 6, 0]]
KIRSCH_HOLE = [[HOLE_RADIUS * math.cos(-2 * math.pi * i / 40),
                HOLE_RADIUS * math.sin(-2 * math.pi * i / 40), 0.0]
               for i in range(40)]


def principal_angle(x, y):
    """Major principal stress direction of the Kirsch solution at ``(x, y)``.

    The closed form for a plate with a circular hole under uniaxial tension.
    Returned as an angle because that is all a cross field can use -- and note
    that ONE family is enough: the minor principal direction is everywhere
    perpendicular to the major one, which is exactly the other arm of the cross.
    Feeding both would be feeding the same constraint twice.
    """
    r2 = x * x + y * y
    if r2 < (HOLE_RADIUS * 1.001) ** 2:
        return None
    th = math.atan2(y, x)
    a2, a4, r4 = HOLE_RADIUS ** 2, HOLE_RADIUS ** 4, r2 * r2
    c2, s2 = math.cos(2 * th), math.sin(2 * th)
    srr = 0.5 * (1 - a2 / r2) + 0.5 * (1 + 3 * a4 / r4 - 4 * a2 / r2) * c2
    stt = 0.5 * (1 + a2 / r2) - 0.5 * (1 + 3 * a4 / r4) * c2
    srt = -0.5 * (1 - 3 * a4 / r4 + 2 * a2 / r2) * s2
    c, s = math.cos(th), math.sin(th)
    sxx = srr * c * c + stt * s * s - 2 * srt * s * c
    syy = srr * s * s + stt * c * c + 2 * srt * s * c
    sxy = (srr - stt) * s * c + srt * (c * c - s * s)
    return 0.5 * math.atan2(2 * sxy, sxx - syy)


def trajectory(y0, step=0.2, limit=400):
    """One major-principal trajectory, traced left to right by midpoint RK2.

    The direction field has no orientation, so each step matches the arm to the
    previous one -- the same trick ``trace.py`` uses on the cross field, and for
    the same reason: without it the streamline flips 180 degrees at an arbitrary
    step and folds back on itself.
    """
    points = [[-9.6, y0, 0.0]]
    previous = 0.0
    for _ in range(limit):
        x, y, _ = points[-1]
        angle = principal_angle(x, y)
        if angle is None:
            break
        if math.cos(angle - previous) < 0:
            angle += math.pi
        mid = principal_angle(x + step * 0.5 * math.cos(angle),
                              y + step * 0.5 * math.sin(angle))
        if mid is None:
            break
        if math.cos(mid - angle) < 0:
            mid += math.pi
        nx, ny = x + step * math.cos(mid), y + step * math.sin(mid)
        if not (-9.7 <= nx <= 9.7 and -5.7 <= ny <= 5.7):
            break
        if nx * nx + ny * ny < (HOLE_RADIUS * 1.05) ** 2:
            break
        points.append([nx, ny, 0.0])
        previous = mid
    return points


def load_paths():
    return [trajectory(y0) for y0 in (-5.0, -3.5, -2.4, 2.4, 3.5, 5.0)]


def part6():
    print('--- 6. a REAL force-line pattern: load paths round a hole ---------')
    print()
    print('    Principal stress trajectories of the Kirsch solution, traced and')
    print('    fed in as guides on a 20x12 plate with a radius-{} hole. Not a'
          .format(HOLE_RADIUS))
    print('    pattern chosen to be easy or hard -- the one a designer would')
    print('    actually want the mesh aligned to.')
    print()
    paths = load_paths()
    print('    {} trajectories, {} points each.'.format(
        len(paths), len(paths[0])))
    print()
    print('    Mean mesh-to-load-path angle, split by distance from the hole --')
    print('    the far field is where the walls already imply the answer, and')
    print('    near the hole is where the guides had something to add.')
    print()
    print('    {:26s} {:>10s} {:>10s} {:>14s} {:>7s} {:>5s}'.format(
        'case', 'near hole', 'far field', 'route', 'min ang', 'floor'))

    for label, guides, kw in [
            ('walls only', None, {}),
            ('load paths, soft w=1', paths, dict(guide_weight=1.0)),
            ('load paths, hard', paths, dict(guide_weight=None, guide_band=1.0)),
            ('load paths, soft, relax', paths, dict(guide_weight=1.0, relax=True)),
    ]:
        try:
            options = dict(guide_weight=None, guide_band=None)
            options.update(kw)
            relax = options.pop('relax', False)
            d = FieldDecomposition.from_boundary(
                KIRSCH_PLATE, inner_boundaries=[KIRSCH_HOLE], guides=guides,
                target_length=0.6, relax=relax, **options)
            dense = d.quad_mesh(target_length=QUAD)
            q = d.quality()
            ok, _ = hard_floor(q)
            near, far = [], []
            for path in paths:
                for (px, py), value in curve_alignment_profile(dense, path):
                    bucket = near if px * px + py * py < 25.0 else far
                    bucket.append(value)
            print('    {:26s} {:10.1f} {:10.1f} {:>14s} {:7.1f} {:>5s}'.format(
                label, statistics.mean(near), statistics.mean(far),
                q['route'], q['min_angle'], 'ok' if ok else 'BAD'))
        except Exception as exc:                                  # noqa: BLE001
            print('    {:26s} RAISED {}: {}'.format(
                label, type(exc).__name__, str(exc)[:50]))

    print()
    print('    The walls-only row is already within a few degrees of the load')
    print('    paths, because the plate edges and the hole between them imply')
    print('    most of the answer. Against that, the guides under the default')
    print('    solver make the mesh WORSE at the very curves they asked for --')
    print('    and the soft run breaches the hard floor doing it.')
    print()
    print('    Only relax=True beats the baseline. That is the honest headline')
    print('    of this file: real force lines are worth feeding in, and today')
    print('    they are only worth feeding in with relax=True.')


# ---------------------------------------------------------------------------
# scene
# ---------------------------------------------------------------------------

def view(rows):
    """One row per pattern, four panels: guides, field, layout, mesh.

    Same convention as ``18_pipeline`` -- a group per pattern with the stages
    nested inside, so a stage can be switched off across every row at once.
    The comparison worth making is column 3 against column 4: which patterns
    have a layout that looks like the cables, and which have one that ignored
    them.
    """
    from compas.geometry import Polyline                          # noqa: E402
    from compas_viewer import Viewer                              # noqa: E402

    from compas_singular.framefield.viz import Grid               # noqa: E402
    from compas_singular.framefield.viz import WALL               # noqa: E402
    from compas_singular.framefield.viz import add_boundary       # noqa: E402
    from compas_singular.framefield.viz import add_dense          # noqa: E402
    from compas_singular.framefield.viz import add_field          # noqa: E402
    from compas_singular.framefield.viz import add_guides         # noqa: E402
    from compas_singular.framefield.viz import add_layout         # noqa: E402
    from compas_singular.framefield.viz import add_singularities  # noqa: E402

    stages = ['1 force lines', '2 cross field', '3 coarse layout', '4 quad mesh']
    pitch = 16.0
    viewer = Viewer()
    grid = Grid(pitch=pitch, cols=len(stages))

    for row, (label, guides, d, dense, exc) in enumerate(rows):
        if exc is not None:
            continue
        q = d.quality()
        group = viewer.scene.add_group(name='{}  ({}, min angle {:.1f} deg)'.format(
            label, q['route'], q['min_angle']))

        def panel(step):
            dx, dy = grid.cell(row * len(stages) + step)
            return dx + 0.35 * pitch - 5.0, dy + 0.35 * pitch - 5.0

        dx, dy = panel(0)
        sub = viewer.scene.add_group(name=stages[0], parent=group)
        loop = [[p[0] + dx, p[1] + dy, 0.0] for p in PLATE]
        sub.add(Polyline(loop + loop[:1]), linecolor=WALL, linewidth=3,
                name='plate')
        add_guides(sub, guides, dx, dy)

        dx, dy = panel(1)
        sub = viewer.scene.add_group(name=stages[1], parent=group)
        add_field(sub, d, dx, dy, budget=260)
        add_singularities(sub, d, dx, dy)
        add_boundary(sub, d, dx, dy)
        add_guides(sub, guides, dx, dy)

        dx, dy = panel(2)
        sub = viewer.scene.add_group(name=stages[2], parent=group)
        add_layout(sub, d, dx, dy, faces=True)
        add_guides(sub, guides, dx, dy)

        dx, dy = panel(3)
        sub = viewer.scene.add_group(name=stages[3], parent=group)
        add_dense(sub, dense, dx, dy)
        add_boundary(sub, d, dx, dy)
        add_guides(sub, guides, dx, dy)

    # A wider margin than Grid.frame's default. Ten rows against four columns is
    # the tallest grid in the suite, and at the default 1.1 the first and last
    # rows sit outside the vertical field of view.
    grid.frame(viewer, margin=1.35)
    viewer.show()


def main(argv):
    show = '--no-view' not in argv
    wanted = [a for a in argv[1:] if not a.startswith('-')]

    part1()
    print()
    rows = run_patterns(wanted)
    part2(rows)
    part2b(wanted)
    print()
    part3()
    print()
    part4()
    print()
    part5()
    print()
    part6()
    print()
    print('IN SHORT')
    print('  The constraint input is not the bottleneck -- the field takes every')
    print('  pattern in this file. The LAYOUT is the bottleneck, and it fails in')
    print('  one specific way: a guide pushes singularities into the walls, and')
    print('  a singularity nearer a wall than {} * target_length loses every'.format(
        PROBE_MIN_SCALE))
    print('  separatrix it should have launched. Turn relax=True on before')
    print('  feeding real force lines in; it is the only setting under which')
    print('  they beat the walls alone.')

    if show and rows:
        view(rows)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
