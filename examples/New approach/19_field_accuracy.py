"""CHECK -- is the FIELD right? The one external number this module has.

Every other check here is self-consistency. ``poincare_hopf`` compares the field
against itself; ``15_baseline.py`` compares today's numbers against yesterday's.
Neither can notice that the field has been solving a slightly different problem
than the one intended, because both would go on agreeing with themselves while
it did.

This script uses a number from outside the repository. For a boundary-aligned
cross field on a disc, the four index ``+1`` singularities sit at

    r ~ 0.85    (fraction of the disc radius)

reached two independent ways in the literature: Beaufort et al. 2017 characterise
the minimisers of the renormalized energy on circular domains in terms of the
singularity radius, and Dai/Qiao/Wang 2026 section 6.2 recover the same value
both from an MBO reference solution and from their own scheme as the
regularisation is refined. They also report the failure signature when the
discretisation takes over instead of the model: the radii COLLAPSE, to r ~ 0.2.

What this found
---------------

``CrossField.solve(relax=False)`` -- the original solver -- puts them at
0.34 to 0.57, and the radius drifts steadily INWARD as the background is
refined, with no sign of a limit. On a shape with four-fold symmetry the four
radii are not even equal to each other. Doubling the boundary sampling changes
nothing, so this is the model rather than the discretisation.

``CrossField.solve(relax=True)`` puts them at 0.85 at every resolution.

Read the ``relax=False`` column as a standing defect, not as a failure of this
script: it is the default solver, and the whole committed baseline was measured
with it. See ``paper-diffusion-generated-crossfields.md`` section 3.

Run:
    python 19_field_accuracy.py              # the disc radius check
    python 19_field_accuracy.py --winding    # Poincare-Hopf, both solvers, 16 domains
    python 19_field_accuracy.py --timing     # cost of the relaxation

No viewer -- numbers only.
"""
import sys
import time
from math import cos
from math import hypot
from math import pi
from math import sin

from compas_singular.framefield import BackgroundMesh
from compas_singular.framefield import CrossField


# --- the reference -----------------------------------------------------------

#: Published singularity radius for a boundary-aligned cross field on a disc,
#: as a fraction of the disc radius. See the module docstring for the two
#: independent derivations.
REFERENCE_RADIUS = 0.85

#: How far from it we accept. The singularity is localised to a background
#: TRIANGLE (``trace.py`` takes its centroid), so on a disc of radius 5 at
#: target_length 0.5 one triangle is already 0.1 of the radius. Anything inside
#: 0.10 is agreement; the defect this script exists to catch is a factor of two.
RADIUS_TOLERANCE = 0.10

CENTRE = (5.0, 5.0)
RADIUS = 5.0
RESOLUTIONS = (0.50, 0.35, 0.25)
SAMPLINGS = (48, 96)


def disc(n):
    return [[CENTRE[0] + RADIUS * cos(2 * pi * i / n),
             CENTRE[1] + RADIUS * sin(2 * pi * i / n), 0.0] for i in range(n)]


def singularity_radii(field):
    """Normalised radius and index of every singular face, sorted by radius."""
    mesh = field.background.mesh
    out = []
    for fkey, k in field.singularities():
        p = mesh.face_centroid(fkey)
        out.append((hypot(p[0] - CENTRE[0], p[1] - CENTRE[1]) / RADIUS, k))
    out.sort()
    return out


# --- the check ---------------------------------------------------------------

def radius_check():
    print('=' * 78)
    print('DISC SINGULARITY RADIUS -- measured against the published {:.2f}'.format(
        REFERENCE_RADIUS))
    print('=' * 78)
    print()
    print('{:>6s} {:>6s} {:>8s} {:>6s} {:>5s}  {:s}'.format(
        'pts', 'tl', 'solver', 'sing', 'err', 'radii'))
    print('-' * 78)

    failures = []
    for n in SAMPLINGS:
        outline = disc(n)
        for tl in RESOLUTIONS:
            background = BackgroundMesh.from_boundary(outline, target_length=tl)
            for relax in (False, True):
                field = CrossField.solve(background, relax=relax)
                radii = singularity_radii(field)
                worst = max((abs(r - REFERENCE_RADIUS) for r, _ in radii),
                            default=float('nan'))
                print('{:>6d} {:>6.2f} {:>8s} {:>6d} {:>5.2f}  {:s}'.format(
                    n, tl, 'relax' if relax else 'plain', len(radii), worst,
                    ' '.join('{:.3f}({:+d})'.format(r, k) for r, k in radii)))

                if relax:
                    if len(radii) != 4:
                        failures.append('{} pts, tl {}: {} singularities, expected 4'
                                        .format(n, tl, len(radii)))
                    elif worst > RADIUS_TOLERANCE:
                        failures.append('{} pts, tl {}: worst radius error {:.3f} > {:.2f}'
                                        .format(n, tl, worst, RADIUS_TOLERANCE))
            print()

    print('-' * 78)
    print('The `plain` rows are the DEFAULT solver and are the standing defect:')
    print('radius ~0.45 against 0.85, drifting inward under refinement, and the')
    print('four radii unequal on a four-fold symmetric shape. Only the `relax`')
    print('rows are asserted -- see paper-diffusion-generated-crossfields.md.')
    print()

    if failures:
        print('FAIL')
        for line in failures:
            print('  ' + line)
    else:
        print('PASS -- relax=True lands within {:.2f} of {:.2f} at every resolution'
              .format(RADIUS_TOLERANCE, REFERENCE_RADIUS))
    return not failures


# --- the winding check, both solvers -----------------------------------------

def winding_check():
    """Poincare-Hopf on the whole domain suite, under both solvers.

    The relaxation must not change the TOPOLOGY it starts from -- it moves
    singularities, it does not invent or destroy them. Index sum against
    boundary winding is the statement of that, and it is the check that a broken
    relaxation would fail loudly.
    """
    domains = _suite()

    print('=' * 78)
    print('POINCARE-HOPF -- interior index sum vs boundary winding, both solvers')
    print('=' * 78)
    print()
    print('{:18s} {:>10s} {:>10s} {:>10s} {:>10s} {:>7s}'.format(
        'domain', 'plain sum', 'winding', 'relax sum', 'winding', 'faces'))
    print('-' * 78)

    failures = []
    for label, outer, holes in domains:
        background = BackgroundMesh.from_boundary(outer, holes, target_length=0.5)
        reports = {}
        for relax in (False, True):
            field = CrossField.solve(background, relax=relax)
            reports[relax] = field.report()
        a, b = reports[False], reports[True]
        print('{:18s} {:>10d} {:>10d} {:>10d} {:>10d} {:>3d}/{:<3d}'.format(
            label, a['interior_index_sum'], a['boundary_winding'],
            b['interior_index_sum'], b['boundary_winding'],
            a['singular_faces'], b['singular_faces']))
        for name, report in (('plain', a), ('relax', b)):
            if not report['ok']:
                failures.append('{} ({}): index sum {} != winding {}'.format(
                    label, name, report['interior_index_sum'],
                    report['boundary_winding']))
        if a['boundary_winding'] != b['boundary_winding']:
            failures.append('{}: relaxation changed the boundary winding {} -> {}'
                            .format(label, a['boundary_winding'], b['boundary_winding']))

    print()
    if failures:
        print('FAIL')
        for line in failures:
            print('  ' + line)
    else:
        print('PASS -- every domain balances under both solvers, and the')
        print('relaxation leaves the boundary winding alone.')
    return not failures


# --- cost --------------------------------------------------------------------

def timing_check():
    print('=' * 78)
    print('COST -- the relaxation against the single solve it wraps')
    print('=' * 78)
    print()
    print('{:>6s} {:>8s} {:>10s} {:>10s} {:>7s} {:>10s}'.format(
        'tl', 'vertices', 'plain (s)', 'relax (s)', 'iters', 'residual'))
    print('-' * 78)
    for tl in RESOLUTIONS:
        background = BackgroundMesh.from_boundary(disc(48), target_length=tl)
        t0 = time.time()
        CrossField.solve(background, relax=False)
        plain = time.time() - t0
        t0 = time.time()
        field = CrossField.solve(background, relax=True)
        relaxed = time.time() - t0
        print('{:>6.2f} {:>8d} {:>10.3f} {:>10.3f} {:>7d} {:>10.1e}'.format(
            tl, background.mesh.number_of_vertices(), plain, relaxed,
            field.iterations, field.residual))
    print()
    print('The matrix is factorised once, so an iteration is a pair of')
    print('triangular solves. Against trace + arrangement + repair this is free.')
    return True


def _suite():
    """The 15_baseline domain suite, imported rather than retyped.

    Retyping a domain definition is how the ellipse trap happened -- the suite
    discretises it with 60 points and the obvious ad-hoc transcription uses 48,
    which is a different shape that measures plausibly and is not in the suite.
    """
    import importlib
    module = importlib.import_module('15_baseline')
    return [(label, outer, holes) for label, outer, holes, _guides, _kw in module.DOMAINS]


def main():
    ok = True
    if '--winding' in sys.argv:
        ok = winding_check() and ok
    elif '--timing' in sys.argv:
        ok = timing_check() and ok
    else:
        ok = radius_check() and ok
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
