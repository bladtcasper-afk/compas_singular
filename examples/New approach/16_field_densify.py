"""EXAMPLE 7 -- field-aware densification, measured.

``12_cables`` part 1 measured half the objective working: a hard-constrained
cable puts the FIELD 0.0 degrees off itself. The other half -- the field
reaching the MESH -- is what ``framefield/densify.py`` adds and what this
example measures.

Three checks, in the order they are worth running:

  --constant-field   TRAP 3, and the cheapest possible test of the integrator.
                     On a square with boundary constraints only, the field is
                     EXACTLY constant. Integrating a constant field across a
                     rectangular patch must give back the same axis-aligned
                     grid ``discrete_coons_patch`` gives -- not nearly, exactly.
                     Compares the two meshes vertex by vertex. If this moves,
                     the integrator is wrong and nothing further is meaningful.

  --alignment        The deliverable. Angle between dense mesh edges and the
                     cable, sampled along the cable, guided against unguided.
                     Replaces ``12_cables`` part 2's median-angle column, which
                     read 0.0 for every row only because the meshes it compared
                     were byte-identical.

  --patches          What the densifier did per patch class over the whole
                     domain suite: relaxed, guarded, or skipped as a pole.

Run:
    python 16_field_densify.py               # all three
    python 16_field_densify.py --constant-field
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402
from compas_singular.framefield.quality import curve_alignment           # noqa: E402
from compas_singular.framefield.quality import curve_alignment_profile   # noqa: E402


#: RULE ZERO -- outlines are IMPORTED, never retyped. The suite's ellipse is
#: sampled with 60 points and the obvious ad-hoc transcription uses 48; the
#: result is deterministic, plausible, and about a shape that is not under test.
_SUITE = __import__('15_baseline')

SQUARE = _SUITE.SQUARE
CABLE = _SUITE.CABLE
RING = _SUITE.RING

QUAD_TARGET = _SUITE.QUAD_TARGET
RESOLUTIONS = _SUITE.RESOLUTIONS

#: the kwargs 15_baseline's ``square+cable`` row is measured with
CABLE_KWARGS = dict(_SUITE.DOMAINS[14][4])


def build(outer, holes=None, guides=None, target_length=0.5, field_aware=True,
          **kwargs):
    d = FieldDecomposition.from_boundary(outer, holes, guides=guides,
                                         target_length=target_length, **kwargs)
    d.field_aware = field_aware
    return d, d.quad_mesh(target_length=QUAD_TARGET)


# --- TRAP 3 ------------------------------------------------------------------

def constant_field():
    """A constant field must reproduce the grid. Exactly, not nearly."""
    print('--- TRAP 3: a constant field reproduces the grid -----------------')
    print()
    print('    Square, boundary constraints only. The field is exactly')
    print('    constant, so every alignment residual is already zero on the')
    print('    Coons grid and the solver must have nothing to do.')
    print()
    print('    {:5s} {:>9s} {:>9s} {:>12s} {:>8s} {:>8s} {:>7s}'.format(
        'bg', 'coons', 'field', 'max move', 'min', 'max', 'ARmax'))

    worst = 0.0
    for tl in RESOLUTIONS:
        da, mesh_a = build(SQUARE, target_length=tl, field_aware=False)
        db, mesh_b = build(SQUARE, target_length=tl, field_aware=True)
        moved = max(
            max(abs(mesh_a.vertex_coordinates(v)[k]
                    - mesh_b.vertex_coordinates(v)[k]) for k in range(3))
            for v in mesh_a.vertices())
        worst = max(worst, moved)
        q = db.quality(mesh_b)
        print('    {:5.2f} {:9d} {:9d} {:12.3e} {:8.2f} {:8.2f} {:7.2f}'.format(
            tl, mesh_a.number_of_faces(), mesh_b.number_of_faces(), moved,
            q['min_angle'], q['max_angle'], q['aspect_max']))

    print()
    ok = worst <= 1e-12
    print('    largest vertex movement over all four backgrounds: {:.3e}'.format(worst))
    print('    {}'.format(
        'PASS -- the integrator leaves a constant field alone.' if ok else
        'FAIL -- the integrator moves a constant field. It is wrong.'))
    return 0 if ok else 1


# --- the deliverable ---------------------------------------------------------

def _would_accept(d, coarse, mu):
    """Whether the guard takes this stiffness -- asked without the guard on."""
    d.densify(coarse, stiffness=(mu,), guard=True)
    return bool(d.densify_stats['relaxed'])


def tradeoff():
    """The alignment / element-quality trade, measured rather than asserted.

    ``square+cable`` is the suite's hardest case for this and the one the
    stiffness schedule was calibrated on. Regenerates the table quoted in
    ``densify.STIFFNESS``.
    """
    from compas_singular.framefield import densify as D

    print('--- alignment vs element quality, square+cable @ bg 0.50 ---------')
    print()
    print('    Weakest stiffness first. Weaker = more alignment, worse')
    print('    elements. The accepted row is the first that both settles and')
    print('    passes densify._accepts.')
    print()
    d = FieldDecomposition.from_boundary(SQUARE, guides=[CABLE],
                                         target_length=0.5, **CABLE_KWARGS)
    coarse = d.decomposition_mesh()
    coarse.collect_strips()
    coarse.set_strips_density_target(QUAD_TARGET)
    print('    {:>9s} {:>8s} {:>8s} {:>7s} {:>8s} {:>8s} {:>10s}'.format(
        'stiff', 'min', 'max', 'folds', 'ARmax', 'cable', 'verdict'))
    for mu in D.STIFFNESS:
        # guard off: show what each stiffness PRODUCES, then say separately
        # whether the guard would have taken it. Guarded, every rejected row
        # would print the Coons numbers and the table would say nothing.
        mesh = d.densify(coarse, stiffness=(mu,), guard=False)
        q = d.quality(mesh)
        s = d.densify_stats
        verdict = 'unconverged' if not s['relaxed'] else (
            'ACCEPTED' if _would_accept(d, coarse, mu) else 'rejected')
        print('    {:9.2f} {:8.2f} {:8.2f} {:7d} {:8.2f} {:8.1f} {:>10s}'.format(
            mu, q['min_angle'], q['max_angle'], s['folds'], q['aspect_max'],
            curve_alignment(mesh, CABLE), verdict))
    print()
    print('    Read the accepted row against the one above it: 0.50 buys 3.6')
    print('    degrees of alignment and pays an aspect ratio of 22.3 and a')
    print('    163.75 degree corner for it. _accepts refuses that trade. The')
    print('    ceiling is not the solver -- it is that the patch boundary IS')
    print('    the four walls, so the mesh must turn to meet them however well')
    print('    the field is integrated. See --alignment for where it turns.')


def alignment():
    """TIER 3. Does the cable reach the mesh?"""
    print('--- the deliverable: does the cable reach the MESH? --------------')
    print()
    print('    Angle between dense mesh edges and the cable, sampled along it.')
    print('    Two reference points, both measured in 12_cables part 1:')
    print('      a hard-constrained cable puts the FIELD    0.0 deg off')
    print('      no cable at all leaves it                 36.9 deg off')
    print('    The guided mesh should move decisively toward the first. The')
    print('    unguided mesh must stay near the second.')
    print()
    print('    "mid cable" is the middle half of the samples -- the part clear')
    print('    of the walls. The patch boundary is fixed, so the mesh cannot')
    print('    turn where the cable runs into one, and the whole-cable mean')
    print('    averages that limit in.')
    print()
    print('    {:14s} {:5s} {:>10s} {:>10s} {:>8s} {:>10s}'.format(
        'mode', 'bg', 'unguided', 'guided', 'gain', 'mid cable'))

    rows = []
    for mode in ('tangent', 'perpendicular'):
        for tl in RESOLUTIONS:
            _, plain = build(SQUARE, target_length=tl)
            _, guided = build(SQUARE, guides=[CABLE], target_length=tl,
                              guide_band=3.0, guide_weight=None, mode=mode)
            a = curve_alignment(plain, CABLE)
            b = curve_alignment(guided, CABLE)
            profile = [v for _, v in curve_alignment_profile(guided, CABLE)]
            n = len(profile)
            middle = profile[n // 4: n - n // 4] or profile
            mid = sum(middle) / len(middle)
            rows.append((mode, tl, a, b, mid))
            print('    {:14s} {:5.2f} {:10.1f} {:10.1f} {:8.1f} {:10.1f}'.format(
                mode, tl, a, b, a - b, mid))
    print()
    print('    A cross is 4-fold symmetric, so tangent and perpendicular are')
    print('    one constraint and give the same number -- same reason part 1')
    print('    gives. They separate only for a non-orthogonal frame.')
    return rows


# --- what the densifier did --------------------------------------------------

def patches():
    # RULE ZERO -- the outlines come from 15_baseline, never retyped here.
    mod = __import__('15_baseline')

    print('--- per-patch verdicts over the domain suite ---------------------')
    print()
    print('    relaxed = interior integrated from the field')
    print('    flat    = no interior node to move, or none that moved -- the')
    print('              rectilinear plates, whose field is exactly constant,')
    print('              and every patch of a triangulation backstop, which is')
    print('              densified at 1 and so is 2x2 with no interior at all')
    print('    pole    = pseudo-quad, kept its Coons interior by design')
    print('    guarded = relaxed, then rejected; Coons kept')
    print()
    print('    A "triangulation" row is the backstop being re-densified, not a')
    print('    coarse layout: its patch count is the backstop\'s face count.')
    print()
    print('    {:18s} {:5s} {:>8s} {:>8s} {:>6s} {:>6s} {:>8s} {:>14s}'.format(
        'domain', 'bg', 'patches', 'relaxed', 'flat', 'pole', 'guarded', 'route'))
    for label, outer, holes, guides, kwargs in mod.DOMAINS:
        for tl in RESOLUTIONS:
            try:
                d = FieldDecomposition.from_boundary(
                    outer, holes, guides=guides, target_length=tl, **kwargs)
                d.quad_mesh(target_length=QUAD_TARGET)
            except Exception as exc:
                print('    {:18s} {:5.2f}  {}'.format(
                    label, tl, type(exc).__name__))
                continue
            s = d.densify_stats
            if not s:
                print('    {:18s} {:5.2f} {:>8s}'.format(
                    label, tl, 'route ' + d.route()))
                continue
            print('    {:18s} {:5.2f} {:8d} {:8d} {:6d} {:6d} {:8d} {:>14s}'.format(
                label, tl, s['patches'], s['relaxed'],
                s['flat'] + s.get('coons', 0), s['poles'], s['guarded'],
                d.route()))


#: this example never opens a viewer, but the suite is driven with --no-view
#: across the board, so it must not be mistaken for a section selector
SECTIONS = ('--constant-field', '--alignment', '--tradeoff', '--patches')


def main():
    only = [a for a in sys.argv[1:] if a in SECTIONS]
    status = 0
    if not only or '--constant-field' in only:
        status |= constant_field()
        print()
    if not only or '--alignment' in only:
        alignment()
        print()
    if not only or '--tradeoff' in only:
        tradeoff()
        print()
    if not only or '--patches' in only:
        patches()
    return status


if __name__ == '__main__':
    sys.exit(main())
