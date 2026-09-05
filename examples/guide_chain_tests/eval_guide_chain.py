"""Measure the automatic guide selection: does it pick clean, long polyedge lines?

Runs :func:`compas_singular.editing.guide_chain` over three meshes and fourteen guides,
in all three boundary modes, and prints what each chain is. Nothing is asserted here --
:mod:`test_guide_chain` does that. This is the table you read when you want to know
whether a rule change was an improvement.

How to read it
--------------

``bnd_i``
    Boundary vertices in the MIDDLE of the chain. Must be 0. Anything else is the outline
    being selected, and the outline gets dragged off the wall when the chain is projected.
``bnd_e``
    Boundary vertices at an END. These are anchors: a cable drawn wall to wall reaching
    the wall. Allowed, and free -- an attached boundary vertex is never MOVED onto the
    guide, it is only allowed to slide along its own outline.
``pur``
    The share of the chain that continues topologically straight. 1.00 means the chain is
    a piece of one polyedge -- true by construction now, measured as a check on that.
``align``
    How closely the chain's own edges run along the guide. Low means the guide does not
    follow any course of this mesh, whatever was selected.
``cover``
    The share of the guide's length the chain spans. Low with high ``align`` means the
    chain was cut short; low with low ``align`` means the guide is at odds with the mesh.
``minang`` / ``fold``
    The worst face angle and the number of inverted faces AFTER attaching the chain,
    before any smoothing. This is what the old radius-based selection destroyed; a course
    must not. It tracks ``off_max``: an interior vertex is yanked as far as it sat off the
    guide, so the distance tolerance sets this, and the angle gate is what drops the
    vertices whose line had already stopped following.
``drift``
    How far an attached boundary vertex ends up off its own outline, in edge lengths.
    Must be 0.00: a boundary vertex slides, it does not move.

What this does NOT measure, deliberately: the quality of the mesh after the constrained
smoothing that follows in the command. The control at the bottom shows why -- on these
meshes that smoothing degrades on its own, with no guide attached at all, so a number
taken after it says more about the smoother than about the selection.

Run:

    python examples/guide_chain_tests/eval_guide_chain.py
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys
from math import degrees

from compas.geometry import angle_vectors

from compas_singular.datastructures import automated_boundary_constraints
from compas_singular.datastructures import constrained_smoothing
from compas_singular.editing import GuideCurve
from compas_singular.editing import chain_quality
from compas_singular.editing import collect_polyedges
from compas_singular.editing import guide_chain
from compas_singular.editing import mean_edge_length

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_meshes import cases  # noqa: E402


MODES = ['anchor', 'exclude', 'free']
TOLERANCES = [0.4, 0.6, 0.8, 1.0, 1.5]
ANGLES = [15.0, 20.0, 30.0, 45.0, 90.0]
DAMPING = 0.5

HEADER = ('guide', 'mode', 'n', 'bnd_i', 'bnd_e', 'path', 'pur', 'runs', 'align',
          'cover', 'off_max', 'minang', 'fold', 'drift')
WIDTHS = (10, 8, 4, 6, 6, 5, 5, 5, 6, 6, 8, 7, 5, 6)


def row(values):
    return '  '.join(str(value).rjust(width) for value, width in zip(values, WIDTHS))


def min_face_angle(mesh):
    """The worst face corner angle in the whole mesh, in degrees."""
    worst = None
    for face in mesh.faces():
        vertices = mesh.face_vertices(face)
        count = len(vertices)
        for index in range(count):
            a = mesh.vertex_coordinates(vertices[index - 1])
            b = mesh.vertex_coordinates(vertices[index])
            c = mesh.vertex_coordinates(vertices[(index + 1) % count])
            u = [a[i] - b[i] for i in range(3)]
            v = [c[i] - b[i] for i in range(3)]
            if not any(u) or not any(v):
                return 0.0
            angle = degrees(angle_vectors(u, v))
            if worst is None or angle < worst:
                worst = angle
    return worst


def smoothed(mesh, kmax):
    """The mesh smoothed with a sliding boundary and nothing else constrained."""
    copy = mesh.copy()
    constrained_smoothing(copy, kmax=kmax, damping=DAMPING,
                          constraints=automated_boundary_constraints(copy),
                          algorithm='area')
    return copy


def tolerance_sweep(mesh, guides, polyedges, average):
    """What the tolerance buys and what it costs, on the guides that follow a course.

    The tolerance decides how far a vertex is allowed to sit off the guide, and therefore
    how far it is YANKED when it is projected onto it -- which is the whole of the damage
    the projection does. Too tight and the chain is empty or in pieces; too loose and it
    grabs the next course over.
    """
    print('-' * 116)
    print('tolerance sweep (in mean edge lengths of {:.3f}) -- n / coverage / off_max / minang'.format(average))
    print('-' * 116)
    print('  {:<12}'.format('guide') + ''.join('{:>24}'.format('tol {:.1f}'.format(t)) for t in TOLERANCES))
    for guide_name, points in guides:
        guide = GuideCurve(points)
        cells = []
        for factor in TOLERANCES:
            chain, _ = guide_chain(mesh, guide, tolerance=factor * average,
                                   polyedges=polyedges)
            if not chain:
                cells.append('{:>24}'.format('-'))
                continue
            quality = chain_quality(mesh, chain, guide)
            cells.append('{:>24}'.format('{} / {:.0f}% / {:.2f} / {:.0f}'.format(
                quality['n'], 100 * quality['coverage'], quality['off_max'],
                quality['min_face_angle'])))
        print('  {:<12}'.format(guide_name) + ''.join(cells))
    print()


def angle_sweep(mesh, guides, polyedges, average):
    """What the ANGLE gate does, with the distance gate held wide open at 1.5.

    Distance asks whether a vertex is near the guide; angle asks whether its polyedge is
    still going the guide's way. They disagree exactly where a line runs beside the guide
    and then veers off -- the vertex stays near for a step or two after it has stopped
    following, and that vertex is the one whose projection bends the mesh. 90 is the gate
    turned off.
    """
    print('-' * 116)
    print('angle sweep (degrees), distance held at 1.5 edge lengths -- n / coverage / off_max / minang')
    print('-' * 116)
    print('  {:<12}'.format('guide') + ''.join('{:>24}'.format('{:.0f} deg'.format(a)) for a in ANGLES))
    for guide_name, points in guides:
        guide = GuideCurve(points)
        cells = []
        for angle in ANGLES:
            chain, _ = guide_chain(mesh, guide, tolerance=1.5 * average, max_angle=angle,
                                   polyedges=polyedges)
            if not chain:
                cells.append('{:>24}'.format('-'))
                continue
            quality = chain_quality(mesh, chain, guide)
            cells.append('{:>24}'.format('{} / {:.0f}% / {:.2f} / {:.0f}'.format(
                quality['n'], 100 * quality['coverage'], quality['off_max'],
                quality['min_face_angle'])))
        print('  {:<12}'.format(guide_name) + ''.join(cells))
    print()


def main():
    print(__doc__.split('Run:')[0].strip())
    print()

    for name, mesh, guides in cases():
        polyedges = collect_polyedges(mesh)
        average = mean_edge_length(mesh)
        print('=' * 116)
        print('{}: {} vertices, {} faces, {} singularities, {} polyedges, mean edge {:.3f}'.format(
            name, mesh.number_of_vertices(), mesh.number_of_faces(),
            len(mesh.singularities()), len(polyedges), average))
        print('=' * 116)
        print(row(HEADER))

        for guide_name, points in guides:
            guide = GuideCurve(points)
            for mode in MODES:
                chain, info = guide_chain(mesh, guide, boundary=mode, polyedges=polyedges)
                quality = chain_quality(mesh, chain, guide)
                if not chain:
                    print(row((guide_name, mode, 0, '-', '-', '-', '-', '-', '-', '-',
                               '-', '-', '-', '-')) + '   ' + info['reason'])
                    continue
                print(row((
                    guide_name, mode, quality['n'], quality['boundary_interior'],
                    quality['boundary_ends'], 'ok' if quality['valid_path'] else 'BAD',
                    '{:.2f}'.format(quality['purity']), quality['runs'],
                    '{:.2f}'.format(quality['alignment']),
                    '{:.0f}%'.format(100 * quality['coverage']),
                    '{:.2f}'.format(quality['off_max']),
                    '{:.1f}'.format(quality['min_face_angle']),
                    quality['folded_faces'],
                    '{:.2f}'.format(quality['end_drift']))))
            print()

        tolerance_sweep(mesh, guides, polyedges, average)
        angle_sweep(mesh, guides, polyedges, average)

    smoothing_control()


def smoothing_control():
    """The smoother's own behaviour, with no guide attached to anything.

    Here so that nobody reads a post-smoothing quality number as a verdict on the
    selection. Constrained smoothing with a sliding boundary does not settle on these
    meshes: it improves for a few iterations and then gives the improvement back.
    """
    print('=' * 116)
    print('control: worst face angle after constrained smoothing, boundary sliding, NO guide')
    print('=' * 116)
    print('  {:<14}{:>10}{:>10}{:>10}{:>10}'.format('mesh', 'raw', 'k=10', 'k=50', 'k=100'))
    for name, mesh, _ in cases():
        angles = [min_face_angle(smoothed(mesh, kmax)) for kmax in (10, 50, 100)]
        print('  {:<14}{:>10.1f}'.format(name, min_face_angle(mesh))
              + ''.join('{:>10.1f}'.format(angle) for angle in angles))
    print()
    print('  A guide-attached run has to be read against this row, not against 90 degrees.')


if __name__ == '__main__':
    main()
