"""Milestone 1b -- separatrix tracing.

Checks step 5 of the plan in ``frame-field-front-end-module-layout.md``:

* every singularity launches ``4 - index`` separatrices (3 arms for valence-3,
  5 for valence-5) -- an independent check on the index computation;
* every trace terminates for a legitimate reason (boundary or another
  singularity), not by hitting the arc-length cap or spiralling;
* separatrices leave a singularity at roughly equal angles.

The square, L-shape and annulus have no singularities at all, so they have no
separatrices either -- for those the domain is a single patch and the layout is
trivially the boundary. The interesting cases are the disc and the pentagon.

Opens a compas_viewer scene with the traces drawn on the field that produced them.

Run:

    python 02_separatrices.py
    python 02_separatrices.py --no-view    # checks only
"""
import os
import sys
from math import degrees

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from compas_singular.framefield import BackgroundMesh  # noqa: E402
from compas_singular.framefield import CrossField      # noqa: E402
from compas_singular.framefield import from_boundary   # noqa: E402
from compas_singular.framefield import from_curves     # noqa: E402
from compas_singular.framefield.trace import Tracer    # noqa: E402

SQUARE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]
PENTAGON = [[0, 0, 0], [10, 0, 0], [13, 7, 0], [5, 12, 0], [-2, 7, 0]]


def circle(radius=5.0, n=48, centre=(5.0, 5.0)):
    from math import cos, pi, sin
    return [[centre[0] + radius * cos(2 * pi * i / n),
             centre[1] + radius * sin(2 * pi * i / n), 0.0] for i in range(n)]


def run(name, boundary, holes=None, guides=None, target_length=0.5,
        mode='perpendicular'):
    bg = BackgroundMesh.from_boundary(boundary, holes, target_length=target_length)
    constraints = from_boundary(bg)
    if guides:
        constraints += from_curves(bg, guides, mode=mode)
    field = CrossField.solve(bg, constraints)

    tracer = Tracer(field)
    seps, report = tracer.separatrices()

    print('--- {} ---'.format(name))
    print('  singularities : {}'.format(
        [(f, '{:+d}'.format(k)) for f, k in field.singularities()] or 'none'))
    print('  separatrices  : {}  reasons {}'.format(report['count'], report['reasons']))
    if report['arm_mismatch']:
        print('  ARM MISMATCH  : {}'.format(report['arm_mismatch']))
    lengths = [sum(((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5
                   for a, b in zip(s.points, s.points[1:])) for s in seps]
    if lengths:
        print('  lengths       : min {:.2f}  max {:.2f}  (domain diagonal {:.2f})'.format(
            min(lengths), max(lengths), tracer.diagonal))
    print('  {}'.format('OK' if report['ok'] else 'FAIL'))
    return bg, field, tracer, seps, report


def launch_angles(field, tracer):
    """Angles between consecutive separatrix launches at each singularity.

    They should be roughly 360/(4-index) apart. A pair much closer than that
    means the radial-crossing search found a spurious root, which would send two
    separatrices down the same corridor and collapse a patch.
    """
    from math import atan2, pi
    out = []
    for fkey, k in field.singularities():
        centre = tracer.mesh.face_centroid(fkey)
        dirs = tracer.launch_directions(centre, k)
        angles = sorted(atan2(d[1], d[0]) for d in dirs)
        if len(angles) < 2:
            continue
        gaps = [degrees((b - a) % (2 * pi))
                for a, b in zip(angles, angles[1:] + angles[:1])]
        out.append((fkey, k, len(dirs), round(min(gaps), 1), round(max(gaps), 1)))
    return out


def main():
    results = {}
    results['square'] = run('square (no singularities expected)', SQUARE)
    results['disc'] = run('disc', circle())
    results['pentagon'] = run('pentagon', PENTAGON)

    guide = [[1.0, 2.0, 0.0], [5.0, 5.0, 0.0], [9.0, 8.0, 0.0]]
    results['guided'] = run('square + diagonal guide', SQUARE, guides=[guide])

    print()
    for name in ('disc', 'pentagon', 'guided'):
        _, field, tracer, _, _ = results[name]
        rows = launch_angles(field, tracer)
        if rows:
            print('launch gaps {:9s}: {}'.format(
                name, [(f, '{:+d}'.format(k), n, mn, mx) for f, k, n, mn, mx in rows]))

    print()
    bad = [k for k, r in results.items() if not r[4]['ok']]
    if bad:
        print('FAILED: {}'.format(', '.join(bad)))
    else:
        print('all checks passed')

    if '--no-view' not in sys.argv:
        show(results)

    return 1 if bad else 0


def show(results):
    """Traced separatrices fanning out of each singularity."""
    from compas_viewer import Viewer

    from compas_singular.framefield.viz import (Grid, add_background, add_boundary, add_field,
                                                add_singularities)
    from compas.geometry import Polyline
    from compas.colors import Color

    viewer = Viewer()
    grid = Grid(pitch=15.0, cols=2)

    for i, (name, (bg, field, tracer, seps, report)) in enumerate(results.items()):
        dx, dy = grid.cell(i)
        shim = _Shim(bg, field)
        group = viewer.scene.add_group(name='{} -- {} separatrices, {}'.format(
            name, report['count'], report['reasons'] or 'none'))
        add_background(group, shim, dx, dy)
        add_field(group, shim, dx, dy, budget=110)
        add_boundary(group, shim, dx, dy)
        # drawn straight from the TRACES, before repair snaps or decimates them,
        # so what the tracer actually produced is what is on screen
        for j, s in enumerate(seps):
            if len(s.points) < 2:
                continue
            group.add(Polyline([[p[0] + dx, p[1] + dy, 0.04] for p in s.points]),
                      linecolor=Color.orange(), linewidth=4,
                      name='separatrix {} (ended: {})'.format(j, s.reason))
        add_singularities(group, shim, dx, dy)

    grid.frame(viewer)
    viewer.show()


class _Shim(object):
    """Two-attribute stand-in -- see 01_field.py."""

    def __init__(self, background, field):
        self.background = background
        self.field = field


if __name__ == '__main__':
    sys.exit(main())
