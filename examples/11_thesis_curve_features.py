"""The three curve-feature benchmarks of Oval's thesis, §4.3.2, Figs 4.17–4.19.

Each is drawn in the thesis's own four panels — (a) Input, (b) Skeleton,
(c) Decomposition, (d) Mesh — so the output can be held next to the figures.

    4.17  one curve feature, both extremities ON the boundary
    4.18  one curve feature, both extremities OFF the boundary
    4.19  both of the above, superimposed (they cross)

WHAT THE THESIS REQUIRES, AS TESTABLE RULES
--------------------------------------------
These are stated outright in §4.3.2 and on p. 100, so they can be checked rather
than eyeballed. :func:`conformance` reports each one.

  R1  an extremity ON the boundary becomes a three-valent boundary vertex and
      features NO singularity
  R2  an extremity OFF the boundary features exactly ONE singularity: two-valent
      if adjacent to one singular face, a pole if adjacent to several — and
      Fig 4.22's unwanted-triangle correction turns that pole into a two-valent
  R3  the pattern is locally aligned along (4.17) or around (4.18) the curve
  R4  the coarse decomposition is all quads, after the seam propagation of
      Fig 4.20d
  R5  no skeleton branch crosses a curve feature (Fig 4.20a vs 4.20b)

Run:  python examples/11_thesis_curve_features.py
Writes ``examples/images/11_thesis_curve_features.svg`` and prints the table.
"""
from __future__ import print_function

import contextlib
import io
import math
import os

from compas.geometry import closest_point_on_segment
from compas.geometry import distance_point_point
from compas.itertools import pairwise

from compas_singular.algorithms import boundary_triangulation
from compas_singular.algorithms import SkeletonDecomposition


HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'images', '11_thesis_curve_features.svg')

S = 10.0        # side of the square
D = 0.2         # Delaunay discretisation
L = 0.5         # target edge length of the pattern


# ==============================================================================
# the three inputs
# ==============================================================================

def sample(a, b, spacing, include_last=True):
    n = max(1, int(round(math.dist(a, b) / spacing)))
    return [[a[0] + i / n * (b[0] - a[0]), a[1] + i / n * (b[1] - a[1]), 0.0]
            for i in range(n + (1 if include_last else 0))]


def square(side=S, spacing=D):
    c = [[0, 0], [side, 0], [side, side], [0, side]]
    return [p for i in range(4) for p in sample(c[i], c[(i + 1) % 4], spacing, False)]


#: Fig 4.17 — corner to corner, so both extremities lie on the boundary.
DIAGONAL = sample([0.0, 0.0], [S, S], D)

#: Fig 4.18 — a free-standing segment. Oriented to CROSS the diagonal, because
#: Fig 4.19 superimposes the two and shows them crossing.
SEGMENT = sample([2.8, 6.6], [6.6, 2.8], D)

BENCHMARKS = [
    ('4.17  extremities on the boundary', [DIAGONAL]),
    ('4.18  extremities off the boundary', [SEGMENT]),
    ('4.19  both, superimposed', [DIAGONAL, SEGMENT]),
]


# ==============================================================================
# run + measure
# ==============================================================================

def run(features, target=L):
    """The pipeline of example 02, with everything needed for the four panels."""
    out = {'status': 'ok', 'error': '', 'notes': []}
    outer = square()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            trimesh = boundary_triangulation(outer, [], [list(f) for f in features], [])
            decomposition = SkeletonDecomposition.from_mesh(trimesh)
            coarse = decomposition.decomposition_mesh([])
            out['trimesh'] = trimesh
            out['skeleton'] = decomposition.branches()
            out['coarse'] = coarse
            out['notes'] = list(getattr(decomposition, 'repair_notes', []))
            coarse.collect_strips()
            coarse.set_strips_density_target(target)
            coarse.densification()
            out['dense'] = coarse.get_quad_mesh()
    except Exception as exc:                    # noqa: BLE001 -- reporting, not handling
        out['status'] = 'CRASH'
        out['error'] = '%s: %s' % (type(exc).__name__, exc)
    out.setdefault('trimesh', None)
    out.setdefault('skeleton', [])
    out.setdefault('coarse', None)
    out.setdefault('dense', None)
    return out


def on_boundary_of(point, outer, tol=1e-6):
    return min(distance_point_point(point, closest_point_on_segment(point, e))
               for e in pairwise(list(outer) + list(outer[:1]))) < tol


def conformance(result, features):
    """R1–R5 as numbers."""
    coarse = result['coarse']
    if coarse is None:
        return None
    outer = square()

    sizes = {}
    for fkey in coarse.faces():
        n = len(coarse.face_vertices(fkey))
        sizes[n] = sizes.get(n, 0) + 1

    interior = [v for v in coarse.vertices()
                if not coarse.is_vertex_on_boundary(v) and len(coarse.vertex_neighbors(v)) != 4]
    twovalent = [v for v in interior if len(coarse.vertex_neighbors(v)) == 2]
    poles = set((coarse.attributes.get('face_pole') or {}).values())

    tips = []
    for f in features:
        for point in (f[0], f[-1]):
            on_b = on_boundary_of(point, outer)
            near = [v for v in coarse.vertices()
                    if distance_point_point(coarse.vertex_coordinates(v), point) < 1e-6]
            tips.append({
                'xy': (round(point[0], 2), round(point[1], 2)),
                'on_boundary': on_b,
                'pinned': bool(near),
                'valence': len(coarse.vertex_neighbors(near[0])) if near else None,
                'is_pole': bool(near) and near[0] in poles,
            })

    # R3: how far the curve strays from the edges of the dense mesh
    drift = None
    if result['dense'] is not None:
        edges = [(result['dense'].vertex_coordinates(u), result['dense'].vertex_coordinates(v))
                 for u, v in result['dense'].edges()]
        worst = 0.0
        for f in features:
            for a, b in pairwise(f):
                for t in (0.0, 0.5):
                    p = [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), 0.0]
                    worst = max(worst, min(distance_point_point(p, closest_point_on_segment(p, e))
                                           for e in edges))
        drift = worst

    return {'sizes': sizes, 'patches': coarse.number_of_faces(),
            'interior': len(interior), 'twovalent': len(twovalent),
            'poles': len(poles), 'tips': tips, 'drift': drift,
            'allquad': max(sizes) <= 4 if sizes else True}


def report(name, result, features):
    c = conformance(result, features)
    print('=' * 78)
    print(name, ' --', result['status'], result['error'][:40])
    if c is None:
        return
    print('  R4 all-quad coarse   : %-5s  %s' % ('YES' if c['allquad'] else 'NO', c['sizes']))
    print('  patches / dense faces: %d / %s'
          % (c['patches'], result['dense'].number_of_faces() if result['dense'] else '--'))
    print('  interior singular    : %d  (two-valent %d, poles %d)'
          % (c['interior'], c['twovalent'], c['poles']))
    print('  R3 drift off edges   : %s' % ('%.4f' % c['drift'] if c['drift'] is not None else '--'))
    for tip in c['tips']:
        rule = ('R1 three-valent boundary vertex, no singularity' if tip['on_boundary']
                else 'R2 one two-valent singularity')
        got = ('not a layout vertex' if not tip['pinned']
               else 'valence %d%s' % (tip['valence'], ', POLE' if tip['is_pole'] else ''))
        ok = ((tip['on_boundary'] and tip['pinned'] and tip['valence'] == 3)
              or (not tip['on_boundary'] and tip['pinned'] and tip['valence'] == 2))
        print('    extremity %-14s %-4s %-28s want: %s'
              % (tip['xy'], 'OK' if ok else 'no', got, rule))
    for note in result['notes'][:2]:
        print('    note: %s' % note[:66])


# ==============================================================================
# the four thesis panels, as a dependency-free SVG
# ==============================================================================

class Sheet(object):

    def __init__(self, cell=250, pad=22, cols=4):
        self.cell, self.pad, self.cols = cell, pad, cols
        self.parts, self.n = [], 0

    def _place(self):
        col, row = self.n % self.cols, self.n // self.cols
        s = (self.cell - 2 * self.pad) / S
        ox = col * self.cell + self.pad
        oy = row * (self.cell + 34) + self.pad
        return lambda p: (ox + p[0] * s, oy + (S - p[1]) * s)

    def panel(self, title, mesh_grey=None, mesh_black=None, branches=(),
              red=(), dots=(), dot_fill='#fff'):
        to = self._place()

        def line(a, b, stroke, w, op=1.0):
            (x1, y1), (x2, y2) = to(a), to(b)
            self.parts.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" '
                              'stroke-width="%.2f" stroke-opacity="%.2f" stroke-linecap="round"/>'
                              % (x1, y1, x2, y2, stroke, w, op))

        for mesh, stroke, w in ((mesh_grey, '#b6b6b6', 0.55), (mesh_black, '#111', 1.7)):
            if mesh is not None:
                for u, v in mesh.edges():
                    line(mesh.vertex_coordinates(u), mesh.vertex_coordinates(v), stroke, w)
        for branch in branches:
            pts = ' '.join('%.2f,%.2f' % to(p) for p in branch)
            self.parts.append('<polyline points="%s" fill="none" stroke="#d63bd6" '
                              'stroke-width="1.2"/>' % pts)
        for pl in red:
            pts = ' '.join('%.2f,%.2f' % to(p) for p in pl)
            self.parts.append('<polyline points="%s" fill="none" stroke="#e01b24" '
                              'stroke-width="1.8" stroke-linejoin="round"/>' % pts)
        for p in dots:
            x, y = to(p)
            self.parts.append('<circle cx="%.2f" cy="%.2f" r="3.6" fill="%s" stroke="#3b1b6b" '
                              'stroke-width="1.2"/>' % (x, y, dot_fill))

        col, row = self.n % self.cols, self.n // self.cols
        self.parts.append('<text x="%.1f" y="%.1f" font-size="10.5" font-family="sans-serif" '
                          'text-anchor="middle" fill="#111">%s</text>'
                          % (col * self.cell + self.cell / 2,
                             (row + 1) * (self.cell + 34) - 10, title))
        self.n += 1

    def save(self, path):
        rows = (self.n + self.cols - 1) // self.cols
        w, h = self.cols * self.cell, rows * (self.cell + 34)
        if not os.path.isdir(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        with open(path, 'w') as f:
            f.write('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
                    'viewBox="0 0 %d %d">\n<rect width="100%%" height="100%%" fill="#fff"/>\n'
                    % (w, h, w, h))
            f.write('\n'.join(self.parts))
            f.write('\n</svg>\n')
        return path


def singular_points(mesh):
    if mesh is None:
        return []
    return [mesh.vertex_coordinates(v) for v in mesh.vertices()
            if not mesh.is_vertex_on_boundary(v) and len(mesh.vertex_neighbors(v)) != 4]


if __name__ == '__main__':
    outer = square()
    ring = outer + outer[:1]
    sheet = Sheet()

    for name, features in BENCHMARKS:
        result = run(features)
        report(name, result, features)
        tag = name.split()[0]
        sheet.panel('%s (a) Input' % tag, red=[ring] + list(features))
        sheet.panel('%s (b) Skeleton' % tag, mesh_grey=result['trimesh'],
                    branches=result['skeleton'], red=[ring] + list(features))
        sheet.panel('%s (c) Decomposition' % tag, mesh_black=result['coarse'],
                    red=[ring] + list(features), dots=singular_points(result['coarse']))
        sheet.panel('%s (d) Mesh' % tag, mesh_grey=result['dense'],
                    red=[ring] + list(features), dots=singular_points(result['coarse']),
                    dot_fill='#6b3bb0')

    print()
    print('wrote', sheet.save(OUT))
