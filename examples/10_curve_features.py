"""Curve ("line") features in the skeleton workflow -- Oval, thesis 4.3.2.

Reproduces the four feature types on one square, headless, with the pipeline of
``examples/02_decomposition_discrete_planar.py`` and nothing else:

    boundary_triangulation
        -> SkeletonDecomposition.from_mesh
        -> decomposition_mesh
        -> collect_strips / set_strips_density_target / densification

All four work. Curve features did not, between March 2021 and September 2026.

WHAT WAS BROKEN
---------------
Fig. 4.20 names an extra step: the topological cut along a curve feature leaves
pentagonal and higher-valency faces, which "become quad faces by propagating the
seams of the discrepancies on the curve features". That step is
``SkeletonDecomposition.quadrangulate_polygonal_faces``, and its call site was
commented out.

It had been live since Robin Oval wrote it in February 2019. Commit ``7d69ca0b``
(16 March 2021) turned it off, and its own message says why: *"fix consequences
of new method mesh.vertices_on_boundaries"*. COMPAS had renamed
``vertices_on_boundary()`` and changed its return from a flat list of vertices to
a list of boundary loops. There are three call sites in ``decomposition.py``. Two
were migrated with the same flattening idiom; the third, inside
``quadrangulate_polygonal_faces``, was missed, and the call was commented out
instead. The 2026 compas-2 port then renamed ``mesh_explode`` and
``geometric_key`` INSIDE the dead method without ever running it, so the
un-migrated line survived that pass too.

A second, older bug sat underneath: ``self.mesh = mesh_weld(mesh)`` was inside
the ``for mesh in supermesh.exploded()`` loop, so ``self.mesh`` ended up as the
last component only. That is the defect behind the method's ``# WIP`` marker.

Two further defects, independent of that one:

* ``boundary_triangulation`` did not weld its feature input, which
  ``RhinoSurface.discrete_mapping`` has always done. Two crossing polylines with
  no shared point gave an 11-gon.
* ``store_pole_data`` recorded a pole only for a triangle with a POINT FEATURE at
  a corner. A triangle is a legitimate pseudo-quad, but only if its pole is
  recorded; without one, ``collect_strips`` raises ``KeyError``. Curve extremities
  and medial-axis degeneracies produce triangles that no point feature covers.

Still not implemented: the constrained Delaunay of Fig. 4.20 (Chew 1989). It is
currently harmless -- across seven feature configurations, 0% of feature segments
were missing from the unconstrained triangulation, because a feature is sampled
at the same density as everything else.

Still not implemented: Fig. 4.22's "unwanted triangle" operation, which turns the
pole at a free curve extremity into a two-valent singularity. Extremities here
come out as poles, and a guide can drift from the mesh edges within one patch of
a free tip.

Run:  python examples/10_curve_features.py
Writes ``examples/images/10_curve_features.svg``.
"""
from __future__ import print_function

import contextlib
import io
import math
import os
import traceback

from compas.datastructures.graph.operations.join import graph_polylines
from compas.itertools import pairwise

from compas_singular.algorithms import boundary_triangulation
from compas_singular.algorithms import SkeletonDecomposition
from compas_singular.datastructures import Network


HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'images', '10_curve_features.svg')


# ==============================================================================
# input geometry
# ==============================================================================

def sample(a, b, spacing, include_last=False):
    """Sample the segment a->b at about `spacing`."""
    n = max(1, int(round(math.dist(a, b) / spacing)))
    return [[a[0] + i / n * (b[0] - a[0]), a[1] + i / n * (b[1] - a[1]), 0.0]
            for i in range(n + (1 if include_last else 0))]


def square(side, spacing):
    corners = [[0, 0], [side, 0], [side, side], [0, side]]
    return [p for i in range(4) for p in sample(corners[i], corners[(i + 1) % 4], spacing)]


def circle(cx, cy, r, spacing, close=False):
    n = max(8, int(round(2 * math.pi * r / spacing)))
    pts = [[cx + r * math.cos(2 * math.pi * i / n), cy + r * math.sin(2 * math.pi * i / n), 0.0]
           for i in range(n)]
    return pts + pts[:1] if close else pts


def weld(*chains):
    """Kept for reference: ``boundary_triangulation`` now does this itself.

    See :func:`compas_singular.algorithms.weld_polyline_features`. Curves that
    share a point become one welded node in the network, and ``graph_polylines``
    splits the result at every node that is not 2-valent. Skipping it is why a
    crossing used to come back as an 11-gon.
    """
    return graph_polylines(Network.from_lines([(u, v) for c in chains for u, v in pairwise(c)]))


# ==============================================================================
# the original pipeline, verbatim
# ==============================================================================

def decompose(outer, inner, features, points, target=0.5):
    """Run the pipeline of example 02 and report how far it gets."""
    result = {'trimesh': None, 'coarse': None, 'dense': None, 'faces': {},
              'status': 'ok', 'error': ''}
    try:
        with contextlib.redirect_stdout(io.StringIO()):   # store_pole_data is chatty
            trimesh = boundary_triangulation(outer, inner, features, points)
            result['trimesh'] = trimesh

            decomposition = SkeletonDecomposition.from_mesh(trimesh)
            coarse = decomposition.decomposition_mesh(points)
            result['coarse'] = coarse
            for fkey in coarse.faces():
                n = len(coarse.face_vertices(fkey))
                result['faces'][n] = result['faces'].get(n, 0) + 1

            coarse.collect_strips()
            coarse.set_strips_density_target(target)
            coarse.densification()
            result['dense'] = coarse.get_quad_mesh()
    except Exception as exc:                     # noqa: BLE001 -- reporting, not handling
        result['status'] = 'crash'
        result['error'] = '%s: %s' % (type(exc).__name__, exc)
        result['trace'] = traceback.format_exc()
    return result


def singularities(mesh):
    return [mesh.vertex_coordinates(v)[:2] for v in mesh.vertices()
            if not mesh.is_vertex_on_boundary(v) and len(mesh.vertex_neighbors(v)) != 4]


# ==============================================================================
# a dependency-free SVG of the result, laid out like the thesis figure
# ==============================================================================

class Sheet(object):

    def __init__(self, cell=260, pad=26, cols=4):
        self.cell, self.pad, self.cols = cell, pad, cols
        self.parts, self.n = [], 0

    def _place(self, bbox):
        col, row = self.n % self.cols, self.n // self.cols
        (xmin, ymin), (xmax, ymax) = bbox
        s = (self.cell - 2 * self.pad) / max(xmax - xmin, ymax - ymin)
        ox = col * self.cell + self.pad
        oy = row * (self.cell + 30) + self.pad
        return lambda p: (ox + (p[0] - xmin) * s, oy + (ymax - p[1]) * s)

    def panel(self, title, dense, coarse, red, red_points=()):
        meshes = [m for m in (dense, coarse) if m is not None]
        xs = [m.vertex_coordinates(v)[0] for m in meshes for v in m.vertices()]
        ys = [m.vertex_coordinates(v)[1] for m in meshes for v in m.vertices()]
        for pl in red:
            xs += [p[0] for p in pl]
            ys += [p[1] for p in pl]
        to = self._place(((min(xs), min(ys)), (max(xs), max(ys))))

        def line(a, b, stroke, w):
            (x1, y1), (x2, y2) = to(a), to(b)
            self.parts.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" '
                              'stroke-width="%.2f" stroke-linecap="round"/>'
                              % (x1, y1, x2, y2, stroke, w))

        if dense is not None:
            for u, v in dense.edges():
                line(dense.vertex_coordinates(u), dense.vertex_coordinates(v), '#b3b3b3', 0.6)
        if coarse is not None:
            for u, v in coarse.edges():
                line(coarse.vertex_coordinates(u), coarse.vertex_coordinates(v), '#111', 2.0)
        for pl in red:
            pts = ' '.join('%.2f,%.2f' % to(p) for p in pl)
            self.parts.append('<polyline points="%s" fill="none" stroke="#e01b24" '
                              'stroke-width="2" stroke-linejoin="round"/>' % pts)
        if coarse is not None:
            for p in singularities(coarse):
                x, y = to(p)
                self.parts.append('<circle cx="%.2f" cy="%.2f" r="4.2" fill="#fff" stroke="#000" '
                                  'stroke-width="1.3"/>' % (x, y))
        for p in red_points:
            x, y = to(p)
            self.parts.append('<circle cx="%.2f" cy="%.2f" r="4.2" fill="#e01b24" '
                              'stroke="#e01b24"/>' % (x, y))

        col, row = self.n % self.cols, self.n // self.cols
        self.parts.append('<text x="%.1f" y="%.1f" font-size="11" font-family="sans-serif" '
                          'text-anchor="middle" fill="#111">%s</text>'
                          % (col * self.cell + self.cell / 2,
                             (row + 1) * (self.cell + 30) - 8, title))
        self.n += 1

    def save(self, path):
        rows = (self.n + self.cols - 1) // self.cols
        w, h = self.cols * self.cell, rows * (self.cell + 30)
        if not os.path.isdir(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        with open(path, 'w') as f:
            f.write('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
                    'viewBox="0 0 %d %d">\n<rect width="100%%" height="100%%" fill="#fff"/>\n'
                    % (w, h, w, h))
            f.write('\n'.join(self.parts))
            f.write('\n</svg>\n')
        return path


# ==============================================================================
# the four feature types, on one square
# ==============================================================================

if __name__ == '__main__':

    S = 10.0            # side of the square
    D = 0.2             # Delaunay discretisation (Rhino uses 0.01 * bbox diagonal)
    C = S / 2.0
    L = 0.5             # target edge length of the pattern

    outer = square(S, D)

    cases = [
        ('(a) outer boundary',
         dict(inner=[], features=[], points=[])),

        ('(b) inner boundary',
         dict(inner=[circle(C, C, 1.8, D)], features=[], points=[])),

        ('(c) point feature',
         dict(inner=[], features=[], points=[[C, C, 0.0]])),

        ('(d) Fig 4.17: extremities on the boundary',
         dict(inner=[], features=[sample([0.0, C], [S, C], D, True)], points=[])),

        ('(d) Fig 4.18: extremities off the boundary',
         dict(inner=[], features=[sample([2.0, C], [8.0, C], D, True)], points=[])),

        ('(d) Fig 4.19: several curve features',
         dict(inner=[], features=[sample([C, C], [x, y], D, True)
                                  for x, y in ((0, 0), (S, 0), (S, S), (0, S))], points=[])),

        ('(d) cross, mid-edge to mid-edge',
         dict(inner=[], features=[sample([C, C], [x, y], D, True)
                                  for x, y in ((0, C), (S, C), (C, 0), (C, S))], points=[])),

        ('(d) two parallel curve features',
         dict(inner=[], features=[sample([2.0, 3.5], [8.0, 3.5], D, True),
                                  sample([2.0, 6.5], [8.0, 6.5], D, True)], points=[])),

        ('(d) closed curve feature',
         dict(inner=[], features=[circle(C, C, 2.5, D, close=True)], points=[])),
    ]

    sheet = Sheet()
    print('%-48s %-7s %-30s %s' % ('case', 'result', 'coarse faces by size', 'error'))
    print('-' * 120)
    for title, kw in cases:
        r = decompose(outer, kw['inner'], kw['features'], kw['points'], target=L)
        print('%-48s %-7s %-30s %s' % (title, r['status'], r['faces'], r['error'][:40]))
        sheet.panel('%s -- %s' % (title, 'densified' if r['status'] == 'ok'
                                  else 'CRASH, coarse mesh only'),
                    r['dense'], r['coarse'],
                    [outer + outer[:1]] + [b + b[:1] for b in kw['inner']] + list(kw['features']),
                    kw['points'])
    print()
    print('wrote', sheet.save(OUT))
