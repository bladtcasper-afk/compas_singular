"""Stage 1 of the overview page: domain outlines + baseline figures -> domains.json.

Reads the domain definitions from ``15_baseline.py`` and the measurements from
``baseline.json``. Nothing here is typed by hand -- an outline transcribed from
memory rather than imported once cost two rounds of false corrections (the
suite's ellipse is sampled with 60 points; the obvious guess is 48).

Run ``build_page.py``; it calls this first.
"""
import importlib.util
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)

#: ``15_baseline`` is not a legal module name, so it is loaded by path. It
#: guards its CLI behind ``if __name__ == '__main__'``, so importing it is safe.
_spec = importlib.util.spec_from_file_location(
    'baseline_suite', os.path.join(PARENT, '15_baseline.py'))
suite = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(suite)

#: SVG user units per panel, and the margin kept clear inside it.
BOX = 100.0
PAD = 8.0


def _fit(loops):
    """A transform placing every loop inside the box, aspect preserved.

    SVG's y axis grows downwards, so the y coordinate is flipped -- without it
    every plate is drawn upside down, which on the symmetric ones is invisible.
    """
    points = [p for loop in loops for p in loop]
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    scale = (BOX - 2 * PAD) / max(x1 - x0, y1 - y0)
    ox = PAD + ((BOX - 2 * PAD) - (x1 - x0) * scale) / 2.0
    oy = PAD + ((BOX - 2 * PAD) - (y1 - y0) * scale) / 2.0

    def transform(p):
        return ox + (p[0] - x0) * scale, BOX - (oy + (p[1] - y0) * scale)

    return transform


def _path(loop, transform, close=True):
    out = []
    for i, p in enumerate(loop):
        x, y = transform(p)
        out.append('%s%.2f %.2f' % ('M' if i == 0 else 'L', x, y))
    if close:
        out.append('Z')
    return ' '.join(out)


def build():
    """domains.json: one entry per domain, outlines plus baseline figures."""
    with open(os.path.join(PARENT, 'baseline.json'), 'r') as f:
        rows = json.load(f)['rows']

    out = []
    for label, outer, holes, guides, _kwargs in suite.DOMAINS:
        transform = _fit([outer] + list(holes or []) + list(guides or []))
        mine = [r for r in rows if r['domain'] == label]
        mine.sort(key=lambda r: -r['target_length'])
        out.append({
            'name': label,
            'outer': _path(outer, transform),
            'holes': [_path(h, transform) for h in (holes or [])],
            # a guide is an open curve: closing it would draw a chord back
            'guides': [_path(g, transform, close=False) for g in (guides or [])],
            'routes': sorted({r['route'] for r in mine}),
            'patches': sorted({r['coarse_faces'] for r in mine}),
            'poles': sorted({r['poles'] for r in mine}),
            'minang': [min(r['min_angle'] for r in mine),
                       max(r['min_angle'] for r in mine)],
            'maxang': [min(r['max_angle'] for r in mine),
                       max(r['max_angle'] for r in mine)],
            'aspect': max(r['aspect_max'] for r in mine),
            'cover': [min(r['coverage'] for r in mine),
                      max(r['coverage'] for r in mine)],
            'floor': 'fail' if any(r['hard_floor'] != 'pass' for r in mine) else 'pass',
        })

    with open(os.path.join(HERE, 'domains.json'), 'w') as f:
        json.dump(out, f, indent=1)
    return out


if __name__ == '__main__':
    domains = build()
    print('%d domains from %s' % (len(domains), os.path.join(PARENT, 'baseline.json')))
    for d in domains:
        print('  %-18s %-14s patches %-10s min %.1f-%.1f  floor %s' % (
            d['name'], '/'.join(d['routes']), d['patches'],
            d['minang'][0], d['minang'][1], d['floor']))
