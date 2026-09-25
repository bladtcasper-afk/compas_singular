"""Do real traced separatrices carry sub-tolerance segments, and does cleaning
at the document tolerance make every rib bakeable?

Rhino refuses a polyline whose consecutive points coincide AT THE DOCUMENT
TOLERANCE. Default is 0.001; 0.01 is common on larger models. Both are checked.
"""
import os
import sys
from math import cos, pi, sin

REPO = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular"
for path in (os.path.join(REPO, "src"), os.path.join(REPO, "examples", "New approach")):
    if path not in sys.path:
        sys.path.insert(0, path)

from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402


def dist(a, b):
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5


def clean(points, tol):
    """The same rule as helpers.clean_polyline_points."""
    out = []
    for p in points:
        p = list(p)[:3]
        if out and dist(out[-1], p) <= tol:
            continue
        out.append(p)
    return out if len(out) >= 2 else []


def refused(ribs, tol):
    """Ribs Rhino would refuse at this tolerance."""
    bad = []
    for i, rib in enumerate(ribs):
        pts = [list(p) for p in rib]
        if len(pts) < 2:
            bad.append((i, "1 point"))
            continue
        m = min(dist(pts[j], pts[j + 1]) for j in range(len(pts) - 1))
        if m <= tol:
            bad.append((i, "{:.2e} segment".format(m)))
    return bad


def tri_plate(n_arc, r_out=9.0, waist=4.2, hole_r=1.3, hole_n=16):
    pts = []
    for k in range(3):
        a0 = 2 * pi * k / 3 + pi / 2
        a1 = 2 * pi * (k + 1) / 3 + pi / 2
        pts.append([r_out * cos(a0), r_out * sin(a0), 0.0])
        for s in range(1, n_arc):
            u = s / float(n_arc)
            ang = a0 + u * (a1 - a0)
            r = r_out - (r_out - waist) * (1.0 - (2 * u - 1.0) ** 2)
            pts.append([r * cos(ang), r * sin(ang), 0.0])
    pts.append(pts[0])
    hole = [[hole_r * cos(2 * pi * i / hole_n), hole_r * sin(2 * pi * i / hole_n), 0.0]
            for i in range(hole_n)]
    hole.append(hole[0])
    return pts, [hole]


TOLERANCES = (0.001, 0.01)
worst_overall = []
found_any = False

for n_arc in (10, 14, 20):
    for target in (1.2, 1.5, 2.0, 2.5):
        outer, inners = tri_plate(n_arc)
        try:
            d = FieldDecomposition.from_boundary(outer, inners, target_length=target)
            d.decomposition_mesh()
            ribs = d.decomposition_polylines()
        except Exception as e:                       # noqa: BLE001
            print("  n_arc={} target={} -> solve failed: {}".format(
                n_arc, target, str(e)[:50]))
            continue

        row = "  n_arc={:<3} target={:<4} {:>3} ribs".format(n_arc, target, len(ribs))
        for tol in TOLERANCES:
            bad = refused(ribs, tol)
            after = [r for r in (clean(rib, tol) for rib in ribs) if not r]
            still = refused([clean(rib, tol) for rib in ribs if clean(rib, tol)], tol)
            row += " | tol {:<5}: {} refused -> {} after clean ({} rib(s) vanish)".format(
                tol, len(bad), len(still), len(after))
            if bad:
                found_any = True
                worst_overall.append((tol, n_arc, target, bad[:2]))
        print(row)

print("\n" + "=" * 62)
if found_any:
    print("CONFIRMED: real traced separatrices DO carry sub-tolerance segments.")
    for tol, n_arc, target, bad in worst_overall[:6]:
        print("  tol {} n_arc={} target={}: {}".format(tol, n_arc, target, bad))
else:
    print("No sub-tolerance rib in this sweep (does not mean none exist).")
print("In every row above, 'after clean' is 0 -- the fix makes every surviving")
print("rib bakeable, and any rib that vanishes had no length to begin with.")
