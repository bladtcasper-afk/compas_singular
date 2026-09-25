"""Every strip of several layouts: what does deleting it actually produce?

The question the gate depends on: is 'non-manifold after delete' a pole thing,
a boundary thing, or everywhere?
"""
import os
import sys
from math import cos, pi, sin

REPO = r"C:\Users\Casper\libraries\carbcomn\compas_singular\compas_singular"
for path in (os.path.join(REPO, "src"), os.path.join(REPO, "examples", "New approach")):
    if path not in sys.path:
        sys.path.insert(0, path)

from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402
from compas_singular.datastructures.mesh_quad.grammar_pattern import (   # noqa: E402
    collateral_strip_deletions, delete_strip, total_boundary_deletions)


def circle(cx, cy, r, n=24):
    pts = [[cx + r * cos(2 * pi * i / n), cy + r * sin(2 * pi * i / n), 0.0]
           for i in range(n)]
    return pts + pts[:1]


def box(cx, cy, r):
    return [[cx - r, cy - r, 0.0], [cx + r, cy - r, 0.0], [cx + r, cy + r, 0.0],
            [cx - r, cy + r, 0.0], [cx - r, cy - r, 0.0]]


SQUARE = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0],
          [0.0, 0.0, 0.0]]
L = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 4.0, 0.0], [4.0, 4.0, 0.0],
     [4.0, 10.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 0.0]]


def sides(mesh):
    out = {}
    for f in mesh.faces():
        n = len(mesh.face_vertices(f))
        out[n] = out.get(n, 0) + 1
    return out


def touches_pole(mesh, skey):
    """Does any FACE of the strip carry a pole -- not just its end edges."""
    fp = mesh.attributes.get("face_pole") or {}
    return any(f in fp for f in mesh.strip_faces(skey))


def sweep(label, outer, inners):
    d = FieldDecomposition.from_boundary(outer, inners, target_length=1.0)
    d.decomposition_mesh()
    base = d.mesh
    base.collect_strips()
    print("\n{}".format(label))
    print("  {} patches, {} strips, poles {}, sides {}".format(
        base.number_of_faces(), len(list(base.strips())),
        len(base.poles()), sides(base)))
    print("  {:>4} {:>6} {:>6} {:>7} {:>6} {:>9} {:>9}  {}".format(
        "skey", "faces", "poleE", "poleF", "bnd", "manifold", "F in->out", "note"))

    rows = []
    for skey in sorted(base.strips()):
        work = base.copy()
        work.collect_strips()
        if skey not in list(work.strips()):
            continue
        pole_edges = work.has_strip_poles(skey)
        pole_faces = touches_pole(work, skey)
        bnd = len(total_boundary_deletions(work, [skey]))
        coll = len(collateral_strip_deletions(work, [skey]))
        f_in = work.number_of_faces()
        note = ""
        manifold = None
        f_out = None
        try:
            delete_strip(work, skey)
            f_out = work.number_of_faces()
            manifold = work.is_manifold()
        except Exception as e:
            note = "RAISED {}: {}".format(type(e).__name__, str(e)[:40])
        rows.append((skey, f_in, f_out, pole_edges, pole_faces, bnd, coll, manifold, note))
        print("  {:>4} {:>6} {:>6} {:>7} {:>6} {:>9} {:>9}  {}".format(
            skey, len(base.strip_faces(skey)),
            "yes" if pole_edges else "-", "yes" if pole_faces else "-",
            bnd, str(manifold), "{}->{}".format(f_in, f_out), note))

    good = [r for r in rows if r[7] is True]
    bad = [r for r in rows if r[7] is False]
    err = [r for r in rows if r[8]]
    print("  --> manifold {} / non-manifold {} / raised {}".format(
        len(good), len(bad), len(err)))
    if bad:
        pe = sum(1 for r in bad if r[3])
        pf = sum(1 for r in bad if r[4])
        print("      of the non-manifold: {} had pole EDGES, {} had pole FACES".format(pe, pf))
    if good:
        pf = sum(1 for r in good if r[4])
        print("      of the manifold    : {} had pole FACES".format(pf))
    return rows


sweep("square + square hole (no poles)", SQUARE, [box(5.0, 5.0, 2.0)])
sweep("square + circular hole (poles)", SQUARE, [circle(5.0, 5.0, 1.5)])
sweep("L-plate", L, None)
sweep("square + two circular holes", SQUARE,
      [circle(3.0, 3.0, 1.0), circle(7.0, 7.0, 1.0)])
