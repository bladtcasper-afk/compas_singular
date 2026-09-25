"""EXAMPLE 2 -- the field layout is an ordinary compas_singular coarse mesh.

This is the point of the whole exercise. ``FieldDecomposition`` returns a
``CoarsePseudoQuadMesh``, so everything compas_singular already does to a coarse
mesh keeps working, unchanged and un-ported:

* per-strip and global density control
* ``densification`` with ``edges_to_curves``, so curved separatrices stay curved
* the strip grammar (``add_strip`` / ``delete_strip``)
* strip and polyedge queries
* saving to JSON and reloading

Nothing in ``src/compas_singular`` was modified to make any of this work.

Opens a compas_viewer scene showing one layout driving several densities.

Run:
    python 11_downstream.py
    python 11_downstream.py --no-view      # numbers only
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas_singular.datastructures.mesh_quad.grammar.add_strip import add_strip  # noqa: E402

from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402


PENTAGON = [[0, 0, 0], [10, 0, 0], [13, 7, 0], [5, 12, 0], [-2, 7, 0]]
SQUARE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]
SQUARE_HOLE = [[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]]


def fresh(boundary=PENTAGON, holes=None, target_length=0.5):
    d = FieldDecomposition.from_boundary(boundary, holes, target_length=target_length)
    coarse = d.decomposition_mesh()
    coarse.collect_strips()
    return d, coarse


def demo_density():
    """Density control, THROUGH the field-aware densifier.

    This is the check that field integration did not cost the strip structure.
    ``densify`` moves patch INTERIORS only and leaves every boundary polyline
    exactly where ``discrete_coons_patch`` would have put it, so opposite edges
    of a patch keep matching densities and the patches still weld -- but that is
    an argument, and this section is the measurement. Each call goes through
    ``d.densify`` rather than ``coarse.densification`` so a regression in the
    new code shows up here rather than only in the baseline's angle columns.
    """
    print('--- 1. density control -------------------------------------------')
    d, coarse = fresh()
    print('    {} patches, {} strips'.format(coarse.number_of_faces(), coarse.number_of_strips()))

    def densified(label):
        mesh = d.densify(coarse)
        quads = all(len(mesh.face_vertices(f)) == 4 for f in mesh.faces())
        ok = mesh.is_manifold() and quads
        print('    {:24s} -> {:4d} faces  manifold {}  all-quad {}'.format(
            label, mesh.number_of_faces(), mesh.is_manifold(), quads))
        return mesh, ok

    good = True
    coarse.set_strips_density(4)
    _, ok = densified('uniform density 4')
    good = good and ok

    coarse.set_strips_density_target(0.6)
    _, ok = densified('target edge length 0.6')
    good = good and ok

    coarse.set_mesh_density_face_target(400)
    _, ok = densified('target ~400 faces')
    good = good and ok

    # one strip graded independently of the rest -- this is the reason to keep a
    # coarse layout at all rather than emit a dense mesh directly
    coarse.set_strips_density(3)
    first = sorted(coarse.strips())[0]
    coarse.set_strip_density(first, 12)
    _, ok = densified('strip {} at 12, rest at 3'.format(first))
    good = good and ok

    # the strips must still be collectable FROM the coarse layout afterwards --
    # densification must not have left it in a state the grammar cannot read
    before = coarse.number_of_strips()
    coarse.collect_strips()
    print('    strips re-collected: {} -> {}  {}'.format(
        before, coarse.number_of_strips(),
        'unchanged' if before == coarse.number_of_strips() else 'CHANGED'))
    good = good and before == coarse.number_of_strips()
    print('    {}'.format('OK' if good else 'FAILED'))
    return good


def demo_edges_to_curves():
    print('--- 2. edges_to_curves keeps separatrices curved ------------------')
    d, coarse = fresh()
    coarse.set_strips_density(6)

    coarse.densification()
    straight = coarse.get_quad_mesh()
    coarse.densification(edges_to_curves=d.edges_to_curves())
    curved = coarse.get_quad_mesh()

    # how far the two dense meshes disagree: with straight chording, every
    # interior separatrix is replaced by the line between its endpoints
    def spread(mesh):
        pts = [mesh.vertex_coordinates(v) for v in mesh.vertices()]
        return sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)

    print('    without: {} faces, centroid {:.4f} {:.4f}'.format(
        straight.number_of_faces(), *spread(straight)))
    print('    with   : {} faces, centroid {:.4f} {:.4f}'.format(
        curved.number_of_faces(), *spread(curved)))
    print('    mapped curves: {} of {} coarse edges'.format(
        len(d.edges_to_curves()), coarse.number_of_edges()))


def demo_strip_grammar():
    print('--- 3. the strip grammar still applies ----------------------------')
    d, coarse = fresh(SQUARE, [SQUARE_HOLE], target_length=0.4)
    coarse.collect_polyedges()
    before = coarse.number_of_strips(), coarse.number_of_faces()

    slot = None
    for _pkey, polyedge in coarse.collect_polyedges():
        polyedge = list(polyedge)
        if len(polyedge) < 3 or polyedge[0] == polyedge[-1]:
            continue
        if not (coarse.is_vertex_on_boundary(polyedge[0])
                and coarse.is_vertex_on_boundary(polyedge[-1])):
            continue
        if all(coarse.is_vertex_on_boundary(v) for v in polyedge):
            continue
        slot = polyedge
        break

    if slot is None:
        print('    no interior wall-to-wall polyedge in this layout, so no slot '
              'for add_strip. Not a failure of the field front end -- '
              'guide_lines.py documents the same scarcity on skeleton layouts, '
              'and the fix there is refinement.')
        return

    add_strip(coarse, slot)
    coarse.collect_strips()
    print('    add_strip on polyedge {}'.format(slot))
    print('    strips {} -> {}, faces {} -> {}'.format(
        before[0], coarse.number_of_strips(), before[1], coarse.number_of_faces()))
    print('    still all quads: {}'.format(
        all(len(coarse.face_vertices(f)) == 4 for f in coarse.faces())))


def demo_roundtrip():
    print('--- 4. it is just a mesh ------------------------------------------')
    import json
    import tempfile

    d, coarse = fresh()
    path = os.path.join(tempfile.gettempdir(), 'framefield_coarse.json')
    coarse.to_json(path)

    from compas_singular.datastructures import CoarsePseudoQuadMesh
    reloaded = CoarsePseudoQuadMesh.from_json(path)
    reloaded.collect_strips()
    reloaded.set_strips_density(5)
    reloaded.densification()
    print('    saved and reloaded: {} patches, {} dense faces'.format(
        reloaded.number_of_faces(), reloaded.get_quad_mesh().number_of_faces()))
    dense = reloaded.get_quad_mesh()
    interior = [v for v in dense.singularities() if not dense.is_vertex_on_boundary(v)]
    print('    singularities in the DENSE mesh: {} ({} interior, {} on the wall)'.format(
        len(dense.singularities()), len(interior),
        len(dense.singularities()) - len(interior)))
    print('    the wall ones are the pentagon\'s own corners -- a valence-2 boundary '
          'vertex is irregular by compas_singular\'s definition. Only the interior '
          'count should match the field\'s singularities.')
    del json


def show():
    """One layout, several densities -- the reason to keep a coarse mesh.

    Every panel is the SAME 5 patches; only the strip densities differ. The
    coarse layout is drawn over each so the patch boundaries stay visible.
    """
    from compas_viewer import Viewer

    from compas_singular.framefield.viz import Grid, add_dense, add_layout, add_singularities

    viewer = Viewer()
    grid = Grid(pitch=17.0, cols=3)

    d, coarse = fresh()
    dx, dy = grid.cell(0)
    g = viewer.scene.add_group(name='1. the layout -- {} patches, {} strips'.format(
        coarse.number_of_faces(), coarse.number_of_strips()))
    add_layout(g, d, dx, dy)
    add_singularities(g, d, dx, dy)

    for i, target in enumerate((2.0, 1.0, 0.5), start=1):
        coarse.set_strips_density_target(target)
        coarse.densification(edges_to_curves=d.edges_to_curves())
        dx, dy = grid.cell(i)
        g = viewer.scene.add_group(name='{}. quad size {:.1f} -- {} quads'.format(
            i + 1, target, coarse.get_quad_mesh().number_of_faces()))
        add_dense(g, coarse.get_quad_mesh(), dx, dy)
        add_layout(g, d, dx, dy, faces=False)

    # one strip graded against the rest
    coarse.set_strips_density(3)
    first = sorted(coarse.strips())[0]
    coarse.set_strip_density(first, 12)
    coarse.densification(edges_to_curves=d.edges_to_curves())
    dx, dy = grid.cell(4)
    g = viewer.scene.add_group(
        name='5. strip {} at 12, the rest at 3'.format(first))
    add_dense(g, coarse.get_quad_mesh(), dx, dy)
    add_layout(g, d, dx, dy, faces=False)

    # a hole changes nothing downstream
    dh, ch = fresh(SQUARE, [SQUARE_HOLE], target_length=0.4)
    ch.set_strips_density_target(0.5)
    ch.densification(edges_to_curves=dh.edges_to_curves())
    dx, dy = grid.cell(5)
    g = viewer.scene.add_group(name='6. a hole is no different -- {} quads'.format(
        ch.get_quad_mesh().number_of_faces()))
    add_dense(g, ch.get_quad_mesh(), dx, dy)
    add_layout(g, dh, dx, dy, faces=False)

    grid.frame(viewer)
    viewer.show()


def main():
    demo_density()
    print()
    demo_edges_to_curves()
    print()
    demo_strip_grammar()
    print()
    demo_roundtrip()

    if '--no-view' not in sys.argv:
        show()
    return 0


if __name__ == '__main__':
    sys.exit(main())
