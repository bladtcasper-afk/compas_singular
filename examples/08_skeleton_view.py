"""Visualise the topological skeleton (medial axis) of any baseline domain.

Four layers per domain:

* the Delaunay triangulation the skeleton is derived from (grey, faint),
* the skeleton itself -- lines between the circumcentres of adjacent triangles,
  bundled into branch polylines (red),
* the singular points -- circumcentres of the triangles with three neighbours,
  which become the singularities of the coarse quad mesh (blue),
* the full decomposition polylines: skeleton + branches out to the boundary +
  the split boundary arcs (green). This is what ``decomposition_mesh`` actually
  builds the patches from.

Two things about the order matter, and both are easy to get wrong.

**Draw the skeleton BEFORE calling ``decomposition_mesh()``.** That call runs
``decomposition_polylines()``, whose fix-up passes (``branches_splitting_*``)
insert vertices into the triangulation, so the skeleton it reports afterwards is
not the one that produced the singularities you are looking at.

**Densify the boundary first.** The medial axis is built from a Delaunay
triangulation of the boundary VERTICES, and the suite's outlines are coarse --
``SQUARE`` is four points. Handed in raw, a square gives the skeleton almost
nothing to work with and it returns a 2-patch decomposition with a 3 degree
corner. ``densify_loop`` is reused from ``17_vs_skeleton.py``, at the same
spacing the field front end uses for its background: fairness correction 2 there,
and the difference between a picture of the method and a picture of a strawman.

The domains come from ``New approach/15_baseline.py`` (``suite.DOMAINS``), so
adding a shape there adds it here with no edit.

Run:
    python 08_skeleton_view.py                     # all baseline domains
    python 08_skeleton_view.py disc L-shape        # just these
    python 08_skeleton_view.py --list              # names, no geometry built
    python 08_skeleton_view.py disc --field        # + the field's own analogue
    python 08_skeleton_view.py --no-view           # table only
"""
import importlib.util
import os
import sys

from compas.colors import Color
from compas.geometry import Point
from compas.geometry import Polyline

from compas_singular.algorithms import SkeletonDecomposition
from compas_singular.algorithms import boundary_triangulation

HERE = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.join(HERE, 'New approach')

#: ``17_vs_skeleton.py`` already owns the domain suite loader and the boundary
#: resampler, and neither is worth a second copy that can drift from it. Neither
#: file is a legal module name, so it is loaded by path; both guard their CLI
#: behind ``__main__``, so importing runs nothing. (``framefield`` no longer
#: needs that path -- it is ``compas_singular.framefield`` now.)
_spec = importlib.util.spec_from_file_location(
    'vs_skeleton', os.path.join(NEW, '17_vs_skeleton.py'))
vs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vs)

suite = vs.suite
densify_loop = vs.densify_loop

#: Boundary sampling. ``17_vs_skeleton.SPACING``, i.e. the same number the field
#: front end uses for its background triangulation.
SPACING = vs.SPACING

TRI = Color(0.6, 0.6, 0.6)
BRANCH = Color.red()
SINGULAR = Color.blue()
POLYLINE = Color.green()


def skeleton_of(outer, holes, spacing=SPACING, poles=()):
    """Every skeleton layer of one domain, in the order they must be computed.

    Plain ``SkeletonDecomposition`` is enough here -- the interior-triangle pole
    bug that ``RobustSkeletonDecomposition`` fixes lives in ``store_pole_data``,
    which only ``decomposition_mesh`` reaches, and nothing below calls it.

    Returns
    -------
    dict
        ``trimesh``, ``lines``, ``branches``, ``singular_points``,
        ``polylines``, and ``error`` (a string if the decomposition polylines
        raised -- the skeleton layers are still valid when they do).
    """
    trimesh = boundary_triangulation(
        densify_loop(outer, spacing),
        [densify_loop(hole, spacing) for hole in (holes or [])], [], list(poles))
    decomposition = SkeletonDecomposition.from_mesh(trimesh)

    layers = {
        'trimesh': trimesh,
        'lines': decomposition.lines(),
        'branches': decomposition.branches(),
        'singular_points': decomposition.singular_points(),
        'polylines': [],
        'error': None,
    }
    # last, and only after the layers above have been read off: this mutates
    # the triangulation
    try:
        layers['polylines'] = decomposition.decomposition_polylines()
    except Exception as exc:                                      # noqa: BLE001
        layers['error'] = '%s: %s' % (type(exc).__name__, exc)
    return layers


def bbox(points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return min(xs), min(ys), max(xs), max(ys)


def main(argv):
    show = '--no-view' not in argv
    with_field = '--field' in argv
    wanted = [a for a in argv[1:] if not a.startswith('-')]

    if '--list' in argv:
        for label, _o, holes, guides, _k in suite.DOMAINS:
            print('%-18s %s%s' % (label,
                                  '%d hole(s) ' % len(holes) if holes else '',
                                  'guided' if guides else ''))
        return 0

    results = []
    for label, outer, holes, guides, _kwargs in suite.DOMAINS:
        if wanted and label not in wanted:
            continue
        layers = skeleton_of(outer, holes)
        results.append((label, outer, holes, guides, layers))

    if not results:
        print('no domain matched %s -- try --list' % wanted)
        return 1

    print()
    print('%-18s %6s %6s %8s %6s %8s   %s'
          % ('domain', 'tri', 'lines', 'branches', 'sing.', 'polylns', 'note'))
    print('-' * 78)
    for label, _o, _h, guides, k in results:
        print('%-18s %6d %6d %8d %6d %8d   %s'
              % (label, k['trimesh'].number_of_faces(), len(k['lines']),
                 len(k['branches']), len(k['singular_points']),
                 len(k['polylines']),
                 k['error'] or ('guide ignored by the skeleton' if guides else '')))

    if show:
        view(results, with_field)
    return 0


def view(results, with_field=False):
    from compas_viewer import Viewer

    from compas_singular.framefield.viz import Grid

    # widest domain in the selection sets the cell size, so RV (a metre-scale
    # curve OBJ) and the square (10 units) can share one screen
    extents = []
    for _l, outer, _h, _g, _k in results:
        x0, y0, x1, y1 = bbox(outer)
        extents.append(max(x1 - x0, y1 - y0))
    pitch = max(extents) * (2.6 if with_field else 1.35)
    cols = max(1, int(len(results) ** 0.5 + 0.5))

    # point size and line width are in SCREEN pixels, not model units, so the
    # 20px dot that reads as a singularity on one domain covers a whole panel of
    # twenty-one. They scale with how far the camera had to pull back.
    dot = 20 if len(results) == 1 else 6
    wide, thin = (4, 2) if len(results) == 1 else (2, 1)

    viewer = Viewer()
    # Grid also owns the camera: the default one frames a single object, so a
    # multi-panel scene opens zoomed into one corner of itself.
    grid = Grid(pitch=pitch, cols=cols)
    for i, (label, outer, _holes, guides, k) in enumerate(results):
        x0, y0, _x1, _y1 = bbox(outer)
        cell_x, cell_y = grid.cell(i)
        # each domain carries its own origin (RV is nowhere near 0,0), so the
        # bbox corner is pulled onto the cell rather than the raw coordinates
        dx = cell_x - x0
        dy = cell_y - y0
        group = viewer.scene.add_group(name=label)

        group.add(translated(k['trimesh'], dx, dy), show_points=False,
                  show_lines=True, show_faces=True, opacity=0.25,
                  linecolor=TRI, name='Delaunay (%d faces)' % k['trimesh'].number_of_faces())
        for polyline in k['polylines']:
            group.add(shifted(polyline, dx, dy, 0.02), linewidth=thin,
                      linecolor=POLYLINE, name='decomposition polyline')
        for branch in k['branches']:
            group.add(shifted(branch, dx, dy, 0.04), linewidth=wide,
                      linecolor=BRANCH, name='skeleton branch')
        for point in k['singular_points']:
            group.add(Point(point[0] + dx, point[1] + dy, 0.06), pointsize=dot,
                      pointcolor=SINGULAR, name='singular point')

        if with_field:
            add_field_panel(viewer, label, outer, _holes, guides,
                            dx + pitch * 0.5, dy)

    grid.frame(viewer)
    viewer.show()


def add_field_panel(viewer, label, outer, holes, guides, dx, dy):
    """The field front end's analogue of the skeleton, one panel to the right.

    There is no medial axis anywhere in ``framefield`` -- the field front end
    replaces it. What plays the same role is the SEPARATRIX GRAPH: field lines
    traced out of the singularities of a smoothed cross field, which cut the
    domain into patches exactly as the skeleton branches do. Same job, different
    object: branches are equidistant from the walls, separatrices follow a
    direction field, so they are the ones that can be steered by a guide curve.
    """
    from compas_singular.framefield.decomposition import FieldDecomposition
    from compas_singular.framefield.viz import add_separatrices
    from compas_singular.framefield.viz import add_singularities
    from compas_singular.framefield.viz import add_layout

    group = viewer.scene.add_group(name=label + ' (field)')
    decomposition = FieldDecomposition.from_boundary(
        outer, inner_boundaries=holes, guides=guides, target_length=SPACING)
    decomposition.decomposition_mesh()
    add_layout(group, decomposition, dx, dy, faces=False)
    add_separatrices(group, decomposition, dx, dy)
    add_singularities(group, decomposition, dx, dy)


def translated(mesh, dx, dy):
    """Copy of ``mesh`` moved by ``(dx, dy)`` -- ``framefield.viz.translated``,
    inlined so this file works with or without the field package present."""
    cloned = mesh.copy()
    for vkey in cloned.vertices():
        x, y, z = cloned.vertex_coordinates(vkey)
        cloned.vertex_attributes(vkey, 'xyz', [x + dx, y + dy, z])
    return cloned


def shifted(points, dx, dy, dz=0.0):
    return Polyline([[p[0] + dx, p[1] + dy, p[2] + dz] for p in points])


if __name__ == '__main__':
    sys.exit(main(sys.argv))
