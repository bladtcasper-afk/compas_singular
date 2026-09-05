"""EXAMPLE 5 -- the same shapes through BOTH front ends.

The field front end exists to replace ``SkeletonDecomposition``'s medial-axis
decomposition. Every measurement so far has compared it against ITSELF -- the
locked baseline says whether a change made things better or worse, never whether
the thing being changed was worth having. This runs both front ends over the
same 18 outlines and measures them with the same function.

    field     FieldDecomposition.from_boundary(...).quad_mesh(...)
    skeleton  boundary_triangulation -> SkeletonDecomposition -> densification

THREE FAIRNESS CORRECTIONS, all of which favour the skeleton
------------------------------------------------------------
Run naively, the skeleton looks far worse than it is. None of these is a real
property of the method, so all three are corrected before measuring:

1. **This correction is now obsolete, and kept only as a record.** It used to
   read: plain ``SkeletonDecomposition`` raises on 15 of the 18 outlines --
   ``pole missing`` and then ``KeyError`` out of ``collect_strips`` -- because
   ``solve_triangular_faces`` only resolves triangles that touch the boundary
   and a fully-interior triangular patch survived with no pole.
   ``RobustSkeletonDecomposition``, below, gave every such triangle a pole.

   Since 2026-09-01 the BASE class does that, as part of the curve-feature
   restoration (``HOW_IT_WORKS.md`` section 5): ``store_pole_data`` records a
   pole for every triangular face, choosing the corner whose two edges are
   closest to equal. Re-measured over all 21 domains of the suite:

       raw outlines        plain raises on 5, robust shim raises on the same 5
       sampled at SPACING  plain raises on 0, robust shim raises on 0

   So the shim now changes nothing -- it is retained because this file is meant
   to stand alone, and because the numbers below were produced with it. The five
   raw-outline failures (L-shape, T-plate, U-plate, comb-plate, square+sq hole)
   are a DIFFERENT defect: the shim never fixed them either. They do not affect
   this comparison, which samples the boundary first (correction 2).

2. **The suite's outlines are coarse polylines** -- a square is four points.
   The medial axis is computed from a Delaunay triangulation of the boundary
   VERTICES, so a four-point square gives it almost nothing to work with, and it
   returns a 2-patch decomposition with a 3 degree corner. That is a strawman.
   Its own examples sample the boundary first (see ``sample_rectangle`` in
   ``06_interior_triangle_pole_fix.py``), so the boundary is densified here to
   the same spacing the field front end uses for its background. On the square
   that alone moves the skeleton from 3.01 degrees to a perfect 90.

3. **The skeleton was densifying every coarse edge as a straight chord.**
   ``quad_mesh`` hands the field's own separatrices and boundary arcs to
   ``densification(edges_to_curves=...)``, so the field front end has always had
   its curves back. Calling the skeleton's ``densification()`` bare gave it
   chords for the same walls, which on a round boundary cuts every corner. Its
   ``decomposition_polylines()`` carries exactly the curves needed; ``_curve_map``
   below matches them to coarse edges by geometric key. Worth an enormous amount
   on the annulus -- minimum angle 34.76 -> 70.39 degrees -- and it is what makes
   the round-holed disc a clear win for the skeleton rather than a near thing.

THE ONE ROW THE SKELETON WINS
-----------------------------
With all three corrections the tally is field 15, skeleton 1, not comparable 2
(the two cables). The single skeleton win is ``square+circle`` -- a round hole in
a square plate -- and it is worth more attention than a 1 deserves, because it is
the only place in the suite the medial axis produces a BETTER HANDLE than the
field:

    square+circle   field 20 patches, min 64.14   skeleton 12 patches, min 63.72
    circle+square   field 20 patches, min 72.44   skeleton 22 patches, min 49.02

The angles are a tie on the first row; the patch count is not. The medial axis of
a round hole in a square is close to what a person would draw, and the skeleton
gets a compact layout straight out of it. The field spends 20 patches on both
mixed-curvature annuli, against 4 for ``disc+round hole`` and 8 for
``square+sq hole`` -- so the cost is not the hole and not the curvature, it is
MIXING them: a wall with corners facing a wall without any. That is the case
where the field's separatrices have nothing natural to terminate on, and it is
recorded here as an open question, not a diagnosis.

Outlines are IMPORTED from ``15_baseline.py``, never retyped -- the suite samples
the ellipse with 60 points and the obvious guess is 48, which is a different
polygon that measures differently and convincingly.

Opens a compas_viewer scene: one group per domain, field on the left and
skeleton on the right, coarse patch edges over each dense mesh.

Run:
    python 17_vs_skeleton.py
    python 17_vs_skeleton.py --no-view
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas.tolerance import TOL                                   # noqa: E402
from compas_singular.algorithms import SkeletonDecomposition        # noqa: E402
from compas_singular.algorithms import boundary_triangulation       # noqa: E402

from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402
from compas_singular.framefield.quality import mesh_quality              # noqa: E402


#: The domain suite, imported so there is exactly one definition of each shape.
#: ``15_baseline`` is not a legal module name and is loaded by path; it guards
#: its CLI behind ``__main__``, so importing it runs nothing.
_spec = importlib.util.spec_from_file_location(
    'baseline_suite', os.path.join(HERE, '15_baseline.py'))
suite = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(suite)

#: Background spacing for the field, and the boundary sampling handed to the
#: skeleton. The same number for both: see fairness correction 2.
SPACING = 0.5

#: Quad size, held equal for both front ends so the only difference measured is
#: the layout each produced.
QUAD = 1.0

#: Domains that get a THIRD panel: the skeleton again, with a point feature at
#: this position, which ``SkeletonDecomposition`` turns into a pseudo-quad pole.
#:
#: The disc is the case that motivates it. Without a pole the skeleton finds a
#: 2-patch decomposition whose silhouette is a perfectly good circle and whose
#: interior hides a 0.09 degree sliver -- a plausible mesh that is not a logical
#: one. A pole at the centre gives the radial fan a disc actually wants.
#:
#: **This is a skeleton-only capability.** ``FieldDecomposition`` accepts a
#: ``poles`` argument but cannot act on one here: it only chooses WHICH corner of
#: an existing triangular patch carries the collapsed side, and the field's disc
#: layout is five clean quads with no triangle to act on. Passing a pole to the
#: field returns a byte-identical mesh. Asking the field for this shape needs a
#: prescribed SINGULARITY -- the disc's total index of +4, currently spread over
#: four +1 singularities in a ring, concentrated at one point -- which is a
#: constraint, not a pole, and does not exist yet.
CENTRE_POLES = {
    'disc': [[5.0, 5.0, 0.0]],
}


class RobustSkeletonDecomposition(SkeletonDecomposition):
    """``SkeletonDecomposition`` that also resolves fully-interior triangles.

    Verbatim from ``06_interior_triangle_pole_fix.py``, and **redundant since
    2026-09-01**: the base class now records a pole for every triangular face,
    and picks a better corner than ``face_vertices(fkey)[0]`` does. Measured
    over the suite's 21 domains, this subclass and the base class now succeed
    and fail on exactly the same ones. Kept so the file stands alone and so the
    numbers below stay reproducible. See fairness correction 1 above.
    """

    def store_pole_data(self, poles):
        mesh = self.mesh
        pole_map = tuple(TOL.geometric_key(pole) for pole in poles)
        face_poles = {}
        for fkey in mesh.faces():
            if len(mesh.face_vertices(fkey)) == 3:
                chosen = None
                for vkey in mesh.face_vertices(fkey):
                    if TOL.geometric_key(mesh.vertex_coordinates(vkey)) in pole_map:
                        chosen = vkey
                        break
                if chosen is None:
                    # interior triangle with no pole among its corners: make one
                    chosen = mesh.face_vertices(fkey)[0]
                face_poles[fkey] = chosen
        mesh.attributes['face_pole'] = face_poles


def densify_loop(loop, spacing):
    """``loop`` resampled at roughly ``spacing``, corners preserved.

    Every original vertex is kept and each segment is subdivided, so a corner of
    the domain is never rounded off by the resampling.
    """
    out = []
    for a, b in zip(loop, list(loop[1:]) + [loop[0]]):
        length = ((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2) ** 0.5
        n = max(1, int(round(length / spacing)))
        for i in range(n):
            t = i / float(n)
            out.append([a[0] + (b[0] - a[0]) * t,
                        a[1] + (b[1] - a[1]) * t, 0.0])
    return out


def run_field(outer, holes, guides, kwargs):
    """The field front end. Returns (coarse, dense, metrics, route)."""
    d = FieldDecomposition.from_boundary(
        outer, inner_boundaries=holes, guides=guides,
        target_length=SPACING, **(kwargs or {}))
    dense = d.quad_mesh(target_length=QUAD)
    return d.mesh, dense, mesh_quality(dense), d.route()


def _curve_map(decomposition, coarse):
    """``{(u, v): polyline}`` for :meth:`densification`, from the skeleton's own
    branch polylines. Fairness correction 3.

    ``decomposition_polylines`` returns the medial-axis branches and the boundary
    arcs between splits -- the true curves the coarse edges are chords of. They
    are matched to coarse vertices by geometric key, the same way
    ``03_decomposition_continuous_curved.py`` does it for Rhino.

    **The mapping must be TOTAL.** Given a mapping at all, ``densification`` looks
    up every edge in it; there is no per-edge straight-chord branch, only a
    per-call one, so a partial dict raises ``KeyError`` several frames away. Any
    coarse edge no polyline matched gets an explicit two-point chord, which is
    exactly what the no-mapping branch would have drawn.

    Returns
    -------
    (dict, int, int)
        The mapping, how many edges a real curve was found for, and the total.
    """
    gkey_vertex = {TOL.geometric_key(coarse.vertex_coordinates(vkey)): vkey
                   for vkey in coarse.vertices()}
    mapping = {}
    for polyline in decomposition.decomposition_polylines():
        u = gkey_vertex.get(TOL.geometric_key(polyline[0]))
        v = gkey_vertex.get(TOL.geometric_key(polyline[-1]))
        if u is not None and v is not None and u != v:
            mapping[u, v] = [list(point) for point in polyline]

    matched = sum(1 for u, v in coarse.edges()
                  if (u, v) in mapping or (v, u) in mapping)
    for u, v in coarse.edges():
        if (u, v) not in mapping and (v, u) not in mapping:
            mapping[u, v] = [coarse.vertex_coordinates(u),
                             coarse.vertex_coordinates(v)]
    return mapping, matched, coarse.number_of_edges()


def run_skeleton(outer, holes, poles=()):
    """The medial-axis front end. Returns (coarse, dense, metrics, note).

    Guides are not passed: the skeleton has nowhere to put them. That is not an
    oversight in this script, it is the finding -- see the cable rows.

    ``poles`` ARE passed, as ``point_features``, because the skeleton can take
    them -- see ``CENTRE_POLES``. Note that a pole makes the result sensitive to
    the boundary sampling in a way the unpoled run is not: the Delaunay of a
    perfectly symmetric outline is cocircular-prone, and on the disc a pole
    survives cleanly at ``SPACING`` 0.5, 0.8 and 1.0 but shatters the layout into
    144 patches at 0.4 and 0.3. ``07_square_centre_pole.py`` documents the same
    effect and works around it by retrying at several resolutions.
    """
    trimesh = boundary_triangulation(
        densify_loop(outer, SPACING),
        [densify_loop(h, SPACING) for h in (holes or [])], [], list(poles))
    decomposition = RobustSkeletonDecomposition.from_mesh(trimesh)
    coarse = decomposition.decomposition_mesh(list(poles))
    coarse.collect_strips()
    coarse.set_strips_density_target(QUAD)
    mapping, matched, total = _curve_map(decomposition, coarse)
    coarse.densification(edges_to_curves=mapping)
    dense = coarse.get_quad_mesh()
    return coarse, dense, mesh_quality(dense), '%d/%d curved' % (matched, total)


def verdict(field_q, skel_q, field_route, guides):
    """Which front end produced the better mesh of this domain, and why.

    Minimum angle first -- it is what makes a mesh unusable -- then the coarse
    patch count, because a coarse layout is the thing compas_singular's strip
    editing operates on and more patches is a worse handle, not a better one.

    **A guided domain is scored ``n/a``, never as a skeleton win.** The skeleton
    returns a beautiful 90-degree grid there and would top any angle comparison,
    but it does so by ignoring the cable completely -- its mesh of
    ``square+cable`` is byte-identical to its mesh of the plain square. Scoring
    that as a win would rank "answered a different question" above "answered
    this one imperfectly". The field pays real quality to follow the guide; that
    is the trade being made, not a defeat.
    """
    if skel_q is None:
        return 'field', 'skeleton raised'
    if guides:
        return 'n/a', 'skeleton cannot accept a guide'
    if field_route != 'field':
        return 'skeleton', 'field fell back to ' + field_route
    lo_f, lo_s = field_q['min_angle'], skel_q['min_angle']
    if lo_f > lo_s + 5.0:
        return 'field', 'min angle %+.0f deg' % (lo_f - lo_s)
    if lo_s > lo_f + 5.0:
        return 'skeleton', 'min angle %+.0f deg' % (lo_s - lo_f)
    pf, ps = field_q['patches'], skel_q['patches']
    if pf < ps:
        return 'field', 'same angles, %d fewer patches' % (ps - pf)
    if ps < pf:
        return 'skeleton', 'same angles, %d fewer patches' % (pf - ps)
    return 'tie', 'equivalent'


def main(argv):
    show = '--no-view' not in argv
    wanted = [a for a in argv[1:] if not a.startswith('-')]

    results = []
    for label, outer, holes, guides, kwargs in suite.DOMAINS:
        if wanted and label not in wanted:
            continue
        f_coarse, f_dense, f_q, f_route = run_field(outer, holes, guides, kwargs)
        f_q['patches'] = f_coarse.number_of_faces()
        try:
            s_coarse, s_dense, s_q, _ = run_skeleton(outer, holes)
            s_q['patches'] = s_coarse.number_of_faces()
        except Exception as exc:                      # noqa: BLE001
            s_coarse = s_dense = s_q = None
            print('  %s: skeleton raised %s' % (label, type(exc).__name__))
        results.append((label, guides, f_coarse, f_dense, f_q, f_route,
                        s_coarse, s_dense, s_q))

    # ---- table ------------------------------------------------------------
    print()
    print('%-18s %26s %26s   %s' % ('', '--- FIELD ---', '--- SKELETON ---', ''))
    print('%-18s %5s %5s %6s %6s  %5s %5s %6s %6s   %s'
          % ('domain', 'patch', 'faces', 'min', 'max',
             'patch', 'faces', 'min', 'max', 'better'))
    print('-' * 104)
    wins = {'field': 0, 'skeleton': 0, 'tie': 0, 'n/a': 0}
    for (label, guides, _fc, _fd, f_q, f_route, _sc, _sd, s_q) in results:
        who, why = verdict(f_q, s_q, f_route, guides)
        wins[who] += 1
        if s_q is None:
            right = '%5s %5s %6s %6s' % ('--', '--', 'raised', '--')
        else:
            right = '%5d %5d %6.2f %6.2f' % (s_q['patches'], s_q['faces'],
                                             s_q['min_angle'], s_q['max_angle'])
        print('%-18s %5d %5d %6.2f %6.2f  %s   %-9s %s'
              % (label, f_q['patches'], f_q['faces'], f_q['min_angle'],
                 f_q['max_angle'], right, who, why))

    print()
    print('better mesh: field %d, skeleton %d, tie %d, not comparable %d'
          % (wins['field'], wins['skeleton'], wins['tie'], wins['n/a']))

    # ---- the cable rows, which are the point -------------------------------
    cabled = [r for r in results if r[1]]
    if cabled:
        print()
        print('THE CABLES -- the skeleton is not merely worse on these, it is blind')
        plain = [r for r in results if r[0] == 'square'][0]
        for (label, _g, _fc, _fd, f_q, _fr, _sc, _sd, s_q) in cabled:
            if s_q is None:
                continue
            same = (plain[8] is not None
                    and abs(s_q['min_angle'] - plain[8]['min_angle']) < 1e-9
                    and s_q['faces'] == plain[8]['faces'])
            print('    %-18s skeleton %d faces, min %.2f deg  %s'
                  % (label, s_q['faces'], s_q['min_angle'],
                     '== its mesh of the PLAIN square' if same else ''))
        print('    The guide curve is not an input the medial axis can accept, so the')
        print('    cable changes nothing. The field front end bends the mesh onto it.')

    if show:
        view(results)
    return 0


def view(results):
    from compas_viewer import Viewer                                # noqa: E402

    from compas_singular.framefield.viz import Grid        # noqa: E402
    from compas_singular.framefield.viz import add_dense   # noqa: E402
    from compas_singular.framefield.viz import add_guides  # noqa: E402
    from compas_singular.framefield.viz import PATCH_EDGE  # noqa: E402
    from compas_singular.framefield.viz import translated  # noqa: E402

    viewer = Viewer()
    grid = Grid(pitch=80, cols=4)
    for i, (label, guides, f_coarse, f_dense, _fq, f_route,
            s_coarse, s_dense, _sq) in enumerate(results):
        dx, dy = grid.cell(i)
        group = viewer.scene.add_group(name=label)
        # left: the field front end
        add_dense(group, f_dense, dx, dy, name='field dense (%s)' % f_route)
        group.add(translated(f_coarse, dx, dy), show_faces=False, show_points=True,
                  linecolor=PATCH_EDGE, linewidth=3, pointcolor=PATCH_EDGE,
                  name='field patches')
        # right: the medial-axis front end, one panel over
        if s_dense is not None:
            add_dense(group, s_dense, dx + 40.0, dy, name='skeleton dense')
            group.add(translated(s_coarse, dx + 40.0, dy), show_faces=False,
                      show_points=True, linecolor=PATCH_EDGE, linewidth=3,
                      pointcolor=PATCH_EDGE, name='skeleton patches')
        for guide in (guides or []):
            add_guides(group, [guide], dx, dy)
            add_guides(group, [guide], dx + 40.0, dy)
    grid.frame(viewer)
    viewer.show()


if __name__ == '__main__':
    sys.exit(main(sys.argv))
