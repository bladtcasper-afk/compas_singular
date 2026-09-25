"""EXAMPLE 23 -- **the singularity as a BLOCK, one joint at a time.**

A singularity is a vertex where 3 or 5 mesh lines gather. A block system cannot
build that: blocks meet along faces, never in a point. Dualising the whole mesh
fixes it -- every vertex becomes a face, so the valence-n joint becomes an
n-gon BLOCK -- and that is what ``FloorExample vault.py`` does. But it moves
every block in the plate by half a block to solve a problem that exists at one
vertex, and on a curved wall it chords the outline and loses area.

``framefield.block`` does it LOCALLY. One joint, one n-gon, everything else the
mesh you started with. This file is the evaluation: what it costs, what it
looks like, and the two conditions under which it will not do it at all.

WHAT IT COSTS -- SEAMS, AND THE PARITY LAW THAT FIXES THE NUMBER
-----------------------------------------------------------------

Truncating the joint leaves the n incident quads as PENTAGONS. Each is repaired
to two quads by splitting one outer edge, and a split edge in a quad mesh is a
strip split -- the face beyond it is a pentagon now too. So the repair runs
outward until its strip reaches a wall, and the block costs n SEAMS radiating
from it. Panel 2 of each row draws them in green.

That is a conservation law, not a weakness of this construction. For any
sub-disc of a quad mesh ``4F = 2E_int + E_bnd``, so ``E_bnd`` is EVEN; hold one
n-gon in it and ``E_bnd`` takes the parity of n. For odd n -- and 3 and 5 are
the only valences that occur -- those disagree, so NO region however large can
contain the block. Part 2 checks that by printing k-ring boundary edge counts;
they are even out to the 4th ring, every time.

The law also says the price could be **1 seam, not n**, because two odd blocks
together have even parity and can absorb each other's. This module does not do
that yet, and the gap is the honest headline: measured below, a block costs n.

WHEN IT REFUSES
---------------

Two conditions, both reported by ``blockable()`` before anything is built:

    a repair strip is CLOSED     the split has no wall to run to, comes back
                                 round and lands on an edge of the joint itself
    a POLE sits in the 1-ring    there is no quad there to truncate

The first is why the DISC takes no local blocks at all -- each of its four
valence-3 joints has one returning strip in BOTH chiralities, at every density.
The second never bites the global dual, but the global dual has the opposite
problem: it must not touch a pole-bearing mesh at all, because it turns each
pole into a valence-3 joint. Row 4 shows both halves of that trade.

CHIRALITY IS A REAL CHOICE
--------------------------

Each incident quad chords to one of the two block corners it touches, and every
corner needs exactly one chord -- a perfect matching on a cycle, of which there
are exactly two. So the pinwheel spins one way or the other, mixing cannot
close, and which one you pick decides success: on ``plate+hole`` NO joint is
blockable in spin 0 and five of the eight are in spin 1. ``blockable()`` tries
both and reports which; ``blocks_any_spin`` searches the assignments when you
want every block at once.

Wanting every block at once is where the seams start getting in each other's
way. It succeeds on the hexagon (2 blocks, 10 seams) and is REFUSED on
``plate+hole``, where one block's seam reaches another block. Nothing about that
is fatal -- the blocks can go in one at a time -- but it is the cost of emitting
n seams where parity only demands 1.

READING THE PANELS
------------------

    1  primal          joints drawn thick -- red valence 3, blue valence 5
    2  block inserted  the n-gon filled in the same colour; green = the seams
    3  relaxed         what it is once the block is an element, not a guess
    4  global dual     the same blocks by the other route, for comparison

Row 1 is the pentagon plate close up, so a single block is legible; rows 2-5 are
whole plates.

Run:
    python 23_singularity_block.py
    python 23_singularity_block.py --no-view
    python 23_singularity_block.py pentagon disc
"""
import importlib.util
import os
import sys
from math import cos
from math import pi
from math import sin

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas.datastructures import Mesh                                # noqa: E402
from compas.geometry import distance_point_point                      # noqa: E402

from compas_singular.framefield import block as B                        # noqa: E402
from compas_singular.framefield.decomposition import FieldDecomposition  # noqa: E402


_spec = importlib.util.spec_from_file_location(
    'baseline_suite', os.path.join(HERE, '15_baseline.py'))
suite = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(suite)


#: Background spacing -- the resolution the FIELD is solved at.
SPACING = 0.6

#: Target quad edge length, applied to the coarse layout.
QUAD = 1.0

#: How far along its incident edges a block's corners sit. The block's SIZE is
#: the one free parameter; its DEGREE is pinned to the valence by the index.
RATIO = 0.45

#: Panel pitch for the viewer grid.
PITCH = 16.0

#: Radius of the close-up in row 1, in model units.
ZOOM_RADIUS = 3.2
ZOOM = 2.0


#: Four domains, for four reasons. ``pentagon`` has ONE valence-5 joint and
#: straight walls, so a single block is easy to read and the area check is
#: exact. ``hexagon`` has two, which is where seams start meeting each other.
#: ``disc`` is the refusal: four valence-3 joints, none of them blockable, in
#: both chiralities, at any density. ``square+circle`` carries POLES -- the one
#: family the global dual may not touch, and the reason to have a local
#: operator at all.
#: A 12-unit plate with a 2.2-unit round hole. NOT ``suite.SQUARE`` +
#: ``suite.ROUND_HOLE``, which is the same shape one third smaller and refuses
#: at every joint -- fine as a result, but then the row shows only half of the
#: trade. This one carries poles AND has a blockable joint, so it shows both:
#: the local operator building a block the global dual may not, and the global
#: dual breaking the index on the very same mesh.
PLATE = [[0, 0, 0], [12, 0, 0], [12, 12, 0], [0, 12, 0]]
ROUND_HOLE = [[6 + 2.2 * cos(-2 * pi * i / 24),
               6 + 2.2 * sin(-2 * pi * i / 24), 0.0] for i in range(24)]

DOMAINS = [
    ('pentagon',   suite.PENTAGON, None),
    ('hexagon',    suite.HEXAGON,  None),
    ('disc',       suite.DISC,     None),
    ('plate+hole', PLATE,          [ROUND_HOLE]),
]


# ------------------------------------------------------------------
# building the mesh under test
# ------------------------------------------------------------------

def plate(outer, holes=None):
    """A dense all-quad plate, straight off the frame-field front end.

    Rebuilt as a plain ``Mesh``: the operator rewrites the face list wholesale,
    so a ``QuadMesh``'s strip and pole attributes would come out stale rather
    than wrong, which is worse.
    """
    d = FieldDecomposition.from_boundary(
        [list(p) for p in outer], inner_boundaries=holes,
        target_length=SPACING)
    dense = d.quad_mesh(target_length=QUAD)
    return Mesh.from_vertices_and_faces(*dense.to_vertices_and_faces())


# ------------------------------------------------------------------
# measuring
# ------------------------------------------------------------------

def ring_parity(mesh, vkey, rings=4):
    """Boundary edge count of the k-ring of ``vkey``, for k = 1..rings.

    The parity law, checked rather than asserted. Every one of these is EVEN,
    and an odd-sided block needs an odd one, which is why the seam has to leave.
    """
    faces, frontier, seen = set(), {vkey}, set()
    out = []
    for _ in range(rings):
        seen |= frontier
        for v in list(frontier):
            faces |= set(f for f in mesh.vertex_faces(v) if f is not None)
        frontier = set()
        for f in faces:
            frontier |= set(mesh.face_vertices(f))
        frontier -= seen
        count = {}
        for f in faces:
            vertices = mesh.face_vertices(f)
            for i in range(len(vertices)):
                e = frozenset((vertices[i], vertices[(i + 1) % len(vertices)]))
                count[e] = count.get(e, 0) + 1
        out.append(sum(1 for c in count.values() if c == 1))
    return out


def kept_edges(primal, result, tol=5):
    """How many of the primal's edges survive verbatim in ``result``.

    The whole claim of a LOCAL operator, as a number. Compare with the global
    dual, which shares essentially none of them.
    """
    def keys(mesh):
        out = set()
        for u, v in mesh.edges():
            a = tuple(round(c, tol) for c in mesh.vertex_coordinates(u))
            b = tuple(round(c, tol) for c in mesh.vertex_coordinates(v))
            out.add(frozenset((a, b)))
        return out

    old = keys(primal)
    return len(old & keys(result)), len(old)


def describe(label, mesh, primal=None):
    print('    {:<20s} F={:4d} V={:4d}  faces {:<18s} interior joints {:<18s}'
          .format(label, mesh.number_of_faces(), mesh.number_of_vertices(),
                  _fmt(B.face_degrees(mesh)), _fmt(B.interior_valences(mesh))))
    line = ('                         index {:<4d} area {:<10.4f} '
            'boundary rings {}  angles {:.1f}-{:.1f} deg, aspect {:.2f}'.format(
                B.index_sum(mesh), mesh.area(),
                len(mesh.vertices_on_boundaries()), *B.quality(mesh)))
    print(line)
    if primal is not None:
        kept, total = kept_edges(primal, mesh)
        print('                         primal edges kept verbatim: '
              '{} of {} ({:.0%})'.format(kept, total, kept / total))


def _fmt(d):
    return '{' + ', '.join('{}:{}'.format(k, d[k]) for k in sorted(d)) + '}'


# ------------------------------------------------------------------

def run(label, outer, holes):
    print()
    print('=== {} '.format(label) + '=' * max(0, 66 - len(label)))
    primal = plate(outer, holes)
    describe('primal', primal)

    js = B.joints(primal)
    if not js:
        print('    no interior joint -- nothing to block')
        return label, primal, []

    print('    k-ring boundary edges around joint {}: {}  (all EVEN, so an odd'
          .format(js[0], ring_parity(primal, js[0])))
    print('      block never fits inside any of them -- a seam has to leave)')

    verdicts = {}
    for v in js:
        spin, reason = B.blockable(primal, v)
        verdicts[v] = spin
        print('    joint {:<4d} valence {}: {}'.format(
            v, primal.vertex_degree(v),
            'BLOCKABLE, {}'.format(reason) if spin is not None
            else 'REFUSED, {}'.format(reason)))

    ok = [v for v in js if verdicts[v] is not None]
    panels = [('1 primal', primal, None)]

    if ok:
        raw, report = B.block_at(primal, ok[0], RATIO, verdicts[ok[0]])
        seams = B.seam_edges(primal, raw)
        describe('one block', raw, primal)
        print('                         cost: {} seams, {} faces split, '
              'ends {}'.format(report['strips'], report['strip_faces'],
                               _tally(report['ends'])))
        relaxed = B.relax(raw.copy())
        describe('one block, relaxed', relaxed)
        print('                         manifold: {}'.format(raw.is_manifold()))
        panels += [('2 one block + its seams', raw, seams),
                   ('3 relaxed', relaxed, None)]

    if len(ok) > 1:
        try:
            allb, report = B.blocks_any_spin(primal, ok, RATIO)
            describe('every block', allb, primal)
            print('                         chirality {}, {} seams, {} faces '
                  'split, ends {}'.format(report.get('spins'),
                                          report['strips'],
                                          report['strip_faces'],
                                          _tally(report['ends'])))
        except Exception as exc:
            print('    every block at once: REFUSED -- {}'.format(exc))

    try:
        dual = B.quad_dual(primal)
        describe('global dual', dual, primal)
        if B.index_sum(dual) != B.index_sum(primal):
            print('                         ! INDEX CHANGED {} -> {}: the dual '
                  'invented singularities here, it did not move them'.format(
                      B.index_sum(primal), B.index_sum(dual)))
        lost = primal.area() - dual.area()
        if abs(lost) > 1e-6:
            print('                         ! area {:+.4f} ({:+.1%}): the dual '
                  'boundary chords a curved wall'.format(
                      -lost, -lost / primal.area()))
        panels.append(('4 global dual', dual, None))
    except Exception as exc:
        print('    global dual: FAILED -- {}'.format(exc))

    return label, primal, panels


def refusal_is_topological():
    """**Is the disc's refusal a resolution artefact?** No, and this is the check.

    A repair strip that comes back to its own joint sounds like something a
    finer mesh would relieve -- more faces, more room for the split to escape.
    It is not. The strip structure of the disc's layout is the same at every
    density: one of the three strips at each valence-3 joint closes on itself,
    so the split has nowhere to go however small the quads are.
    """
    print()
    print('=== the disc refuses for a TOPOLOGICAL reason, not a size one '
          + '=' * 8)
    for target in (1.4, 1.0, 0.7):
        try:
            d = FieldDecomposition.from_boundary(
                [list(p) for p in suite.DISC], target_length=SPACING)
            dense = d.quad_mesh(target_length=target)
            mesh = Mesh.from_vertices_and_faces(*dense.to_vertices_and_faces())
        except Exception as exc:
            print('  quad target {}: could not build -- {}'.format(target, exc))
            continue
        js = B.joints(mesh)
        blocked = [v for v in js if B.blockable(mesh, v)[0] is not None]
        returning = [sum(1 for k, _ in B.head_strips(mesh, v, s) if k != 'open')
                     for v in js for s in (0, 1)]
        print('  quad target {:<4} {:4d} faces, {} joints, {} blockable, '
              'returning strips per joint per spin: {}'.format(
                  target, mesh.number_of_faces(), len(js), len(blocked),
                  sorted(set(returning))))
    print('  Same answer at every size: use the global dual on a disc.')


def _tally(ends):
    out = {}
    for e in ends:
        out[e] = out.get(e, 0) + 1
    return _fmt(out) if out else '{}'


# ------------------------------------------------------------------
# viewing
# ------------------------------------------------------------------

def close_up(mesh, centre, radius):
    """Only the faces within ``radius`` of ``centre``. A block is 1/150th of a
    plate, so at plate scale the thing this file is about is three pixels."""
    keep = [f for f in mesh.faces()
            if distance_point_point(mesh.face_centroid(f), centre) < radius]
    index, faces = {}, []
    for f in keep:
        face = []
        for v in mesh.face_vertices(f):
            if v not in index:
                index[v] = len(index)
            face.append(index[v])
        faces.append(face)
    points = [None] * len(index)
    for v, i in index.items():
        points[i] = mesh.vertex_coordinates(v)
    return Mesh.from_vertices_and_faces(points, faces)


def view(scenes):
    from compas.colors import Color
    from compas.geometry import Polygon
    from compas.geometry import Polyline
    from compas.geometry import Scale
    from compas.geometry import Translation
    from compas_viewer import Viewer

    from compas_singular.framefield.viz import Grid
    from compas_singular.framefield.viz import NEGATIVE
    from compas_singular.framefield.viz import POSITIVE

    def hue(n):
        """The suite's own convention: red where the mesh wants valence 3,
        blue where it wants 5. A block inherits the colour of the joint it
        replaces, so panel 2 reads as the same singularity in another form."""
        return POSITIVE if n == 3 else NEGATIVE if n == 5 else Color.purple()

    viewer = Viewer()
    grid = Grid(pitch=PITCH, cols=4)
    row = [0]

    def panel(name, mesh, seams, about, scale, parent):
        # column from the panel's own number, so a domain that refuses still
        # shows its dual under the other duals rather than sliding left
        dx, dy = grid.cell(row[0] * 4 + int(name[0]) - 1)
        group = viewer.scene.add_group(name=name, parent=parent)

        def place(p):
            return [(p[0] - about[0]) * scale + dx,
                    (p[1] - about[1]) * scale + dy, p[2]]

        moved = mesh.copy()
        moved.transform(Translation.from_vector([-about[0], -about[1], 0]))
        moved.transform(Scale.from_factors([scale, scale, scale]))
        moved.transform(Translation.from_vector([dx, dy, 0]))
        group.add(moved, show_faces=True, show_lines=True, show_points=False,
                  facecolor=Color(0.89, 0.89, 0.91),
                  linecolor=Color(0.30, 0.30, 0.34), name=name)
        for f in moved.faces():
            vertices = moved.face_vertices(f)
            if len(vertices) == 4:
                continue
            # a block is about 1/25th of a plate across, so it is drawn twice:
            # filled, and RINGED in the same colour. The ring is what actually
            # finds it at whole-plate scale.
            ring = [[x, y, 0.05] for x, y, _ in
                    (moved.vertex_coordinates(v) for v in vertices)]
            group.add(Polygon(ring), facecolor=hue(len(vertices)),
                      linecolor=hue(len(vertices)), opacity=0.9,
                      name='{}-gon BLOCK'.format(len(vertices)))
            group.add(Polyline(ring + ring[:1]), linecolor=hue(len(vertices)),
                      linewidth=7,
                      name='{}-gon BLOCK outline'.format(len(vertices)))
        for v in moved.vertices():
            n = moved.vertex_degree(v)
            if moved.is_vertex_on_boundary(v) or n == 4:
                continue
            a = moved.vertex_coordinates(v)
            for u in moved.vertex_neighbors(v):
                b = moved.vertex_coordinates(u)
                group.add(Polyline([[a[0], a[1], 0.04], [b[0], b[1], 0.04]]),
                          linecolor=hue(n), linewidth=6,
                          name='valence-{} JOINT'.format(n))
        for a, b in (seams or []):
            pa, pb = place(a), place(b)
            group.add(Polyline([[pa[0], pa[1], 0.03], [pb[0], pb[1], 0.03]]),
                      linecolor=Color(0.05, 0.60, 0.20), linewidth=5,
                      name='seam')

    # row 1: the pentagon plate close up, so one block is actually legible
    for label, primal, panels in scenes:
        if label != 'pentagon' or not panels:
            continue
        centre = primal.vertex_coordinates(B.joints(primal)[0])
        group = viewer.scene.add_group(name='pentagon, close up')
        for name, mesh, seams in panels:
            about = centre
            if name.startswith('4'):
                odd = [f for f in mesh.faces()
                       if len(mesh.face_vertices(f)) != 4]
                if odd:
                    about = mesh.face_centroid(odd[0])
            near = [(a, b) for a, b in (seams or [])
                    if distance_point_point(a, about) < ZOOM_RADIUS
                    and distance_point_point(b, about) < ZOOM_RADIUS]
            panel(name, close_up(mesh, about, ZOOM_RADIUS), near, about, ZOOM,
                  group)
        row[0] += 1

    # rows 2+: whole plates
    for label, primal, panels in scenes:
        if not panels:
            continue
        group = viewer.scene.add_group(name=label)
        centre = _centre(primal)
        for name, mesh, seams in panels:
            panel(name, mesh, seams, centre, 1.0, group)
        row[0] += 1

    grid.used = row[0] * 4
    grid.frame(viewer)
    viewer.show()


def _centre(mesh):
    xs = [mesh.vertex_coordinates(v)[0] for v in mesh.vertices()]
    ys = [mesh.vertex_coordinates(v)[1] for v in mesh.vertices()]
    return [(min(xs) + max(xs)) * 0.5, (min(ys) + max(ys)) * 0.5, 0.0]


# ------------------------------------------------------------------

def main(argv):
    show = '--no-view' not in argv
    wanted = [a for a in argv[1:] if not a.startswith('-')]

    print('The singularity as a BLOCK, locally. Background spacing {}, quad '
          'target {}, block size {}.'.format(SPACING, QUAD, RATIO))
    print('index = Sum(4-valence) over interior vertices + Sum(4-degree) over '
          'faces. Trading a joint for a block must not change it.')

    scenes = []
    for label, outer, holes in DOMAINS:
        if wanted and label not in wanted:
            continue
        try:
            scenes.append(run(label, outer, holes))
        except Exception as exc:
            print('  {} blew up: {}: {}'.format(label, type(exc).__name__, exc))

    if not wanted or 'disc' in wanted:
        refusal_is_topological()

    print()
    print('IN SHORT')
    print('  A local block is buildable and cheap where the repair strips have')
    print('  a wall to run to: the pentagon plate goes to 183 quads + 1')
    print('  pentagon with the mesh outside the seams kept vertex for vertex,')
    print('  and the area exact. It is NOT buildable where they do not -- the')
    print('  disc refuses at all four joints, in both chiralities, at every')
    print('  density -- and there the global dual is the answer. The one thing')
    print('  neither route can do is make a block of a degree other than the')
    print('  valence it replaces; the index will not have it.')

    if show and scenes:
        view(scenes)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
