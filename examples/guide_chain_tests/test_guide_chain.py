"""Checks on the automatic guide selection. No Rhino, no pytest, no viewer.

Everything asserted here was measured first with :mod:`eval_guide_chain`; these are the
properties that must not regress, not a wish list. Run:

    python examples/guide_chain_tests/test_guide_chain.py
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys
import traceback

from compas.geometry import Point

from compas_singular.datastructures import QuadMesh
from compas_singular.editing import GuideCurve
from compas_singular.editing import attach_chain
from compas_singular.editing import chain_quality
from compas_singular.editing import collect_polyedges
from compas_singular.editing import guide_chain
from compas_singular.editing import mean_edge_length

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_meshes import cases  # noqa: E402


PASSED, FAILED = [], []


def check(condition, message):
    (PASSED if condition else FAILED).append(message)
    print('  {} {}'.format('ok  ' if condition else 'FAIL', message))


def grid(n=6):
    """A plain n x n grid: the one mesh whose right answer can be written down."""
    vertices = [[i, j, 0] for j in range(n + 1) for i in range(n + 1)]
    faces = [[j * (n + 1) + i, j * (n + 1) + i + 1,
              (j + 1) * (n + 1) + i + 1, (j + 1) * (n + 1) + i]
             for j in range(n) for i in range(n)]
    return QuadMesh.from_vertices_and_faces(vertices, faces)


# ---------------------------------------------------------------------------
# On a grid, where the answer is known by hand
# ---------------------------------------------------------------------------

def test_grid():
    print('grid 6x6')
    mesh = grid()
    polyedges = collect_polyedges(mesh)

    chain, info = guide_chain(mesh, GuideCurve([[0, 3, 0], [6, 3, 0]]), polyedges=polyedges)
    check(chain == [21, 22, 23, 24, 25, 26, 27], 'a guide down a row selects that whole row')
    check(info['coverage'] == 1.0, 'and covers the guide end to end')
    quality = chain_quality(mesh, chain, GuideCurve([[0, 3, 0], [6, 3, 0]]))
    check(quality['purity'] == 1.0 and quality['runs'] == 1, 'as one straight run')
    check(quality['off_max'] == 0.0, 'sitting exactly on the guide')
    check(quality['min_face_angle'] == 90.0, 'and projecting it changes nothing')

    across = GuideCurve([[3, 0, 0], [3, 6, 0]])
    chain, _ = guide_chain(mesh, across, polyedges=polyedges)
    check(chain == [3, 10, 17, 24, 31, 38, 45],
          'a guide across the other family selects the column, not the row')
    check(chain_quality(mesh, chain, across)['alignment'] > 0.99,
          'both families are equally selectable -- alignment is measured, not assumed')

    short = GuideCurve([[1.4, 3.0, 0], [4.6, 3.0, 0]])
    chain, _ = guide_chain(mesh, short, polyedges=polyedges)
    check(chain == [22, 23, 24, 25, 26], 'a short guide selects a TRIMMED run, not the row')

    far = GuideCurve([[0, 50, 0], [6, 50, 0]])
    chain, info = guide_chain(mesh, far, polyedges=polyedges)
    check(chain == [] and info['reason'], 'a guide nowhere near the mesh returns a reason')

    diagonal = GuideCurve([[0, 0, 0], [6, 6, 0]])
    chain, info = guide_chain(mesh, diagonal, polyedges=polyedges)
    check(chain_quality(mesh, chain, diagonal)['alignment'] < 0.8 if chain else True,
          'a 45 degree guide is reported as poorly aligned, not forced into a staircase')


# ---------------------------------------------------------------------------
# The boundary rule
# ---------------------------------------------------------------------------

def test_boundary():
    print('boundary rule')
    mesh = grid()
    polyedges = collect_polyedges(mesh)
    average = mean_edge_length(mesh)
    hugging = GuideCurve([[0.2, 0.4, 0], [5.8, 0.4, 0]])   # nearer the wall than the row

    chain, _ = guide_chain(mesh, hugging, tolerance=average, boundary='free',
                           polyedges=polyedges)
    quality = chain_quality(mesh, chain, hugging)
    check(quality['boundary_interior'] > 0,
          "without the rule the outline itself is selected ('free' is the old behaviour)")

    chain, _ = guide_chain(mesh, hugging, tolerance=average, boundary='anchor',
                           polyedges=polyedges)
    quality = chain_quality(mesh, chain, hugging)
    check(quality['boundary_interior'] == 0, 'with it, no boundary vertex is in the middle')
    check(quality['boundary_ends'] == 2, 'but the two ends may anchor on the wall')

    chain, _ = guide_chain(mesh, hugging, tolerance=average, boundary='exclude',
                           polyedges=polyedges)
    quality = chain_quality(mesh, chain, hugging)
    check(quality['boundary_interior'] == 0 and quality['boundary_ends'] == 0,
          "'exclude' drops the anchors as well")


# ---------------------------------------------------------------------------
# Every guide on every mesh
# ---------------------------------------------------------------------------

def _sweep(tolerance_factor=None):
    """Every mesh, every guide, both boundary modes -> the quality of each chain."""
    out = []
    for _, mesh, guides in cases():
        polyedges = collect_polyedges(mesh)
        tolerance = (None if tolerance_factor is None
                     else tolerance_factor * mean_edge_length(mesh))
        for _, points in guides:
            guide = GuideCurve(points)
            for mode in ('anchor', 'exclude'):
                chain, _ = guide_chain(mesh, guide, tolerance=tolerance, boundary=mode,
                                       polyedges=polyedges)
                if chain:
                    out.append(chain_quality(mesh, chain, guide))
    return out


def test_all_cases():
    print('all meshes, all guides, all modes')
    sweep = _sweep()
    check(len(sweep) >= 20, 'the sweep actually ran ({} chains)'.format(len(sweep)))
    check(all(q['boundary_interior'] == 0 for q in sweep),
          'no chain has a boundary vertex anywhere but at an end')
    check(all(q['valid_path'] for q in sweep),
          'every chain is a connected path of real edges, no vertex twice')
    check(all(q['purity'] == 1.0 and q['runs'] == 1 for q in sweep),
          'every chain is exactly one run of one polyedge')


def test_folding_is_the_price_of_reach():
    """Where the folds come from, pinned so nobody has to rediscover it.

    A vertex is pulled onto the guide by exactly as far as it sat off it, so a wide
    distance tolerance inverts the faces around the vertices it reaches for. The
    SELECTION is not what breaks -- every chain above is still one clean run -- so this
    is a property of the reach that was asked for, not a defect to be fixed in the
    selection. Casper set the default at 2.0 knowing this.
    """
    print('folding is the price of reach, not a selection fault')
    tight = _sweep(0.8)
    wide = _sweep(2.0)
    check(all(q['folded_faces'] == 0 for q in tight),
          'at 0.8 edge lengths nothing folds')
    folded = [q for q in wide if q['folded_faces']]
    check(bool(folded), 'at 2.0 some chains do ({} of {})'.format(len(folded), len(wide)))
    check(all(q['off_max'] > 1.0 for q in folded),
          'and every one of them was reaching more than a whole edge length to do it')
    check(all(q['purity'] == 1.0 and q['valid_path'] for q in folded),
          'the chains themselves are still clean -- it is the projection that folds')


def test_boundary_never_moves():
    """A boundary vertex may be selected, but it may not be taken off its wall."""
    print('a selected boundary vertex slides, it does not move')
    mesh = grid()
    polyedges = collect_polyedges(mesh)
    hugging = GuideCurve([[0.2, 0.4, 0], [5.8, 0.4, 0]])

    chain, _ = guide_chain(mesh, hugging, tolerance=mean_edge_length(mesh),
                           polyedges=polyedges)
    boundary = [vertex for vertex in chain if mesh.is_vertex_on_boundary(vertex)]
    check(len(boundary) == 2, 'the chain has two boundary ends to hold')

    moves, constraints = attach_chain(mesh, chain, hugging)
    check(all(vertex not in moves for vertex in boundary),
          'no boundary vertex is given a new position')
    check(all(vertex in moves for vertex in chain if vertex not in boundary),
          'every interior vertex is')
    check(all(not isinstance(constraints[vertex], Point) for vertex in boundary),
          'a boundary vertex is constrained to a curve -- it slides, it is not pinned')
    check(all(isinstance(constraints[vertex], Point)
              for vertex in chain if vertex not in boundary),
          "and an interior one is pinned under hold='fixed'")

    _, sliding = attach_chain(mesh, chain, hugging, hold='sliding')
    check(all(not isinstance(sliding[vertex], Point) for vertex in chain),
          "under hold='sliding' nothing is pinned")

    drift = 0.0
    for _, mesh, guides in cases():
        polyedges = collect_polyedges(mesh)
        for _, points in guides:
            guide = GuideCurve(points)
            chain, _ = guide_chain(mesh, guide, polyedges=polyedges)
            if chain:
                drift = max(drift, chain_quality(mesh, chain, guide)['end_drift'])
    check(drift == 0.0,
          'over every mesh and guide, no attached vertex ends up off its outline')


def test_angle_gate():
    """The angle gate drops what the distance gate keeps: a line that has veered off."""
    print('the angle metric')
    mesh = grid()
    polyedges = collect_polyedges(mesh)

    diagonal = GuideCurve([[0, 0, 0], [6, 6, 0]])
    chain, info = guide_chain(mesh, diagonal, tolerance=3.0, max_angle=30.0,
                              polyedges=polyedges)
    check(not chain and 'degrees' in info['reason'],
          'a 45 degree guide is refused by angle, however close the mesh is')
    chain, _ = guide_chain(mesh, diagonal, tolerance=3.0, max_angle=90.0,
                           polyedges=polyedges)
    check(bool(chain), 'and selected again with the gate off')

    across = GuideCurve([[3, 0, 0], [3, 6, 0]])
    chain, _ = guide_chain(mesh, across, polyedges=polyedges)
    check(chain == [3, 10, 17, 24, 31, 38, 45],
          'the gate is not cross-symmetric: a perpendicular course is still selected '
          'for a guide that runs along it')

    # the real case: on the L-plate the elbow guide keeps a vertex whose polyedge has
    # already turned away, and that vertex is the one that squashes a face
    plate = [case for case in cases() if case[0] == 'L-plate'][0]
    _, plate_mesh, plate_guides = plate
    plate_polyedges = collect_polyedges(plate_mesh)
    elbow = GuideCurve(dict(plate_guides)['elbow'])
    average = mean_edge_length(plate_mesh)
    loose, _ = guide_chain(plate_mesh, elbow, tolerance=1.0 * average, max_angle=90.0,
                           polyedges=plate_polyedges)
    gated, _ = guide_chain(plate_mesh, elbow, tolerance=1.0 * average, max_angle=30.0,
                           polyedges=plate_polyedges)
    check(len(gated) < len(loose), 'the gate shortens a chain that veers off at the end')
    check(chain_quality(plate_mesh, gated, elbow)['min_face_angle']
          > chain_quality(plate_mesh, loose, elbow)['min_face_angle'],
          'and the vertices it drops are the ones that were squashing a face')


def test_no_collapsed_edge():
    print('a guide shorter than the polyedge it follows')
    # everything past the end of a guide projects onto the SAME endpoint, so two such
    # vertices would land on one point and the edge between them would vanish
    mesh = grid()
    polyedges = collect_polyedges(mesh)
    guide = GuideCurve([[1.5, 3.0, 0], [4.5, 3.0, 0]])
    chain, _ = guide_chain(mesh, guide, tolerance=10.0, polyedges=polyedges)
    projected = [tuple(round(value, 6) for value in guide.closest_point(
        mesh.vertex_coordinates(vertex))) for vertex in chain]
    check(len(set(projected)) == len(projected),
          'no two selected vertices project onto the same point')
    check(chain_quality(mesh, chain, guide)['min_face_angle'] > 0.0,
          'and no face is left with a zero angle')


if __name__ == '__main__':
    for test in (test_grid, test_boundary, test_boundary_never_moves, test_angle_gate,
                 test_all_cases, test_folding_is_the_price_of_reach,
                 test_no_collapsed_edge):
        try:
            test()
        except Exception:                                   # noqa: BLE001
            FAILED.append(test.__name__)
            traceback.print_exc()
        print()
    print('{} passed, {} failed'.format(len(PASSED), len(FAILED)))
    sys.exit(1 if FAILED else 0)
