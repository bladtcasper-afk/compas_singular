"""Tests for the QuadMesh polyedge/strip machinery.

Covers the batch 1-2 port of compas_quad: the null guards, the accessors, and the
opt-in flags (``strict``, ``oriented``, ``both_sides``, ``legacy``) that carry the
new behaviour without moving any default.
"""
import pytest

from compas_singular.datastructures import QuadMesh


GRID_VERTICES = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0],
                 [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0],
                 [0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0]]
GRID_FACES = [[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]]


@pytest.fixture
def grid():
    """Standard 2 x 2 quad grid."""
    return QuadMesh.from_vertices_and_faces(GRID_VERTICES, GRID_FACES)


@pytest.fixture
def annulus():
    """Four quads around a square hole -- two nested boundaries."""
    vertices = [[0, 0, 0], [3, 0, 0], [3, 3, 0], [0, 3, 0],
                [1, 1, 0], [2, 1, 0], [2, 2, 0], [1, 2, 0]]
    faces = [[0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]
    return QuadMesh.from_vertices_and_faces(vertices, faces)


@pytest.fixture
def kinked_grid():
    """2 x 2 grid with an extra quad hanging off vertex 1.

    Vertex 1 is then a degree-4 vertex *on the boundary*, i.e. singular -- the
    case a naive ``len(nbrs) == 4`` rule walks straight through.
    """
    vertices = GRID_VERTICES + [[1.0, -1.0, 0.0], [2.0, -1.0, 0.0]]
    faces = GRID_FACES + [[1, 9, 10, 2]]
    return QuadMesh.from_vertices_and_faces(vertices, faces)


# ----------------------------------------------------------------------------
# null guards
# ----------------------------------------------------------------------------

def test_face_opposite_edge_returns_none_on_a_boundary_halfedge(grid):
    u, v = 1, 0  # halfedge on the outside of the boundary edge (0, 1)
    assert grid.halfedge[u][v] is None
    assert grid.face_opposite_edge(u, v) is None


def test_face_opposite_edge_still_works_inside_a_face(grid):
    assert grid.face_opposite_edge(0, 1) == (4, 3)


def test_collect_strip_survives_a_boundary_halfedge_seed(grid):
    # would previously raise while unpacking the None
    assert grid.collect_strip(1, 0)


# ----------------------------------------------------------------------------
# accessors
# ----------------------------------------------------------------------------

def test_polyedge_accessors(grid):
    grid.collect_polyedges()
    pkey = next(iter(grid.polyedges()))

    assert grid.number_of_polyedges() == len(list(grid.polyedges()))
    assert grid.polyedge_vertices(pkey) == grid.attributes['polyedges'][pkey]
    assert grid.polyline(pkey) == [grid.vertex_coordinates(v) for v in grid.polyedge_vertices(pkey)]

    edges = grid.polyedge_edges(pkey)
    assert isinstance(edges, list)      # not a one-shot iterator
    assert len(edges) == len(grid.polyedge_vertices(pkey)) - 1
    assert list(edges) == list(edges)   # re-iterable

    assert grid.polyedge_length(pkey) == pytest.approx(2.0)
    assert len(grid.polyedge_midpoint(pkey)) == 3


# ----------------------------------------------------------------------------
# collect_polyedge: both_sides / oriented
# ----------------------------------------------------------------------------

def test_collect_polyedge_default_loses_the_seed_direction(grid):
    # the walk hits an extremity on the first side and returns reversed
    assert grid.collect_polyedge(2, 1) == [0, 1, 2]


def test_collect_polyedge_oriented_keeps_the_seed_direction(grid):
    assert grid.collect_polyedge(2, 1, oriented=True) == [2, 1, 0]


def test_collect_polyedge_one_sided_stops_at_the_first_extremity(grid):
    assert grid.collect_polyedge(2, 1, both_sides=False) == [2, 1, 0]


def test_collect_strip_one_sided_returns_a_list_of_edges(grid):
    strip = grid.collect_strip(1, 0, both_sides=False)
    assert strip == [(1, 0)]


# ----------------------------------------------------------------------------
# vertex_opposite_vertex: strict
# ----------------------------------------------------------------------------

def test_strict_does_not_turn_an_interior_polyedge_onto_a_boundary(annulus):
    # vertex 0 is on the outer boundary, 4 on the inner one, and (0, 4) is an
    # interior edge; the default rule only tests the vertices and turns onto the
    # hole boundary
    assert annulus.vertex_opposite_vertex(0, 4) == 5
    assert annulus.vertex_opposite_vertex(0, 4, strict=True) is None


def test_strict_separates_the_two_boundary_loops(annulus):
    annulus.collect_polyedges(strict=True)
    # each ring closes on itself, joined by four clean radial polyedges
    assert sorted(annulus.attributes['polyedges'].values()) == [
        [2, 3, 0, 1, 2],        # outer ring
        [4, 0], [5, 1], [6, 2],  # radials
        [6, 7, 4, 5, 6],        # inner ring
        [7, 3],                 # radial
    ]


def test_default_walk_mixes_the_two_boundary_loops(annulus):
    """What strict=True fixes -- kept as the record of the default behaviour."""
    annulus.collect_polyedges()
    polyedges = list(annulus.attributes['polyedges'].values())
    # every polyedge drags vertices of both rings along
    inner, outer = {4, 5, 6, 7}, {0, 1, 2, 3}
    assert any(set(p) & inner and set(p) & outer for p in polyedges)


def test_strict_does_not_walk_through_a_four_valent_boundary_vertex(kinked_grid):
    # vertex 1 has four neighbours but sits on the boundary, so it is singular;
    # an unguarded `len(nbrs) == 4` rule would return 2 here, and 9 from vertex 4
    assert kinked_grid.vertex_neighbors(1, ordered=True) == [0, 4, 2, 9]
    assert kinked_grid.vertex_opposite_vertex(0, 1, strict=True) is None
    assert kinked_grid.vertex_opposite_vertex(4, 1, strict=True) is None


def test_strict_defaults_off(annulus):
    """The default walk must keep its historical (wrong) answer."""
    assert annulus.collect_polyedge(0, 4) == [0, 4, 5, 1, 0]


# ----------------------------------------------------------------------------
# polyedge_graph: legacy
# ----------------------------------------------------------------------------

def test_legacy_graph_is_the_default_and_emits_self_loops(grid):
    grid.collect_polyedges()
    _, edges = grid.polyedge_graph()
    assert [e for e in edges if e[0] == e[1]]


def test_new_graph_has_no_self_loops_and_fewer_edges(grid):
    grid.collect_polyedges()
    _, legacy = grid.polyedge_graph(legacy=True)
    _, new = grid.polyedge_graph(legacy=False)
    assert [e for e in new if e[0] == e[1]] == []
    assert len(new) < len(legacy)


def test_new_graph_pairs_the_two_polyedges_crossing_at_a_vertex(grid):
    grid.collect_polyedges()
    nodes, edges = grid.polyedge_graph(legacy=False)
    vertex_to_polyedges = {}
    for pkey, polyedge in grid.polyedges(data=True):
        for vkey in polyedge:
            vertex_to_polyedges.setdefault(vkey, set()).add(pkey)
    expected = sum(1 for vkey, pkeys in vertex_to_polyedges.items()
                   if len(pkeys) == 2 and not grid.is_vertex_singular(vkey))
    assert len(edges) == expected
    assert set(nodes) == set(grid.polyedges())
