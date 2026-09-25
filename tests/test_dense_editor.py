"""Strip deletion on :class:`DenseMeshEditor`: one gated path, one plan shape."""
import pytest

from compas_singular.datastructures import QuadMesh
from compas_singular.editing import DenseMeshEditor


def _grid(n):
    vertices = [[float(x), float(y), 0.0] for y in range(n + 1) for x in range(n + 1)]
    faces = [[y * (n + 1) + x, y * (n + 1) + x + 1,
              (y + 1) * (n + 1) + x + 1, (y + 1) * (n + 1) + x]
             for y in range(n) for x in range(n)]
    return QuadMesh.from_vertices_and_faces(vertices, faces)


@pytest.fixture
def editor():
    return DenseMeshEditor(_grid(3))


def test_a_refused_plan_has_the_same_keys_as_an_accepted_one(editor):
    ok_plan = editor.plan_strip_deletion((0, 1))
    bad_plan = editor.plan_strip_deletion((0, 99))
    assert ok_plan['ok'] and not bad_plan['ok']
    assert set(ok_plan) == set(bad_plan)


def test_a_non_edge_is_reported_as_a_non_edge(editor):
    plan = editor.plan_strip_deletion((0, 99))
    assert 'not an edge' in plan['reason']


def test_remove_strip_runs_the_same_checks_as_remove_line(editor):
    # a pentagon on the strip through (0, 1) -- remove_line refuses it, so must remove_strip
    editor.mesh.delete_face(0)
    editor.mesh.delete_face(1)
    editor.mesh.add_face([0, 1, 2, 6, 5, 4])
    ok_line, notes_line = editor.remove_line((0, 4))
    ok_strip, notes_strip = editor.remove_strip((0, 4))
    assert not ok_line and not ok_strip
    assert notes_strip['error'] == notes_line['error']


def test_a_strip_can_be_planned_and_removed_by_key(editor):
    work = editor.mesh.copy()
    work.collect_strips()
    skey = work.edge_strip((0, 1))
    plan = editor.plan_strip_deletion(skey=skey)
    assert plan['ok'], plan['reason']
    ok, notes = editor.remove_strip(skey=skey)
    assert ok, notes
    assert notes['faces_after'] == plan['faces_after']


def test_pre_split_is_always_a_dict(editor):
    ok, notes = editor.remove_line((0, 1))
    assert ok, notes
    assert notes['pre_split'] == {}
