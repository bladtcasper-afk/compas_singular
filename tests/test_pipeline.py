"""Headless regression tests for the core compas_singular pipelines.

These exercise the decomposition / densification / grammar logic against the
current COMPAS 2.x API (no visualisation, so they run anywhere).
"""
import os
import json

import pytest

from compas_singular.datastructures import CoarseQuadMesh
from compas_singular.datastructures import CoarsePseudoQuadMesh
from compas_singular.datastructures import QuadMesh
from compas_singular.datastructures.lizard import Lizard
from compas_singular.algorithms import boundary_triangulation
from compas_singular.algorithms import SkeletonDecomposition

HERE = os.path.dirname(__file__)
DATA = os.path.abspath(os.path.join(HERE, '..', 'examples', 'data'))


def test_coarse_quad_densification():
    mesh = CoarseQuadMesh.from_json(os.path.join(DATA, 'coarse_quad_mesh_british_museum.json'))
    assert mesh.number_of_faces() == 12
    mesh.collect_strips()
    mesh.set_strips_density(3)
    mesh.densification()
    assert mesh.get_quad_mesh().number_of_faces() > mesh.number_of_faces()


def test_pseudo_quad_densification_preserves_poles():
    mesh = CoarsePseudoQuadMesh.from_json(os.path.join(DATA, 'coarse_quad_mesh_british_museum_poles.json'))
    # the pole map must survive JSON (de)serialisation with integer face keys
    assert all(isinstance(fkey, int) for fkey in mesh.attributes['face_pole'])
    assert sorted(mesh.poles()) == [0, 3, 16, 19]
    mesh.collect_strips()
    mesh.set_strips_density_target(t=.5)
    mesh.densification()
    assert mesh.get_quad_mesh().number_of_faces() > 0


def test_skeleton_decomposition_planar():
    with open(os.path.join(DATA, '01_decomposition.json')) as f:
        outer_boundary, inner_boundaries, polyline_features, point_features = json.load(f)
    trimesh = boundary_triangulation(outer_boundary, inner_boundaries, polyline_features, point_features)
    decomposition = SkeletonDecomposition.from_mesh(trimesh)
    coarsemesh = decomposition.decomposition_mesh(point_features)
    coarsemesh.collect_strips()
    coarsemesh.set_strips_density_target(0.5)
    coarsemesh.densification()
    densemesh = coarsemesh.get_quad_mesh()
    assert densemesh.number_of_faces() > coarsemesh.number_of_faces()


@pytest.mark.parametrize('string', ['ata', 'atta', 'attta', 'attpptta'])
def test_lizard_grammar_produces_manifold_mesh(string):
    vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0],
                [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0],
                [0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0]]
    faces = [[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]]
    mesh = QuadMesh.from_vertices_and_faces(vertices, faces)
    mesh.collect_strips()
    lizard = Lizard(mesh)
    lizard.initiate()
    lizard.from_string_to_rules(string)
    assert mesh.is_manifold()
    assert all(len(mesh.halfedge[vkey]) > 0 for vkey in mesh.vertices())
