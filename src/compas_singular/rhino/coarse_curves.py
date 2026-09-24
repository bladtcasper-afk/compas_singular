"""Moved to :mod:`compas_singular.datastructures.mesh_quad_coarse.coarse_curves`.

Nothing in it ever touched Rhino -- it rebuilds a coarse layout's edge curvature
from point lists alone -- and living under ``compas_singular.rhino`` meant
``algorithms`` could not reach it without importing the CAD package.
``SkeletonDecomposition.edges_to_curves`` needs exactly that, so it moved.

This module stays so the ``CMD_`` commands keep importing from where they always
have. Import from ``compas_singular.datastructures`` in anything new.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import annotations

from compas_singular.datastructures.mesh_quad_coarse.coarse_curves import BoundaryLoop  # noqa: F401
from compas_singular.datastructures.mesh_quad_coarse.coarse_curves import coarse_edges_to_curves  # noqa: F401
from compas_singular.datastructures.mesh_quad_coarse.coarse_curves import mean_edge_length  # noqa: F401
from compas_singular.datastructures.mesh_quad_coarse.coarse_curves import snap_corners_to_walls  # noqa: F401


__all__ = [
    'coarse_edges_to_curves',
    'snap_corners_to_walls',
    'mean_edge_length',
    'BoundaryLoop',
]
