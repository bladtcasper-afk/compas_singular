"""
********************************************************************************
compas_singular.symmetry
********************************************************************************

**Symmetry for any meshing route: find it, mesh one unit, expand it.**

.. code-block:: python

    report = decomposition.find_symmetry()
    coarse_unit = decomposition.symmetry_unit(keys=('M0', 'M45'))
    coarse_unit.collect_strips()
    coarse_unit.set_strips_density_target(0.5)
    quad_unit = coarse_unit.quad_mesh()
    quad_mesh = quad_unit.expand_symmetrically()

The route (skeleton, field, or anything that turns a domain into a coarse quad
layout) only ever meshes the UNIT -- one fundamental region of the chosen group,
cut out along seams. The global mesh is the unit's images under the group,
welded along the seams, so it is symmetric by construction rather than repaired
into symmetry afterwards.

.. currentmodule:: compas_singular.symmetry

.. autosummary::
    :toctree: generated/

    find_symmetry
    SymmetryReport
    SymmetryGroup
    Domain
    SymmetricUnit
    SymmetricQuadUnit
"""
from __future__ import absolute_import
from __future__ import annotations

from typing import Any

from compas_singular.symmetry.group import SymmetryGroup  # noqa: F401
from compas_singular.symmetry.domain import Domain  # noqa: F401
from compas_singular.symmetry.report import SymmetryReport  # noqa: F401
from compas_singular.symmetry.detect import find_symmetry  # noqa: F401


def __getattr__(name: str) -> Any:
    # The unit classes import compas_singular.datastructures, which is heavier and
    # would make ``compas_singular.symmetry`` impossible to import from inside
    # ``compas_singular.datastructures``' own initialisation. Resolve them lazily.
    if name in ('SymmetricUnit', 'SymmetricQuadUnit', 'build_unit'):
        from compas_singular.symmetry import unit
        return getattr(unit, name)
    raise AttributeError(name)


__all__ = ['find_symmetry', 'SymmetryReport', 'SymmetryGroup', 'Domain',
           'SymmetricUnit', 'SymmetricQuadUnit', 'build_unit']
