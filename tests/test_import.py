import importlib
import inspect
import types

import pytest

import compas_singular
import compas_singular.datastructures as ds


def test_trivial():
    print(compas_singular.__version__)
    assert True


# ----------------------------------------------------------------------------
# package namespace
#
# ``from .foo import *`` binds ``foo`` as an attribute of the package. When the
# package __all__ is computed with ``dir()`` that module object gets re-exported
# and shadows any function or subpackage of the same name, which used to make
# ``compas_singular.datastructures.add_strip`` a module rather than a function.
# ----------------------------------------------------------------------------

PACKAGES = [
    'compas_singular.algorithms',
    'compas_singular.datastructures',
    'compas_singular.datastructures.lizard',
    'compas_singular.datastructures.mesh',
    'compas_singular.datastructures.mesh_quad',
    'compas_singular.datastructures.mesh_quad.grammar',
    'compas_singular.datastructures.mesh_quad_coarse',
    'compas_singular.datastructures.mesh_quad_pseudo',
    'compas_singular.datastructures.mesh_quad_pseudo_coarse',
    'compas_singular.datastructures.network',
    'compas_singular.datastructures.skeleton',
    'compas_singular.geometry',
    'compas_singular.topology',
    'compas_singular.utilities',
]


@pytest.mark.parametrize('name', PACKAGES)
def test_subpackage_is_importable(name):
    importlib.import_module(name)


@pytest.mark.parametrize('name', PACKAGES)
def test_package_does_not_re_export_submodules(name):
    package = importlib.import_module(name)
    leaked = [n for n in getattr(package, '__all__', [])
              if isinstance(getattr(package, n, None), types.ModuleType)]
    assert leaked == []


@pytest.mark.parametrize('name', ['add_strip', 'add_strips', 'delete_strip', 'delete_strips'])
def test_grammar_functions_are_callable(name):
    assert callable(getattr(ds, name))


@pytest.mark.parametrize('name', ['delete_strip', 'delete_strips'])
def test_deletion_keeps_the_preserve_boundaries_signature(name):
    # algorithms/twocoloring.py and algorithms/mapping.py both rely on this
    assert 'preserve_boundaries' in inspect.signature(getattr(ds, name)).parameters


def test_datastructures_mesh_quad_is_the_package():
    assert hasattr(ds.mesh_quad, '__path__')
