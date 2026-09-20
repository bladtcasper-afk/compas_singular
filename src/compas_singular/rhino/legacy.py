"""**Reading a document written before the session existed.**

Until 2026-09-18 the commands kept their state in the document itself and beside
it: the settings as document user text under ``"settings"``, the layout baked on
``Skeleton::Mesh`` + ``Poles`` with its separatrices on ``Skeleton::Polylines``,
and the layout's attributes and the field as JSON side-cars in a hidden folder
beside the ``.3dm`` (``.<name>.compas_singular/``). The first time a command
runs on such a document, :func:`read_legacy` fills a new session from them, and
that command's :meth:`~RhinoSession.record` moves them into the document for good.

The layout, best first: the side-car, if it still describes the drawn mesh (it
carries strips, densities and patterns); otherwise the drawn mesh itself, whose
strips and densities are then re-derived. Only a SAVED document's side-cars are
read: an unsaved one used a single temp folder shared by every unsaved document,
so what is there may belong to any of them.

Delete this module when no ``.3dm`` from before that date matters any more.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import json
import os

from compas_singular.settings import Settings


__all__ = ['read_legacy']


def read_legacy(session):
    """Fill ``session`` from its document's old settings, drawing and side-cars."""
    doc = session.doc
    text = doc.Strings.GetValue('settings')
    if text:
        session.settings = Settings.model_validate(json.loads(text))

    drawn = _drawn_layout()
    coarse = _side_car(doc, 'coarse.json', 'CoarsePseudoQuadMesh')
    if coarse is None or drawn is None or not _same_corners(coarse, drawn):
        coarse = drawn
    if coarse is not None and not coarse.shape_polylines():
        coarse.set_shape_polylines(_drawn_polylines())
    session.coarse = coarse
    session.field = _side_car(doc, 'field.json', 'CrossField')


def _side_car(doc, filename, kind):
    """A side-car beside a SAVED document, or ``None``."""
    if not doc.Path:
        return None
    path = os.path.join(os.path.dirname(doc.Path),
                        '.' + os.path.splitext(doc.Name)[0] + '.compas_singular', filename)
    if kind == 'CrossField':
        from compas_singular.framefield.field import CrossField
        return CrossField.load_from_json(path, default=None)
    from compas_singular.datastructures import CoarsePseudoQuadMesh
    return CoarsePseudoQuadMesh.load_from_json(path, default=None)


def _drawn_layout():
    """The layout baked on ``Skeleton::Mesh`` with its poles, or ``None``."""
    from .helpers import read_coarse
    try:
        coarse, _poles = read_coarse()
    except RuntimeError:                 # nothing drawn
        return None
    return coarse


def _drawn_polylines():
    """The separatrices baked on ``Skeleton::Polylines``."""
    import rhinoscriptsyntax as rs

    from .helpers import read_polylines
    from .project import layer_path

    layer = layer_path('Polylines')
    return read_polylines(layer) if rs.IsLayer(layer) else []


def _same_corners(one, other):
    """Whether two layouts have the same ROUNDED corners -- a bake renumbers and
    stores single precision, so neither order nor exact coordinates prove anything."""
    from compas.tolerance import TOL

    def keys(mesh):
        return set(TOL.geometric_key(mesh.vertex_coordinates(v)) for v in mesh.vertices())

    return keys(one) == keys(other)
