#! python3

# r: compas

"""**Strip densities, and the one identity they are safe to be stored against.**

Shared by step 5 (``CMD_densities``, where the user sets them) and step 6
(``CMD_quad_mesh``, where they are applied). Both files carried a byte-identical
copy of every function below; a fix to one silently did not reach the other, and
a reader had no way to tell which copy was authoritative.

Not named ``CMD_`` on purpose: it is not a command, and Rhino picks commands up
by filename.

**Densities are stored against geometry, never against the strip index, and that
is the one trap in this whole design.** ``collect_strips`` builds its list by
popping ``self.edges()``, so indices follow vertex insertion order -- and a
layout that has been baked to Rhino and read back has whatever order the
geometry came back in. Measured: bake a 4-strip layout, read it back, **0 of 4
indices still point at the same strip**, permuted 0->2 1->3 2->0 3->1. A density
map keyed by index does not raise; it silently applies the wrong numbers to the
wrong bands. So each density is stored against the rounded midpoint of the edge
that was picked, and the strip owning that edge is looked up again on the other
side. Editing the layout in step 4 can move that midpoint, and an override that
no longer finds its strip is reported as lost rather than dropped in silence.

The library has the same idea under a different key --
``compas_singular.agent.core.rebuild.apply_densities`` over
``address.strip_address``. It is not used here because the key format below is
what is already written into open documents, and changing it would silently
orphan every density a user has set.
"""

#Temporary import of compas_singular development library
from CMD_start import import_compas_singular
import_compas_singular()

import json

import rhinoscriptsyntax as rs

from compas.tolerance import TOL

from CMD_start import get_settings


def edge_key(a, b):
    """The identity a density is stored against: the edge's rounded midpoint.

    Order-independent by construction, and rounded the way the rest of the
    codebase rounds -- ``TOL.geometric_key``, 3 decimals, the same resolution
    ``from_polylines`` matches endpoints at.
    """
    mid = [(a[i] + b[i]) / 2.0 for i in range(3)]
    return TOL.geometric_key(mid)


def get_densities():
    """The overrides stored on the document. ``{}`` when none were set."""
    blob = rs.GetDocumentUserText(get_settings()["density_key"])
    return json.loads(blob) if blob else {}


def set_densities(densities):
    """Write the overrides back to the document."""
    rs.SetDocumentUserText(get_settings()["density_key"], json.dumps(densities))
    return densities


def strip_edge_map(coarse):
    """``{edge midpoint key: skey}`` for every edge of every strip.

    Built once and shared: this is the lookup that replaces the strip index
    everywhere it would otherwise be used across a bake.
    """
    out = {}
    for skey in coarse.strips():
        for u, v in coarse.strip_edges(skey):
            if u == v:
                continue
            out[edge_key(coarse.vertex_coordinates(u),
                         coarse.vertex_coordinates(v))] = skey
    return out


def strip_of_edge(coarse, a, b, lookup=None):
    """Which strip owns the edge between ``a`` and ``b``. ``None`` if no match."""
    lookup = strip_edge_map(coarse) if lookup is None else lookup
    return lookup.get(edge_key(a, b))


def apply_densities(coarse, target_length, densities=None, verbose=True):
    """Set every strip's density: the target length, then the stored overrides.

    Returns ``(applied, lost)`` -- how many overrides found their strip, and how
    many did not because the layout changed under them. A lost override is worth
    saying out loud: quietly meshing at the default is how a density the user set
    goes missing.
    """
    if densities is None:
        densities = get_densities()

    coarse.collect_strips()
    coarse.set_strips_density_target(target_length)

    lookup = strip_edge_map(coarse)
    applied = lost = 0
    for key, density in densities.items():
        skey = lookup.get(key)
        if skey is None:
            lost += 1
            continue
        coarse.set_strip_density(skey, int(density))
        applied += 1

    if verbose and (applied or lost):
        print("densities: {} override(s) applied, {} lost to a changed layout"
              .format(applied, lost))
    return applied, lost
