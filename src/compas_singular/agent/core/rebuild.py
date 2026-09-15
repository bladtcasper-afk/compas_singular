"""**Densifying without throwing the caller's work away.**

:meth:`FieldDecomposition.quad_mesh` is the right entry point for a person: it
owns the whole chain and every fallback in it, and it always returns a mesh.
Both of those are wrong for a session that is holding state.

**Its densities are keyed by STRIP KEY.** ``quad_mesh`` used to reset every
strip on every call; it no longer does -- it keeps whatever is already set when
no size argument is given, and takes a ``densities={skey: d}`` map otherwise.
That is enough for a script that builds a layout and meshes it in the same
breath, and not enough here: a strip key is an index into
``attributes['strips']``, rebuilt from the current edge order by
``collect_strips``, so it names a DIFFERENT strip after any edit. A session that
sets densities on turn 4 and re-densifies on turn 30 needs a name that survives
the layout changing under it, which is what ``strip_address`` gives.

**And when densification fails it replaces the layout.** The triangulation
backstop branch ends with ``self.mesh = mesh`` and ``self._edited = False``, so
a coarse layout somebody hand-edited is *gone*, silently, because a density was
too coarse to densify. For a person clicking a button that is the right
trade -- a mesh of the right shape beats an exception. For a session that has
been editing for twenty steps it is the single most destructive thing in the
API.

:func:`densify_preserving` is ``quad_mesh``'s field-route branch with both of
those answered: densities come from an ADDRESS-keyed map the caller owns, and a
failure RAISES :class:`DensifyRefused` with the reason instead of falling back.
Nothing here ever writes ``decomposition.mesh``. Choosing the backstop is then a
decision the caller makes explicitly, having been told what it costs.

**Densities are keyed geometrically**, by :func:`~compas_singular.agent.core.address.strip_address`,
because strip keys renumber on every edit -- see ``address.py``. A map survives a
density change untouched (the coarse geometry does not move) and legitimately
loses entries after a layout edit, which is reported rather than hidden.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ...editing.repair import densifiable
from .address import strip_address


__all__ = ['DensifyRefused', 'apply_densities', 'read_densities',
           'densify_preserving']


class DensifyRefused(Exception):
    """Densification did not produce a usable mesh, and nothing was replaced.

    Attributes
    ----------
    reason : str
        Why, in the words of ``densifiable`` or ``_acceptable``.
    metrics : dict
        Element quality of the rejected mesh, when it got far enough to measure.
    stage : str
        ``'layout'`` when the coarse mesh could not be densified at all,
        ``'result'`` when a mesh was produced and rejected.
    """

    def __init__(self, reason, metrics=None, stage='result'):
        super(DensifyRefused, self).__init__(reason)
        self.reason = reason
        self.metrics = metrics or {}
        self.stage = stage


def read_densities(coarse):
    """``{strip address: density}`` for a coarse mesh with strips collected.

    The inverse of :func:`apply_densities`. Use it to capture what a uniform
    pass produced before overriding individual strips.
    """
    out = {}
    for skey in coarse.strips():
        address = strip_address(coarse, skey)
        if address is None:
            continue
        out[address] = coarse.get_strip_density(skey)
    return out


def apply_densities(coarse, densities=None, target_length=None, density=None):
    """Set every strip's density. Strips must already be collected.

    A base pass first -- uniform ``density``, or ``target_length`` converted per
    strip -- and then ``densities`` overrides the strips it names. That order is
    deliberate: a map from before a layout edit will not name every strip, and a
    strip with no density at all raises ``KeyError`` from ``get_strip_density``
    several frames away inside the densification.

    Parameters
    ----------
    coarse : CoarsePseudoQuadMesh
    densities : dict, optional
        ``{strip address: density}``. Entries naming a strip that no longer
        exists are counted and returned, not raised.
    target_length : float, optional
        Base pass, converted per strip by ``set_strips_density_target``.
    density : int, optional
        Base pass, uniform. Takes precedence over ``target_length``.

    Returns
    -------
    dict
        ``applied`` -- how many strips took a value from ``densities``;
        ``unmatched`` -- addresses in ``densities`` that named no strip;
        ``base`` -- what the base pass was.
    """
    if density is not None:
        coarse.set_strips_density(density)
        base = {'density': density}
    elif target_length is not None:
        coarse.set_strips_density_target(target_length)
        base = {'target_length': target_length}
    else:
        raise ValueError('apply_densities needs a base: target_length or density')

    report = {'applied': 0, 'unmatched': [], 'base': base}
    if not densities:
        return report

    by_address = {}
    for skey in coarse.strips():
        address = strip_address(coarse, skey)
        if address is not None:
            by_address[address] = skey

    for address, value in densities.items():
        skey = by_address.get(address)
        if skey is None:
            report['unmatched'].append(address)
            continue
        coarse.set_strip_density(skey, int(value))
        report['applied'] += 1

    report['unmatched'].sort()
    return report


def densify_preserving(decomposition, densities=None, target_length=None,
                       density=None, coarse=None):
    """**Densify the layout, keeping per-strip densities. Never falls back.**

    The field route of :meth:`FieldDecomposition.quad_mesh`, with the density
    clobber and the triangulation fallback taken out. See the module docstring
    for why both had to go.

    Parameters
    ----------
    decomposition : FieldDecomposition
        Its ``field_aware`` flag is read by ``densify`` and is NOT touched here
        -- set it on the decomposition before calling if the caller wants the
        Coons interior instead.
    densities : dict, optional
        ``{strip address: density}``, applied over the base pass.
    target_length : float, optional
        Base pass. Defaults to the background spacing -- the same default
        ``quad_mesh`` uses. NOT the final element size when ``densities``
        overrides it.
    density : int, optional
        Uniform base pass. Takes precedence over ``target_length``.
    coarse : CoarsePseudoQuadMesh, optional
        The layout to densify. Defaults to ``decomposition.mesh``.

    Returns
    -------
    (QuadMesh, dict)
        The dense mesh, and a report: ``metrics`` (element quality),
        ``densities`` (what every strip ended up with, by address),
        ``density_report`` (:func:`apply_densities`' own), ``route``,
        ``coverage``, ``densify_stats`` and ``notes`` -- the repair notes this
        call added, and nothing from before it.

    Raises
    ------
    DensifyRefused
        The layout will not densify, or the result was rejected. **Nothing on
        the decomposition has been replaced** -- ``mesh`` and ``_edited`` are as
        they were, so the caller's layout survives a bad density.
    """
    if coarse is None:
        coarse = decomposition.mesh
    if coarse is None:
        raise ValueError(
            'no coarse layout -- call decomposition_mesh() or edit_coarse() first')

    ok, why = densifiable(coarse)
    if not ok:
        raise DensifyRefused(why, stage='layout')

    if density is None and target_length is None:
        target_length = decomposition.background.target_length

    mark = len(decomposition.repair_notes)

    coarse.collect_strips()
    density_report = apply_densities(coarse, densities=densities,
                                     target_length=target_length,
                                     density=density)

    try:
        dense = decomposition.densify(coarse)
    except Exception as exc:
        raise DensifyRefused(
            '{} during densification: {}'.format(type(exc).__name__, exc),
            stage='result')

    # ``_acceptable`` and ``_note_quality`` are private to FieldDecomposition,
    # and used here on purpose: they are the SAME structural gate and the same
    # quality note ``quad_mesh`` applies, and reimplementing either would let
    # this route drift from the one the rest of the pipeline is measured on.
    good, why, metrics = decomposition._acceptable(dense)
    if not good:
        raise DensifyRefused('densified mesh rejected: {}'.format(why),
                             metrics=metrics, stage='result')

    # Quality never reroutes -- see ``_acceptable`` -- but a degenerate element
    # is said out loud, by name and value.
    decomposition._note_quality(dense, metrics)
    decomposition.dense = dense

    return dense, {
        'metrics': metrics,
        'densities': read_densities(coarse),
        'density_report': density_report,
        'route': decomposition.route(),
        'coverage': decomposition._coverage(dense),
        'densify_stats': dict(decomposition.densify_stats or {}),
        'notes': list(decomposition.repair_notes[mark:]),
        'faces': dense.number_of_faces(),
    }
