"""Step 7b -- **global relaxation** of the dense mesh, after densification.

WHY THIS EXISTS AND WHAT IT IS NOT
----------------------------------

``densify.relax_patch`` solves ``alignment + stiffness * Laplacian`` **per patch, with
the four patch boundaries held fixed** -- the load-bearing decision that keeps patches
welding into a manifold and keeps opposite sides at matching densities. Two things are
structurally out of its reach because of it, and both were measured:

* **Seams.** The worst element in a mesh is routinely ON a patch boundary, which no
  stiffness setting can move: ``square+cable`` min 42.88 on the seams against 66.43
  inside the patches.
* **Anything Coons already got wrong.** ``_accepts`` clamps its floor to
  ``min(PATCH_MIN_ANGLE, ref_lo)``, so it can stop a patch being spent below its Coons
  reference but can never make one BETTER than it. Raising the floors 30 -> 45 -> 55
  moves ``square+ring cable`` not at all: min 13.73 at every setting.

WHAT IS CONSTRAINED, AND WHY NOT THE GUIDES
-------------------------------------------

Every constraint is chosen **topologically**, never by proximity, and that is measured
rather than stylistic. Pinning the interior vertices *nearest* a guide takes
``square+ring cable`` from 13.73/135.00 to **11.56/167.37**, worse than no post-pass at
all; constraining those same vertices ONTO the guide curve **folds faces** (max 180.00).
A radius around a curve collects a ragged band from both sides, and projecting that band
onto a line collapses it.

======================  ==================================================
 domain boundary         slides along the boundary curve, corners pinned
 patch seams             slide along their own ``edges_to_curves`` polyline
 singularities, poles    pinned
 everything else         free
 guides                  **nothing** -- see below
======================  ==================================================

Guides get no constraint. They are already in the field and in the separatrix routing,
and there is nothing well defined to constrain anyway: the field route never snaps a
polyedge onto a guide, so no chain of vertices *is* the guide -- ``square+ring cable``
carries **2** vertices within 0.05 of its ring in a 207-face mesh. What protects the
guide through this pass is that the layout carrying it is topological.

TWO GUARDS, BOTH LEARNED THE HARD WAY
-------------------------------------

**Order along the curve.** A vertex free to slide anywhere on its polyline can slide
past its neighbour, collapsing the edge between them. Unguarded, this pass took the ring
cable to **min 0.00 / max 180.00 / aspect inf** -- a hard-floor breach, far worse than
the mesh it started from. So every chain of vertices sharing a curve keeps its original
order along it, clamped by arc-length parameter with a minimum gap. This is the same
guard the coarse relaxation in ``guide_lines`` needed, for the same reason.

**Improvement only, unless there is something to buy.** On an already-good domain this
pass makes things worse: ``disc+round hole`` starts at 80.63/99.37 and a free relaxation
takes it to 74.26/108.00, because an area-based smoother redistributes boundary vertices
that were already right. So the result is accepted only if it does not worsen the
minimum angle, the maximum angle or the aspect ratio -- exactly the rule
``densify._accepts`` applies to an unguided patch. A domain with nothing to gain is left
untouched rather than nudged.
"""
from __future__ import annotations

from math import cos
from math import radians
from typing import Any
from typing import Sequence
from typing import TYPE_CHECKING

from compas.tolerance import TOL

from compas_singular.geometry.polyline import closest_on_polyline

if TYPE_CHECKING:
    from compas_singular.datastructures import Mesh


__all__ = ['relax_mesh', 'relaxation_constraints']

#: Smoothing schedule, gentlest first. The first setting whose result is accepted
#: wins, so a mesh takes the least relaxation that helps it and no more -- the same
#: shape of rule as ``densify.STIFFNESS``.
SCHEDULE = ((20, 0.3), (50, 0.5), (100, 0.5))

#: Minimum spacing between two vertices sliding on the same curve, as a fraction of
#: their original spacing. Zero would let an edge collapse; 1 would freeze the chain.
MIN_GAP = 0.25


# ------------------------------------------------------------------
# curves and parameters
# ------------------------------------------------------------------

class _Curve(object):
    """A polyline that reports arc-length parameters as well as closest points."""

    def __init__(self, points: list[list[float]], closed: bool = False) -> None:
        self.points = [list(p)[:3] for p in points]
        self.closed = closed
        self.cumulative = [0.0]
        for a, b in zip(self.points, self.points[1:]):
            self.cumulative.append(self.cumulative[-1] + _distance(a, b))
        self.length = self.cumulative[-1]

    def project(self, xyz: list[float]) -> tuple[list[float] | None, float]:
        """``(point, parameter)`` of the closest point, parameter as arc length.

        The arc length is an O(1) lookup rather than a second pass: the search
        reports WHICH segment won, and :attr:`cumulative` already holds the
        length up to it.
        """
        index, t, q, d = closest_on_polyline(xyz, self.points)
        if d == float('inf'):
            return None, 0.0
        a, b = self.points[index], self.points[index + 1]
        # not _distance: it squares with ``** 2``, which is not always exactly
        # ``x * x`` here, and the cumulative table this adds to was built the
        # other way round
        abx, aby, abz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        span = (abx * abx + aby * aby + abz * abz) ** 0.5
        return q, self.cumulative[index] + t * span

    def point_at(self, parameter: float) -> list[float]:
        """The point at an arc-length parameter, clamped to the curve."""
        parameter = max(0.0, min(self.length, parameter))
        for i, (lo, hi) in enumerate(zip(self.cumulative, self.cumulative[1:])):
            if parameter <= hi:
                span = hi - lo
                t = 0.0 if span <= 0.0 else (parameter - lo) / span
                a, b = self.points[i], self.points[i + 1]
                return [a[k] + (b[k] - a[k]) * t for k in range(3)]
        return list(self.points[-1])


def _distance(a: list[float], b: list[float]) -> float:
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5


# ------------------------------------------------------------------
# constraints
# ------------------------------------------------------------------

def _singularities(mesh: Mesh) -> set[int]:
    """Interior vertices whose valency is not 4, plus every pole.

    A pole is a collapsed side, so it is a singularity a valency count cannot see --
    ``face_pole`` is where it is recorded.

    **Off by default**, and the reasoning that said otherwise was wrong. Pinning these
    looks principled -- they are the layout's irregular nodes -- but moving one does
    not change its valency, so it is not a topology change at all, just a vertex
    finding a better position. Measured on ``square+ring cable``, pinning them is what
    stopped this pass working: singularities+poles pinned, rejected; poles only,
    rejected; nothing pinned, **accepted at min 13.73 -> 17.63, max 135.00 -> 129.78,
    aspect 4.69 -> 3.32**. Which is the same mistake as pinning the vertices nearest a
    guide, wearing better clothes.
    """
    out = set()
    for vertex in mesh.vertices():
        if mesh.is_vertex_on_boundary(vertex):
            continue
        if len(mesh.vertex_neighbors(vertex)) != 4:
            out.add(vertex)
    out.update(mesh.attributes.get('face_pole', {}).values())
    return out


def _boundary_curves(mesh: Mesh) -> list[tuple[list[int], _Curve]]:
    """The mesh's own boundary loops, as closed :class:`_Curve` objects."""
    curves = []
    for loop in mesh.vertices_on_boundaries():
        keys = list(loop)
        while len(keys) > 1 and keys[-1] == keys[0]:
            keys.pop()
        if len(keys) < 3:
            continue
        points = [mesh.vertex_coordinates(v) for v in keys]
        curves.append((keys, _Curve(points + points[:1], closed=True)))
    return curves


def relaxation_constraints(
    mesh: Mesh,
    seam_edge: dict[str, tuple[int, int, int]] | None = None,
    edges_to_curves: dict[tuple[int, int], list[list[float]]] | None = None,
    corner_angle: float = 30.0,
    seams: str = 'free',
    pin_singularities: bool = False,
) -> tuple[dict[_Curve, list[int]], set[int]]:
    """Work out what each vertex may do.

    Returns
    -------
    (dict, set)
        ``chains`` maps a :class:`_Curve` to the ordered list of vertices sliding on
        it, and ``pinned`` is the set of vertices that may not move at all.
    """
    pinned = set(_singularities(mesh)) if pin_singularities else set()
    chains = {}

    # -- the domain outline: slide along it, pin the corners ------------------
    limit = cos(radians(180.0 - corner_angle))
    for keys, curve in _boundary_curves(mesh):
        count = len(keys)
        for i, vertex in enumerate(keys):
            a = mesh.vertex_coordinates(keys[i - 1])
            b = mesh.vertex_coordinates(vertex)
            c = mesh.vertex_coordinates(keys[(i + 1) % count])
            u = [b[k] - a[k] for k in range(3)]
            v = [c[k] - b[k] for k in range(3)]
            lu, lv = _distance(a, b), _distance(b, c)
            if lu <= 0.0 or lv <= 0.0:
                continue
            cosine = sum(u[k] * v[k] for k in range(3)) / (lu * lv)
            if cosine < limit:                       # a kink: this is a domain corner
                pinned.add(vertex)
        chains[curve] = [v for v in keys if v not in pinned]

    # -- patch seams ----------------------------------------------------------
    if seams == 'free':
        return chains, pinned

    if seam_edge and edges_to_curves:
        by_edge = {}
        for vertex in mesh.vertices():
            if mesh.is_vertex_on_boundary(vertex):
                continue                              # the outline already claimed it
            found = seam_edge.get(TOL.geometric_key(mesh.vertex_coordinates(vertex)))
            if found is None:
                continue
            u, v, index = found
            by_edge.setdefault((u, v), []).append((index, vertex))

        for edge, members in by_edge.items():
            points = edges_to_curves.get(edge) or edges_to_curves.get((edge[1], edge[0]))
            if not points or len(points) < 2:
                continue
            members.sort()
            ordered = [vertex for _, vertex in members if vertex not in pinned]
            if not ordered:
                continue
            if seams == 'fixed':
                pinned.update(ordered)
            else:
                chains[_Curve(points)] = ordered

    return chains, pinned


# ------------------------------------------------------------------
# the pass
# ------------------------------------------------------------------

def _quality(mesh: Mesh) -> tuple[float, float, float]:
    """``(min angle, max angle, max aspect)`` -- the three the gate compares."""
    from compas_singular.framefield.quality import mesh_quality
    q = mesh_quality(mesh)
    return q['min_angle'], q['max_angle'], q['aspect_max']


def _positions(mesh: Mesh) -> dict[int, list[float]]:
    return {v: mesh.vertex_coordinates(v) for v in mesh.vertices()}


def _restore(mesh: Mesh, positions: dict[int, list[float]]) -> None:
    for vertex, xyz in positions.items():
        mesh.vertex_attributes(vertex, 'xyz', xyz)


def _reproject(mesh: Mesh, chains: dict[_Curve, list[int]]) -> None:
    """Put every sliding vertex back on its curve, keeping its order along it."""
    for curve, ordered in chains.items():
        if not ordered:
            continue
        wanted = []
        for vertex in ordered:
            _, parameter = curve.project(mesh.vertex_coordinates(vertex))
            wanted.append(parameter)

        # the spacing the chain started with, as the floor each gap may not go below
        floor = MIN_GAP * (curve.length / (len(ordered) + 1))

        # one forward sweep then one backward sweep leaves the sequence increasing
        # and inside the curve, which is what stops a vertex passing its neighbour
        for i in range(1, len(wanted)):
            wanted[i] = max(wanted[i], wanted[i - 1] + floor)
        upper = curve.length if not curve.closed else curve.length - floor
        for i in range(len(wanted) - 1, -1, -1):
            cap = upper if i == len(wanted) - 1 else wanted[i + 1] - floor
            wanted[i] = min(wanted[i], cap)

        for vertex, parameter in zip(ordered, wanted):
            mesh.vertex_attributes(vertex, 'xyz', curve.point_at(parameter))


def _smooth(
    mesh: Mesh,
    chains: dict[_Curve, list[int]],
    pinned: set[int],
    kmax: int,
    damping: float,
) -> None:
    """Area-weighted smoothing, with the sliding vertices reprojected every round."""
    fixed = set(pinned)
    for k in range(kmax):
        centroid = {}
        for vertex in mesh.vertices():
            if vertex in fixed:
                continue
            areas, points = [], []
            for face in mesh.vertex_faces(vertex, ordered=True):
                if face is None:
                    continue
                areas.append(mesh.face_area(face))
                points.append(mesh.face_centroid(face))
            if not areas or not sum(areas):
                continue
            total = sum(areas)
            centroid[vertex] = [
                sum(a * p[i] for a, p in zip(areas, points)) / total for i in range(3)]

        for vertex, target in centroid.items():
            xyz = mesh.vertex_coordinates(vertex)
            mesh.vertex_attributes(vertex, 'xyz', [
                xyz[i] + damping * (target[i] - xyz[i]) for i in range(3)])

        _reproject(mesh, chains)


def relax_mesh(
    mesh: Mesh,
    seam_edge: dict[str, tuple[int, int, int]] | None = None,
    edges_to_curves: dict[tuple[int, int], list[list[float]]] | None = None,
    corner_angle: float = 30.0,
    schedule: Sequence[tuple[int, float]] = SCHEDULE,
    seams: str = 'free',
    pin_singularities: bool = False,
) -> dict[str, Any]:
    """**Relax a dense mesh globally**, seams and outline sliding, singularities pinned.

    Modifies ``mesh`` in place, and only if the result is an improvement -- see the
    module docstring on why a domain with nothing to gain is left alone.

    Parameters
    ----------
    mesh : Mesh
        The dense mesh.
    seam_edge : dict, optional
        ``{geometric_key: (u, v, index)}`` from ``field_densification``'s stats, i.e.
        ``FieldDecomposition.densify_stats['seam_edge']``. Without it the patch seams
        are not recognised and only the outline is constrained.
    edges_to_curves : dict, optional
        ``{(u, v): [point, ...]}`` from :meth:`FieldDecomposition.edges_to_curves`.
    corner_angle : float, optional
        Kink angle, in degrees, past which a boundary vertex is a domain corner.
    schedule : sequence[(int, float)], optional
        ``(kmax, damping)`` settings, gentlest first.
    seams : {'free', 'slide', 'fixed'}, optional
        What the patch-seam vertices may do. ``'free'`` -- the default, and the only
        one measured to help -- lets them leave the separatrix. ``'slide'`` keeps them
        on their own curve, which never hurts but never helps either: sliding cannot
        repair a 13.73 degree angle that the seam's own shape causes. ``'fixed'``
        reproduces today's behaviour, patch interiors only.
    pin_singularities : bool, optional
        Hold the irregular vertices and poles in place. **Off by default** -- see
        :func:`_singularities` for the measurement that says pinning them is what
        stops this pass working.

    Returns
    -------
    dict
        ``accepted`` -- the setting used, or None if the mesh was left untouched.
        ``before`` / ``after`` -- ``(min angle, max angle, max aspect)``.
        ``sliding`` / ``pinned`` -- how many vertices were in each category.
    """
    chains, pinned = relaxation_constraints(
        mesh, seam_edge=seam_edge, edges_to_curves=edges_to_curves,
        corner_angle=corner_angle, seams=seams, pin_singularities=pin_singularities)

    before = _quality(mesh)
    original = _positions(mesh)
    sliding = sum(len(v) for v in chains.values())

    report = {'accepted': None, 'before': before, 'after': before,
              'sliding': sliding, 'pinned': len(pinned)}

    for kmax, damping in schedule:
        _restore(mesh, original)
        _smooth(mesh, chains, pinned, kmax, damping)
        after = _quality(mesh)
        # improvement only: no worse on any of the three
        if after[0] >= before[0] and after[1] <= before[1] and after[2] <= before[2]:
            report['accepted'] = (kmax, damping)
            report['after'] = after
            return report

    _restore(mesh, original)
    return report
