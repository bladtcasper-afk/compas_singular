"""Step 6b -- ``FieldDecomposition``, the drop-in for ``SkeletonDecomposition``.

Deliberately mirrors that class's method names, so an existing script changes by
two lines:

    decomposition = SkeletonDecomposition.from_mesh(trimesh)
    coarsemesh = decomposition.decomposition_mesh(point_features)

becomes

    decomposition = FieldDecomposition.from_boundary(outer, guides=cables)
    coarsemesh = decomposition.decomposition_mesh()

Everything after that -- ``collect_strips``, ``set_strips_density_target``,
``densification``, ``add_strip``, the guide_lines helpers -- is untouched.
"""
from compas.geometry import distance_point_point
from compas.itertools import pairwise
from compas_singular.datastructures import CoarsePseudoQuadMesh
from compas_singular.geometry.polyline import distance_to_polyline
from compas_singular.geometry.polyline import project_on_polyline

from compas_singular.framefield.arrangement import faces_from_arrangement
from compas_singular.framefield.arrangement import faces_with_repeated_vertices
from compas_singular.framefield.arrangement import planar_arrangement
from compas_singular.framefield.background import BackgroundMesh
from compas_singular.framefield.constraints import as_curve_list
from compas_singular.framefield.constraints import from_boundary
from compas_singular.framefield.constraints import from_curves
from compas_singular.framefield.densify import field_densification
# These moved to ``editing`` -- they weld, snap and repair a layout from plain
# geometry and never needed the field. The import direction is one-way:
# ``framefield`` may use ``editing``, never the reverse.
from compas_singular.editing.rebuild import coarse_from_skeleton
from compas_singular.editing.rebuild import mesh_from_faces
from compas_singular.editing.rebuild import warp_polyline
from compas_singular.framefield.field import CrossField
from compas_singular.framefield.quality import hard_floor
from compas_singular.framefield.quality import mesh_quality
from compas_singular.framefield.separatrix_network import SHARP_TURN
from compas_singular.framefield.separatrix_network import build_network
from compas_singular.editing.repair import densifiable
from compas_singular.editing.repair import solve_non_quad_faces
from compas_singular.editing.repair import topological_quad_split
from compas_singular.framefield.trace import Tracer


__all__ = ['FieldDecomposition']


def _on_loop(point, loop, tol=1e-6):
    """Whether a point sits on a boundary loop (used to pick the outer arcs)."""
    from compas.geometry import distance_point_point
    from compas.itertools import pairwise
    for a, b in pairwise(list(loop) + list(loop[:1])):
        abx, aby = b[0] - a[0], b[1] - a[1]
        length2 = abx * abx + aby * aby
        if length2 == 0.0:
            continue
        t = ((point[0] - a[0]) * abx + (point[1] - a[1]) * aby) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + abx * t, a[1] + aby * t, 0.0]
        if distance_point_point(point, q) < tol:
            return True
    return False


def _loop_chord_tolerance(loop):
    """How far off a loop the CURVE it samples can lie. ``1e-6`` if it is a
    polygon.

    A round hole reaches the pipeline as a polygon -- 48 points, say -- and its
    loop is the chords between them, which run INSIDE the circle. So a corner
    that is genuinely on the circle is off the loop by a sagitta, and
    :func:`_on_loop`'s 1e-6 says it is not on the wall at all. That matters as
    soon as a corner comes from anywhere but the loop's own point list: a Rhino
    edit snaps it to the document CURVE, and a coarse mesh that round-tripped
    through ``rs.AddMesh`` is single precision besides. Measured on a radius-0.8
    hole sampled at 48 points, a corner slid 40 degrees along the true circle:
    ``_on_loop`` False, both its ring edges chorded, the hole 0.150 inside its
    own circle -- 19% of the radius -- while every other hole stayed round.

    The tolerance is taken from the loop itself rather than from the domain,
    because that is where the answer is: the distance from each vertex to the
    chord joining its neighbours IS the local sampling error, and it needs no
    knowledge of the radius. Vertices turning more than ``repair.SHARP_TURN``
    are skipped -- **a 90-degree corner of a plate is not a badly sampled
    curve**, and measuring it would hand a square a tolerance the size of its
    own patches. That is the same 45-degree question, and the same constant,
    that :func:`repair.boundary_corners` settles for boundary corners.

    Halved because that vertex-to-chord distance spans TWO segments, and what a
    point on the curve is bounded by is ONE segment's sagitta -- a quarter of it
    on a circle. So this still carries a factor of two in hand, and no more:
    measured on a radius-0.7 circle sampled at 9 points, the loosest sampling
    this workflow produces, 0.082 against a true sagitta of 0.042.

    A polygon keeps 1e-6 exactly, and only a genuinely sampled curve widens.
    """
    from math import acos

    pts = list(loop)
    n = len(pts)
    if n < 3:
        return 1e-6

    worst = 0.0
    for i in range(n):
        a, b, c = pts[i - 1], pts[i], pts[(i + 1) % n]
        ux, uy = b[0] - a[0], b[1] - a[1]
        wx, wy = c[0] - b[0], c[1] - b[1]
        lu = (ux * ux + uy * uy) ** 0.5
        lw = (wx * wx + wy * wy) ** 0.5
        if lu < 1e-12 or lw < 1e-12:
            continue
        cosv = max(-1.0, min(1.0, (ux * wx + uy * wy) / (lu * lw)))
        if acos(cosv) > SHARP_TURN:
            continue
        worst = max(worst, distance_to_polyline(b, [a, c]))
    return max(1e-6, 0.5 * worst)


def _arcs_between(loop, pa, pb):
    """BOTH ways round a closed loop between two points ON it, shorter first.

    Two, not one, because **length does not decide which way round is the
    edge**. The usual case says it does -- a coarse edge spans a fraction of the
    boundary, so the shorter arc is the edge -- and that holds for the outer
    wall of a plate, where one edge is a small part of a long loop. It fails on
    a HOLE, whose loop is short and whose ring may be as few as three corners:
    an edge spanning more than half the circle then has its shorter arc running
    the wrong way, through the ring's other corners.

    Measured on a radius-0.7 hole with ring corners at 0, 200 and 280 degrees:
    the 0 -> 200 edge was handed the 160-degree complement,
    :meth:`FieldDecomposition._arc_is_one_edge` rejected it -- correctly, it runs
    through the corner at 280 -- and the edge fell to a chord while its two
    neighbours came out round. Several identical circles in one domain then mesh
    round or faceted depending on where the layout happened to put their
    corners, which is the symptom to recognise. Nothing is wrong with the
    geometry: the right arc was the one never offered.

    Returning both and letting the caller keep the first that no other corner
    lies on puts the decision where the answer actually is. Shorter stays first,
    so any edge that already resolved resolves the same way.

    Both ways are ``None`` for a degenerate edge -- the two ends projecting to
    the same place on the loop -- which is what :func:`_arc_between` already
    returned there.
    """
    from compas.geometry import distance_point_point
    from compas.itertools import pairwise

    pts = list(loop)
    ring = pts + pts[:1]
    cum = [0.0]
    for a, b in pairwise(ring):
        cum.append(cum[-1] + distance_point_point(a, b))
    total = cum[-1]
    if total == 0.0:
        return []

    def param(point):
        best_d, best_s = float('inf'), 0.0
        for i, (a, b) in enumerate(pairwise(ring)):
            abx, aby = b[0] - a[0], b[1] - a[1]
            length2 = abx * abx + aby * aby
            if length2 == 0.0:
                continue
            t = ((point[0] - a[0]) * abx + (point[1] - a[1]) * aby) / length2
            t = max(0.0, min(1.0, t))
            q = [a[0] + abx * t, a[1] + aby * t, 0.0]
            d = distance_point_point(point, q)
            if d < best_d:
                best_d, best_s = d, cum[i] + t * (cum[i + 1] - cum[i])
        return best_s

    def point_at(s):
        s = s % total
        for i, (c0, c1) in enumerate(zip(cum, cum[1:])):
            if c0 - 1e-9 <= s <= c1 + 1e-9:
                a, b = ring[i], ring[i + 1]
                length = c1 - c0
                t = 0.0 if length == 0 else (s - c0) / length
                return [a[k] + (b[k] - a[k]) * t for k in range(3)]
        return list(ring[-1])

    def span_arc(s0, span, reverse):
        arc = [list(pa if not reverse else pb)]
        for i in range(1, len(pts) + 1):
            arc.append(point_at(s0 + span * i / float(len(pts) + 1)))
        arc.append(list(pb if not reverse else pa))
        if reverse:
            arc.reverse()
        arc[0], arc[-1] = list(pa), list(pb)
        return arc

    sa, sb = param(pa), param(pb)
    forward = (sb - sa) % total
    backward = total - forward
    if forward <= 1e-12 or backward <= 1e-12:
        return []

    ahead = (sa, forward, False)
    behind = (sb, backward, True)
    order = (ahead, behind) if forward <= backward else (behind, ahead)
    return [span_arc(*args) for args in order]


def _arc_between(loop, pa, pb):
    """The shorter arc of a closed loop between two points ON it, pa -> pb."""
    found = _arcs_between(loop, pa, pb)
    return found[0] if found else None


class FieldDecomposition(object):
    """Coarse quad layout from a cross field's separatrices.

    Attributes
    ----------
    background : BackgroundMesh
    field : CrossField
    separatrices : list[Separatrix]
    mesh : CoarsePseudoQuadMesh or None
        Set by :meth:`decomposition_mesh`.
    """

    def __init__(self, background, field, separatrices, tracer):
        self.background = background
        self.field = field
        #: Every argument :meth:`from_boundary` was called with, so the same route
        #: can be rebuilt on another domain -- a symmetric unit, for one.
        self.inputs = {}
        #: What :meth:`find_symmetry` found, or ``None``. Deliberately not
        #: :attr:`symmetry`, which is the in-solver group of ``framefield.symmetry``.
        self.symmetry_report = None
        self.separatrices = separatrices
        self.tracer = tracer
        self.mesh = None
        #: The last mesh :meth:`quad_mesh` returned. Kept separately because on
        #: the field route ``self.mesh`` is the COARSE layout -- ``quad_mesh``
        #: returns a densification of it without overwriting it -- so
        #: ``self.mesh`` is the wrong thing to measure element quality on.
        self.dense = None
        self.polylines = None
        #: The guide curves the field was solved with, if any. Kept because
        #: densification asks: element quality is spent to satisfy a guide and
        #: never otherwise -- see ``densify._accepts``.
        self.guides = []
        self.trace_report = {}
        self.repair_notes = []
        #: The symmetry group of the input, or ``None``. Set by
        #: :meth:`from_boundary`; read by :meth:`_build`, which has to pass it to
        #: ``build_network`` -- the layout stage makes DISCRETE choices with
        #: tolerance tie-breaks, so it is the one place downstream of the field
        #: that cannot be made symmetric by giving it symmetric input alone.
        self.symmetry = None
        #: What ``symmetry.snap_singularities`` did, including its Poincare-Hopf
        #: re-check. Empty when no symmetry was applied.
        self.snap_report = {}
        self._network = None
        #: Whether patch INTERIORS are integrated from the field (``densify``)
        #: or left to ``discrete_coons_patch``. Only ever set to ``False`` to
        #: measure the difference -- see ``16_field_densify.py``.
        self.field_aware = True
        #: What :meth:`densify` did, per patch class. Set by every
        #: :meth:`quad_mesh` call that reaches the field route.
        self.densify_stats = {}
        #: Whether :attr:`mesh` is a HAND-EDITED layout rather than a generated
        #: one. Changes two things and nothing else: :meth:`quad_mesh` stops
        #: regenerating the layout, and the coverage check stops rerouting on it
        #: -- a user who deleted a patch changed the layout on purpose.
        self._edited = False
        #: The ``poles`` :attr:`mesh` was built with, as a rounded tuple. The
        #: cache key of :meth:`decomposition_mesh`: the layout is a pure
        #: function of the (cached) network and these, so asking twice must give
        #: the same OBJECT back, not an equal one -- see that method.
        self._mesh_poles = None
        #: What :meth:`edit_coarse` and :meth:`edges_to_curves` made of the
        #: edit. Empty until :meth:`edit_coarse` is called.
        self.edit_notes = {}
        #: Curves the USER drew on the layout -- a line added in Rhino that cut
        #: a patch in two. They belong to the layout and not to the traced
        #: network, so ``_build`` neither produces them nor caches them, and
        #: :meth:`edit_coarse` has to put them back into :attr:`polylines` on
        #: every call or the added edge densifies as a straight chord. Set by
        #: ``rhino.edit_coarse.CoarseLayoutObject.commit``; see
        #: :meth:`edges_to_curves`, whose exact-match branch is what finds them.
        self.user_curves = []

    def warnings(self):
        """Anything that will silently degrade the layout.

        ``arm_mismatch`` is the one that matters. If the probe circles around a
        singularity never produce exactly ``4 - index`` launch directions, the
        tracer falls back to the best partial set and the layout comes out with
        the wrong number of patches -- a pentagon at one background resolution
        gives 5 clean quads, and at a slightly finer one gives 4 arms, one
        5-sided patch, and a densification failure. Nothing about the coarse mesh
        says which happened, so check this.

        ``QUALITY FAILURE`` is the other one to look for. :meth:`quad_mesh`
        always returns a mesh, so a degenerate element -- an angle at or near 0
        or 180 degrees, a zero-length edge -- reaches the caller as a warning
        rather than as an exception. It names the metric and its value. See
        :meth:`quality` for the numbers behind it.
        """
        out = []
        mismatch = self.trace_report.get('arm_mismatch') or []
        for fkey, got, expected in mismatch:
            out.append('singularity {}: found {} launch directions, expected {}'.format(
                fkey, got, expected))
        out.extend(self.repair_notes)
        for reason in ('cycle', 'lost', 'length'):
            n = (self.trace_report.get('reasons') or {}).get(reason)
            if n:
                out.append('{} separatrix/ices ended with reason {!r}'.format(n, reason))
        return out

    # ------------------------------------------------------------------
    # constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_boundary(cls, outer_boundary, inner_boundaries=None, guides=None,
                      mode='perpendicular', target_length=None, orthogonal=None,
                      guide_weight=1.0, guide_band=None, relax=False, field_tau=None,
                      symmetry='auto', solve=True):
        """Solve the field and trace its separatrices.

        Parameters
        ----------
        outer_boundary : list[[x, y, z]]
        inner_boundaries : list[list[[x, y, z]]], optional
        guides : list[list[[x, y, z]]], optional
            Guide curves -- cables, force lines. A LIST of curves, like
            ``inner_boundaries`` and unlike ``outer_boundary``; wrap a lone
            guide in a list. Each curve may be a list of points or a compas
            ``Polyline``, and its points lists, tuples or compas ``Point``s.
        mode : {'perpendicular', 'tangent'}, optional
        target_length : float, optional
            Background triangulation spacing. NOT the final quad size -- that is
            set later by ``set_strips_density_target`` on the coarse mesh.
        orthogonal : bool, optional
            ``None`` (default) means auto. Currently always solves a cross field;
            the frame-field fallback arrives with ``warp.py`` (milestone 3), and
            this parameter exists so callers do not have to change then.
        relax : bool, optional
            Solve the field by diffusion + normalisation rather than by a single
            Dirichlet solve. Off by default. See ``field.py``'s module docstring:
            the default solver is measurably wrong about singularity PLACEMENT
            (disc radius 0.45 against a published 0.85), and turning this on
            moves the layout on every curved domain. ``19_field_accuracy.py`` is
            the check.
        field_tau : float, optional
            Diffusion time per relaxation step. ``relax`` only.
        symmetry : Symmetry or 'auto' or None, optional
            ``'auto'`` (the default) detects the symmetry group of the input --
            outer boundary, holes and guides together -- and keeps it all the way
            to the layout. ``None`` reproduces the behaviour from before
            2026-08-14 on any domain.

            **A symmetric problem does not otherwise give a symmetric layout.**
            The solve is innocent: on a symmetric mesh it is exact to 1e-15,
            being a convex problem with a unique minimiser. But the background
            triangulation is asymmetric by construction, and the layout stage
            then amplifies sub-triangle drift into whole extra patches, because
            it snaps at ``0.8 * target_length`` and a corner either has a twin
            inside that tolerance or does not. Measured on a square with four
            arcs that are exact mirror images: 18 patches with 2 pole triangles
            and 15% of corners having a twin, against 16 clean quads at 100%
            under all eight symmetries, deviation 0.0000. ``symmetry.py``
            documents all five places it had to be told, and
            ``24_symmetry.py`` measures them one at a time.

            An ASYMMETRIC domain is unaffected: the detected group is trivial,
            every step switches itself off, and the layout is byte-identical to
            what it was.
        solve : bool, optional
            ``False`` stores the inputs and returns without solving or tracing.
            Such a decomposition supports :meth:`find_symmetry`,
            :meth:`symmetry_unit` and :meth:`solve`; anything that needs the field
            raises until :meth:`solve` is called.
        """
        inputs = {'outer_boundary': outer_boundary, 'inner_boundaries': inner_boundaries,
                  'guides': guides, 'mode': mode, 'target_length': target_length,
                  'orthogonal': orthogonal, 'guide_weight': guide_weight,
                  'guide_band': guide_band, 'relax': relax, 'field_tau': field_tau,
                  'symmetry': symmetry}
        if not solve:
            # Inputs only. Enough for find_symmetry() and symmetry_unit(), which
            # mesh a unit of the domain and never need the whole-domain field.
            out = cls(None, None, [], None)
            out.inputs = inputs
            return out

        if orthogonal is False:
            raise NotImplementedError(
                'non-orthogonal frame fields need warp.py -- milestone 3. '
                'Cases A and B (orthogonal stress field, single cable family) '
                'are handled by the cross field; see the evaluation doc.')

        # Background, constraints, solve and symmetrisation are the FIELD's own
        # business -- nothing in them is specific to decomposing. They live on
        # ``CrossField.from_boundary``, which a caller can use directly to get a
        # field with no separatrix tracing behind it, and then densify any layout
        # with it. Everything below this line is the tracing half.
        field = CrossField.from_boundary(
            outer_boundary, inner_boundaries=inner_boundaries, guides=guides,
            mode=mode, target_length=target_length, guide_weight=guide_weight,
            guide_band=guide_band, relax=relax, tau=field_tau, symmetry=symmetry)

        tracer = Tracer(field, singularity_points=field.singularity_points)
        separatrices, trace_report = tracer.separatrices()
        out = cls(field.background, field, separatrices, tracer)
        out.guides = list(field.guides)
        out.trace_report = trace_report
        out.symmetry = field.symmetry
        out.snap_report = field.snap_report
        out.inputs = inputs
        return out

    @property
    def is_solved(self):
        return self.field is not None

    def _require_solved(self):
        if self.field is None:
            raise ValueError('this FieldDecomposition was built with solve=False and has no field '
                             'yet -- call solve() first, or use symmetry_unit()')

    def solve(self):
        """**The solved decomposition** for the stored inputs. Returns a NEW object.

        For a decomposition built with ``solve=False``. A solved one returns itself.
        """
        if self.field is not None:
            return self
        return type(self).from_boundary(**self.inputs)

    # ------------------------------------------------------------------
    # symmetry
    # ------------------------------------------------------------------

    def find_symmetry(self, tol=None, include=('walls', 'holes', 'guides', 'poles'), max_order=12):
        """**Detect the symmetry of the domain this decomposition was built from.**

        Walls, holes and guides are tested. Works on an unsolved decomposition
        (``solve=False``), which is the cheap way in: the unit is meshed on its
        own, so a whole-domain solve would be thrown away. Kept on
        :attr:`symmetry_report`.

        Returns
        -------
        :class:`compas_singular.symmetry.SymmetryReport`
        """
        from compas_singular.symmetry import find_symmetry
        from compas_singular.symmetry.routes import domain_of
        self.symmetry_report = find_symmetry(domain=domain_of(self), tol=tol, include=include,
                                             max_order=max_order)
        return self.symmetry_report

    def symmetry_unit(self, keys=None, centre='route', seam=None, report=None):
        """**Solve and trace one symmetric unit of the domain.**

        The unit is cut out along the seams of the symmetries in ``keys`` and
        decomposed with the settings this decomposition was built with, and its
        own field -- which ``quad_mesh()`` then uses for the patch interiors.
        ``symmetry`` is not applied inside the unit: the unit is already one
        fundamental region, and the symmetry comes from expanding it.

        Parameters are those of ``SkeletonDecomposition.symmetry_unit``.

        Returns
        -------
        :class:`compas_singular.symmetry.SymmetricUnit`
        """
        from compas_singular.symmetry import build_unit
        from compas_singular.symmetry.routes import mesher_for
        report = report or self.symmetry_report or self.find_symmetry()
        return build_unit(report, mesher_for(self), keys=keys, centre=centre, seam=seam)

    @classmethod
    def from_mesh(cls, trimesh, **kwargs):
        """Construct from a triangulation, mirroring ``SkeletonDecomposition.from_mesh``.

        Only the BOUNDARY of ``trimesh`` is used. compas_singular's
        ``boundary_triangulation`` puts no vertices in the interior at all, so
        its triangles cannot carry a field -- see ``background.py``. The
        triangulation is rebuilt at a spacing a field can live on.
        """
        loops = trimesh.boundaries()
        loops = sorted(loops, key=len, reverse=True)
        outer = [trimesh.vertex_coordinates(v) for v in loops[0]]
        inners = [[trimesh.vertex_coordinates(v) for v in loop] for loop in loops[1:]]
        return cls.from_boundary(outer, inners, **kwargs)

    # ------------------------------------------------------------------
    # the SkeletonDecomposition-shaped API
    # ------------------------------------------------------------------

    def decomposition_polylines(self):
        """All the polylines forming the decomposition -- separatrices and walls."""
        boundary, others, _ = self._build()
        self.polylines = boundary + others
        return self.polylines

    def decomposition_mesh(self, poles=(), force=False):
        """The coarse quad mesh. **The same object every time you ask.**

        The layout is a pure function of the polyline network -- which
        :meth:`_build` caches -- and ``poles``. It used to be rebuilt from
        scratch on every call anyway, which made this the one method in the
        pipeline you could not call twice::

            coarse = d.decomposition_mesh()
            coarse.collect_strips()
            coarse.set_strip_density(skey, 9)
            d.quad_mesh()          # calls decomposition_mesh() again ->
                                   # a DIFFERENT mesh; the 9 is on the old one

        Nothing errored and the layout looked identical, so the only symptom was
        a density that had no effect. It now returns :attr:`mesh` unchanged
        whenever the network and ``poles`` are the ones it was built from, so
        holding a layout across calls is safe.

        A HAND-EDITED layout is not cached: an explicit call regenerates and
        discards the edit, which is the behaviour :meth:`edit_coarse` documents
        and :meth:`quad_mesh` relies on. Pass ``force=True`` to regenerate an
        unedited one.

        Parameters
        ----------
        poles : list[[x, y, z]], optional
            Preferred pole positions, as in
            ``SkeletonDecomposition.decomposition_mesh``. A field singularity is
            not a pole and none is needed for the domains in ``14_arrangement``,
            but :func:`repair.solve_non_quad_faces` resolves any triangular
            patch it cannot otherwise repair as a pseudo-quad pole, and a
            position given here decides which of that triangle's corners the
            collapsed side sits at. Different poles miss the cache.
        force : bool, optional
            Rebuild even on a cache hit. For a caller that wants a clean layout
            back after mutating the one it was given -- densities, strips, moved
            corners -- since the cached object is the one it mutated.

        Returns
        -------
        CoarsePseudoQuadMesh
        """
        self._require_solved()
        poles_key = tuple(tuple(round(float(c), 6) for c in point)
                          for point in (poles or ()))
        if (not force and self.mesh is not None and not self._edited
                and self._mesh_poles == poles_key):
            return self.mesh

        boundary, others, _ = self._build()
        self.polylines = boundary + others
        self.repair_notes = []

        if not others and self.background.inners:
            # No interior line AND a hole: there is no single face here to
            # build. The arcs handed in belong to two different loops, and
            # threading them into one ring produces a patch that runs straight
            # across the hole -- measured on the round-holed disc, a two-quad
            # "layout" covering 66% of the domain. That used to be caught by
            # accident, because the 6-gon it made quad-split into inverted
            # faces; it is not something to rely on.
            mesh = self._triangulation_fallback(
                'no interior separatrices and {} boundary loops -- a domain with '
                'a hole has no single-patch layout'.format(1 + len(self.background.inners)))
        elif not others:
            # from_polylines discards any face whose vertices are ALL on the
            # boundary, so a domain with no interior line at all comes back with
            # zero faces. A square is exactly that case: no singularities, no
            # reflex corners, and a perfectly good one-quad layout. Build the
            # single face directly from the boundary arcs' endpoints -- they are
            # already in loop order.
            mesh = self._single_patch(boundary)
        elif self._cut_holes():
            # A hole opened by ``Tracer.hole_launches``. Its cuts run boundary to
            # boundary, so EVERY vertex of every resulting patch lies on a
            # boundary loop -- and that is exactly the face ``from_polylines``
            # throws away (it keeps a face only ``if len(notonboundary)``). It is
            # the same rule that makes a square come back empty above, but here
            # it discards the real patches and leaves artifacts: measured on the
            # round-holed disc, two slivers covering 16% of the domain.
            #
            # So walk the arrangement and recover the faces directly. Restricted
            # to this case on purpose: every other domain has an interior line
            # and is served correctly by ``from_polylines`` today.
            faces = faces_from_arrangement(boundary, others, self._loops())
            mesh = self._mesh_from_faces(faces)
            if mesh is None:
                mesh = self._triangulation_fallback(
                    'a hole was cut but the arrangement yielded no face')
        else:
            mesh = CoarsePseudoQuadMesh.from_polylines(boundary, others)

            # A face that visits one vertex twice is the signature of a network
            # that was not a planar arrangement -- the face walk went through an
            # unresolved crossing, or up a dangling arm and back. ``arrangement``
            # is supposed to have made that impossible; if one appears anyway the
            # layout is already wrong here, several frames before ``densifiable``
            # reports it as an inverted quad with no clue where it came from.
            repeated = faces_with_repeated_vertices(mesh)
            if repeated:
                self.repair_notes.append(
                    'ARRANGEMENT WARNING: {} recovered face(s) repeat a vertex '
                    '({}). The polyline network is not a planar arrangement; the '
                    'layout that follows is unreliable.'.format(
                        len(repeated), repeated[:5]))

        # The layout is only useful if it densifies, and densification needs
        # quads. Whatever the network produced, it leaves here all-quad.
        mesh, note = solve_non_quad_faces(mesh, CoarsePseudoQuadMesh, self._loops(),
                                          poles=poles)
        if note:
            self.repair_notes.append(note)

        ok, reason = densifiable(mesh)
        if not ok:
            # The network did not close. Subdividing it further does not help --
            # a patch that snakes round on itself is self-overlapping, and every
            # quad a split makes from it is too. Discard the layout and mesh the
            # domain as a single polygon instead.
            #
            # This throws away the field's alignment, which is the whole point of
            # the front end, so it is recorded loudly. But a mesh that ignores
            # the field beats no mesh at all, and the domain outline is at least
            # geometrically sound.
            fallback = self._fallback_patch(boundary, reason)
            if fallback is not None:
                mesh = fallback

        mesh.attributes['decomposition_type'] = 'field'
        self.mesh = mesh
        self._mesh_poles = poles_key
        self._edited = False
        return self.mesh

    def coarse_mesh(self, poles=(), force=False):
        """**The coarse quad layout.** The name both front ends answer to.

        ``SkeletonDecomposition.coarse_mesh`` is the same step by the other
        method, so a script can swap one class for the other and change nothing
        else::

            decomposition = FieldDecomposition.from_boundary(...)
            coarse = decomposition.coarse_mesh()
            coarse.set_strips_density_target(t=0.5)
            dense = coarse.densification(field=decomposition.get_field())

        :meth:`decomposition_mesh` is the implementation and keeps working; this
        is its name in the shared workflow.
        """
        return self.decomposition_mesh(poles=poles, force=force)

    def get_field(self):
        """**The cross field, to hand to** ``densification(field=...)``.

        The field is the whole of what densification needs -- it carries its own
        background and builds its own point locator, so it travels on its own::

            field = decomposition.get_field()
            dense = any_coarse_layout.densification(field=field)

        which is how a layout that did NOT come from this decomposition -- a
        skeleton one, a hand-built one, one read back out of a document -- gets
        patch interiors that follow the field. The reverse is just as useful:
        densify this decomposition's own layout WITHOUT passing the field, and
        the interiors fall back to ``discrete_coons_patch``.

        Named to match ``CoarseQuadMesh.get_quad_mesh``. :attr:`field` is the
        same object and stays a plain attribute.
        """
        self._require_solved()
        return self.field

    def edit_coarse(self, geometry, poles=None, snap_tol=None, strict=True):
        """**Take a hand-edited coarse layout in place of the generated one.**

        The one seam in the pipeline. Hand back what came out of
        :meth:`decomposition_mesh` after editing it -- as a mesh, or as one
        closed polyline per patch with ONE POINT PER CORNER (see
        ``edit.py``, and use :func:`edit.face_polylines` to bake it). The layout
        is welded, snapped to the walls, repaired to all-quad and checked; from
        then on :meth:`quad_mesh` densifies THAT rather than regenerating one.

        **The field is not re-solved and not touched.** This sets
        :attr:`mesh` and nothing else. :attr:`field`, :attr:`tracer`,
        :attr:`background`, :attr:`separatrices` and the cached polyline network
        are all left exactly as they were -- an edit to the layout is not
        evidence about the field, and the separatrices are needed intact so
        :meth:`edges_to_curves` can warp them onto the edited corners.

        Vertex moves and topology changes are both supported, because nothing
        carries an index across the round trip: patches are re-welded from their
        corner coordinates, so a deleted patch is simply absent and a split one
        is two.

        **A bad edit raises rather than falling back.** :meth:`quad_mesh`
        promises to always return a mesh, but honouring that here would mean
        silently discarding the user's work and returning a triangulation that
        looks nothing like what they drew. Better to fail while they are still
        in Rhino and can fix it. ``strict=False`` downgrades it to a note.

        Parameters
        ----------
        geometry : mesh or list
            The edited layout. See :func:`edit.coarse_from_skeleton`.
        poles : list[[x, y, z]], optional
            Preferred pole positions. Defaults to the current layout's poles, so
            an untouched pseudo-quad keeps its collapsed corner where it was.
        snap_tol : float, optional
            How far a corner ON THE LAYOUT'S BOUNDARY may be projected onto a
            domain wall. Interior corners are never moved, whatever their
            distance. Defaults to half the background spacing. ``0`` disables
            it -- but read ``edit.py`` on why snapping is load-bearing first.
        strict : bool, optional

        Returns
        -------
        CoarsePseudoQuadMesh
            A **new** mesh: the weld and repair rebuild it, so vertex keys are
            RENUMBERED and no strips are collected on it. Two consequences for a
            caller holding state across the call:

            * anything keyed by vertex or strip INDEX must be re-keyed
              geometrically -- rounded edge midpoints work;
            * :meth:`densify` needs ``collect_strips`` and a density set first,
              or it raises ``KeyError`` from ``get_strip_density`` several
              frames away. :meth:`quad_mesh` does both itself.
        """
        # The separatrices have to exist before the edited layout arrives:
        # ``edges_to_curves`` warps against ``self.polylines``. ``_build`` is
        # cached, so on a decomposition that already produced a layout this is
        # free, and it never recomputes the field either way.
        # ``user_curves`` are re-appended on EVERY call, not once: this line is
        # the only thing keeping a curve the user drew on the layout alive. The
        # network from ``_build`` does not contain it -- nothing traced it -- so
        # without this the added edge finds no polyline in
        # :meth:`edges_to_curves` and densifies as a straight chord, which is
        # exactly the alignment loss the whole edit seam exists to avoid.
        boundary, others, _ = self._build()
        self.polylines = boundary + others + list(self.user_curves)

        if poles is None:
            poles = ([self.mesh.vertex_coordinates(v) for v in self.mesh.poles()]
                     if self.mesh is not None and hasattr(self.mesh, 'poles') else [])
        if snap_tol is None:
            snap_tol = self.background.target_length * 0.5

        mesh, notes = coarse_from_skeleton(
            geometry, loops=self._loops(), poles=poles, snap_tol=snap_tol,
            cls=CoarsePseudoQuadMesh)

        if mesh is None:
            raise ValueError(
                'the edited layout has no usable patch ({} face(s) in, {} '
                'degenerate). A patch must be a closed polyline with one point '
                'per corner.'.format(notes['faces_in'], notes['degenerate']))

        ok, why = densifiable(mesh)
        if not ok and strict:
            raise ValueError(
                'the edited layout will not densify: {}. It is NOT being '
                'replaced by a fallback -- fix the layout and hand it back. '
                '({} patch(es) in, {} out, side counts {})'.format(
                    why, notes['faces_in'], notes['faces_out'], notes['sides']))

        # ``route()`` classifies a run by substring, so this note must not
        # contain 'triangulation' or 'DISCARDED' -- either would relabel a
        # perfectly good edited field layout as a fallback.
        note = ('EDITED coarse layout adopted: {} patch(es) in, {} out, {} '
                'corner(s) snapped to a wall'.format(
                    notes['faces_in'], notes['faces_out'], notes['snapped']))
        if notes['off_wall']:
            note += (', {} corner(s) on the layout edge sit on no wall (further '
                     'than the snap tolerance from every loop) -- expected where '
                     'a patch was deleted, a mistake otherwise'.format(
                         notes['off_wall']))
        if notes['degenerate']:
            note += ', {} degenerate patch(es) dropped'.format(notes['degenerate'])
        if notes['repair']:
            note += '; repair: {}'.format(notes['repair'])
        if not ok:
            note += '; WILL NOT DENSIFY: {}'.format(why)
        self.repair_notes.append(note)

        self.mesh = mesh
        self._edited = True
        self.edit_notes = notes
        return self.mesh

    def set_user_curves(self, curves):
        """Adopt the curves a user drew on the layout, and publish them.

        Assigning :attr:`user_curves` alone is not enough once
        :meth:`edit_coarse` has already run: it appends them to
        :attr:`polylines` as it goes, so a curve handed over afterwards -- which
        is when the Rhino side knows the corrected endpoints, because the weld
        renumbers and the snap can move them -- would sit in ``user_curves``
        while the densifier read a :attr:`polylines` without it, and every added
        edge would quietly come out a straight chord.

        Rebuilt rather than appended so committing twice does not double them.
        ``_build`` is cached, so this costs nothing and never touches the field.
        """
        self.user_curves = [[list(point) for point in curve] for curve in curves]
        if self.polylines is not None:
            boundary, others, _ = self._build()
            self.polylines = boundary + others + list(self.user_curves)

    def quad_mesh(self, target_length=None, density=None, coarse=None,
                  densities=None):
        """**The deliverable: an all-quad mesh of the domain. Always.**

        Use this rather than assembling ``decomposition_mesh`` +
        ``collect_strips`` + ``densification`` by hand. Those three can each fail
        on a domain whose separatrix network did not close, and a front end that
        returns no mesh is useless however good its field was. This method owns
        the whole chain and every fallback in it.

        Three routes, best first, and :meth:`warnings` always says which was
        taken:

        1. **The field layout**, densified by :meth:`densify` -- coarse edges
           follow separatrices and patch interiors are integrated from the
           field. Element flow follows the field on both counts. This is the
           point of the whole front end and it is what every rectilinear plate,
           the pentagon and the square-holed plate get.
        2. **The domain as one polygon**, when the layout did not close but the
           domain is simply connected. Ignores the field; covers the shape.
        3. **A quad-split triangulation** at the requested element size, when
           even that is unavailable -- a multiply-connected domain, say. Ugly,
           high-valence, and correct.

        Route 1 is also what a HAND-EDITED layout takes -- see
        :meth:`edit_coarse`. It is the same route, with the layout supplied
        rather than generated; the field, the separatrices and the interior
        integration are identical.

        Parameters
        ----------
        **Per-strip densities are respected.** The rule is: an explicit size
        argument wins, and with no size argument whatever is already set wins::

            d.quad_mesh()                        # keeps every density already
                                                 # set; fills unset strips from
                                                 # the background spacing
            d.quad_mesh(target_length=0.3)       # 0.3 everywhere -- a size
                                                 # asked for is a size meant
            d.quad_mesh(target_length=0.3,       # 0.3 everywhere except these
                        densities={skey: 9})

        so asking for a coarser mesh still works, and setting a strip by hand
        beforehand is no longer silently discarded. It used to reset every strip
        on every call, unconditionally.

        Parameters
        ----------
        target_length : float, optional
            Target quad edge length, applied to every strip that ``densities``
            does not name. Omit it to keep the densities already on the layout;
            strips with none are filled from the background spacing.
        density : int, optional
            Fixed subdivision per coarse edge instead of a target length. Same
            rule as ``target_length`` and takes precedence over it.
        coarse : mesh or list, optional
            A HAND-EDITED coarse layout to densify instead of generating one --
            the same thing :meth:`edit_coarse` takes, and it is put through it.
            Once a layout has been edited it stays in use, so this only has to
            be passed once.
        densities : dict, optional
            ``{strip key: density}``, applied over the base pass. A key naming
            no strip of this layout is reported in :meth:`warnings` rather than
            raised -- an edit can legitimately remove the strip a density was
            set on. **Strip keys renumber whenever the layout is rebuilt**, so
            anything that outlives one call should key its densities
            geometrically and resolve them just before the call.

        Returns
        -------
        :class:`compas_singular.datastructures.QuadMesh`
        """
        self._require_solved()
        # A size that was asked for is a size that was meant, and it applies to
        # every strip. With neither argument there is nothing to override with,
        # so whatever the layout already carries stands.
        keep_preset = target_length is None and density is None
        target_length = target_length or self.background.target_length
        if coarse is not None:
            coarse = self.edit_coarse(coarse)
        elif self._edited:
            # An edit is not thrown away by asking for a second density.
            coarse = self.mesh
        else:
            coarse = self.decomposition_mesh()

        ok, why = densifiable(coarse)
        if ok:
            try:
                coarse.collect_strips()
                # Strips whose density is NOT the base pass's to set: the ones
                # named here, plus -- when no size was asked for -- the ones
                # that already carry one. ``skeys=`` on both setters is what
                # makes this a filter rather than a second density pass.
                pinned = set()
                for skey, value in (densities or {}).items():
                    if skey in coarse.attributes['strips']:
                        coarse.set_strip_density(skey, int(value))
                        pinned.add(skey)
                    else:
                        self.repair_notes.append(
                            'DENSITY: strip {!r} is not a strip of this layout; '
                            'its density was dropped. Strip keys renumber on '
                            'every rebuild -- key them geometrically.'.format(skey))
                if keep_preset:
                    pinned |= set(coarse.get_strip_densities())
                todo = [skey for skey in coarse.strips() if skey not in pinned]
                if density is not None:
                    coarse.set_strips_density(density, skeys=todo)
                else:
                    coarse.set_strips_density_target(target_length, skeys=todo)
                dense = self.densify(coarse)
                good, why, metrics = self._acceptable(dense)
                if good:
                    # Quality never reroutes -- see _acceptable -- but a
                    # degenerate element is said out loud, by name and value.
                    self._note_quality(dense, metrics)
                    self.dense = dense
                    return dense
                why = 'densified mesh rejected: {}'.format(why)
            except Exception as exc:
                why = '{} during densification'.format(type(exc).__name__)

        # Densifying the backstop is what fails -- ``collect_strips`` cannot make
        # consistent strips of a quad-split triangulation. So do not densify it:
        # build it AT the requested element size and return it directly. A coarse
        # all-quad mesh at the right size is the mesh.
        mesh = self._triangulation_fallback(why, target_length)
        good, bad_why, metrics = self._acceptable(mesh)
        if not good:
            self.repair_notes.append(
                'WARNING: even the triangulation backstop is imperfect ({}). The '
                'mesh is returned anyway -- inspect it before using it.'.format(bad_why))
        self._note_quality(mesh, metrics)
        self.mesh = mesh
        # The backstop replaced whatever layout was here, edited or not, and the
        # note above says so. Do not leave the edited flag set over a mesh that
        # is no longer the user's.
        self._edited = False
        # ...and it is a quad-split triangulation, not the separatrix layout, so
        # it must not be handed back as one. Clearing the key makes the next
        # ``decomposition_mesh()`` rebuild rather than cache-hit on it.
        self._mesh_poles = None
        self.dense = mesh
        return mesh

    def densify(self, coarse=None, **kwargs):
        """**The layout, densified with the field steering patch interiors.**

        The other half of the objective. ``edges_to_curves`` already makes the
        coarse EDGES follow separatrices; this makes the patch INTERIORS follow
        the field too, instead of filling them with a bilinear blend of their
        own boundaries that never consults it. A cable can then steer the
        inside of a patch without having to split the patch first -- which is
        the only way a cable reaches a mesh on a domain, like a plain square,
        where it forces no topology at all.

        Boundaries are held fixed, so the result still welds into a manifold
        mesh whose strips can be collected and edited. See ``densify.py``.

        Parameters
        ----------
        coarse : CoarsePseudoQuadMesh, optional
            Strips collected and densities set. Defaults to :attr:`mesh`.
        **kwargs
            ``stiffness``, ``guard``, ``iterations``, ``tikhonov`` -- passed
            straight to :func:`densify.field_densification`. Only
            ``16_field_densify.py --tradeoff`` uses them.

        Returns
        -------
        QuadMesh
        """
        if coarse is None:
            coarse = self.mesh
        if self.field_aware:
            # One implementation, on the field -- see ``CrossField.densify``.
            # ``spend`` is left to its default there, which is the same rule:
            # quality is spent to satisfy a guide and never otherwise.
            dense, stats = self.field.densify(
                coarse, edges_to_curves=self.edges_to_curves(), **kwargs)
        else:
            # ``field_aware=False`` reproduces ``discrete_coons_patch`` exactly
            # and is only ever set to measure the difference -- 16_field_densify.
            dense, stats = field_densification(
                coarse, self.field, self.tracer,
                edges_to_curves=self.edges_to_curves(), field_aware=False,
                **dict({'spend': bool(self.guides)}, **kwargs))
        self.densify_stats = stats
        if stats.get('guarded') and self.guides:
            # Only worth saying on a GUIDED domain, where a guarded patch means
            # a curve the user asked for did not reach that patch's interior.
            # Without a guide it is the ordinary outcome -- there was no
            # alignment to buy, so the Coons interior was already right.
            self.repair_notes.append(
                'DENSIFY: {} of {} patch(es) rejected the field-integrated '
                'interior and kept the Coons one.'.format(
                    stats['guarded'], stats['patches']))
        return dense

    def route(self):
        """Which of :meth:`quad_mesh`'s three routes produced the last mesh.

        ``'field'``    -- the separatrix layout; element flow follows the field.
        ``'polygon'``  -- the domain as one polygon; ignores the field.
        ``'triangulation'`` -- the backstop; ignores the field and is ugly.

        Worth printing next to any result: a mesh that came out of route 2 or 3
        is a mesh of the right SHAPE but says nothing about the cables or force
        lines that were fed in.
        """
        if any('triangulation' in n for n in self.repair_notes):
            return 'triangulation'
        if any('DISCARDED' in n for n in self.repair_notes):
            return 'polygon'
        return 'field'

    def quality(self, mesh=None, low_angle=None):
        """**The numbers to gate on.** Element quality of a mesh, as a dict.

        Structural correctness -- all quads, manifold, right area -- says
        nothing about whether the elements are usable. Those three properties
        certified a mesh carrying a literal 180 degree angle, and that mesh is
        live in the suite: ``12_cables``' ring cable. So measure the elements
        too, and let callers gate on the numbers rather than on a verdict
        somebody typed out by hand.

        Parameters
        ----------
        mesh : Mesh, optional
            Defaults to the last mesh :meth:`quad_mesh` returned. NOT
            ``self.mesh``, which on the field route is the coarse layout.
        low_angle : float, optional
            Threshold for ``share_below``, in degrees. Defaults to
            ``quality.LOW_ANGLE``.

        Returns
        -------
        dict
            Everything :func:`quality.mesh_quality` returns -- ``min_angle``,
            ``max_angle``, ``aspect_max``, ``share_below``,
            ``irregular_interior``, ``poles``, ``faces``, ``worst_face`` --
            plus ``route``, ``coverage`` and ``coarse_faces``, which are
            properties of the decomposition rather than of the mesh.

        Pseudo-quads are handled throughout; see ``quality.py``'s docstring for
        why a naive reading of one reports 0 degrees and divides by zero.
        """
        if mesh is None:
            mesh = self.dense if self.dense is not None else self.mesh
        if mesh is None:
            raise ValueError('no mesh yet -- call quad_mesh() or '
                             'decomposition_mesh() first')

        kwargs = {} if low_angle is None else {'low_angle': low_angle}
        out = mesh_quality(mesh, **kwargs)
        out['route'] = self.route()
        out['coverage'] = self._coverage(mesh)
        out['coarse_faces'] = (self.mesh.number_of_faces()
                               if self.mesh is not None else None)
        return out

    def _coverage(self, mesh):
        """The mesh's area as a fraction of the domain's. 1.0 is exact."""
        def polygon_area(loop):
            pts = list(loop)
            return abs(0.5 * sum(p[0] * q[1] - q[0] * p[1]
                                 for p, q in zip(pts, pts[1:] + pts[:1])))

        want = polygon_area(self.background.outer)
        want -= sum(polygon_area(loop) for loop in self.background.inners)
        if want <= 0:
            return 1.0
        got = 0.0
        for fkey in mesh.faces():
            pts = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
            got += 0.5 * sum(p[0] * q[1] - q[0] * p[1]
                             for p, q in zip(pts, pts[1:] + pts[:1]))
        return abs(got) / want

    def _acceptable(self, mesh):
        """Does this mesh represent the domain, and are its elements usable?

        Two separate questions, and the return value keeps them separate on
        purpose.

        **Structure** -- manifold, all quads, covering the right AREA -- is what
        ``ok`` reports. The area check is the one that catches a layout which
        parsed and densified but does not match the shape: a round-holed plate
        produced a manifold all-quad mesh covering 108% of its domain, because
        the patches ran across the hole. Nothing structural is wrong with such a
        mesh; it is simply not a mesh of this domain.

        Triangles are allowed where they are the FAN AROUND A POLE, which is
        what a coarse pseudo-quad densifies to and what ``face_pole`` records.
        A pole is a legitimate quad-mesh singularity, not a defect; a triangle
        that is not one still fails.

        **Element quality** is returned as ``metrics`` and deliberately does NOT
        move ``ok``. Structure decides the ROUTE, and rerouting on quality would
        make things worse, not better: the only place left to reroute to is the
        triangulation backstop, whose own minimum angle is 21.5 degrees with 802
        irregular interior vertices on the round-holed disc. Swapping a poor
        field mesh for that is a downgrade on every metric except the one that
        triggered it. So :meth:`quad_mesh` returns the mesh and says loudly what
        is wrong with it -- see :func:`quality.hard_floor`.

        Returns
        -------
        (bool, str, dict)
            Structural verdict, the first structural reason when false, and the
            element-quality metrics. ``metrics`` is ``{}`` when the mesh was
            rejected before it was worth measuring.
        """
        if not mesh.faces():
            return False, 'no faces', {}
        face_pole = mesh.attributes.get('face_pole') or {}
        for fkey in mesh.faces():
            n = len(mesh.face_vertices(fkey))
            if n == 4 or (n == 3 and fkey in face_pole):
                continue
            return False, 'not all quads (face {} has {} sides)'.format(fkey, n), {}
        if not mesh.is_manifold():
            return False, 'not manifold', {}

        metrics = mesh_quality(mesh)

        ratio = self._coverage(mesh)
        if not 0.97 <= ratio <= 1.03:
            if self._edited:
                # A HAND-EDITED layout is not measured against the domain. The
                # area check exists to catch a GENERATED layout that parsed but
                # does not match the shape; a user who deleted a patch to leave
                # a void has changed the layout on purpose, and rerouting them
                # to the triangulation backstop would throw the edit away and
                # return a mesh that looks nothing like what they drew. Say it
                # and keep the mesh.
                self.repair_notes.append(
                    'EDITED layout covers {:.0%} of the domain area. Kept -- the '
                    'area check does not overrule a hand edit. If that was not '
                    'intended, a patch is missing or overlapping.'.format(ratio))
                return True, '', metrics
            return False, 'covers {:.0%} of the domain area'.format(ratio), metrics
        return True, '', metrics

    def _note_quality(self, mesh, metrics):
        """Record a hard-floor breach loudly, without changing the route.

        The note must not contain the words ``triangulation`` or ``DISCARDED``:
        :meth:`route` reads ``repair_notes`` by substring, so either word here
        would silently relabel a perfectly good field layout as a fallback.
        """
        ok, why = hard_floor(metrics or {})
        if ok:
            return
        self.repair_notes.append(
            'QUALITY FAILURE: {}. The mesh is returned anyway -- the '
            '"always returns a mesh" contract stands -- but this element is '
            'degenerate, not merely poor, and nothing downstream should '
            'assume otherwise.'.format(why))

    def _loops(self):
        """Every boundary loop, for boundary-aware repair."""
        return [self.background.outer] + list(self.background.inners)

    def _fallback_patch(self, boundary, reason):
        """Mesh the domain as one polygon, ignoring the separatrices.

        Only for simply-connected domains: with a hole, a single polygon is not
        the right topology and there is nothing honest to fall back to, so the
        broken layout is returned as-is with the reason recorded.
        """
        if not self.background.inners:
            outer_arcs = [arc for arc in boundary
                          if _on_loop(arc[0], self.background.outer)]
            mesh = self._single_patch(outer_arcs or boundary)
            mesh, _ = solve_non_quad_faces(mesh, CoarsePseudoQuadMesh, self._loops())
            ok, why = densifiable(mesh)
            if ok:
                self.repair_notes.append(
                    'FIELD LAYOUT DISCARDED -- not densifiable ({}). Meshed the '
                    'domain as a single {}-gon instead; the result ignores the '
                    'field.'.format(reason, len(outer_arcs or boundary)))
                return mesh
            reason = '{}; single-polygon fallback also failed ({})'.format(reason, why)

        return self._triangulation_fallback(reason)

    def _triangulation_fallback(self, reason, target_length=None):
        """Last resort: quad-split a coarse triangulation of the domain.

        Every triangle becomes three quads. This works for ANY domain topology --
        holes, several holes, whatever -- because it never has to reason about
        the layout at all: the domain is already triangulated, and splitting a
        triangulation is unconditional.

        It is the ugly answer. Each original triangle vertex becomes a
        high-valence singularity and the element flow is the triangulation's, not
        the field's. But the alternative for a multiply-connected domain whose
        separatrix network did not close is NO MESH, and a mesh that covers the
        domain exactly and densifies is worth more than nothing. The note says
        plainly what happened.
        """
        xs = [p[0] for p in self.background.outer]
        ys = [p[1] for p in self.background.outer]
        diagonal = ((max(xs) - min(xs)) ** 2 + (max(ys) - min(ys)) ** 2) ** 0.5

        # A quad split halves the element size, so triangulate at twice the
        # target to land near it -- but never COARSER than the background, which
        # was already validated against this domain. A comb-shaped plate with
        # 2-unit slots triangulated at 1.2 bridges its own slots: 104% coverage
        # and a non-manifold result.
        spacing = (target_length * 2.0) if target_length else (diagonal / 6.0)
        spacing = min(spacing, self.background.target_length)
        coarse_bg = BackgroundMesh.from_boundary(
            self.background.outer, self.background.inners, target_length=spacing)

        index = {v: i for i, v in enumerate(coarse_bg.mesh.vertices())}
        vertices = [coarse_bg.mesh.vertex_coordinates(v) for v in coarse_bg.mesh.vertices()]
        faces = [[index[v] for v in coarse_bg.mesh.face_vertices(f)]
                 for f in coarse_bg.mesh.faces()]
        vertices, faces = topological_quad_split(vertices, faces, self._loops())
        mesh = CoarsePseudoQuadMesh.from_vertices_and_faces(vertices, faces)
        mesh, _ = solve_non_quad_faces(mesh, CoarsePseudoQuadMesh, self._loops())

        ok, why = densifiable(mesh)
        self.repair_notes.append(
            'FIELD LAYOUT DISCARDED -- not densifiable ({}). Fell back to a '
            'quad-split triangulation of the domain: {} patches that cover the '
            'shape but follow the triangulation, not the field.{}'.format(
                reason, mesh.number_of_faces(),
                '' if ok else ' STILL not densifiable: {}'.format(why)))
        return mesh

    def _single_patch(self, boundary):
        """One face spanning a domain with no interior lines.

        Any corner count, not just four: a hexagonal plate has no singularities
        and no reflex corners, so nothing cuts it, yet it is a perfectly good
        6-gon that :func:`solve_non_quad_faces` turns into six quads. Raising
        here instead -- as this did -- means the caller gets no mesh at all for a
        shape that has an obvious one.
        """
        corners = [list(arc[0]) for arc in boundary]
        if len(corners) < 3:
            raise ValueError(
                'domain has neither interior separatrices nor enough boundary '
                'corners to bound a face ({} found)'.format(len(corners)))
        if len(corners) != 4:
            self.repair_notes.append(
                'no interior separatrices and {} boundary corners: built as a '
                'single {}-gon'.format(len(corners), len(corners)))
        return CoarsePseudoQuadMesh.from_vertices_and_faces(
            corners, [list(range(len(corners)))])

    def _cut_holes(self):
        """Inner loops that needed cutting open -- no corner, no singularity.

        The condition ``Tracer.hole_launches`` uses, asked again here because it
        decides how the faces are recovered: a domain whose holes were all cut
        this way has every patch corner on a boundary, which is the one case
        ``from_polylines`` cannot return. See :meth:`decomposition_mesh`.
        """
        from compas_singular.framefield.separatrix_network import boundary_corners

        return [loop for loop in self.background.inners
                if not boundary_corners(loop)]

    def _mesh_from_faces(self, faces):
        """A coarse mesh from faces given as lists of corner points.

        Corners are welded by rounded coordinate, the same resolution
        ``from_polylines`` matches endpoints at, so two patches meeting at a cut
        share that vertex instead of each carrying its own copy.

        The same weld an edited layout comes back through, which is why it lives
        in ``edit.py`` -- an arrangement face and a Rhino patch outline are the
        same kind of thing here, a ring of corner points with no indices.
        """
        mesh, _dropped = mesh_from_faces(faces, CoarsePseudoQuadMesh)
        return mesh

    def edges_to_curves(self):
        """``{(u, v): polyline}`` for densification.

        Separatrices are CURVED. Chording each coarse edge straight would throw
        away exactly the alignment the field was solved for, so the traced
        polyline for every coarse edge is handed to ``densification``. The payoff
        is larger here than for the boundary-curve case it was first built for.

        Four ways an edge finds its curve, best first. The first two are all a
        GENERATED layout ever needs -- its corners are the traced endpoints, so
        the exact match hits. The third exists for a HAND-EDITED layout, where a
        moved corner means no traced polyline ends there any more:

        1. **exact match** -- the traced polyline whose ends are this edge's
           ends, by geometric key. Untouched edges take this and are unchanged;
        2. **boundary arc** -- both ends on one wall, so re-derive the arc from
           the loop. Without it a curved domain loses area outright: measured on
           an ellipse, 65% coverage, and 74% on a disc. **A HOLE is where this
           branch is hard**, and for one reason: its loop is short. One coarse
           edge can be more than half of it, so which way round is not settled by
           length (:func:`_arcs_between` offers both); its ring corners can be a
           tenth of the background spacing apart, so "another corner lies on this
           arc" is not settled by distance (:meth:`_arc_is_one_edge` measures
           along it); and the loop is a POLYGON sampling a circle, so a corner
           really on the wall is a sagitta off the loop
           (:func:`_loop_chord_tolerance`). Each of the three, on its own, turned
           one circle of several into a polygon while its neighbours stayed
           round;
        3. **warped separatrix** -- the traced polyline this edge came from,
           moved onto its new endpoints by :func:`edit.warp_polyline`. This is
           what makes an edit cost the nudge rather than the curvature.
           **ONLY on an edited layout**, and that restriction is not caution,
           it is correctness: on a generated layout the first two branches are
           by construction the complete answer, so anything reaching here is an
           edge that is MEANT to be straight -- a quad-split piece, a fallback
           patch. Warping those changes results nobody asked to change, and it
           did: allowing it cost ``square+circle @ bg 0.60`` a 0.001 degree
           element and moved 3 of the baseline's 84 rows;
        4. **chord** -- a straight line, and the count is reported.

        Which branch each edge took goes into :attr:`edit_notes` under
        ``'edges'``, so a user who moved a corner can see that, say, 3 of 40
        edges ended up straight.
        """
        self._require_solved()
        if self.mesh is None:
            self.decomposition_mesh()

        from compas.tolerance import TOL
        lookup = {}
        for polyline in (self.polylines or []):
            key = (TOL.geometric_key(polyline[0]), TOL.geometric_key(polyline[-1]))
            lookup[key] = polyline
            lookup[(key[1], key[0])] = list(reversed(polyline))

        tally = {'exact': 0, 'boundary': 0, 'warped': 0, 'chord': 0}
        claimed = set()
        # BOUNDARY corners only -- see _arc_is_one_edge. An interior corner
        # cannot end a boundary edge, so it is not evidence that an arc spans
        # several of them.
        corners = [self.mesh.vertex_coordinates(w) for w in self.mesh.vertices()
                   if self.mesh.is_vertex_on_boundary(w)]
        # Once per loop, not once per edge: a loop that samples a CURVE is a
        # chord of it, so a corner really on the wall is a sagitta off the loop.
        # See _loop_chord_tolerance -- a polygon still gets 1e-6.
        walls = [(loop, _loop_chord_tolerance(loop))
                 for loop in [self.background.outer] + list(self.background.inners)]
        out = {}
        for u, v in self.mesh.edges():
            pa = self.mesh.vertex_coordinates(u)
            pb = self.mesh.vertex_coordinates(v)
            curve = lookup.get((TOL.geometric_key(pa), TOL.geometric_key(pb)))
            how = 'exact' if curve is not None else None

            if curve is None:
                # Both ends on a wall, but no polyline runs between them: the
                # edge is a PIECE of a boundary arc, from a quad split, a
                # fallback patch, or a corner an edit snapped to the wall.
                # Same idea as guide_lines.boundary_arc_between.
                for loop, wall_tol in walls:
                    if _on_loop(pa, loop, wall_tol) and _on_loop(pb, loop, wall_tol):
                        # Both ways round, shorter first -- on a hole the
                        # shorter one is not always the edge. See _arcs_between.
                        for arc in _arcs_between(loop, pa, pb):
                            if self._arc_is_one_edge(arc, pa, pb, corners):
                                curve = arc
                                how = 'boundary'
                                break
                        break

            if curve is None and self._edited:
                curve = self._warped_curve(pa, pb, claimed)
                how = 'warped' if curve is not None else None

            # ``densification`` looks EVERY edge up once it is given a mapping at
            # all -- it has no straight-chord branch per edge, only per call --
            # so an incomplete dict raises KeyError. Fill the gaps with the
            # chord, which is what the no-mapping branch would have used anyway.
            if (curve is None
                    or distance_point_point(curve[0], pa) > 1e-6
                    or distance_point_point(curve[-1], pb) > 1e-6):
                curve = [list(pa), list(pb)]
                how = 'chord'
            tally[how] += 1
            out[u, v] = curve

        self.edit_notes['edges'] = tally
        return out

    def _arc_is_one_edge(self, arc, pa, pb, corners):
        """Is this boundary arc ONE edge of the layout, or several?

        The boundary branch of :meth:`edges_to_curves` assumes that an edge with
        both ends on a wall IS a piece of that wall. That is true of a quad-split
        piece or a fallback patch's side, and false of a **chord across a
        corner** -- an edge whose two ends happen to sit on the same loop but
        whose patch lies inside it.

        The two are told apart by asking what the arc runs through. A genuine
        piece of wall runs between two ADJACENT corners of the layout, so no
        other corner lies on it. An arc that passes through one has gone round
        a corner of the domain and come back, and taking it makes the patch
        retrace its own other sides.

        ``corners`` must be the layout's BOUNDARY corners only. What this
        detects is an arc spanning several boundary EDGES, and only a boundary
        vertex can be one of their ends, so an interior vertex is not evidence
        of anything and letting one veto is a category error. Not a rare one: a
        hole near another feature puts interior corners just outside its wall,
        and one 0.076 from a radius-0.7 hole -- against a ``tol`` of 0.119 --
        was what left that hole with two chord edges slicing 21% of the way
        across it while its two identical siblings came out round.

        This is not hypothetical. ``22_force_lines``' arch: the bottom-left
        patch is the triangle (0.00, 2.11) (0.00, 0.00) (1.66, 0.00), whose third
        side has both ends on the outer wall. The arc it was given ran along the
        bottom wall to (0, 0) and then up the left wall -- 46 of its 82 points on
        ``x = 0`` -- so the patch enclosed nothing, densified into two zero-area
        quads (minimum angle 0.00, maximum 180.00), and the whole 195-face field
        mesh was rejected for a 2520-face triangulation that says nothing about
        the cables. Guarding this leaves it on the field route at 195 faces,
        minimum angle 10.75, coverage 1.0000.

        **THROUGH is measured along the arc, not across it.** A corner that is
        near the arc because it is near one of its ENDS is not one the arc runs
        through, and vetoing on distance alone made that mistake constantly on
        small holes: ``tol`` is a quarter of the background spacing, a length set
        by the whole domain, while the thing it has to resolve is the corner
        spacing on one small circle. Measured on a radius-0.7 hole at background
        0.5 (so ``tol`` = 0.125 against a 4.4-long loop): any two ring corners
        closer than 12 degrees vetoed BOTH of their outward neighbours -- two of
        the ring's four edges -- so one circle came out round and an identical
        one beside it a polygon. Requiring the offending corner to sit a margin
        clear of both ends restores those without weakening the guard: in the
        ``22_force_lines`` case the domain corner is in the MIDDLE of the arc.

        Returns ``True`` when the arc may be used; the caller tries the other way
        round the loop, then falls back to the chord.
        """
        if not arc or len(arc) < 3:
            return bool(arc)

        tol = self.background.target_length * 0.25
        xs = [p[0] for p in arc]
        ys = [p[1] for p in arc]
        lo_x, hi_x = min(xs) - tol, max(xs) + tol
        lo_y, hi_y = min(ys) - tol, max(ys) + tol

        for point in corners:
            # the arc's own endpoints are not "another corner"
            if (distance_point_point(point, pa) < 1e-6
                    or distance_point_point(point, pb) < 1e-6):
                continue
            if not (lo_x <= point[0] <= hi_x and lo_y <= point[1] <= hi_y):
                continue
            d, s, total = project_on_polyline(point, arc)
            if d >= tol:
                continue
            # capped so the margin can never swallow the arc it is protecting
            margin = min(tol, 0.25 * total)
            if s <= margin or s >= total - margin:
                continue
            return False
        return True

    def _warped_curve(self, pa, pb, claimed):
        """The traced separatrix this edge came from, moved onto its endpoints.

        Only ever reached on an edited layout: on a generated one every edge
        matches a polyline exactly or is a boundary arc.

        **The match is ANCHORED on an end that did not move.** That is the rule
        that makes this well posed rather than a nearest-neighbour guess, and it
        follows from what an edit actually is. Dragging one corner moves one end
        of each edge around it and leaves the other exactly where the tracer put
        it -- so the right separatrix is the one still ending at the still
        point, and among those, the one whose free end was nearest. No distance
        threshold has to be invented for the end that MOVED, which is the one
        the user may have dragged as far as they liked. An earlier version put a
        one-background-spacing cap on both ends and it rejected every real edit:
        a 1.8-unit drag on a 0.6 background left five edges as straight chords,
        which is exactly the alignment loss this branch exists to prevent.

        Only when NEITHER end is anchored -- both corners of one edge moved --
        does it fall back to nearest-total within one background spacing, where
        the correspondence really is a guess and should be a timid one.

        Two guards remain, because the wrong curve warped is worse than no curve
        at all -- it densifies into a patch that bulges through its neighbour
        rather than one that is merely straight:

        * :func:`edit.warp_polyline` rejects a warp that reaches further than
          the candidate's own length, which would fold it;
        * a candidate whose length is wildly out of proportion to the new chord
          is a different curve, not this one moved.

        A polyline already claimed by another edge is only reused when nothing
        unclaimed fits -- the case where one patch was split in two and both
        halves genuinely lie along one separatrix.
        """
        if not self.polylines:
            return None

        anchor = self.background.target_length * 0.1
        loose = self.background.target_length
        chord = distance_point_point(pa, pb)

        best = None
        for i, polyline in enumerate(self.polylines):
            for points in (polyline, list(reversed(polyline))):
                da = distance_point_point(points[0], pa)
                db = distance_point_point(points[-1], pb)

                if da <= anchor or db <= anchor:
                    # anchored: one end never moved, so trust the other freely
                    tier, free = 0, (db if da <= anchor else da)
                elif da <= loose and db <= loose:
                    tier, free = 1, da + db
                else:
                    continue

                length = sum(distance_point_point(a, b) for a, b in pairwise(points))
                if chord > 1e-9 and not (0.25 * chord <= length <= 4.0 * chord):
                    continue

                score = (i in claimed, tier, free)
                if best is None or score < best[0]:
                    best = (score, i, points)

        if best is None:
            return None
        _score, index, points = best
        curve = warp_polyline(points, pa, pb)
        if curve is None:
            return None
        claimed.add(index)
        return curve

    # ------------------------------------------------------------------

    def _build(self):
        """The polyline network, as a PLANAR ARRANGEMENT.

        ``build_network`` snaps and splits, but it leaves polylines crossing in
        their interiors with no node at the crossing -- and ``from_polylines``
        never computes one either, so the graph its face search walks is not a
        planar embedding and the faces it recovers repeat vertices. That is what
        sent every curved domain to the polygon fallback; see ``arrangement.py``
        for the measurements. Computing the arrangement here means
        ``decomposition_polylines``, ``decomposition_mesh`` and
        ``edges_to_curves`` all see the same, consistent network.
        """
        self._require_solved()
        if self._network is None:
            boundary, others, info = build_network(
                self.field, self.separatrices,
                singularity_points=self.tracer._singularity_points(),
                symmetry=self.symmetry)
            info = dict(info)
            boundary, others = planar_arrangement(
                boundary, others, self.background.target_length * 0.8, info,
                loops=self._loops())
            self._network = (boundary, others, info)
        return self._network

    def report(self):
        boundary, others, info = self._build()
        out = dict(info)
        out.update(self.field.report())
        if self.mesh is not None:
            faces = [len(self.mesh.face_vertices(f)) for f in self.mesh.faces()]
            out['faces'] = len(faces)
            out['quads'] = sum(1 for n in faces if n == 4)
            out['non_quads'] = sorted(set(n for n in faces if n != 4))
        return out
