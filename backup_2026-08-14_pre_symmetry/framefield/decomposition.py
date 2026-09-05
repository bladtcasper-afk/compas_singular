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

from .arrangement import faces_from_arrangement
from .arrangement import faces_with_repeated_vertices
from .arrangement import planar_arrangement
from .background import BackgroundMesh
from .constraints import from_boundary
from .constraints import from_curves
from .densify import field_densification
from .edit import coarse_from_skeleton
from .edit import mesh_from_faces
from .edit import warp_polyline
from .field import CrossField
from .quality import hard_floor
from .quality import mesh_quality
from .repair import build_network
from .repair import densifiable
from .repair import solve_non_quad_faces
from .repair import topological_quad_split
from .trace import Tracer


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


def _distance_to_polyline(point, points):
    """Shortest distance from a point to an OPEN polyline.

    ``edit.closest_on_loop`` is the closed-loop version and is the wrong tool
    here: it joins last to first, which on an arc invents a segment the arc does
    not have.
    """
    best = float('inf')
    for a, b in pairwise(points):
        abx, aby = b[0] - a[0], b[1] - a[1]
        length2 = abx * abx + aby * aby
        if length2 == 0.0:
            continue
        t = ((point[0] - a[0]) * abx + (point[1] - a[1]) * aby) / length2
        t = max(0.0, min(1.0, t))
        q = [a[0] + abx * t, a[1] + aby * t, 0.0]
        best = min(best, distance_point_point(point, q))
    return best


def _arc_between(loop, pa, pb):
    """The shorter arc of a closed loop between two points ON it, pa -> pb.

    The SHORTER of the two ways round: a coarse edge spans a fraction of the
    boundary, never most of it, so the short arc is always the intended one.
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
        return None

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

    sa, sb = param(pa), param(pb)
    forward = (sb - sa) % total
    reverse = forward > total - forward
    s0, span = (sb, total - forward) if reverse else (sa, forward)
    if span <= 1e-12:
        return None

    def point_at(s):
        s = s % total
        for i, (c0, c1) in enumerate(zip(cum, cum[1:])):
            if c0 - 1e-9 <= s <= c1 + 1e-9:
                a, b = ring[i], ring[i + 1]
                length = c1 - c0
                t = 0.0 if length == 0 else (s - c0) / length
                return [a[k] + (b[k] - a[k]) * t for k in range(3)]
        return list(ring[-1])

    arc = [list(pa if not reverse else pb)]
    steps = 1
    while steps * (total / max(1, len(pts))) < span:
        steps += 1
    for i in range(1, len(pts) + 1):
        s = s0 + span * i / float(len(pts) + 1)
        arc.append(point_at(s))
    arc.append(list(pb if not reverse else pa))
    if reverse:
        arc.reverse()
    arc[0], arc[-1] = list(pa), list(pb)
    return arc


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
                      guide_weight=1.0, guide_band=None, relax=False, field_tau=None):
        """Solve the field and trace its separatrices.

        Parameters
        ----------
        outer_boundary : list[[x, y, z]]
        inner_boundaries : list[list[[x, y, z]]], optional
        guides : list[list[[x, y, z]]], optional
            Guide curves -- cables, force lines.
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
        """
        if orthogonal is False:
            raise NotImplementedError(
                'non-orthogonal frame fields need warp.py -- milestone 3. '
                'Cases A and B (orthogonal stress field, single cable family) '
                'are handled by the cross field; see the evaluation doc.')

        background = BackgroundMesh.from_boundary(
            outer_boundary, inner_boundaries, target_length=target_length)

        constraints = from_boundary(background)
        if guides:
            constraints += from_curves(background, guides, mode=mode,
                                       weight=guide_weight, band=guide_band)

        field = CrossField.solve(background, constraints, relax=relax, tau=field_tau)
        tracer = Tracer(field)
        separatrices, trace_report = tracer.separatrices()
        out = cls(background, field, separatrices, tracer)
        out.guides = list(guides or [])
        out.trace_report = trace_report
        return out

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

    def decomposition_mesh(self, poles=()):
        """The coarse quad mesh.

        Parameters
        ----------
        poles : list[[x, y, z]], optional
            Preferred pole positions, as in
            ``SkeletonDecomposition.decomposition_mesh``. A field singularity is
            not a pole and none is needed for the domains in ``14_arrangement``,
            but :func:`repair.solve_non_quad_faces` resolves any triangular
            patch it cannot otherwise repair as a pseudo-quad pole, and a
            position given here decides which of that triangle's corners the
            collapsed side sits at.

        Returns
        -------
        CoarsePseudoQuadMesh
        """
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

        self.mesh = mesh
        self._edited = False
        return self.mesh

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

    def quad_mesh(self, target_length=None, density=None, coarse=None):
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
        target_length : float, optional
            Target quad edge length. Defaults to the background spacing.
        density : int, optional
            Fixed subdivision per coarse edge instead of a target length.
        coarse : mesh or list, optional
            A HAND-EDITED coarse layout to densify instead of generating one --
            the same thing :meth:`edit_coarse` takes, and it is put through it.
            Once a layout has been edited it stays in use, so this only has to
            be passed once.

        Returns
        -------
        :class:`compas_singular.datastructures.QuadMesh`
        """
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
                if density is not None:
                    coarse.set_strips_density(density)
                else:
                    coarse.set_strips_density_target(target_length)
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
        dense, stats = field_densification(
            coarse, self.field, self.tracer,
            edges_to_curves=self.edges_to_curves(),
            field_aware=self.field_aware,
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
        from .repair import boundary_corners

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
           an ellipse, 65% coverage, and 74% on a disc;
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
        corners = [self.mesh.vertex_coordinates(w) for w in self.mesh.vertices()]
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
                for loop in [self.background.outer] + list(self.background.inners):
                    if _on_loop(pa, loop) and _on_loop(pb, loop):
                        arc = _arc_between(loop, pa, pb)
                        if self._arc_is_one_edge(arc, pa, pb, corners):
                            curve = arc
                            how = 'boundary'
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

        This is not hypothetical. ``22_force_lines``' arch: the bottom-left
        patch is the triangle (0.00, 2.11) (0.00, 0.00) (1.66, 0.00), whose third
        side has both ends on the outer wall. The arc it was given ran along the
        bottom wall to (0, 0) and then up the left wall -- 46 of its 82 points on
        ``x = 0`` -- so the patch enclosed nothing, densified into two zero-area
        quads (minimum angle 0.00, maximum 180.00), and the whole 195-face field
        mesh was rejected for a 2520-face triangulation that says nothing about
        the cables. Guarding this leaves it on the field route at 195 faces,
        minimum angle 10.75, coverage 1.0000.

        Returns ``True`` when the arc may be used; the caller falls back to the
        chord otherwise.
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
            if _distance_to_polyline(point, arc) < tol:
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
        if self._network is None:
            boundary, others, info = build_network(self.field, self.separatrices)
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
