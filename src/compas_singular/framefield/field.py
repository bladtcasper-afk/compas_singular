"""Steps 2 and 3 -- the cross field, and its singularities.

A cross (four directions 90 degrees apart) is represented by one complex number
per vertex, ``u = exp(i * 4 * theta)``. Raising to the 4th power is what makes
the 90-degree-symmetric object single-valued, which turns "find the smoothest
cross field" into a plain sparse linear solve -- no integer variables, no period
jumps to track. On a PLANAR domain every tangent space is world XY, so there is
no connection to parallel-transport across either.

TWO SOLVERS, and the default one is known to be inaccurate
----------------------------------------------------------

``solve(relax=False)`` -- the original -- minimises the Dirichlet energy of ``u``
with ``|u|`` left completely FREE. The magnitude collapsing toward 0 in the
interior is not a solver failure, and those zeros sit where the singularities
are. What this docstring used to claim is that the collapse IS the
Ginzburg-Landau relaxation. It is not: with no penalty holding ``|u|`` near 1
there is no relaxation, only its symptom. This is the ``epsilon -> infinity`` end
of the Ginzburg-Landau family, and it puts the singularities in the wrong place.

Measured on the suite's disc, radii normalised by the disc radius, against the
published value of ``r ~ 0.85`` (Dai/Qiao/Wang 2026 section 6.2; Beaufort et al.
2017's renormalized energy for circular domains)::

    target_length   0.50    0.35    0.25
    relax=False     0.52..  0.39..  0.34..     drifting INWARD, no limit
    relax=True      0.85    0.85    0.85       resolution-independent

``solve(relax=True)`` adds the missing step: alternate an implicit diffusion
solve with pointwise normalisation back onto the unit circle, to steady state.
That is the diffusion generated method -- Merriman/Bence/Osher, applied to cross
fields by Viertel & Osting 2019 and analysed by Dai/Qiao/Wang 2026. It is the
whole of their Algorithm 4.1 that we need; their diffuse domain method and FFT
solve address a problem (avoiding a triangulation) that we do not have.

See ``19_field_accuracy.py``, which is the check, and
``paper-diffusion-generated-crossfields.md`` section 3 for the full measurement.

``relax`` DEFAULTS TO FALSE. Turning it on moves the layout on every curved
domain, so it is opt-in until ``baseline.json`` says which way those rows went.
"""
from cmath import phase
from math import cos
from math import sin
from numbers import Number

from compas.data import Data
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import factorized
from scipy.sparse.linalg import spsolve

from compas_singular.framefield.constraints import PERIOD
from compas_singular.framefield.constraints import from_boundary
from compas_singular.framefield.constraints import representation


__all__ = ['CrossField', 'field_provenance', 'wrap_to_period']


#: Bumped when :meth:`CrossField.save_to_json`'s layout changes in a way a
#: reader has to notice. :meth:`CrossField.load_from_json` refuses anything
#: newer than it understands rather than guessing.
JSON_VERSION = 1

#: Written into every file and checked on load, so pointing this at a coarse
#: mesh's JSON fails with a sentence instead of a KeyError.
JSON_TYPE = 'compas_singular.framefield.CrossField'

#: Coordinate rounding for :func:`field_provenance`. Well below anything
#: geometric in the pipeline (``TOL.geometric_key`` works at 3) but coarse enough
#: to absorb the last-bit noise a CAD curve round-trip leaves on a point that has
#: not moved. **Changing it makes every stored field report its outline changed.**
COORDINATE_DIGITS = 9


def _canonical_curve(curve):
    """One curve as a plain, rounded list of ``[x, y, z]``.

    Accepts what ``from_boundary`` accepts -- lists, tuples, compas ``Point``s,
    a compas ``Polyline`` -- and drops a repeated closing point, because a loop
    given closed and the same loop given open are the same input to the solver
    (``background._as_open_loop`` drops it too).
    """
    points = []
    for point in curve:
        points.append([round(float(point[i]), COORDINATE_DIGITS) for i in range(3)])
    if len(points) > 1 and points[0] == points[-1]:
        points = points[:-1]
    return points


def _plain(value):
    """A solve parameter as something JSON can hold and a human can read."""
    if value is None or isinstance(value, bool) or isinstance(value, str):
        return value
    if isinstance(value, Number):
        return float(value)
    # A Symmetry. Its group IS the identity -- centre, elements and the steps it
    # is enabled for -- and none of that is reconstructible from ``repr``.
    if hasattr(value, 'centre') and hasattr(value, 'steps'):
        return {
            'symmetry': [round(float(c), COORDINATE_DIGITS) for c in value.centre],
            'names': sorted(value.names()),
            'steps': sorted(value.steps),
        }
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, dict):
        return dict((str(k), _plain(v)) for k, v in value.items())
    return repr(value)


def _describe_difference(before, after):
    """The first differing entry of two readable dicts, as ``'name a -> b'``."""
    before = before or {}
    after = after or {}
    for name in sorted(set(before) | set(after)):
        old = before.get(name, '<unset>')
        new = after.get(name, '<unset>')
        if old != new:
            return '{} {} -> {}'.format(name, old, new)
    return 'no visible difference'


def field_provenance(outer_boundary, inner_boundaries=None, guides=None,
                     mode='perpendicular', target_length=None, guide_weight=1.0,
                     guide_band=None, relax=False, tau=None, symmetry='auto'):
    """**Everything that determines a field, canonicalised for comparison.**

    The signature is :meth:`CrossField.from_boundary`'s, defaults included, and
    that is load-bearing: :meth:`CrossField.mismatch` compares what a caller
    would have solved against what was solved, so a parameter that defaults
    differently in the two places would report a change that did not happen. **A
    parameter added to** ``from_boundary`` **has to be added here too.**

    ``'auto'`` symmetry is resolved before recording, because the group is what
    the solve actually used -- storing the string would make a field solved with
    ``symmetry=None`` compare equal to one solved under the full group.

    Each loop keeps its ROLE in the record. An earlier solve cache hashed holes
    and guides into one flat list, so the same polyline in either role compared
    equal. Holes and guides keep their given ORDER, because ``Symmetry.detect``
    and ``from_curves`` both consume them as sequences.

    Returns
    -------
    dict
        ``{'geometry': {'outer', 'inners', 'guides'}, 'params': {...}}``, all
        plain JSON types.
    """
    from compas_singular.framefield.constraints import as_curve_list
    from compas_singular.framefield.symmetry import Symmetry

    guides = as_curve_list(guides)
    if symmetry == 'auto':
        symmetry = Symmetry.detect(
            [outer_boundary] + list(inner_boundaries or []) + list(guides or []))

    params = (('mode', mode), ('target_length', target_length),
              ('guide_weight', guide_weight), ('guide_band', guide_band),
              ('relax', relax), ('tau', tau), ('symmetry', symmetry))
    return {
        'geometry': {
            'outer': _canonical_curve(outer_boundary),
            'inners': [_canonical_curve(loop) for loop in (inner_boundaries or [])],
            'guides': [_canonical_curve(curve) for curve in guides],
        },
        'params': dict((name, _plain(value)) for name, value in params),
    }


def wrap_to_period(delta):
    """Fold an angle difference into (-pi/4, +pi/4], the cross's half-period.

    Two crosses can never differ by more than 45 degrees -- past that they are
    better matched by the next arm round. Every angle comparison in this module
    goes through here.
    """
    return delta - PERIOD * round(delta / PERIOD)


class CrossField(Data):
    """A boundary-aligned cross field on a planar triangulation.

    Attributes
    ----------
    background : BackgroundMesh
    u : dict[int, complex]
        ``exp(i*4*theta)`` per vertex. Magnitude NOT normalised when
        ``relax=False``; exactly 1 everywhere when ``relax=True``.
    theta : dict[int, float]
        A representative cross angle per vertex, in (-pi/4, pi/4].
    relaxed : bool
        Which solver produced this field.
    iterations : int or None
        Diffusion/normalisation steps taken. ``None`` when ``relaxed`` is False.
    residual : float or None
        Final ``max|u_k+1 - u_k|``. Above ``tol`` means the iteration ran out of
        steps rather than converging.
    """

    def __init__(self, background, u, relaxed=False, iterations=None, residual=None):
        super(CrossField, self).__init__()
        self.background = background
        self.u = u
        self.theta = {vkey: phase(value) / 4.0 for vkey, value in u.items()}
        self.relaxed = relaxed
        self.iterations = iterations
        self.residual = residual
        self._singularities = None
        #: The guide curves this field was solved with, if any. Read by
        #: :meth:`densify` to decide whether element quality may be SPENT on
        #: alignment -- see ``densify._accepts``. Set by :meth:`from_boundary`.
        self.guides = []
        #: The symmetry group of the input, or ``None``. Set by
        #: :meth:`from_boundary`; carried so a caller that solves a field on its
        #: own is not holding something less than the decomposition would have.
        self.symmetry = None
        #: Where ``symmetry.snap_singularities`` moved each singularity, and its
        #: report. Needed by whoever builds a TRACING tracer -- point location
        #: does not care. Empty when no symmetry was applied.
        self.singularity_points = {}
        self.snap_report = {}
        #: What this field was solved FROM, canonicalised: the outer boundary,
        #: the holes, the guides and every solver parameter. ``None`` on a field
        #: built through :meth:`solve` directly, which never sees the
        #: boundaries. Set by :meth:`from_boundary`, written by
        #: :meth:`save_to_json`, and compared by :meth:`mismatch` -- a field
        #: read back off disk is worthless unless you can tell whether the
        #: domain moved under it.
        self.inputs = None
        #: Lazily built by :meth:`locator`, and deliberately not pickled.
        self._locator = None

    def __getstate__(self):
        """Everything but the locator.

        Applies to ``pickle`` and ``copy.deepcopy`` alike. The locator holds a
        bucket grid with an entry per background face -- pure derived state that
        rebuilds in milliseconds, inflates the copy, and would be silently stale
        if the background ever moved under it.
        """
        state = dict(self.__dict__)
        state['_locator'] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._locator = None

    # ------------------------------------------------------------------
    # the solve
    # ------------------------------------------------------------------

    @classmethod
    def solve(cls, background, constraints=None, relax=False, tau=None,
              tol=1e-9, max_iterations=20000):
        """Smoothest cross field satisfying ``constraints``.

        Minimises the Dirichlet energy of ``u`` over triangulation edges, subject
        to hard constraints (eliminated) and soft ones (least-squares penalty).

        Parameters
        ----------
        background : BackgroundMesh
        constraints : list[Constraint], optional
            Defaults to ``from_boundary(background)``, i.e. tangent to every wall.
        relax : bool, optional
            Enforce ``|u| = 1`` by diffusion + normalisation instead of leaving
            the magnitude free. See the module docstring: ``False`` (the default)
            is the original solver and is measurably wrong about where
            singularities go.
        tau : float, optional
            Diffusion time per step, ``relax`` only. Defaults to
            ``(bounding-box diagonal / 25)**2`` -- a property of the DOMAIN, not
            of the discretisation, which is the point: the steady state should
            not move when the background is refined.
        tol, max_iterations : optional
            Steady-state test on ``max|u_k+1 - u_k|``, ``relax`` only.

        Returns
        -------
        CrossField
        """
        if constraints is None:
            constraints = from_boundary(background)

        mesh = background.mesh
        vkeys = list(mesh.vertices())
        index = {vkey: i for i, vkey in enumerate(vkeys)}
        n = len(vkeys)

        hard = {}
        soft = {}
        for c in constraints:
            if c.vkey not in index:
                continue
            if c.weight is None:
                hard[index[c.vkey]] = representation(c.direction)
            else:
                # several soft constraints on one vertex accumulate, which is the
                # right behaviour where two guides overlap
                w, acc = soft.get(index[c.vkey], (0.0, 0j))
                soft[index[c.vkey]] = (w + c.weight, acc + c.weight * representation(c.direction))

        # uniform-weight graph Laplacian. The plan says uniform first; cotangent
        # weights are the textbook choice but go negative on obtuse triangles,
        # and the background mesh is built well-shaped precisely so this is enough.
        rows, cols, vals = [], [], []
        for u_key, v_key in mesh.edges():
            i, j = index[u_key], index[v_key]
            rows += [i, j, i, j]
            cols += [i, j, j, i]
            vals += [1.0, 1.0, -1.0, -1.0]
        for i, (w, _) in soft.items():
            rows.append(i)
            cols.append(i)
            vals.append(w)
        A = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr().astype(complex)

        b = np.zeros(n, dtype=complex)
        for i, (_, acc) in soft.items():
            b[i] = acc

        free = [i for i in range(n) if i not in hard]
        if not free:
            raise ValueError('every vertex is hard-constrained; nothing to solve')
        fixed = sorted(hard)
        x_fixed = np.array([hard[i] for i in fixed], dtype=complex)

        A_ff = A[free][:, free]
        A_fh = A[free][:, fixed]
        rhs = b[free] - A_fh.dot(x_fixed)

        x_free = spsolve(A_ff.tocsc(), rhs)

        iterations = residual = None
        if relax:
            x_free, iterations, residual = cls._relax(
                A, b, free, fixed, x_fixed, x_free, background, tau,
                tol, max_iterations)

        values = np.zeros(n, dtype=complex)
        values[free] = x_free
        values[fixed] = x_fixed

        return cls(background, {vkey: complex(values[index[vkey]]) for vkey in vkeys},
                   relaxed=bool(relax), iterations=iterations, residual=residual)

    @staticmethod
    def _relax(A, b, free, fixed, x_fixed, x_start, background, tau, tol,
               max_iterations):
        """Diffusion + normalisation to steady state -- Dai/Qiao/Wang Algorithm 4.1.

        One step is implicit Euler on the gradient flow of the same energy ``A``
        and ``b`` encode, followed by projection back onto ``|u| = 1``::

            (I + tau*A)_ff u_f  =  u_f^prev + tau*(b_f - A_fh x_h)
            u_f <- u_f / |u_f|

        Implicit rather than explicit so the step is unconditionally stable and
        ``tau`` can be chosen for convergence speed rather than to satisfy a CFL
        condition -- which is the practical content of the paper's Theorem 4.1
        for us. The matrix does not change between steps, so it is factorised
        ONCE and each iteration is a pair of triangular solves. Measured at
        ``target_length`` 0.35 on the disc: 0.24 s for 2366 iterations against
        0.017 s for the single solve it replaces.

        Started from the unrelaxed solve rather than from the boundary data --
        it already has the right winding, so the iteration only has to move the
        singularities, not find them.

        A NOTE ON ``tau``, because it looks like a model parameter and is not.
        The steady state is very nearly independent of it (disc radii 0.85 for
        tau 0.25, 0.5 and 1.0 at three resolutions); what tau sets is the
        convergence RATE. Running a fixed iteration budget instead of a
        steady-state test makes it look like a strong parameter -- radii 0.37 to
        0.85 -- because total diffusion time is ``iterations * tau`` and the
        small-tau runs simply had not finished. Do not import the paper's
        ``epsilon = tau = alpha*h`` scaling: their alpha ties the diffusion time
        to the diffuse-interface width, and with hard constraints there is no
        interface.
        """
        from scipy.sparse import identity

        if tau is None:
            xs = [p[0] for p in background.outer]
            ys = [p[1] for p in background.outer]
            diagonal = ((max(xs) - min(xs)) ** 2 + (max(ys) - min(ys)) ** 2) ** 0.5
            tau = (diagonal / 25.0) ** 2

        n = A.shape[0]
        step = (identity(n, format='csr', dtype=complex) + tau * A).tocsr()
        solve_ff = factorized(step[free][:, free].tocsc())
        offset = tau * (b[free] - A[free][:, fixed].dot(x_fixed))

        u = np.asarray(x_start, dtype=complex).copy()
        u /= np.maximum(np.abs(u), 1e-12)

        residual = float('inf')
        iterations = 0
        for iterations in range(1, max_iterations + 1):
            nxt = solve_ff(u + offset)
            nxt /= np.maximum(np.abs(nxt), 1e-12)
            residual = float(np.max(np.abs(nxt - u)))
            u = nxt
            if residual < tol:
                break
        return u, iterations, residual

    # ------------------------------------------------------------------
    # sampling
    # ------------------------------------------------------------------

    def directions(self, vkey):
        """The four arms of the cross at a vertex, as unit vectors."""
        t = self.theta[vkey]
        return [[cos(t + k * PERIOD), sin(t + k * PERIOD), 0.0] for k in range(4)]

    def magnitude(self, vkey):
        """``|u|``. Near 0 means the vertex is close to a singularity."""
        return abs(self.u[vkey])

    def angle_in_face(self, fkey, bary, reference=None):
        """Interpolate the cross angle inside a triangle.

        Interpolates ``u`` -- the complex representation -- and takes the angle
        of the result. NOT a blend of the three vertex ANGLES.

        The difference is not cosmetic. Each vertex angle is only defined modulo
        90 degrees, so blending them means first choosing a branch for each, and
        near a singularity the three genuinely differ by more than a quarter
        period: no branch choice is right and the blend jumps. Walking a circle
        around a pentagon's singularity that way, the unwrapped residual advanced
        by 2, 3 or 4 periods depending only on the probe radius, where the index
        says it must always be 5 -- so an arm went missing and the layout came out
        with a 5-sided patch, at some background resolutions and not others.

        ``u`` has no such ambiguity: it is a single-valued smooth complex field,
        linear interpolation is well defined, and ``arg(u)/4`` is continuous
        wherever ``u != 0``. The 4th power is doing exactly the job it was
        introduced for, and the branch matching disappears.

        Parameters
        ----------
        fkey : int
        bary : (float, float, float)
            Barycentric coordinates, in the order of ``face_vertices(fkey)``.
        reference : float, optional
            Arm to prefer, as an angle. Only selects which of the four arms is
            returned; it no longer affects the interpolation itself.

        Returns
        -------
        float
        """
        vkeys = self.background.mesh.face_vertices(fkey)
        value = 0j
        for w, vkey in zip(bary, vkeys):
            value += w * self.u[vkey]
        if value == 0:
            # exactly on the zero: fall back to the nearest vertex
            value = self.u[max(zip(bary, vkeys))[1]]
        theta = phase(value) / 4.0
        if reference is not None:
            theta += PERIOD * round((reference - theta) / PERIOD)
        return theta

    # ------------------------------------------------------------------
    # singularities
    # ------------------------------------------------------------------

    def singularities(self):
        """Singular faces and their indices, in quarter-turns.

        The index of a face is the winding of the cross angle around its three
        edges, each difference folded into the cross's half-period. The sum is
        necessarily a multiple of 90 degrees; dividing by 90 gives an integer
        number of quarter-turns -- ``+1`` where a quad mesh would put a valence-3
        vertex, ``-1`` for valence-5.

        Returns
        -------
        list[(int, int)]
            ``(fkey, index)`` for every face with a non-zero index.
        """
        if self._singularities is None:
            out = []
            mesh = self.background.mesh
            for fkey in mesh.faces():
                a, b, c = mesh.face_vertices(fkey)
                total = (wrap_to_period(self.theta[b] - self.theta[a])
                         + wrap_to_period(self.theta[c] - self.theta[b])
                         + wrap_to_period(self.theta[a] - self.theta[c]))
                k = int(round(total / PERIOD))
                if k:
                    out.append((fkey, k))
            self._singularities = out
        return self._singularities

    # ------------------------------------------------------------------
    # validation
    # ------------------------------------------------------------------

    def _oriented_boundary_loops(self):
        """Boundary loops with the domain on the left: outer CCW, holes CW."""
        mesh = self.background.mesh
        loops = []
        for loop in mesh.boundaries():
            pts = [mesh.vertex_coordinates(v) for v in loop]
            area2 = sum(p[0] * q[1] - q[0] * p[1]
                        for p, q in zip(pts, pts[1:] + pts[:1]))
            loops.append((abs(area2), area2, loop))
        loops.sort(reverse=True)                     # largest is the outer one
        out = []
        for i, (_, area2, loop) in enumerate(loops):
            want_ccw = (i == 0)
            is_ccw = area2 > 0
            out.append(loop if is_ccw == want_ccw else list(reversed(loop)))
        return out

    def poincare_hopf(self):
        """Check the interior indices against the boundary winding.

        The discrete argument principle: the sum of the interior singularity
        indices must equal the total winding of the field along the boundary,
        both counted in quarter-turns.

        This is the version to assert, NOT ``sum == 4 * chi``. That shortcut only
        holds for a SMOOTH boundary. On a polygon the corners absorb the turning
        -- a square's boundary-tangent cross field is the constant field, with a
        boundary winding of 0 and no interior singularities at all, because each
        90-degree corner turns the tangent by exactly the cross's period and so
        contributes nothing. Comparing against the measured winding is correct in
        both cases and is what actually catches a broken angle unwrapping.

        Returns
        -------
        dict
        """
        winding = 0.0
        for loop in self._oriented_boundary_loops():
            n = len(loop)
            for i in range(n):
                a = self.theta[loop[i]]
                b = self.theta[loop[(i + 1) % n]]
                winding += wrap_to_period(b - a)

        interior = sum(k for _, k in self.singularities())
        expected = int(round(winding / PERIOD))
        return {
            'interior_index_sum': interior,
            'boundary_winding': expected,
            'residual_radians': winding - expected * PERIOD,
            'singular_faces': len(self.singularities()),
            'ok': interior == expected,
        }

    def report(self):
        """One-line summary plus the Poincare-Hopf result.

        ``min_magnitude`` and ``mean_magnitude`` are only informative when
        ``relaxed`` is False -- the relaxation pins every magnitude to 1, so
        under ``relax=True`` they report 1.0 and say nothing except which solver
        ran. ``iterations`` and ``residual`` are the diagnostics that replace
        them there.
        """
        mags = [abs(v) for v in self.u.values()]
        info = dict(self.poincare_hopf())
        info['min_magnitude'] = min(mags)
        info['mean_magnitude'] = sum(mags) / len(mags)
        info['relaxed'] = self.relaxed
        info['iterations'] = self.iterations
        info['residual'] = self.residual
        return info

    # ------------------------------------------------------------------
    # standing on its own: solving from a domain, and densifying with it
    # ------------------------------------------------------------------

    @classmethod
    def from_boundary(cls, outer_boundary, inner_boundaries=None, guides=None,
                      mode='perpendicular', target_length=None, guide_weight=1.0,
                      guide_band=None, relax=False, tau=None, symmetry='auto'):
        """**A field over a domain, with no decomposition around it.**

        Everything :meth:`FieldDecomposition.from_boundary` does up to and
        including the solve -- background, constraints, solve, symmetrisation --
        and none of what comes after it. Separatrix tracing and the planar
        arrangement are the expensive half, and a caller who only wants to
        DENSIFY with the field has no use for either::

            field = CrossField.from_boundary(outer, guides=cables, target_length=0.5)
            dense = coarse.densification(field=field)

        which is what lets a layout from anywhere -- a skeleton decomposition, a
        hand-built mesh, one read back out of a document -- be densified with a
        field. ``FieldDecomposition.from_boundary`` calls this and then traces.

        Parameters are :meth:`FieldDecomposition.from_boundary`'s, and mean the
        same things; ``target_length`` is the BACKGROUND spacing, not the quad
        size.

        Returns
        -------
        CrossField
            With :attr:`guides`, :attr:`symmetry`, :attr:`singularity_points`
            and :attr:`snap_report` set, so nothing is left behind for a caller
            that later wants to trace with it.
        """
        # Imported here, not at module scope: ``trace`` imports
        # ``wrap_to_period`` from this module and ``symmetry`` imports this
        # class, so either at the top is a circular import.
        from compas_singular.framefield.background import BackgroundMesh
        from compas_singular.framefield.constraints import as_curve_list
        from compas_singular.framefield.constraints import from_curves
        from compas_singular.framefield.symmetry import Symmetry
        from compas_singular.framefield.symmetry import snap_singularities
        from compas_singular.framefield.symmetry import symmetrise

        guides = as_curve_list(guides)

        if symmetry == 'auto':
            symmetry = Symmetry.detect(
                [outer_boundary] + list(inner_boundaries or []) + list(guides or []))

        background = BackgroundMesh.from_boundary(
            outer_boundary, inner_boundaries, target_length=target_length,
            symmetry=symmetry)

        constraints = from_boundary(background)
        if guides:
            # ``mode`` is a no-op for a cross field and ``guides`` order no
            # longer is either -- see from_curves' notes.
            constraints += from_curves(background, guides, mode=mode,
                                       weight=guide_weight, band=guide_band)

        field = cls.solve(background, constraints, relax=relax, tau=tau)

        points = {}
        snap_report = {}
        if symmetry is not None and symmetry.enabled('field'):
            # Project the field onto the symmetric subspace. Cheap, and exact:
            # a 90 degree rotation acts TRIVIALLY on exp(i*4*theta) and every
            # reflection conjugates it, so the group average is arithmetic.
            field = symmetrise(field, symmetry)
        if symmetry is not None and symmetry.enabled('singularities'):
            # ... and then take each singularity off the arbitrary triangle it
            # was reported on. Only meaningful after the line above: it snaps to
            # the minimum of a field it assumes is already symmetric.
            field, points, snap_report = snap_singularities(field)

        field.guides = list(guides or [])
        field.symmetry = symmetry
        field.singularity_points = points
        field.snap_report = snap_report
        # ``symmetry`` is the RESOLVED group by now, never the string 'auto',
        # so what is recorded is what was actually solved under.
        field.inputs = field_provenance(
            outer_boundary, inner_boundaries, guides, mode=mode,
            target_length=target_length, guide_weight=guide_weight,
            guide_band=guide_band, relax=relax, tau=tau, symmetry=symmetry)
        return field

    # ------------------------------------------------------------------
    # provenance -- has the domain moved under this field?
    # ------------------------------------------------------------------

    def mismatch(self, outer_boundary, inner_boundaries=None, guides=None,
                 mode='perpendicular', target_length=None, guide_weight=1.0,
                 guide_band=None, relax=False, tau=None, symmetry='auto'):
        """**Why this field does not belong to these inputs.** ``None`` if it does.

        Takes what :meth:`from_boundary` takes. Use it on a field read back off
        disk, before densifying with it::

            field = CrossField.load_from_json('field.json')
            why = field.mismatch(outer, holes, guides=cables, target_length=0.5)
            if why:
                field = CrossField.from_boundary(outer, holes, guides=cables,
                                                 target_length=0.5)

        A field is a function of its domain. Densifying a layout with a field
        solved for a DIFFERENT outline is not an error anywhere downstream --
        every patch interior still integrates, the mesh still welds, the result
        still passes the quality gate -- it is simply aligned to a shape that is
        no longer there. Nothing else in the pipeline can catch that, which is
        the whole reason the inputs are stored.

        Returns
        -------
        str or None
            A one-line description of the FIRST difference found, or ``None``
            when everything that determines the field agrees.
        """
        if self.inputs is None:
            return ('this field carries no record of its inputs -- it was built '
                    'through CrossField.solve rather than from_boundary, so '
                    'whether it matches cannot be answered')

        other = field_provenance(
            outer_boundary, inner_boundaries, guides, mode=mode,
            target_length=target_length, guide_weight=guide_weight,
            guide_band=guide_band, relax=relax, tau=tau, symmetry=symmetry)

        mine = self.inputs.get('geometry') or {}
        theirs = other['geometry']
        labels = {'outer': 'outer boundary', 'inners': 'holes',
                  'guides': 'guide curves'}
        for role in ('outer', 'inners', 'guides'):
            if mine.get(role) != theirs.get(role):
                return 'the {} changed'.format(labels[role])

        if (self.inputs.get('params') or {}) != other['params']:
            return 'solver settings changed -- {}'.format(
                _describe_difference(self.inputs.get('params'), other['params']))
        return None

    def matches(self, *args, **kwargs):
        """Whether this field was solved from these inputs. See :meth:`mismatch`.

        ``False`` when the field carries no provenance at all -- unknown is not
        the same as equal, and the safe reading of "cannot tell" is "re-solve".
        """
        return self.mismatch(*args, **kwargs) is None

    # ------------------------------------------------------------------
    # serialisation
    # ------------------------------------------------------------------

    def save_to_json(self, filepath, pretty=False):
        """**Write this field to a JSON file.**

        The coarse layout has been serialisable all along
        (``CoarseQuadMesh.to_json``); this is the other half, so a whole
        field-driven job survives being closed::

            decomposition.get_field().save_to_json('field.json')
            decomposition.decomposition_mesh().to_json('coarse.json')

            # ... a week later ...
            field = CrossField.load_from_json('field.json')
            coarse = CoarsePseudoQuadMesh.from_json('coarse.json')
            coarse.collect_strips()
            coarse.set_strips_density_target(0.5)
            dense = coarse.densification(field=field)

        Written with ``compas.json_dump``, so the background triangulation goes
        through COMPAS's own encoder and comes back as the right ``Mesh``
        subclass rather than as a dict.

        **Derived state is not stored, and that is deliberate.** :attr:`theta`
        is recomputed from :attr:`u` by ``__init__``; the point locator is a
        bucket grid that rebuilds in milliseconds. Both would inflate the file,
        and neither can disagree with what it came from if it is never written.
        :attr:`_singularities` IS stored, because ``snap_singularities``
        overwrites it and it is not recoverable from ``u`` afterwards.

        Parameters
        ----------
        filepath : str
        pretty : bool, optional
            Indent the JSON. Off by default -- ``u`` has one entry per
            background vertex, so a plate at 0.3 spacing is several thousand
            lines nobody reads.

        Returns
        -------
        str
            ``filepath``, so a save can be chained or logged.
        """
        import os

        import compas

        # Same two differences from a bare ``json_dump`` that
        # ``Mesh.save_to_json`` has from ``to_json``: make the parent directory,
        # and hand the path back.
        folder = os.path.dirname(os.path.abspath(filepath))
        if folder and not os.path.isdir(folder):
            os.makedirs(folder)
        compas.json_dump(self.__jsondata__(), filepath, pretty=pretty)
        return filepath

    @property
    def __data__(self):
        """compas ``Data``: the same payload as :meth:`save_to_json`, so a field
        can travel inside a larger JSON document (a session) and come back
        exactly."""
        return self.__jsondata__()

    @classmethod
    def __from_data__(cls, data):
        return cls.__from_jsondata__(data)

    def __jsondata__(self):
        """The payload :meth:`save_to_json` writes: plain types plus one Mesh."""
        symmetry = None
        if self.symmetry is not None:
            symmetry = {'centre': [float(c) for c in self.symmetry.centre],
                        'names': list(self.symmetry.names()),
                        'steps': list(self.symmetry.steps)}

        return {
            'type': JSON_TYPE,
            'version': JSON_VERSION,
            'inputs': self.inputs,
            'background': {
                'mesh': self.background.mesh,
                'outer': [[float(c) for c in p] for p in self.background.outer],
                'inners': [[[float(c) for c in p] for p in loop]
                           for loop in self.background.inners],
                'target_length': float(self.background.target_length),
            },
            # ``complex`` is the one type here JSON has no opinion about. Python
            # writes floats with ``repr``, so re/im round-trip EXACTLY -- which
            # matters, because ``21_edit_coarse.py`` asserts a field is
            # unchanged to exactly zero.
            'u': {str(vkey): [value.real, value.imag]
                  for vkey, value in self.u.items()},
            'relaxed': bool(self.relaxed),
            'iterations': self.iterations,
            'residual': self.residual,
            'singularities': (None if self._singularities is None else
                              [[int(fkey), int(index)]
                               for fkey, index in self._singularities]),
            'guides': [[[float(c) for c in p] for p in curve]
                       for curve in (self.guides or [])],
            'symmetry': symmetry,
            'singularity_points': {
                str(fkey): [float(c) for c in point]
                for fkey, point in (self.singularity_points or {}).items()},
            'snap_report': self.snap_report or {},
        }

    @classmethod
    def load_from_json(cls, filepath, default=None):
        """**Read a field back.** The inverse of :meth:`save_to_json`.

        Mirrors ``Mesh.load_from_json``, ``default`` included, so a caller
        resuming a job can ask for both halves the same way::

            field = CrossField.load_from_json(field_path, default=None)
            coarse = CoarsePseudoQuadMesh.load_from_json(coarse_path)
            if field is None or field.mismatch(outer, holes, guides=cables):
                field = CrossField.from_boundary(outer, holes, guides=cables)

        Does NOT check that the field still belongs to your domain -- that needs
        the domain, which this does not have. Ask :meth:`mismatch` straight
        afterwards, every time.

        Parameters
        ----------
        filepath : str
        default : optional
            Returned when the file does not exist, instead of raising. For a
            caller resuming from a save that may never have been written.

        Returns
        -------
        CrossField

        Raises
        ------
        ValueError
            If the file is not a field, or was written by a newer version than
            this one understands. Both fail here, with a sentence, rather than
            several frames downstream with a ``KeyError`` about a vertex.
        """
        import os

        import compas

        if not os.path.isfile(filepath):
            return default
        return cls.__from_jsondata__(compas.json_load(filepath))

    @classmethod
    def __from_jsondata__(cls, data):
        """Rebuild from :meth:`__jsondata__`'s payload."""
        from compas_singular.framefield.background import BackgroundMesh
        from compas_singular.framefield.symmetry import ELEMENTS
        from compas_singular.framefield.symmetry import Symmetry

        found = data.get('type') if isinstance(data, dict) else type(data).__name__
        if found != JSON_TYPE:
            raise ValueError(
                'not a CrossField file: expected type {!r}, found {!r}'.format(
                    JSON_TYPE, found))
        version = data.get('version')
        if version is None or version > JSON_VERSION:
            raise ValueError(
                'this file is version {!r}; this compas_singular understands up '
                'to {}. Re-solve the field, or update.'.format(
                    version, JSON_VERSION))

        bg = data['background']
        background = BackgroundMesh(bg['mesh'], bg['outer'], bg['inners'],
                                    bg['target_length'])

        # JSON object keys are always strings. This is the trap
        # ``PseudoQuadMesh.__from_data__`` exists for and ``test_pipeline.py``
        # pins: a vertex key left as '17' looks fine until something indexes the
        # background mesh with it.
        u = {int(vkey): complex(value[0], value[1])
             for vkey, value in data['u'].items()}

        field = cls(background, u, relaxed=data.get('relaxed', False),
                    iterations=data.get('iterations'),
                    residual=data.get('residual'))

        singularities = data.get('singularities')
        if singularities is not None:
            field._singularities = [(int(fkey), int(index))
                                    for fkey, index in singularities]
        field.guides = [[list(p) for p in curve]
                        for curve in data.get('guides') or []]
        field.singularity_points = {
            int(fkey): list(point)
            for fkey, point in (data.get('singularity_points') or {}).items()}
        field.snap_report = data.get('snap_report') or {}
        field.inputs = data.get('inputs')

        symmetry = data.get('symmetry')
        if symmetry is not None:
            by_name = {element[5]: element for element in ELEMENTS}
            field.symmetry = Symmetry(
                symmetry['centre'],
                [by_name[name] for name in symmetry['names']],
                symmetry['steps'])
        return field

    def locator(self):
        """**Point location on this field's background.** Built once, cached.

        A ``Tracer`` has two jobs and this is the second one. Its headline job is
        separatrix integration -- launching field lines out of singularities and
        walking them across the domain -- which needs ``singularity_points``,
        ``step`` and ``max_length``. Densification uses none of that. What it
        uses is ``locate(point, hint) -> (fkey, barycentric)``: a patch-interior
        node has no vertex on the background mesh, so before
        :meth:`angle_in_face` can say what the field is doing there, something
        has to find the containing triangle. That is a bucket grid over face
        bounding boxes and nothing else, so it is derivable from the field alone
        -- measured bit-identical against a decomposition's own tracer.

        Named ``locator`` rather than ``tracer`` for exactly that reason: the
        object can trace, but nothing here asks it to.
        """
        if self._locator is None:
            from compas_singular.framefield.trace import Tracer
            self._locator = Tracer(self, singularity_points=self.singularity_points)
        return self._locator

    def densify(self, coarse, edges_to_curves=None, spend=None, **kwargs):
        """**Densify a coarse layout with this field steering patch interiors.**

        Reached as ``coarse.densification(field=field)``; call it directly to get
        the statistics back as well.

        ``edges_to_curves`` already makes the coarse EDGES follow whatever they
        came from; this makes the patch INTERIORS follow the field too, instead
        of a bilinear blend of their own boundaries that never consults it. A
        guide can then steer the inside of a patch without splitting the patch
        first -- the only way a guide reaches a mesh on a domain, like a plain
        square, where it forces no topology at all.

        Parameters
        ----------
        coarse : CoarseQuadMesh or CoarsePseudoQuadMesh
            Densities already set. Strips are collected on demand.
        edges_to_curves : dict, optional
            ``{(u, v): polyline}``, as ``densification`` takes it.
        spend : bool, optional
            Whether element quality may be spent to satisfy a guide. Defaults to
            whether this field HAS guides, which is the rule
            ``FieldDecomposition`` applies: there is no alignment to buy on an
            unguided domain, so the Coons interior was already right.
        **kwargs
            ``stiffness``, ``guard``, ``iterations``, ``tikhonov`` -- passed to
            :func:`densify.field_densification`.

        Returns
        -------
        (QuadMesh, dict)
            The dense mesh, and the per-patch counts. The mesh is also set on
            ``coarse``, so ``coarse.get_quad_mesh()`` keeps working.
        """
        from compas_singular.framefield.densify import field_densification

        if spend is None:
            spend = bool(self.guides)
        dense, stats = field_densification(
            coarse, self, self.locator(), edges_to_curves=edges_to_curves,
            spend=spend, **kwargs)
        return dense, stats
