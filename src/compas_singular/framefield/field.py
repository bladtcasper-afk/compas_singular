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

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import factorized
from scipy.sparse.linalg import spsolve

from .constraints import PERIOD
from .constraints import from_boundary
from .constraints import representation


__all__ = ['CrossField', 'wrap_to_period']


def wrap_to_period(delta):
    """Fold an angle difference into (-pi/4, +pi/4], the cross's half-period.

    Two crosses can never differ by more than 45 degrees -- past that they are
    better matched by the next arm round. Every angle comparison in this module
    goes through here.
    """
    return delta - PERIOD * round(delta / PERIOD)


class CrossField(object):
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
        self.background = background
        self.u = u
        self.theta = {vkey: phase(value) / 4.0 for vkey, value in u.items()}
        self.relaxed = relaxed
        self.iterations = iterations
        self.residual = residual
        self._singularities = None

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
