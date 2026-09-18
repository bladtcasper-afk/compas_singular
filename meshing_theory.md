# Meshing theory — how the two coarse-mesh generation methods work

This document covers the two independent routes this codebase has for turning a planar boundary
(with optional holes, guide curves and point features) into a coarse quad-mesh layout: the
**frame-field method** (`CrossField`, current default for curved/guided domains) and the
**skeleton-decomposition method** (`SkeletonDecomposition`, the original medial-axis-based route).
Both converge on the same downstream object — a `CoarsePseudoQuadMesh` whose faces are Coons-patched
by `densification()` — but they derive singularity placement and patch topology by entirely different
means. Part I explains the frame field in detail; Part II explains the skeleton method and contrasts
the two. See [field vs skeleton](../17_vs_skeleton.py)-style comparisons in the examples folder for a
quantitative A/B between them.

---

# Part I — `CrossField`: the frame-field method

**Source:** [`framefield/field.py`](../../../src/compas_singular/framefield/field.py)
**Depends on:** [`framefield/background.py`](../../../src/compas_singular/framefield/background.py) (the triangulation it solves on),
[`framefield/constraints.py`](../../../src/compas_singular/framefield/constraints.py) (boundary/guide directions)
**Consumed by:** `trace.py` (streamline launch and integration), `arrangement.py`,
`densify.py` (field-aware patch interiors)

This note explains `CrossField` line by line against the theory it is an implementation of: the
4-symmetry-direction (N-RoSy, N=4) representation of Ray et al., the Dirichlet-energy field design of
Bommes/Zimmer/Kobbelt, and the diffusion-generated / Ginzburg–Landau relaxation of
Merriman–Bence–Osher → Viertel & Osting (2019) → Dai/Qiao/Wang (2026). The companion review
[`paper-diffusion-generated-crossfields.md`](../paper-diffusion-generated-crossfields.md) covers *which
parts of the 2026 paper were worth adopting*; this document explains *what the code that resulted
actually does and why it is correct*.

## 1. What problem this class solves

A **cross** at a point is an unordered set of 4 unit directions, 90° apart — the tangent frame you
want a quad mesh's edges to follow. `CrossField` computes one cross per vertex of a triangulated
planar domain such that:

1. every boundary vertex's cross is tangent to the wall (or free, where two walls disagree — see
   §4), and
2. the field is as smooth as possible everywhere else, i.e. neighbouring crosses differ as little as
   the topology of the domain allows.

Condition 2 cannot always be met with *zero* disagreement — a domain with non-zero Euler
characteristic (a disc, say) forces the field to wind around, and where it winds a full turn faster
or slower than the mesh combinatorics can absorb without a defect, a **singularity** appears. Locating
those singularities correctly is the entire point of this module: they are where a quad mesh must put
a valence-3 or valence-5 vertex, and they are what `trace.py` fires separatrices from to cut the
domain into paintable patches.

## 2. The 4-fold symmetric representation

A cross has no preferred one of its four arms — rotating it by any multiple of 90° gives the same
object. That symmetry is the whole reason this is a linear problem instead of a combinatorial one.

`constraints.representation()` ([constraints.py:33](../../../src/compas_singular/framefield/constraints.py#L33)) encodes a
direction `d = (cos θ, sin θ)` as

```
u = exp(i · 4θ)
```

Because `exp(i·4(θ + kπ/2)) = exp(i·4θ)` for any integer `k`, **all four arms of the same cross map to
the same complex number.** This is the standard *N-RoSy* (N-direction rotationally symmetric)
representation from Ray, Vallet, Li & Lévy, *N-Symmetry Direction Fields on Surfaces* (ACM TOG, 2008):
raising the direction's angle to the `N`-th power (here `N = 4`) is exactly the operation that turns
"a direction defined up to a `2π/N` rotation" into "a single-valued complex number," which is what
converts a combinatorial branch-matching problem (which of the 4 arms of vertex `i`'s cross lines up
with which arm of vertex `j`'s?) into a plain smoothness functional on `u ∈ ℂ`.

`wrap_to_period()` ([field.py:59](../../../src/compas_singular/framefield/field.py#L59)) is the inverse operation used whenever
the code needs a genuine *angle* difference rather than a complex one: it folds a difference into
`(−π/4, π/4]`, i.e. picks whichever of the two crosses' four arms are nearest each other. Every angular
comparison in the module — the winding sum in `singularities()`, the boundary winding in
`poincare_hopf()` — goes through this function, which is precisely what makes those sums well defined
despite the underlying value only existing modulo 90°.

`CrossField.__init__` stores both representations: `u` (the complex field, what everything is computed
from) and `theta = phase(u) / 4` (one representative arm, for sampling — see §5).

## 3. Field design as a linear solve — the Dirichlet energy

### 3.1 The energy

The standard formulation for a smooth N-RoSy field (Ray et al. 2008; Bommes, Zimmer & Kobbelt,
*Mixed-Integer Quadrangulation*, SIGGRAPH 2009 for the quad-meshing application) is: minimise the
Dirichlet energy of `u` over the mesh,

```
E(u) = Σ_{edges (i,j)} w_ij · |u_i − u_j|²
```

subject to hard directional constraints at some vertices. On a curved surface this requires a
*connection* — parallel-transporting `u_j` into `i`'s tangent frame before differencing, because
neighbouring tangent planes are different vector spaces. **This module is planar only**, so every
tangent space is the same copy of world `XY`, the connection is the identity, and the energy above is
literally what `CrossField.solve` assembles — the docstring's opening lines say as much explicitly.
This is also why `background.face_basis` is a stub that returns the identity: it exists purely so that
a future surface case would not need an API change, without doing any work today.

### 3.2 Weights: uniform, not cotangent

[field.py:153–166](../../../src/compas_singular/framefield/field.py#L153-L166) builds `w_ij` as a **uniform** graph Laplacian
(every mesh edge gets weight 1), not the cotangent weights that are the textbook choice for a
discrete Dirichlet energy (Pinkall & Polthier, 1993, in the general discrete-differential-geometry
literature this field design descends from). Cotangent weights go negative on obtuse triangles, which
turns the energy indefinite; the code's stated position is that `background.py` builds a well-shaped
mesh (densified boundary plus a jittered near-regular interior grid) specifically so that uniform
weights are enough and the more fragile weighting is not needed.

### 3.3 Hard vs. soft constraints

Each `Constraint(vkey, direction, weight)` from `constraints.py` is either:

- **hard** (`weight is None`): the vertex's `u` is fixed to `representation(direction)` and removed
  from the unknowns entirely (`hard` dict, [field.py:146](../../../src/compas_singular/framefield/field.py#L146)).
- **soft** (`weight` a float): contributes a least-squares penalty term `w·|u_i − target|²`, which adds
  `w` to the diagonal of the Laplacian at that vertex and `w·target` to the right-hand side
  ([field.py:148–151](../../../src/compas_singular/framefield/field.py#L148-L151)). Several soft constraints on the same vertex
  simply accumulate — this is standard least-squares combination, and it is what lets two guide curves
  disagree near each other without either one breaking the solve; the field reports their compromise
  rather than failing.

### 3.4 Elimination and solve

The assembled system is partitioned into free and fixed (hard-constrained) unknowns
([field.py:172–182](../../../src/compas_singular/framefield/field.py#L172-L182)):

```
A_ff · u_free = b_free − A_fh · u_fixed
```

This is textbook **Dirichlet boundary condition elimination** for a discrete Laplace problem: fixing
`u_fixed` and moving its contribution to the right-hand side reduces the system to one over the free
vertices only, which `scipy.sparse.linalg.spsolve` solves directly (`A` is complex, sparse, and — with
uniform positive weights plus at least one hard constraint — positive definite on the free block).

This whole solve (`relax=False`, the default) is minimising `E(u)` with **`|u|` completely
unconstrained**. That single fact is the subject of the rest of this document.

## 4. `from_boundary` and the 45° corner problem

Before the solve can run it needs constraints, and the default set —
`constraints.from_boundary()` ([constraints.py:44](../../../src/compas_singular/framefield/constraints.py#L44)) — pins every
boundary vertex tangent to the wall. A boundary vertex sits between two wall segments, hence has two
candidate tangents; they are combined by **averaging the two `representation()` values**, not the raw
directions.

This is not an arbitrary choice of averaging space — it is forced by the geometry of a right-angle
corner. Represented as raw unit vectors, two edges 90° apart average to something 45° off from either,
which is wrong: a square's corner cross should be axis-aligned, matching both walls, not diagonal.
Represented in the 4th-power complex space, a right-angle corner's two edge directions map to **the
same complex number** (since a 90° rotation is exactly the cross's period), so their average is that
same number, exactly — the corner is pinned to the one cross both walls actually agree on.

Where the two walls disagree by something other than a multiple of 90° — most sharply, an interior
angle near 45° or 135°, where the two tangents are exactly antipodal in the 4-fold representation —
the average magnitude shrinks toward zero and the vertex is left **unconstrained**
(`corner_tolerance`, [constraints.py:78](../../../src/compas_singular/framefield/constraints.py#L78)) rather than pinned to an
arbitrary compromise direction. The field then decides that vertex for itself, which typically means a
singularity nearby — the honest outcome, and the mechanism behind the `±1/4` boundary-defect behaviour
the Dai/Qiao/Wang curvature catalogue documents at acute and right-angle corners (see
[paper-diffusion-generated-crossfields.md §4.3](../paper-diffusion-generated-crossfields.md)).

## 5. Sampling the field

### 5.1 `directions(vkey)`

Reconstructs the 4 literal arms from `theta` by adding the 4 multiples of `PERIOD = π/2`
([field.py:262](../../../src/compas_singular/framefield/field.py#L262)) — the inverse of the representation map, used only for
drawing/consumption, never for computation.

### 5.2 `angle_in_face` — why it interpolates `u`, not `θ`

This is the one piece of the class that looks like it could be simplified and specifically must not
be. [field.py:271–314](../../../src/compas_singular/framefield/field.py#L271-L314) interpolates the three vertices' **complex**
`u` values with barycentric weights and only then takes `phase(value) / 4`:

```
value = Σ w_k · u_k        (linear interpolation of a single-valued field)
theta = phase(value) / 4
```

The alternative — average the three vertex *angles* — requires first choosing, for each vertex, which
of its four equivalent arms to use, i.e. picking a branch. Near a singularity the three vertices'
crosses genuinely differ by more than a quarter turn from each other, so **no consistent branch choice
exists**, and any per-triangle branch pick makes the interpolated field disagree with the true winding
around a probe loop. The documented failure mode: walking a circle around a pentagon's singularity
with the branch-averaged interpolation, the unwrapped residual advanced by 2, 3 or 4 periods depending
only on the probe radius, where the topological index says it must always be 5 — an arm went missing
and the recovered patch layout came out 5-sided at some background resolutions and not others.

`u` has no such ambiguity because it is smooth and single-valued everywhere `u ≠ 0` (this is precisely
what the 4th-power map was introduced to buy in §2), so **linear interpolation of `u` is unambiguous
and correct**, and `arg(u)/4` is continuous on the interior of any triangle not containing a zero. The
`reference` parameter only chooses which of the four equivalent arms of the *result* to report — it
plays no part in the interpolation itself, and changing it cannot reintroduce the branch bug.

The one degenerate case, `value == 0` (the sample point sits exactly on a zero of the field, i.e. at a
singularity), falls back to the nearest vertex's own `u`, since there is no well-defined direction at
an actual zero.

## 6. Singularities and the discrete argument principle

### 6.1 The winding index

`singularities()` ([field.py:320](../../../src/compas_singular/framefield/field.py#L320)) computes, for every triangular face,
the sum of the three folded angle differences around its boundary:

```
total = wrap(θ_b − θ_a) + wrap(θ_c − θ_b) + wrap(θ_a − θ_c)
```

Because the three differences are each folded into `(−π/4, π/4]`, their sum can only ever land on a
multiple of `π/2` (this is a discrete instance of the argument principle: the total turning of a
locally-defined angular field around a closed loop is quantised). `k = round(total / PERIOD)` is the
face's index in quarter-turns: `+1` is where a quad mesh needs a valence-3 vertex, `−1` a valence-5
vertex — the sign convention used throughout the rest of this module (`trace.py`'s
`target = 4 − index` for how many separatrix arms to launch, and the site's red/blue dot convention
for the two cases).

This is deliberately **not** the paper's approach (locating singularities as local maxima of a 3×3
discrete Dirichlet-energy neighbourhood, Dai/Qiao/Wang §5.2) and the review in
`paper-diffusion-generated-crossfields.md` §5.2 argues why: the winding computation gives an *exact
integer index* directly from the combinatorics, with no threshold and no neighbourhood search, whereas
an energy-maximum heuristic gives a location only and needs the index inferred separately. The exact
version is strictly stronger and it is what a grid-free field (this one) can afford that a grid-based
field cannot.

### 6.2 Poincaré–Hopf, correctly stated for a polygon

`poincare_hopf()` ([field.py:369](../../../src/compas_singular/framefield/field.py#L369)) is the sanity check on all of the
above: **the sum of interior singularity indices must equal the total winding of the field measured
along the boundary**, both in quarter-turns. This is the discrete argument principle applied globally
instead of per-face.

The docstring is explicit that the check must be phrased this way and *not* as `Σk = 4·χ` (Euler
characteristic), which is the version usually quoted for a **smooth** boundary. On a polygon the
corners themselves absorb turning: a square's boundary-tangent cross field is the exactly-constant
field (zero interior singularities, zero boundary winding), because each 90° corner turns the wall
tangent by exactly one cross period and therefore contributes nothing to the winding sum — even though
a smooth convex curve of the same total turning (2π) would force `Σk = 4`. Comparing the interior index
sum against the *measured* boundary winding is correct in both the smooth and polygonal cases and is
what actually catches a broken angle-unwrapping bug, where the naive `4·χ` shortcut would not.

`_oriented_boundary_loops()` ([field.py:352](../../../src/compas_singular/framefield/field.py#L352)) supplies the sign
convention this needs: outer loop CCW, hole loops CW (domain always on the left), determined from the
shoelace-formula signed area of each loop — necessary because the winding sum's sign depends on
consistent orientation and a mesh's raw boundary-loop traversal order does not guarantee it.

## 7. The two solvers, and the theory that separates them

This is the part of the module the docstring is most insistent about, because the original
implementation's own comments used to overclaim what it was doing.

### 7.1 `relax=False` (default) — Dirichlet energy, `|u|` free

Just §3's linear solve, with no constraint at all on `|u|`. Nothing stops the solver from paying for a
singularity by letting `|u|` sag smoothly across a broad interior region instead of concentrating a
sharp defect at a point — and nothing in the energy rewards concentration over spreading. Measured
directly: `min|u|` on the test suite's disc is `0.003`–`0.010`, i.e. the field is close to collapsing
outright, and the winding computed from a wide near-zero blob has no reason to centre where the true
singularity is.

This is the `ε → ∞` limit of the **Ginzburg–Landau functional**

```
E_ε(u) = ∫ |∇u|² + (1/4ε²)(1 − |u|²)² dΩ
```

— the standard regularisation used to give a vector/cross field well-localised point defects (see
Beaufort, Lambrechts, Henrotte, Geuzaine & Remacle, *Computing cross fields — A PDE approach based on
the Ginzburg–Landau theory*, Procedia Engineering 2017, for exactly this functional applied to cross
fields for quad meshing). Dropping the second term (`ε → ∞`) leaves only the Dirichlet term, which is
what `relax=False` solves. The module's docstring used to describe the resulting `|u| → 0` collapse
*as* the Ginzburg–Landau relaxation; it is not — it is what happens in that functional's **absence**,
which is the symptom of a relaxation that was never implemented, not the relaxation itself.

Measured against the one external number available — Beaufort et al.'s and Dai/Qiao/Wang's (2026,
§6.2) published radius for the four `+1/4` singularities of a boundary-tangent field on the unit disc,
`r ≈ 0.85` — this solver gives roughly half that (`r ≈ 0.45`–`0.52`) and, worse, **drifts further
inward as the background mesh is refined**, with no sign of a limit. Refining a discretisation is
supposed to converge a PDE solution toward the true answer, not move it monotonically away; this is
the clearest evidence that `relax=False` is not solving the constrained problem the theory specifies,
only a related but different one.

### 7.2 `relax=True` — diffusion-generated / MBO relaxation

`_relax()` ([field.py:197](../../../src/compas_singular/framefield/field.py#L197)) adds back exactly the missing ingredient: it
enforces `|u| = 1` explicitly, by alternating

```
(I + τA)_ff u_f = u_f^prev + τ·(b_f − A_fh x_h)     (implicit diffusion step)
u_f ← u_f / |u_f|                                    (pointwise projection back to the unit circle)
```

to steady state (`max|u_{k+1} − u_k| < tol`). This is the **Merriman–Bence–Osher (MBO) threshold
dynamics scheme** — originally devised for mean curvature flow by alternating a heat/diffusion step
with a pointwise thresholding/projection step — generalised to `S¹`-valued (cross/vector) fields by
Viertel & Osting, *An approach to quad meshing based on harmonic cross-valued maps and the Ginzburg–
Landau theory*, SIAM J. Sci. Comput. 2019, and specialised into the exact iteration used here (implicit
Euler diffusion + projection, run to steady state) as Algorithm 4.1 of Dai, Qiao & Wang, *An efficient
and stable diffusion generated method for quadrilateral mesh generation in general domains*, 2026.
Diffusion-then-project is the **splitting-scheme realisation** of the Ginzburg–Landau gradient flow in
the `ε → 0` limit: the diffusion step descends the Dirichlet term, the projection step is the limiting
behaviour of the double-well penalty term as it becomes infinitely stiff. This module deliberately
takes only that iteration from the 2026 paper — not its diffuse-domain/FFT machinery, which solves a
different problem (avoiding a triangulation) this module does not have; see
[paper-diffusion-generated-crossfields.md §1](../paper-diffusion-generated-crossfields.md) for that
scoping decision.

Three implementation details, each with a specific reason:

- **Implicit, not explicit, Euler.** Implicit Euler is unconditionally stable, so `τ` can be chosen
  for convergence *speed* rather than to satisfy a CFL-type stability bound — stated by the code as
  "the practical content of Theorem 4.1 for us" (the paper's unconditional-energy-decay result). The
  system matrix `(I + τA)_ff` does not change between iterations, so it is LU-factorised once
  (`scipy.sparse.linalg.factorized`) and every subsequent step is a pair of triangular solves —
  measured at 0.24 s for 2366 iterations at `target_length = 0.35`, against 0.017 s for the one-shot
  solve it replaces, i.e. still negligible next to `trace`/`arrangement`/`repair`.
- **Started from the unrelaxed solve, not from raw boundary data.** `x_start` already has the correct
  topological winding (`relax=False`'s solve got that part right — only the interior *placement* of the
  defect is wrong), so the iteration only has to relocate existing singularities rather than discover
  their number and location from scratch.
- **`τ` is a convergence rate, not a model parameter, and this took a real measurement to establish.**
  Run to a fixed iteration budget instead of a steady-state test, the disc's singularity radius appears
  to depend strongly on `τ` (0.37 to 0.85) — but this is an artefact: total diffusion time is
  `iterations × τ`, and the small-`τ` runs simply had not finished within the budget. Run to an actual
  steady state (`tol = 1e-9`), the radius is `0.85 ± 0.05` at `τ ∈ {0.25, 0.5, 1.0}` and at three mesh
  resolutions — i.e. essentially `τ`-independent, which is what a genuine steady-state relaxation
  should be. The default `τ = (bounding-box diagonal / 25)²` is deliberately a property of the
  **domain's size**, not of the discretisation, precisely so the steady state does not drift as
  `background` is refined — the opposite of `relax=False`'s behaviour in §7.1. The module explicitly
  warns against importing the paper's own `ε = τ = α·h` scaling, since that couples diffusion time to a
  diffuse-interface width this hard-constrained formulation does not have.

Measured result: `r ≈ 0.85` at three resolutions and three `τ` values — matching the published value
and, unlike `relax=False`, **not drifting under refinement**. `singularities()` is unaffected by which
solver ran, because it only ever reads the winding of `θ`, never `|u|`; `angle_in_face` stays exactly
as well defined, since `u` never revisits zero away from an actual singularity once `|u|` is pinned to
1 everywhere else.

### 7.3 Why the default is still `relax=False`

This is a genuine, measured trade-off, not an oversight. Fixing the field moves every singularity on a
curved or guided domain **outward, toward the wall** — correct per the theory, and directly observable
on the disc, where patch count becomes resolution-independent (5/5/9/9 → 5/5/5/5). But
`trace.py`'s separatrix-launch probe radii and `repair.py`'s landing-cluster tolerances were tuned
against the *old*, systematically-too-far-inward field. Moving singularities close to the wall breaks
those downstream assumptions on some domains before it helps them — the suite's own A/B is 18 rows
improved, 14 regressed, all of the movement confined to the 8 curved/guided domains (the 32 purely
rectilinear rows are byte-identical either way, since their field is exactly constant and the
relaxation is a provable no-op there). Turning `relax=True` on by default would move every curved-wall
layout in the committed baseline, so it stays opt-in until the downstream tolerances are widened to
match a correctly-placed field — see `paper-diffusion-generated-crossfields.md` §3.3 and §7 for the
full measurement and the planned next step (a soft boundary-weight sweep, §4.4 of that document).

## 8. Summary — theory-to-code map

| Theory | Where in this codebase |
|---|---|
| N-RoSy / 4-fold symmetric complex representation (Ray et al. 2008) | `constraints.representation()`, `CrossField.u` |
| Dirichlet-energy field design with hard/soft directional constraints (Bommes/Zimmer/Kobbelt 2009 lineage) | `CrossField.solve`, §3 |
| Uniform graph Laplacian in place of cotangent weights | [field.py:153–166](../../../src/compas_singular/framefield/field.py#L153-L166), §3.2 |
| Discrete argument principle / winding index per face | `CrossField.singularities`, §6.1 |
| Poincaré–Hopf, polygon-correct form | `CrossField.poincare_hopf`, §6.2 |
| Ginzburg–Landau functional, `ε → ∞` (no penalty) | `solve(relax=False)`, §7.1 — the default, and known inaccurate |
| Ginzburg–Landau functional, `ε → 0` via MBO threshold dynamics (Viertel & Osting 2019; Dai/Qiao/Wang 2026 Alg. 4.1) | `solve(relax=True)` / `_relax`, §7.2 |
| Branch-free interpolation via the single-valued 4th-power field | `angle_in_face`, §5.2 |

**Primary references**

- Ray, N., Vallet, B., Li, W.-C., Lévy, B. — *N-Symmetry Direction Fields on Surfaces*, ACM TOG 27(2), 2008.
- Bommes, D., Zimmer, H., Kobbelt, L. — *Mixed-Integer Quadrangulation*, ACM SIGGRAPH 2009.
- Beaufort, P.-A., Lambrechts, J., Henrotte, F., Geuzaine, C., Remacle, J.-F. — *Computing cross fields
  — A PDE approach based on the Ginzburg–Landau theory*, Procedia Engineering 203, 2017.
- Viertel, R., Osting, B. — *An approach to quad meshing based on harmonic cross-valued maps and the
  Ginzburg–Landau theory*, SIAM J. Sci. Comput. 41(1), 2019.
- Dai, J., Qiao, Z., Wang, D. — *An efficient and stable diffusion generated method for quadrilateral
  mesh generation in general domains*, arXiv:2605.27854v1, 2026 — reviewed in full in
  [paper-diffusion-generated-crossfields.md](../paper-diffusion-generated-crossfields.md).
- Merriman, B., Bence, J., Osher, S. — *Diffusion generated motion by mean curvature*, 1992 — origin of
  the alternate-diffuse-then-project (MBO) scheme this module's `_relax` specialises.

---

# Part II — `SkeletonDecomposition`: the medial-axis method

**Source:** [`algorithms/decomposition.py`](../../../src/compas_singular/algorithms/decomposition.py),
[`datastructures/skeleton/skeleton.py`](../../../src/compas_singular/datastructures/skeleton/skeleton.py)
**Depends on:** [`algorithms/triangulation.py`](../../../src/compas_singular/algorithms/triangulation.py)
(`boundary_triangulation`, the Delaunay mesh it reads structure from)
**Consumed by:** `CoarsePseudoQuadMesh.densification()` — the same Coons-patch densifier `CrossField`'s
output eventually reaches, via `arrangement.py` on that route.

This is the codebase's original coarse-mesh route, and it predates the frame field entirely. It never
solves a field or minimises an energy: singularity placement and patch topology both fall directly out
of the combinatorics of a Delaunay triangulation of the boundary. Where Part I is a PDE discretised on
a mesh, this is discrete geometry with no continuous relaxation anywhere in the pipeline.

## 1. What problem this class solves — and how it differs from the field

Same target as `CrossField`: a coarse quad-mesh layout — vertices, edges, singular faces — for a
planar domain, ready for Coons-patch densification. The difference is entirely in *how* the
singularities and patch boundaries are found:

- `CrossField` **places** singularities by solving for a smooth tangent field and reading off its
  winding (Part I §6).
- `SkeletonDecomposition` **discovers** singularities as a structural feature of the domain's Delaunay
  triangulation — specifically, points equidistant from three boundary features (the medial axis, Blum
  1967) — with no solve involved.

Both are legitimate: the medial axis is the classical shape-skeleton construction, and a
boundary-tangent smooth cross field's singularities are known to sit near medial-axis branch points on
simple domains. The methods diverge sharply once the domain has smooth curvature or guide curves,
which is exactly where §5 below draws the comparison out.

## 2. Stage 1 — triangulate the boundary only

`boundary_triangulation()`
([triangulation.py:23](../../../src/compas_singular/algorithms/triangulation.py#L23)) builds a Delaunay
triangulation using **only the boundary, hole, and feature points as vertices** — critically, no
interior points are inserted, unlike `CrossField`'s `background.py`, which populates the interior with
a jittered near-regular grid specifically so a *field* has somewhere to be smooth over.

```python
vertices = [pt for boundary in [outer_boundary] + inner_boundaries + polyline_features for pt in boundary] + point_features
faces = delaunay(vertices)
```

Three cleanup passes follow: degenerate (zero-area, collinear) faces are deleted; faces whose
circumcentre falls outside the outer boundary or inside a hole are deleted (this is what keeps the
triangulation conforming to a non-convex or multiply-connected domain, since a raw Delaunay
triangulation of a point set is always convex); and the mesh is topologically cut — unwelded, not
geometrically split — along any `polyline_features`, so that a feature curve can act as an internal
seam the skeleton will route along.

With every vertex on a boundary, the circumcentres of adjacent triangles trace a discrete
approximation of the domain's **medial axis** — the set of points with more than one closest boundary
point. This is the load-bearing geometric fact the entire rest of the class depends on, and it is the
reason no interior points are added: adding one would pull nearby circumcentres off the true medial
axis and toward the inserted point instead.

## 3. Stage 2 — read singularities and skeleton branches off the triangulation

All of `skeleton.py` is pure combinatorics on the Delaunay mesh — no coordinates are solved for,
only classified:

- **singular faces** = triangles with exactly 3 neighbours, i.e. interior triangles with no boundary
  edge ([skeleton.py:44](../../../src/compas_singular/datastructures/skeleton/skeleton.py#L44)). Each
  one becomes a singularity of the coarse mesh; its position is the triangle's circumcentre
  (`trimesh_face_circle`), not an arbitrary choice — it is the medial-axis point equidistant from that
  triangle's three defining boundary features.
- **corner faces** = triangles with exactly 1 neighbour
  ([decomposition.py:81](../../../src/compas_singular/algorithms/decomposition.py#L81)) — these sit at
  domain corners and are excluded from the singularity set.
- **skeleton lines** = segments joining the circumcentres of every pair of adjacent triangles
  ([skeleton.py:68](../../../src/compas_singular/datastructures/skeleton/skeleton.py#L68)), chained into
  polyline **branches** by `network_polylines`, which walks a valence-2 chain until it hits a node of
  different valency or an explicit split point.

Because singularity count and location are read directly off how many triangles happen to end up
3-valent, they are sensitive to **boundary sampling density** in the same way `CrossField`'s
singularities are sensitive to the 45° turn-limit heuristic on curved arcs (see
[framefield-45-degree-turn-limit](../../../../memory/framefield-45-degree-turn-limit.md) in project
memory) — a different mechanism arriving at the same class of fragility: both methods can spray extra
singularities on a domain whose curvature the discretisation samples coarsely.

## 4. Stage 3 — assemble a patch layout, then patch the patcher

`decomposition_polylines()`
([decomposition.py:163](../../../src/compas_singular/algorithms/decomposition.py#L163)) collects branch
polylines from three structural sources and three corrective ones, then re-chains everything through
`network_polylines` once more (splitting additionally at domain corners) to get the final polyline
arrangement:

| source | what it contributes |
|---|---|
| singularity → singularity | skeleton branches whose neither end is a corner circumcentre |
| singularity → boundary | circumcentre of each singular face to each of its 3 triangle vertices ([decomposition.py:132](../../../src/compas_singular/algorithms/decomposition.py#L132)) — these are new branches, not part of the topological skeleton itself, and are what split a singular triangle's neighbourhood into paintable quads |
| boundary | each boundary loop, split at corner vertices and at split vertices (the triangle vertices touched by a singular face) ([decomposition.py:145](../../../src/compas_singular/algorithms/decomposition.py#L145)) |
| *fix* — boundary kinks | a boundary curvature spike the skeleton's discretisation missed at low density gets a forced split ([decomposition.py:305](../../../src/compas_singular/algorithms/decomposition.py#L305), threshold `relative_kink_angle_limit = π/8`) |
| *fix* — collapsed boundaries | a boundary loop with fewer than 3 splits would degenerate to a 1- or 2-sided "polygon"; extra splits are injected at ⅓ and ⅔ (or symmetric offsets) ([decomposition.py:225](../../../src/compas_singular/algorithms/decomposition.py#L225)) |
| *fix* — flipped faces | a singularity-to-singularity branch that turns more than `flip_angle_limit = π/2` in total would fold a patch back on itself; it is subdivided at intermediate singular-face edges instead ([decomposition.py:265](../../../src/compas_singular/algorithms/decomposition.py#L265)) |

These three correction passes exist purely to guarantee that the polyline set forms a clean planar
subdivision — they are the fragile, threshold-tuned part of the method, directly analogous to
`CrossField`'s `corner_tolerance` and separatrix-launch tolerances on the field route: both methods need
a set of hand-tuned thresholds to convert an idealised construction (medial axis; smooth field) into
something robust enough for a real, possibly noisy, boundary discretisation.

## 5. Stage 4 — polylines to a coarse quad mesh, with pole handling

`decomposition_mesh(poles)`
([decomposition.py:200](../../../src/compas_singular/algorithms/decomposition.py#L200)) does the actual
face-finding by delegating to COMPAS core's `CoarsePseudoQuadMesh.from_polylines` — the enclosed
regions of the polyline arrangement become faces — then repairs what that construction cannot itself
guarantee is quad-clean:

- **`solve_triangular_faces()`**
  ([decomposition.py:343](../../../src/compas_singular/algorithms/decomposition.py#L343)) converts
  every triangular face (a sign the polyline arrangement produced a degenerate quad) into a proper quad
  or removes it, by cause: a triangle with **1** boundary vertex means a singularity sits at that
  vertex's position too, so the vertex is duplicated and the triangle becomes a real quad; a triangle
  with **2** boundary vertices means two singularities have coincided at the same point, so the face is
  deleted and the two vertices are welded together at the midpoint of the connecting decomposition
  polyline. Any resulting zero-length edge is then nudged 10% toward its neighbourhood centroid so
  `densification()` never has to divide by an edge length of zero.
- **`split_quads_with_poles(poles)`**
  ([decomposition.py:433](../../../src/compas_singular/algorithms/decomposition.py#L433)) is this
  method's analogue of `CrossField`'s point-feature handling: each `point_feature` you pass in is a
  **pole** — a quad that touches it is split across its diagonal at that vertex into two pseudo-quad
  (triangle-with-a-doubled-vertex) faces. The face→pole map is written to
  `mesh.attributes['face_pole']`, read straight back by `densification()`
  ([mesh_quad_pseudo_coarse.py:70–73](../../../src/compas_singular/datastructures/mesh_quad_pseudo_coarse/mesh_quad_pseudo_coarse.py#L70-L73))
  to insert the `None` edge that turns a Coons-patch quad into a fan around the pole.

## 6. Contrast with `CrossField` — where each one wins

Both routes bottom out at the same `CoarsePseudoQuadMesh` → `densification()` handoff (Part I is
consumed via `arrangement.py`'s patch tracing; this route builds the coarse mesh directly), so the
comparison is entirely about how faithfully each one reads the *shape* of the domain:

- **No solve, no convergence question.** The skeleton method has nothing analogous to Part I §7's
  `relax=False` vs `relax=True` trade-off — there is no energy to minimise and therefore no question of
  whether a discretisation converges to the "true" singularity location. This is a strength (nothing to
  mistune) and a weakness (nothing to relax toward a better answer, either): the medial axis *is* the
  answer, for better or worse.
- **Sensitive to boundary sampling, not to a turn-angle heuristic.** Where `CrossField` sprays
  singularities on a smoothly curved arc because a sampled turn exceeds its 45° corner tolerance (see
  [framefield-45-degree-turn-limit](../../../../memory/framefield-45-degree-turn-limit.md)), the
  skeleton method sprays singularities on the same arc because a finely sampled circular boundary
  produces many near-equidistant Delaunay triangles along it — both are discretisation artefacts of an
  underlying smooth shape, arrived at through unrelated mechanisms.
- **No notion of a guide curve.** `CrossField`'s soft constraints (Part I §3.3) let a guide curve steer
  the field without hard-pinning it — there is no equivalent lever here; the only way to steer the
  skeleton method is to add a `polyline_feature`, which is unwelded into the triangulation as a hard
  topological seam, not a soft directional hint.
- **Branching / limbed domains are the skeleton method's natural strength.** The medial axis of a
  domain with narrow limbs or T-junctions (a floor plan, a branching vault) produces exactly the
  skeleton branches you would sketch by hand; the frame field has to discover the equivalent topology
  through winding, which needs the field to actually resolve a singularity at each branch point rather
  than reading it off structure directly. This is consistent with
  [field vs skeleton](../../../../memory/framefield-vs-skeleton.md) in project memory: field wins on 15
  of 18 compared domains, and the skeleton method's one win is on exactly this kind of branching case.

## 7. Summary — theory-to-code map

| Concept | Where in this codebase |
|---|---|
| Medial axis via Delaunay circumcentres (Blum 1967) | `Skeleton.lines`, `Skeleton.singular_points` |
| Boundary-only Delaunay triangulation, non-convex domain via circumcentre culling | `boundary_triangulation`, §2 |
| Singularity = interior (3-neighbour) Delaunay triangle | `Skeleton.singular_faces`, §3 |
| Skeleton branch chaining (valence-2 walk with splits) | `graph_polylines` (`compas.datastructures`) |
| Polyline arrangement → coarse quad faces | `CoarsePseudoQuadMesh.from_polylines`, §4 |
| Degenerate-triangle repair (coincident singularities/vertices) | `solve_triangular_faces`, §5 |
| Point feature as a pole (fan face around a marked vertex) | `split_quads_with_poles`, §5 — shares the `face_pole` attribute contract with the field route |

**Primary reference**

- Blum, H. — *A transformation for extracting new descriptors of shape*, in *Models for Perception of
  Speech and Visual Forms*, 1967 — origin of the medial-axis / topological-skeleton construction this
  module implements via Delaunay circumcentres.