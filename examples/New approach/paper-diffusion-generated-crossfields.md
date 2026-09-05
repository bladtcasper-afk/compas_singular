# Diffusion-generated cross fields — what is worth taking

**Paper:** J. Dai, Z. Qiao, D. Wang, *An efficient and stable diffusion generated method for
quadrilateral mesh generation in general domains*, arXiv:2605.27854v1, 27 May 2026.
**Date of this review:** 2026-08-05
**Question:** which parts of it are worth implementing in `framefield/`, and do they fix anything
that is actually broken or inaccurate here?

---

## 1. Verdict

The paper's *headline* contribution — moving the whole solve onto a regular background grid via a
diffuse domain method so the field can be computed by FFT — is **not** worth taking. It buys
efficiency we do not need and costs a rewrite of every module downstream of `background.py`.

The paper's *supporting* contribution — the **unit-norm relaxation**, i.e. alternating diffusion
with pointwise normalisation instead of a single Dirichlet-energy solve — **is** worth taking, and
it is not a refinement. Measured on a disc, our current field puts its four singularities in the
wrong place, and the error grows as the background is refined. Adding the normalisation step fixes
it, agrees with the published value, and costs 0.2 s. §3 is the measurement.

That split matters because the two are separable. The paper couples them (both live inside their
Algorithm 4.1), but the normalisation is what corrects the field; the FFT/diffuse-domain machinery
only changes where the field is *stored*. We can take the first without the second, and §3 verifies
that on our own triangulation.

| Paper element | § | Our code | Verdict | Effort |
|---|---|---|---|---|
| Unit-norm relaxation: diffusion **+ normalisation** | 4, Alg 4.1 | `framefield/field.py` | **IMPLEMENTED** — `solve(relax=True)`, opt-in. Field fixed; mesh 18 rows better / 14 worse (§3.3) | done |
| Sub-cell singularity localisation (zero of the interpolated representation) | 5.2 | `trace.py:_singularity_points` | **Take — cheap and exact** | ~20 lines |
| Curvature / corner behaviour catalogue | 6.3, Fig 11 | `15_baseline.py`, `constraints.from_boundary` | **Take as a test oracle** | 6 new domains |
| Boundary alignment as a *penalty*, not a hard BC | 2.3, 3.1 | `constraints.from_boundary(weight=…)` | **Try — the parameter already exists** | 1 line + a sweep |
| Harmonic map per block instead of Coons | 5.3 | `densify.py` | **Try on unguided patches** | ~60 lines |
| Boundary-data extension by closest point + mollification | 5.1, Alg 5.1 | `constraints.from_boundary` | Only if the penalty route is taken | — |
| Diffuse domain + FFT on a regular grid | 3.1, 5.1 | `background.py` + everything downstream | **Do not take** | rewrite |
| Singularity detection by 3×3 energy local max | 5.2 | `field.singularities` | **Do not take — ours is better** | — |
| Unconditional energy decay theorem | 4, Thm 4.1 | — | Not a motivation for us (§5.3) | — |
| Streamline tracing, branch selection by angular deviation | 2.2, 5.2 | `trace.direction_at` | Already identical — confirmation | — |

---

## 2. Stage-by-stage comparison

| Stage | Paper | `framefield/` |
|---|---|---|
| Domain representation | Phase field `φ = G_ε * χ` on a periodic square `Ω₂ ⊃ Ω₁`; no mesh | Delaunay triangulation of `Ω₁` itself, boundary densified + jittered interior grid (`background.py`) |
| Boundary condition | Diffuse penalty `ε⁻³B(φ)|u − g|²`, `g` extended off the wall by closest-point projection + Gaussian mollification | **Hard** constraint, eliminated from the system, one per boundary vertex (`constraints.from_boundary`, `weight=None`) |
| Representation | `u ∈ R²`, `|u| = 1`, i.e. `exp(4iθ)` | Identical: `u = exp(4iθ)`, one complex unknown per vertex |
| Solve | Iterate `u ← normalise(e^{τΔ}u + (τ/ε³)B(φ)g)` to steady state; FFT per step | **One** sparse solve of the graph Laplacian, `|u|` left free (`CrossField.solve`) |
| Singularities | Local maxima of a discrete Dirichlet energy on 3×3 cell neighbourhoods, then the zero of the bilinearly interpolated representation inside the cell | Exact winding index per triangle, in quarter-turns (`CrossField.singularities`); position = **face centroid** |
| Separatrix launch | Points on the cell boundary where a cross branch aligns with the radial direction from the singularity | Level crossings of the unwrapped residual on probe circles, 18 (radius, sample-count) attempts (`Tracer.launch_directions`) |
| Tracing | Branch chosen to minimise angular deviation from the incoming direction; stop on boundary or near a singularity | Same rule, same stopping conditions (`Tracer.direction_at`) |
| Layout repair | Not discussed | `arrangement.py` + `repair.py` — the bulk of this module |
| Block filling | Discrete harmonic map: two Laplace solves per block, finite differences on the logical square | `discrete_coons_patch`, plus field-aware relaxation on guided patches (`densify.py`) |
| Interior directional constraints (cables, force lines) | **None. The paper has no mechanism for them.** | `constraints.from_curves`, and `densify.py` carries them into the mesh |

The last row is the one to keep in view: the paper is a *domain* mesher. Our objective — cables and
force lines as field inputs — is outside its scope, so nothing in it can be adopted wholesale.

---

## 3. The measurement that decides the ranking

The paper reports (§6.2, and consistent with Beaufort et al. [3]'s renormalized-energy analysis for
circular domains) that a boundary-aligned cross field on the unit disc puts its four `+1/4`
singularities at radius **r ≈ 0.85**. They reach it two independent ways: an MBO reference solution
(Fig 9) and their own method as `ε, τ → 0` (Fig 10a). They also record the failure signature when
the discretisation takes over: below `α ≈ 1/150` the radii **collapse to r ≈ 0.2**.

That is an external number we can check ourselves against, and until now we had none — `poincare_hopf`
only verifies the field against *itself*.

### 3.1 Our field, measured

`CrossField.solve` on the suite's `DISC` (radius 5, so radii below are normalised by 5):

| boundary pts | `target_length` | singularities | radii |
|---|---|---|---|
| 48 | 0.50 | 4 × `+1` | 0.520 0.520 0.561 0.567 |
| 48 | 0.35 | 4 × `+1` | 0.386 0.495 0.495 0.504 |
| 48 | 0.25 | 4 × `+1` | 0.340 0.424 0.435 0.435 |
| 96 | 0.50 | 4 × `+1` | 0.520 0.520 0.561 0.567 |
| 96 | 0.35 | 4 × `+1` | 0.386 0.495 0.495 0.504 |
| 96 | 0.25 | 4 × `+1` | 0.281 0.365 0.412 0.412 |

Three things are wrong here, and only the first is obvious:

1. **The radius is ~0.45, not 0.85.** The singularities sit at roughly half the distance from the
   centre that the model says they should.
2. **It drifts inward under refinement** — 0.52 → 0.50 → 0.42 → and further at 96 points. There is
   no sign of a limit. Refining the background is supposed to converge the field, not move it.
3. **The four radii are not equal**, on a shape with 4-fold symmetry (0.386 against 0.504), and the
   spread does not shrink with refinement either.

Doubling the boundary sampling changes nothing at `tl` 0.50 and 0.35 — the numbers are byte-identical
— so this is not boundary discretisation. It is the field model.

**Why.** `CrossField.solve` minimises `∫|∇u|²` with `|u|` completely unconstrained. That is neither
the constrained problem (1.1) nor the Ginzburg–Landau relaxation (1.2); it is the `ε → ∞` end, where
the penalty term is absent altogether. With no term holding `|u|` near 1, the cheapest way to pay for
a singularity is to let `|u|` sag across a broad interior region rather than to concentrate a defect —
so what we detect is the winding around a diffuse low-magnitude blob, not a point, and the blob's
centroid has no reason to sit where the theory puts the singularity. `field.py`'s docstring says
`|u| → 0` "is the Ginzburg–Landau relaxation". It is the *symptom* of a relaxation we never
implemented; the measured `min |u|` is 0.003–0.010, i.e. the field is close to collapsing.

### 3.2 The paper's fix, on our triangulation

Prototype: keep `background.py`, keep the same graph Laplacian, keep the **hard** boundary
constraints, and replace the single solve with the paper's Algorithm 4.1 iteration — implicit
diffusion `(I + τL)u_{k+1} = u_k`, then pointwise normalisation — run to steady state
(`max|u_{k+1} − u_k| < 1e-10`). No diffuse domain, no FFT, no phase field.

| `target_length` | vertices | τ | iterations | radii |
|---|---|---|---|---|
| 0.50 | 378 | 0.25 | 1716 | 0.847 0.847 0.853 0.867 |
| 0.50 | 378 | 0.50 | 2456 | 0.847 0.847 0.853 0.856 |
| 0.50 | 378 | 1.00 | 1673 | 0.776 0.781 0.801 0.847 |
| 0.35 | 699 | 0.25 | 3508 | 0.845 0.901 0.909 0.909 |
| 0.35 | 699 | 0.50 | 2366 | 0.845 0.877 0.877 0.901 |
| 0.35 | 699 | 1.00 | 13638 | 0.845 0.877 0.877 0.901 |
| 0.25 | 1328 | 0.25 | 5170 | 0.808 0.821 0.848 0.849 |
| 0.25 | 1328 | 0.50 | 3046 | 0.808 0.821 0.848 0.849 |
| 0.25 | 1328 | 1.00 | 2545 | 0.808 0.848 0.849 0.853 |

**r = 0.85 ± 0.05, at every resolution and every τ.** The published value, reproduced, and
resolution-independent — which is exactly what the one-shot solve is not.

Two caveats, both real:

- **τ is not free, it is a convergence rate.** An earlier run with `τ = αh²` and a fixed 400-iteration
  budget appeared to show a strong τ-dependence (radii 0.37 to 0.85). It was under-convergence: total
  diffusion time is `iterations × τ`, and the small-τ rows had not finished. Run to steady state the
  τ-dependence almost disappears. Do not read the paper's `ε = τ = αh` scaling as applying here —
  their α couples the diffusion time to the *diffuse-interface width*, and with hard boundary
  constraints there is no interface.
- **It does not converge quickly.** 1.7k–13.6k iterations. But the matrix is prefactorised once
  (`scipy.sparse.linalg.factorized`), so each iteration is a triangular solve: **0.24 s against
  0.017 s** for the one-shot solve at `tl = 0.35`, and 0.12 s at a `1e-6` tolerance. Against the
  cost of `trace` + `arrangement` + `repair` this is free.

### 3.3 Implemented, and what the mesh did

Both of the above are now in the repo: `CrossField.solve(relax=True)` in `field.py`, checked by
[19_field_accuracy.py](19_field_accuracy.py) (`--winding`, `--timing`). **`relax` defaults to False.**

The field is fixed. The **mesh is not uniformly better**, which is what §7's step 2 was there to find
out. `15_baseline.py --relax` runs the same 64 rows under the new solver; A/B against the same code
state with it off:

- **32 rows identical** — every rectilinear plate (square, L, T, U, plus, comb, rect-slot,
  square+sq hole). Their field is exactly constant, so the relaxation is a no-op. That is the
  cheapest available proof it is not perturbing anything it should not.
- **32 rows moved** — the 8 curved or guided domains. **18 improved on minimum angle, 14 regressed.**

| domain | verdict | evidence |
|---|---|---|
| `disc` | **better** | patch count becomes resolution-**independent**: 5/5/9/9 → 5/5/5/5. Max angle 139→132, 142→126 at bg 0.6/0.5. Min angle mixed (+10, +5, −5, −5) |
| `square+ring cable` | **better** — the suite's worst domain | bg 0.5: min 7.4→18.1, max **180.0→149.7**, AR 18.1→4.9, poles 14→3, patches 43→21. bg 0.4: route `polygon`→`field`, max 180.0→137.3. Two of its three hard-floor breaches clear |
| `pentagon`, `stadium`, `square+cable` | neutral | wiggles of a few degrees either way |
| `hexagon` | **worse** | poles 0/0/0/0 → 1/2/3/0, min angle 67.9→50.4, 65.3→24.9, 62.4→17.8 |
| `ellipse` | **worse at bg 0.5** | patches 10→4, poles 0→7, min 65.3→14.3 |
| `disc+round hole` | mixed, and already bad on both | bg 0.6 min 1.05→7.05; bg 0.5 7.88→2.37; bg 0.4 falls back to `triangulation` |

**The regressions have one mechanism, and it is the predicted one.** Every domain that got worse
gained **poles** — the marker that `solve_non_quad_faces` could not make the layout all-quad without
collapsing a side. Measuring how far the singularities sit from the wall, in units of the background
spacing `h`:

| domain | bg | `relax=False` | `relax=True` |
|---|---|---|---|
| hexagon | 0.5 | 6.9h, 6.2h | **1.8h**, 3.1h |
| hexagon | 0.4 | 8.7h, 8.3h | 3.1h, 3.1h |
| ellipse | 0.5 | 0.9h, 0.8h, 4.0h, 4.0h | **0.4h**, **0.4h**, 2.8h, 3.6h |
| ellipse | 0.4 | 1.6h, 2.1h, 4.4h, 4.8h | 1.6h, 1.1h, 1.4h, 1.8h |

The relaxation moves singularities **outward**, toward the wall — on the disc that is the whole point,
and it is what the literature says should happen. At 0.4h from the wall, `launch_directions`' probe
circles (radii 0.45h to 2.0h) partly leave the domain, and a separatrix arm is short enough for
`repair._cluster` to merge its landing with a neighbour's. Those tolerances were tuned against fields
whose singularities sat 4–8h inside the domain.

So the honest reading: **the old field's error was protecting the repair stage.** Singularities
placed too far inside the domain make easier layouts. Fixing the field moves the bottleneck
downstream to `trace` and `repair`, where the wall-proximity assumptions are now the limit. The paper
says as much in its own terms — §2.3 notes the renormalized energy has no intrinsic mechanism to
balance "stay away from each other" against "stay away from the boundary", and §6.3 shows curvature
concentration driving singularities onto the wall.

**This is why `relax` ships off by default**, and why the next move is §4.4 (a boundary weight, which
pushes singularities back *inward* in a controlled way) rather than flipping the default.

---

## 4. Worth implementing

### 4.1 Unit-norm relaxation in `CrossField.solve` — **highest value**

*Paper §4, Algorithm 4.1.* Replace the single solve with diffusion + normalisation, as prototyped in
§3.2. Keep hard boundary constraints; the diffuse-domain half of the paper is not needed for this.

- **Fixes:** singularity placement error (~2× on the disc), its drift under refinement, and the loss
  of symmetry on symmetric domains.
- **Touches:** `field.py` only. `singularities()` uses the winding index, not `|u|`, so it is
  unaffected. `angle_in_face` interpolates `u`, which stays well defined at unit modulus.
  `magnitude()` and `report()`'s `min_magnitude` / `mean_magnitude` become meaningless (all 1.0) —
  they are used only by `01_field.py`'s printout, so replace them with the iteration count and
  residual.
- **Keep the old path.** Add it as an option (`CrossField.solve(..., relax=True)`) and default it off
  until the baseline says otherwise, so the 64 committed rows stay interpretable while it is evaluated.
- **Verify:** disc radius → 0.85 at three resolutions; `poincare_hopf` still `ok` on all 16 domains;
  then `15_baseline.py` diff, and expect movement — the layouts *will* change on curved domains.
- **Risk:** medium, and it is downstream, not in the solve. See §3.3.

### 4.2 Sub-cell singularity localisation

*Paper §5.2.* They locate the singularity as the zero of the interpolated representation inside the
cell. We use `mesh.face_centroid(fkey)` (`trace.py:160`).

On a triangle with `u` interpolated linearly this is a closed form, not a search: find barycentric
`(b₁,b₂,b₃)` with `Σbᵢ = 1` and `Σbᵢuᵢ = 0` — two real equations plus the constraint, a 3×3 solve.
Roughly 20 lines.

- **Fixes:** `launch_directions` walks probe circles centred on the singularity. If the centre is off
  by a fraction of a triangle, the circles are eccentric and the unwrapped residual is noisier than it
  needs to be — which is what the 18-attempt `(scale, samples)` ladder at `trace.py:200` exists to
  survive, and what the docstring's "5 arms at 0.50 and 0.35 but only 4 at 0.60, 0.45, 0.40" records.
  It also improves the coarse patch corner position, which `repair`'s clustering tolerances act on.
- **Verify:** `02_separatrices.py` (arm counts = `4 − index`) should pass on the first `(scale, samples)`
  attempt more often; instrument which attempt succeeds, before and after. `15_baseline.py` unchanged
  or improved.
- **Risk:** low. Falls back to the centroid when the system is singular.
- **Note:** worth doing **after** 4.1, not before — normalising `u` changes where the zero is.

### 4.3 The curvature / corner catalogue as a test oracle

*Paper §6.3, Fig 11.* They take `Ω₁ = [−1.5, 1.5] × [−0.5, 1]` and replace the top edge with six
profiles, three smooth and three with a `C⁰` corner:

| profile | reported behaviour |
|---|---|
| `y = −√(1 − x²)` | two interior `−1/4` singularities |
| `y = 4x²` | the same two, closer to the top boundary |
| `y = 16x²` | closer still — curvature concentration draws them to the wall |
| `y = 4\|x\|` (acute) | the two merge **at the corner** into an effective boundary singularity of index `−1/2` |
| `y = \|x\|/√3` | — |
| `y = \|x\|` (right angle) | a boundary singularity persists, index `+1/4` |
| obtuse corner | an additional **interior** `+1/4` appears |

This is free to adopt and it tests machinery we have repeatedly got wrong. Specifically it exercises
`constraints.from_boundary`'s `corner_tolerance` rule — the one that leaves a vertex free when the two
wall tangents cancel in the 4th-power representation, i.e. near 45° and 135° interior angles — against
published behaviour, at exactly the angles where that rule fires. Our suite has no acute wedge and no
concentrated-curvature case at all; the hexagon's 120° corners are the closest thing.

- **Index convention when writing the assertions:** paper `±1/4` = our `±1` quarter-turn
  (`+1` ↔ valence 3). Their merged `−1/2` boundary defect would be our `−2`, and it sits *on* the
  boundary, where `CrossField.singularities` — which scans faces — will not see it as such. Expect to
  compare against `poincare_hopf`'s `boundary_winding` rather than the interior sum.
- **Also worth noting from their §6.3:** the paper reports the smooth→`C⁰` transition is *not* captured
  by standard Ginzburg–Landau asymptotics. So this is an empirical catalogue, not a theorem — treat a
  disagreement as a question, not automatically as our bug.
- **Where:** add the six outlines to `10_domains.py` and rows to `15_baseline.py`, plus a singularity
  signature assertion alongside the existing `--corners` check.

### 4.4 Boundary alignment as a penalty rather than a hard constraint

*Paper §2.3 closing, §3.1.* Their explicit criticism of model (1.1): "imposing a Dirichlet boundary
condition is often overly restrictive for mesh generation, since strict boundary constraints are not
essential for determining interior singularity locations and separatrix structures." Their `ε⁻³B(φ)`
penalty makes wall alignment tunable, and §6.2 shows the knob moving the disc's singularity radius
monotonically from 0.85 down to ~0.48.

**We already have this parameter and have never used it.** `constraints.from_boundary(background,
weight=w)` takes a float; `field.py` routes finite weights to the least-squares path. Default is
`None` = hard, everywhere.

- **Why it might matter here:** the memory record says pole count drives quality and is *not*
  resolution-independent (ellipse min angle spread 41.5° across backgrounds, entirely the one row with
  poles). A knob that moves singularities relative to the wall is a lever on exactly the failure class
  that `repair`'s tolerances keep being tuned against — and it is a lever in the *model*, which is the
  right place, rather than another tolerance.
- **Cost:** a sweep, not an implementation. Run `15_baseline.py` over boundary weights and record
  patches / poles / min angle.
- **Caveat:** with a soft weight the field is no longer exactly tangent at the wall, so separatrices
  can approach the boundary at an angle — the grazing-hit case `trace.py` guards. Coarse nodes are
  projected onto the boundary in `repair`/`arrangement`, so exact tangency is less load-bearing than it
  looks, but this needs checking rather than assuming.
- **Do this after 4.1.** Under the current one-shot model, softening the boundary just adds a second
  unconstrained direction to a field that is already not solving the right problem.

### 4.5 Harmonic block filling for unguided patches — **demoted, see §8**

*Measured after this was written: our rounded domains do not lose their quality in patch interiors,
so this is not the fix for them. Kept because the argument still stands on its own; dropped down the
order in §7.*


*Paper §5.3.* Each block is filled by a discrete harmonic map: solve `Δx₁ = Δx₂ = 0` on the logical
square `[0,1]²` with Dirichlet data from the four boundary curves, finite differences, then image of a
uniform grid.

Today `densify.py` relaxes patch interiors toward the field **only when the domain carries a guide**
(`spend=bool(d.guides)`); unguided patches keep their `discrete_coons_patch` interiors, and so do all
pole patches. That gating is deliberate and well-argued — spending 12° of minimum angle to chase
discretisation noise was measured and rejected — but it leaves every unguided curved domain (disc,
ellipse, stadium, hexagon) on a plain bilinear blend.

A harmonic map is the classical fix for Coons distortion on curved and non-convex patches, it holds
patch boundaries fixed exactly as we require, and it does not consult the field at all — so it is not
a trade of quality for alignment and does not need `_accepts` to arbitrate.

- **Note the difference from what `densify.py` already does:** our regulariser is the Laplacian of the
  *correction* `x − x_coons`, which keeps the result near Coons. The paper's is the Laplacian of `x`
  itself, which is a different target. On a curved patch these are not the same map.
- **Verify:** `15_baseline.py`. The disc / ellipse / stadium / hexagon min-angle rows are the ones to
  watch; everything rectilinear should be unchanged, since on a rectangular patch with uniformly
  sampled straight sides the harmonic map *is* the Coons grid.
- **Caveat:** harmonic maps are guaranteed bijective onto **convex** targets (Radó–Kneser–Choquet). Our
  patches are not always convex, so folds are possible — the existing `hard_floor` check catches them,
  but this must be a per-patch accept/reject like `relax_patch`, not an unconditional replacement.
- **Not for pole patches**, for the same reason they are excluded from field relaxation: a collapsed
  side has no valid logical square to map from.

---

## 5. Not worth implementing

### 5.1 The diffuse domain method and the FFT solve

*Paper §3.1, §5.1 — their headline contribution.*

- **The efficiency argument does not apply.** Our solve is 0.017 s on 699 vertices, and 0.24 s with
  the §4.1 relaxation added. The field is nowhere near the bottleneck; `trace`, `arrangement`,
  `repair` and `densify` are.
- **The "no triangular mesh needed" argument does not apply either.** `background.py` exists, is
  validated (`validate()` checks orientation, Euler characteristic, boundary loop count), and its
  quality is not what is limiting us — §3 shows the *model* was.
- **The cost is a rewrite, not a swap.** `trace.py`, `arrangement.py`, `repair.py` and `densify.py`
  all sample the field through `angle_in_face(fkey, bary)` — face location and barycentric
  interpolation on a triangulation. A regular grid replaces that with bilinear interpolation and O(1)
  location, which is genuinely nicer, but it is every call site.
- **And it degrades what we care most about.** The diffuse interface makes the field near the wall
  accurate only to `O(ε)`. Wall behaviour is where our layout quality is decided — landings, corners,
  grazing hits, `edges_to_curves` recovery. Trading exact boundary tangency for FFT speed is the wrong
  trade for this project.

The one property genuinely worth envying is that their field is independent of the boundary
discretisation. §3.1 shows ours already is at coarse resolutions (48 vs 96 points, identical), so we
are not paying much for the triangulation.

**Revisit only if** the module ever moves to curved surfaces (§6, LIMIT 6), where the triangulation
stops being incidental — and note the paper's method is planar-only too, so it would not help there.

### 5.2 Singularity detection by local maxima of the discrete Dirichlet energy

*Paper §5.2.* They flag a cell as singular when its energy is a strict local max over the 3×3
neighbourhood.

Explicitly worse for us. That is a thresholded heuristic that gives a *location* and no index; a grid
field needs it because it has no combinatorial structure to count winding on. `CrossField.singularities`
computes the exact winding index per triangle — no threshold, no neighbourhood, and it returns the
index, which `launch_directions` needs (`target = 4 − index`). Keep ours.

The *localisation* half of their §5.2 is a different matter — that is §4.2 above, and worth taking.

### 5.3 The unconditional energy decay theorem

*Paper §4, Theorem 4.1 — their main theoretical result.*

Not a reason for us to change anything. It guarantees their **iteration** is stable for any `(ε, τ)`.
We currently solve a linear system directly, so there is no stability question to answer.

It becomes relevant only *because of* §4.1: adopting the iteration means adopting its convergence
behaviour, and the theorem is the reassurance that it cannot diverge for any τ. Worth citing in the
docstring; not worth restructuring anything for.

### 5.4 An implementation trap in the paper, if any of this is coded from it

The paper is internally inconsistent about the penalty coefficient. Equation (4.1) and Theorem 4.1
both write

```
u_{k+1} = normalise( e^{τΔ}u_k + (τ/ε³)·B(φ)·g )
```

but **Algorithm 4.1, line 5** writes `(2/ε³)·B(φ)·g`, as does the sentence introducing (4.1) in §4.
Theorem 3.2's proof uses `τ/ε³`. Take `τ/ε³` — it is the one the derivation produces and the one that
makes the scheme a supporting-hyperplane step of `E_{τ,ε}`. (A smaller inconsistency: the §3 proof
alternates between `1/τ` and `1/2τ` prefactors.) None of this affects §4.1, which drops the penalty
term entirely in favour of hard constraints.

---

## 6. What the paper does not solve

- **LIMIT 4 — a domain needing a cut from neither a singularity nor a corner.** `disc + round hole` is
  still on the `triangulation` route in `baseline.json` (min angle 17.9–26.9°) and it is this module's
  one clear defeat against `SkeletonDecomposition`. The paper does not help. A smooth annulus has
  outer boundary winding `+4` quarter-turns and inner `−4`, so the total is zero: the boundary-aligned
  field is the singularity-free polar field (concentric circles and radials), under their model as
  much as ours. Their tracing (§5.2) also starts only from singularities, so it has the same gap. The
  annular domains in their Fig 13 all have corners or asymmetry to break it. `Tracer.hole_launches` —
  ours, not theirs — remains the right line of attack.
- **LIMIT 3 — non-orthogonal cable families.** The paper is a pure cross-field method with 4-fold
  symmetry baked in. Nothing there helps case C of the evaluation doc; that still needs a frame field
  and the warp step.
- **LIMIT 2 / the project objective — cables and force lines as field inputs.** The paper has no
  mechanism for interior directional constraints at all. Their `g` is boundary data only. Our
  `constraints.from_curves` and `densify.py` have no counterpart in it.
- **LIMIT 6 — surfaces.** Planar only, and more deeply so than we are: their computational domain is a
  periodic square.
- **Everything in `repair.py` and `arrangement.py`.** The paper devotes about half a page (§5.2, §5.3)
  to going from a field to a mesh and reports no failures. That is roughly 1,900 lines here, all of it
  written against failures that actually occurred. Read the brevity as scope, not as evidence that the
  problems are easy.

---

## 7. Order — steps 1 and 2 are done

Each step ends with `15_baseline.py` and a diff; do not stack two of them before diffing.

1. ~~Commit the §3 measurement as a check script.~~ **Done** — [19_field_accuracy.py](19_field_accuracy.py).
   Default run is the disc radius check; `--winding` runs Poincaré–Hopf over all 16 domains under both
   solvers (all pass, and the relaxation leaves every boundary winding alone); `--timing` is the cost.
2. ~~§4.1 unit-norm relaxation behind a flag.~~ **Done** — `CrossField.solve(relax=...)`,
   `FieldDecomposition.from_boundary(relax=..., field_tau=...)`, `15_baseline.py --relax`.
   **Default stays False**, on the §3.3 numbers: 18 rows better, 14 worse, and the regressions are
   downstream of the field rather than in it.
3. **§4.4 boundary weight sweep — promoted to next.** §3.3 changed its standing: it is no longer a
   speculative knob, it is the lever that pushes singularities back *inward*, which is exactly the
   variable the hexagon and ellipse regressions turn on. Sweep `from_boundary(weight=…)` against
   `--relax` and look for a weight that keeps the disc at r ≈ 0.85 while holding the hexagon's
   singularities clear of the wall.
4. **Widen the wall-proximity assumptions in `trace` and `repair`** — the actual blocker now.
   `launch_directions`' probe radii (0.45h–2.0h) and `repair._cluster`'s landing tolerance both assume
   a singularity several `h` inside the domain. That was true of the old field and is not true of the
   new one. This is the work that would let the default flip.
5. **§4.2 sub-cell localisation** — after 3/4, because normalisation moves the zero. Instrument which
   `(scale, samples)` attempt `launch_directions` succeeds on, before and after.
6. **§4.3 the six curvature/corner domains** as new baseline rows with singularity-signature
   assertions. Independent of everything above, and now more interesting than it was: §6.3 of the paper
   is *about* singularities being driven onto the wall, which is the regime we just entered.
7. **§4.5 harmonic block filling** — **demoted to last** by §8.3: on the rounded domains it was meant
   to help, every sub-40° face is at a patch *corner* and none is in a patch *interior*, which is not
   what a better interior fill fixes.

Steps 3 and 4 are now the same work seen from two directions — §3.3 reached it from the relaxation's
regressions, §8.3 from measuring where rounded domains lose their quality. Both point at separatrix
landings and patch corners on curved walls.

### A note on `baseline.json`

It is **stale as of this work, and not because of it.** Running `15_baseline.py` with `relax=False`
moves 4 rows — all `disc+round hole`, from `triangulation` to `field` with 7 patches. That is
`Tracer.hole_launches` (LIMIT 4) landing after the baseline was last regenerated: the domain has
**zero** singularities and **zero** corner launches, so all 4 of its separatrices come from hole cuts,
and no change to the field solver could have produced them. Regenerate deliberately, as its own
reviewed step — not folded into a field change.

---

## 8. Why the paper's rounded shapes look so good

Their Fig 13 meshes a figure-eight with two round holes, a stepped shape with a round hole and an
airfoil in a rectangle, and all of them come out clean. Ours are the weakest rows in the suite. Four
hypotheses, three of them measured and killed.

### 8.1 It is not the boundary data

Their Algorithm 5.1 mollifies the boundary driving field; we pin every boundary vertex hard, to a
cross built from its two adjacent **chord** directions. That looked like the obvious suspect on a
curve, where a chord is not a tangent.

Measured on the disc — pinning the **analytic** circle tangent instead, at two resolutions, two
boundary samplings, both solvers, and with the boundary deliberately sampled non-uniformly:
**byte-identical results in every case.** Two reasons, and both were already in the code:
`_densify_loop` re-samples the boundary to `target_length` before anything sees it, so an irregular
input sampling does not survive to the solve; and `from_boundary` averages the two chords **in the
4th-power representation**, where the average of a symmetric pair is exactly the tangent.

So Algorithm 5.1 buys us nothing on a smooth boundary. Struck from §4's list.

### 8.2 It is not the interior triangulation

Our interior is a jittered axis-aligned grid, which is not 4-fold symmetric even when the disc and
its sampling are — a plausible source of the residual scatter in the four singularity radii.

Measured by rotating the domain off the 48-gon's symmetry angle (1.0°, 2.5°, 3.7°, 5.0°, 6.2°, so
the grid lands differently against the same shape each time): radii move in the **third decimal**,
spread 0.009 → 0.012. Not it.

*(A first attempt at both tests measured nothing: the rotations were exact multiples of 360/48, so
the polygon mapped onto itself, and on a regular polygon the chord average equals the analytic
tangent identically. Worth recording — the null result looked like a clean answer.)*

### 8.3 The field is now right; the loss is at patch CORNERS

With `relax=True` the disc's field is correct to within a triangle (r = 0.85, spread 0.01). So the
gap is downstream. Splitting the dense mesh's sub-40° faces by where they sit:

| domain | solver | patches | faces <40° at a patch corner | faces <40° in a patch interior |
|---|---|---|---|---|
| disc | plain | 5 | 0 | 0 |
| disc | relax | 5 | 0 | 0 |
| ellipse | plain | 10 | 0 | 0 |
| ellipse | relax | 4 | **13** (worst 8.0°) | **0** |

Every bad element is at a patch corner. None is in an interior. The disc has none at all — its
50–68° minimum is a *uniformly mediocre* mesh, not a mesh with a few bad spots.

That kills the harmonic-map idea (§4.5) as the fix for rounded domains: `discrete_coons_patch` is not
where they lose quality. What is: where separatrices land on a curved wall and where they meet each
other — `repair` and `arrangement`.

### 8.4 What is actually different: their rounded shapes are not uniformly round

This is the part worth internalising, and it is their §6.3 read backwards.

A cross field's singularities are **well-determined where curvature is concentrated and free to
drift where it is uniform**. §6.3 is the demonstration: as the top boundary goes from `y = −√(1−x²)`
to `4x²` to `16x²`, the two `−1/4` defects are drawn progressively and predictably to the wall, and
at a `C⁰` corner they pin exactly to it.

Every rounded shape in their Fig 13 has that property. The airfoil is the clearest case — smooth and
closed, but its leading and trailing edges concentrate curvature, so the field has definite,
well-separated singularities and the separatrix graph is clean. A **perfect circle is the degenerate
opposite**: uniform curvature, nothing to pin to, and in the annular case it admits the
singularity-free polar field outright. The paper does not show a circle in a circle being
decomposed, because there is nothing there for the method to find.

Our suite is unusual in containing exactly that degenerate case (`disc+round hole`, LIMIT 4) and in
reporting it. So the honest comparison is not "they handle round shapes and we do not" — it is that
**uniform curvature is hard for every cross-field method**, ours included, and the response has to be
structural (what `Tracer.hole_launches` does) rather than a better solver.

Two fairness points belong with that, the same ones `17_vs_skeleton.py` needed:

- their figures are **single curated results**; §6.2 says outright that the energy landscape is
  non-convex with multiple locally stable configurations and that `(ε, τ)` selects among them. Our
  numbers are a 64-row matrix that takes the first answer at four resolutions and prints its worst
  rows;
- they show the **decomposition and the mesh**, not a quality metric. We have no minimum angle from
  them to compare against.

### 8.5 What to take from it

1. **Nothing in their field solver or boundary handling** — §8.1 and §8.2 close those off.
2. **Their §6.3 as a diagnostic we do not have.** Curvature concentration predicts where singularities
   go and how firmly they are held. A domain whose boundary curvature is near-uniform over a long run
   is one where the field will not decide the layout for us, and that is knowable **before** tracing.
   Detecting it and routing to a structural cut is a better answer than tuning tolerances.
3. **Aim the next work at patch corners on curved walls**, not at patch interiors — §8.3. This is the
   same conclusion §3.3 reached from the relaxation regressions, arrived at independently.
