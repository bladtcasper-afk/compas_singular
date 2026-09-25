# Frame-field front end for compas_singular — module layout

**Date:** 2026-08-03
**Companion to:** `cross-field-quad-meshing-evaluation.md`
**Scope:** replace `SkeletonDecomposition`'s medial-axis front end with a frame-field →
separatrix → coarse-layout front end, leaving everything downstream untouched.

---

## 0. The integration seam (verified)

`SkeletonDecomposition.decomposition_mesh()` ends at
[decomposition.py:214](../../src/compas_singular/algorithms/decomposition.py#L214):

```python
self.mesh = CoarsePseudoQuadMesh.from_polylines(boundary_polylines, other_polylines)
```

`from_polylines` is inherited from **compas core** `Mesh` (verified in compas 2.15.1,
signature `(boundary_polylines, other_polylines)`).

**So the entire front-end contract is:**

> produce a list of polylines — each a `list[[x, y, z]]` — partitioned into those lying on the
> domain boundary and those that do not.

Nothing else. Everything downstream works unchanged on the result:

| Downstream capability | Where |
|---|---|
| `collect_strips`, `collect_polyedges` | `mesh_quad.py` |
| `set_strips_density_target`, `densification(edges_to_curves=...)` | `mesh_quad_coarse.py`, `mesh_quad_pseudo_coarse.py` |
| `add_strip` / `delete_strip` grammar | `mesh_quad/grammar/` |
| `guide_line_mesh`, `guide_band_mesh` | `guide_lines.py` |
| two-colouring, isomorphism, lizard encoding | `algorithms/`, `datastructures/lizard/` |

This is a much cleaner seam than expected. **The new front end must be polyline-compatible,
not mesh-compatible.**

---

## 1. Two findings that shape the design

### 1.1 `boundary_triangulation` has no interior vertices

[triangulation.py:49](../../src/compas_singular/algorithms/triangulation.py#L49):

```python
vertices = [pt for boundary in [outer_boundary] + inner_boundaries + polyline_features
            for pt in boundary] + point_features
```

Every vertex is a boundary or feature point. The interior is spanned by large, often sliver
Delaunay triangles. That is *by design* — their circumcentres approximate the medial axis, which
is exactly what the skeleton method needs.

**A frame field cannot live on that mesh.** Field smoothness is measured per-edge across a
triangulation; slivers wreck the discretisation and there is no interior resolution to carry a
field at all. `background.py` must refine, or build its own interior triangulation, and keep the
domain outline exactly.

This is the single most important practical consequence of the whole redesign: the two front ends
want *opposite* triangulations from the same boundary.

### 1.2 The existing feature-curve cut is already the §5.2 pre-cutting trick

[triangulation.py:68-70](../../src/compas_singular/algorithms/triangulation.py#L68) does a
topological cut along `polyline_features` by unwelding edges. That is precisely the
"pre-cut the domain along the guide curves" option from the evaluation, already implemented —
so **Option B is reachable with existing code** and is the cheapest first milestone. (Known
caveat from prior work: the headless feature-curve path is fragile and leaves interior triangles.)

---

## 2. Why not just keep using `guide_lines.py`

`guide_lines.py` already places guides on an existing coarse mesh, and with `swing` it turns
the crossing courses too. Its own docstring states the limit
([guide_lines.py:117](../../src/compas_singular/guide_lines.py#L117)):

> *"Refinement subdivides the lines that exist, it does not invent lines in new directions — only
> `add_strip` does that."*

and the failure it causes:

> *"This plate's decomposition has no vertical line within x < 1.2 of the centre, so a guide at
> x = 0.5 is rejected however much you refine."*

That is the exact gap a field-derived layout closes: the layout is *derived from* the guides
rather than fitted to a layout that ignored them. `guide_lines.py` stays useful — it becomes the
**fine-tuning** stage after a field-derived layout, not the primary mechanism.

---

## 3. Module layout

**Status (2026-08-28): promoted.** The package now lives at
`src/compas_singular/framefield/`, imported as `compas_singular.framefield`. The tree below is
the original design-time layout, kept for the milestone history; see the note after it for what
actually shipped.

Originally a self-contained package **inside `examples/New approach/`**, so nothing in `src/`
was touched while the approach was being proven. A script run from that folder put it on
`sys.path[0]`, so `from framefield import ...` worked with no install step.

```
examples/New approach/
    framefield/
        __init__.py      public API
        background.py    BackgroundMesh      triangulation + tangent bases
        constraints.py   Constraint set from curves / boundary / stress tensor
        field.py         CrossField          the solve + singularity extraction
        warp.py          Warp                frame field -> cross field + deformation
        trace.py         separatrix tracing -> raw polyline network
        repair.py        network fix-ups (the analogue of solve_triangular_faces)
        decomposition.py FieldDecomposition  mirrors SkeletonDecomposition's API
    01_field.py          milestone 1a — field + singularities
    02_separatrices.py   milestone 1b — traced network
    03_layout.py         milestone 1c — coarse quad mesh
```

Promotion to `src/compas_singular/framefield/` happened 2026-08-28: a flat folder move, no edits
inside the package (every intra-package import was already relative). The 17 shipped modules
differ from the sketch above -- `viz.py`, `cache.py`, `arrangement.py`, `guides.py`, `quality.py`,
`relax.py`, `symmetry.py` and `edit.py` were added as the approach grew; `warp.py` was folded into
`field.py` and never existed standalone. `block.py` and `edit.py` stayed in `framefield/` rather
than moving to `compas_singular.blocks` / `compas_singular.editing`, since both hold field-layout
specifics those packages don't share. 27 consumer files across `compas_singular` and the sibling
`compas_topology` repo were repointed to `compas_singular.framefield` at the same time.

`decomposition.py` last, and deliberately mirroring `SkeletonDecomposition`'s method names, so
existing example scripts become drop-ins (§5).

---

### `background.py`

```python
class BackgroundMesh:
    @classmethod
    def from_boundary(cls, outer_boundary, inner_boundaries=None,
                      target_length=None, cut_curves=None) -> 'BackgroundMesh'
    mesh: Mesh                      # triangles, well-shaped, interior vertices present
    def face_basis(fkey) -> (e0, e1)        # orthonormal tangent basis per face
    def dual_edges() -> list[(f0, f1)]      # adjacency for the smoothness term
```

Owns finding 1.1. Planar for now, so the tangent basis is trivial (world XY) and no
parallel-transport bookkeeping is needed — a large simplification over the surface case. Keep
`face_basis` in the API anyway so the surface extension does not require a rewrite.

`cut_curves` routes to the existing unweld-based cut (finding 1.2).

---

### `constraints.py`

```python
Constraint = namedtuple('Constraint', 'fkey direction weight hard')

def from_boundary(background, tangent=True) -> list[Constraint]
def from_curves(background, curves, mode='tangent', band=None) -> list[Constraint]
def from_stress(background, tensor_field) -> list[Constraint]
```

`mode` is `'tangent'` or `'perpendicular'`. **Under 4-fold symmetry these are the same constraint
90° apart**, so for a cross field this is a one-line phase shift; it only becomes a genuine
distinction once `warp.py` is in play (see the evaluation doc, §2).

`band` limits a curve's influence to faces within a radius, so a guide that should steer its
neighbourhood does not fight the boundary halfway across the domain.

---

### `field.py`

```python
class FrameField:
    @classmethod
    def solve(cls, background, constraints, orthogonal=True, smoothness=1.0) -> 'FrameField'

    def frame(fkey) -> (u, v)               # the two directions, per face
    def period_jumps() -> dict[(f0, f1), int]
    def singularities() -> list[(fkey, index)]   # index in quarter-turns
    def is_orthogonal(tol=1e-6) -> bool
```

`orthogonal=True` → power-4 cross field, sparse linear solve, one complex unknown per face.
`orthogonal=False` → PolyVector (Diamanti et al.): still a **sparse linear solve with no integer
variables**.

numpy 2.5.1 / scipy 1.18.0 are present in `singular312`, so `scipy.sparse.linalg.spsolve` is
available. This module is genuinely the easy one — see the evaluation doc's difficulty table.

`singularities()` returns indices in quarter-turns (±1 for a valence-3/5 vertex). Sanity check
worth asserting: total index must satisfy Poincaré–Hopf for the domain. A field whose indices do
not sum correctly means the period jumps are wrong, and it is far cheaper to catch here than in
`trace.py`.

---

### `warp.py`

```python
class Warp:
    @classmethod
    def from_frame_field(cls, field) -> 'Warp'
    cross: FrameField                # the warped, now-orthogonal field
    def forward(xyz) -> xyz
    def inverse(xyz) -> xyz
    def is_identity(tol) -> bool
```

Panozzo et al.'s construction: factor the frame field into a cross field plus a per-face
deformation, mesh in the deformed domain, map back.

**Skip this module entirely for milestone 1.** Cases A and B (orthogonal stress field; single
cable family) produce `is_identity() == True`, so building it later costs nothing that milestone 1
has to unlearn — as long as `trace.py` and `decomposition.py` are written to accept an optional
`warp` and call `inverse()` on their output points.

---

### `trace.py`

```python
def separatrices(field, background, warp=None) -> list[list[[x, y, z]]]
def network(field, background, warp=None) -> (boundary_polylines, other_polylines)
```

Integrate the field's directions outward from each singularity, in each of its 4 (or n) sectors,
until hitting the boundary or another separatrix.

**This is the highest-risk module**, and it is where the separatrix route trades its avoidance of
mixed-integer parametrization for a different problem. Known failure modes, all of which need
explicit handling:

- **limit cycles** — a separatrix that spirals without terminating. Needs an arc-length cap plus
  self-proximity detection.
- **near-misses** — two separatrices passing within ε, producing a sliver patch instead of a
  clean junction. Needs a merge tolerance, consistent with the `geometric_key` precision the rest
  of the codebase matches points on.
- **grazing boundary hits** — a separatrix approaching the boundary tangentially.

`network()` returns exactly the `from_polylines` pair. Note the existing partition test at
[decomposition.py:211-213](../../src/compas_singular/algorithms/decomposition.py#L211) uses
`geometric_key` on boundary vertices; match that convention rather than inventing another.

---

### `repair.py`

```python
def split_collapsed_boundaries(polylines, boundary) -> list[polyline]
def split_high_turning_patches(polylines, limit) -> list[polyline]
def solve_non_quad_faces(coarse_mesh) -> None
```

**Budget the most time here.** `SkeletonDecomposition` carries four correction methods —
`branches_splitting_boundary_kinks`, `branches_splitting_collapsed_boundaries`,
`branches_splitting_flipped_faces`, `solve_triangular_faces` — roughly 120 lines of hard-won
special-case handling, and they exist because a raw polyline network is *not* a valid coarse quad
mesh.

They **cannot be reused**: every one of them is written against the Delaunay/skeleton
representation (`trimesh_face_circle`, `self.singular_faces()`, `self.vertex_faces`). The
separatrix analogues address the same underlying problems — patches that are not four-sided,
boundary loops with too few splits — but need reimplementing against the traced network.

The one genuinely reusable idea is `solve_triangular_faces`'s two cases
([decomposition.py:343](../../src/compas_singular/algorithms/decomposition.py#L343)): a triangle
with one boundary vertex is fixed by duplicating that vertex; with two, by merging them. Those
cases recur here.

---

### `decomposition.py`

```python
class FieldDecomposition:
    @classmethod
    def from_mesh(cls, trimesh, guides=None, mode='perpendicular',
                  orthogonal=None, **kw) -> 'FieldDecomposition'

    def decomposition_polylines(self) -> list[polyline]      # same name as Skeleton*
    def decomposition_mesh(self, poles=()) -> CoarsePseudoQuadMesh
    def edges_to_curves(self) -> dict[(u, v), polyline]

    field: FrameField
    warp: Warp | None
```

Deliberately mirrors `SkeletonDecomposition` so existing scripts change by two lines (§5).

`orthogonal=None` means **auto**: solve orthogonal, and fall back to a frame field only if the
guide constraints are mutually inconsistent under 4-fold symmetry. This keeps cases A and B on
the cheap path without the caller having to classify their own input.

`edges_to_curves()` matters more here than for the skeleton front end. Separatrices are *curved*,
so chording each coarse edge straight at densification would throw away exactly the alignment the
field was solved for. The tracer already has the true polyline for every coarse edge — pass it
through. (Prior work established this fix for boundary curves; the same mechanism applies to
interior separatrices, and the payoff is larger.)

---

## 4. Data contracts

| Between | Type |
|---|---|
| `background` → `constraints` | `Mesh` + per-face basis |
| `constraints` → `field` | `list[Constraint]` |
| `field` → `warp` / `trace` | per-face `(u, v)` + period jumps + `list[(fkey, index)]` |
| `trace` → `repair` | `list[list[[x, y, z]]]` |
| `repair` → `decomposition` | `(boundary_polylines, other_polylines)` |
| `decomposition` → **existing code** | `CoarsePseudoQuadMesh` |

Only the last row crosses into existing code, and it is a type that already exists.

---

## 5. The drop-in

Current, [02_decomposition_discrete_planar.py:19-31](../02_decomposition_discrete_planar.py#L19):

```python
trimesh = boundary_triangulation(outer_boundary, inner_boundaries, polyline_features, point_features)
decomposition = SkeletonDecomposition.from_mesh(trimesh)
coarsemesh = decomposition.decomposition_mesh(point_features)
coarsemesh.collect_strips()
coarsemesh.set_strips_density_target(0.5)
coarsemesh.densification()
densemesh = coarsemesh.get_quad_mesh()
```

Target:

```python
trimesh = boundary_triangulation(outer_boundary, inner_boundaries, polyline_features, point_features)
decomposition = FieldDecomposition.from_mesh(trimesh, guides=cables, mode='perpendicular')
coarsemesh = decomposition.decomposition_mesh(point_features)
coarsemesh.collect_strips()
coarsemesh.set_strips_density_target(0.5)
coarsemesh.densification(edges_to_curves=decomposition.edges_to_curves())
densemesh = coarsemesh.get_quad_mesh()
```

Two changed lines. Every existing example, viewer script and strip-editing workflow keeps working.

---

## 6. Build order

| # | Milestone | Modules | Proves |
|---|---|---|---|
| **1** | Boundary-aligned cross field, no guides, reproduce a skeleton-quality layout on the pentagon plate | `background`, `constraints.from_boundary`, `field` (orthogonal), `trace`, `repair`, `decomposition` | The seam works and the tracer is sound. Directly comparable against the existing output on the same input. |
| **2** | Single guide curve, `mode='perpendicular'` | `+ constraints.from_curves` | The thing `guide_lines.py` cannot do — a line where the layout had none. |
| **3** | Multiple guides, non-orthogonal | `+ field(orthogonal=False)`, `warp` | Case C. |
| **4** | Stress tensor input | `+ constraints.from_stress` | Case A end to end. |

Milestone 1 is the honest test and should be built before anything guide-related: if the traced
layout on a domain with *no* guides is worse than what `SkeletonDecomposition` already produces,
the problem is in `trace`/`repair` and no amount of field constraint work will fix it.

---

## 7. Risk register

| Risk | Severity | Note |
|---|---|---|
| `repair.py` is as large as the skeleton's correction layer | **High** | ~120 lines of special cases in the existing code, none reusable. Assume comparable. |
| Separatrix limit cycles / near-misses | **High** | The known weak point of the separatrix route versus MIQ. Needs caps and merge tolerances from day one. |
| Background triangulation quality (§1.1) | Medium | Cheap to fix, easy to overlook, and silently degrades the field if missed. |
| Patches that are not four-sided | Medium | Not automatic from a separatrix partition. `solve_triangular_faces`'s two cases recur. |
| `from_polylines` on a degenerate network | Low | Match the existing `geometric_key` precision convention. |
| Field solve | **Low** | Sparse linear solve; scipy present. Do not let this absorb the schedule — it is the easy module. |

---

## 8. Step-by-step implementation plan

Each step ends in something runnable and checkable. Do not start a step before its predecessor's
check passes — every one of these catches a class of bug that is invisible downstream.

### Step 1 — `background.py`: a triangulation a field can live on

1. Densify the outer boundary (and any inner boundaries) to roughly `target_length`.
2. Generate interior points on a grid of spacing `target_length`, keep those strictly inside the
   outer boundary and outside every inner boundary, with a margin so they do not crowd the wall.
3. `compas.geometry.delaunay_triangulation` over boundary + interior points.
4. Delete degenerate (zero-area) faces and faces whose circumcentre lies outside the domain —
   same test as [triangulation.py:63](../../src/compas_singular/algorithms/triangulation.py#L63).
5. Expose `vertices`, `faces`, `boundary_loops`, and per-boundary-vertex tangent.

> Planar domain ⇒ all tangent spaces are world XY ⇒ **no parallel transport, no per-edge
> connection**. Keep `face_basis()` in the API for the eventual surface case, but return the
> identity for now.

**Check:** every face has positive area; every boundary vertex has exactly 2 boundary neighbours;
Euler characteristic V − E + F matches the domain's.

### Step 2 — `field.py`: the cross field

Represent a cross by one complex number per **vertex**, `u = e^{i·4θ}`. The 4th power is what makes
the 90°-symmetric object single-valued, so a plain linear solve is legitimate.

1. Build the graph Laplacian over triangulation edges (uniform weights first; cotangent later only
   if quality demands it).
2. Constrain boundary vertices: `u = e^{i·4·θ_tangent}`, hard.
3. Harmonic-extend to the interior: solve `L·u_interior = −L_bc·u_boundary`.
4. Recover `θ = arg(u)/4` per vertex. The cross is `{θ, θ+π/2, θ+π, θ+3π/2}`.

Do **not** normalise `u` before extracting singularities — `|u|` dropping toward 0 in the interior
is not solver failure, it is the Ginzburg–Landau relaxation, and the zeros are exactly the
singularities.

**Check:** on a disc-like domain with a smooth boundary the field is near-constant `|u| ≈ 1` away
from a few isolated dips.

### Step 3 — `field.py`: singularity extraction

Per triangle, sum the angle differences around its three edges, each wrapped into `(−π/4, +π/4]`
(the cross's period). The sum is necessarily a multiple of `π/2`; the index in quarter-turns is
`sum / (π/2)`.

**Check — two of them, and you need both.**

*(1) Poincaré–Hopf, as the discrete argument principle.* Compare the interior index sum against
the **measured winding of the field along the boundary**, both in quarter-turns:

```
Σ interior indices  ==  winding of θ around the boundary loops
```

Do **not** assert `Σ == 4·χ`. That shortcut holds only for a *smooth* boundary. On a polygon the
corners absorb the turning: a square's boundary-tangent cross field is the constant field, with
zero boundary winding and no interior singularities at all, because a 90° corner turns the tangent
by exactly the cross's period and therefore costs the field nothing.

*(2) An independent prediction.* Check (1) compares two numbers read off the same angles, so a
**mis-pinned wall passes it happily** — it validates consistency, not correctness. Predict the
answer separately from the quad-mesh index identity:

```
Σ interior (4 − valence)  +  Σ boundary (3 − valence)  =  4·χ
```

counting each convex corner as a valence-2 quad corner (+1) and each reflex corner as valence 4
(−1). Verified values: square 0, L-shape 0, square-with-square-hole 0, smooth disc +4, pentagon −1.

> **Implemented, and it caught a real bug.** Taking the boundary tangent as the chord between a
> vertex's two neighbours — the obvious implementation — is wrong at a corner, and wrong in the
> worst possible way: at a 90° corner the chord bisects to 45°, exactly halfway between the two
> arms of the cross and equally far from both. The field then unwinds that spurious 45° somewhere
> in the interior, and a square comes out with **two** singularities instead of none. Check (1)
> passed throughout. Only check (2) found it.
>
> The fix is to keep both adjacent edge directions and average them **in the 4th-power
> representation**, where a right-angle corner's two edges are literally the same complex number.
> Where they genuinely disagree (interior angle near 45° or 135°) the average cancels to near
> zero; leave that vertex free rather than pinning it to an arbitrary compromise.

### Step 4 — `constraints.py`: guide curves

1. `from_boundary` — already implicit in Step 2; factor it out to a `Constraint` list.
2. `from_curves(background, curves, mode)` — for each curve sample, find the nearest vertex (or
   vertices within `band`), constrain `u = e^{i·4·θ}` where `θ` is the curve tangent
   (`mode='tangent'`) or tangent + 90° (`mode='perpendicular'`).

**Under 4-fold symmetry these two modes are the same constraint** — `e^{i·4·(θ+π/2)} = e^{i·4·θ}`.
Implement `mode` anyway, assert the identity in a test, and leave a comment: it stops being an
identity the moment `warp.py` exists, and a future reader will otherwise "simplify" it away.

**Check:** a guide across a square visibly bends the field; removing it restores the unguided one.

### Step 5 — `trace.py`: separatrices

From each singular triangle, launch one polyline per sector of the cross.

- Step with fixed arc length `h ≈ target_length / 2`.
- At each point: locate the containing triangle, take the three vertex angles, pick from each
  vertex the cross representative **closest to the incoming direction** (this branch matching is
  the whole trick), average, advance.
- Terminate on: boundary crossing, another singularity, arc-length cap, or self-proximity.

Guard limit cycles and near-misses from the first version, not as a later fix — they are the known
weak point of the separatrix route.

**Finding the launch directions is the fiddly part, not the integration.** A separatrix leaves
wherever the field points *radially*, i.e. where the residual `φ − θ` hits a multiple of the
period. Three attempts, and only the third works:

| Attempt | Result on the disc / pentagon |
|---|---|
| Sign changes in the **wrapped** residual | 4 arms on a `+1` singularity (should be 3); two launches 1.5° apart |
| Level crossings of the **unwrapped** residual | 11–13 arms; some 0.2° apart |
| Level crossings of the **running maximum** of the unwrapped residual | Correct: 3 arms per `+1`, 5 per `−1` |

The wrapped residual is a sawtooth, so a sign-change search cannot distinguish a genuine root from
the −45°/+45° jump. Unwrapping fixes that but exposes the real problem: near a singularity
`|u| → 0`, so `θ = arg(u)/4` is ill-conditioned and the residual picks up small oscillations, and
every wiggle across a level counts as another crossing.

The fix inverts the problem. **The arm count is known in advance** — over one lap the residual
advances by exactly `4 − index` periods — so nothing needs to be *discovered*, only *located*.
Searching against the running maximum, which is monotone by construction, consumes each level
exactly once.

**Check:** arm count equals `4 − index` per singularity; every trace terminates at a boundary or
another singularity, never on the length cap or a cycle. Verified — disc: 4 singularities × 3 arms
= 12 traces (8 into singularities, 4 into the boundary); pentagon: 5 arms at 68.9–75.3°, against
72° ideal.

### Step 5b — separatrices also come from CORNERS, not only singularities

**Not in the original design, and without it a whole class of domain has no layout at all.**

An L-shape has no interior singularities — the field is exactly constant — so it produces no
separatrices and no patches. Yet it obviously decomposes into three rectangles, and those cuts
start at the reflex corner.

A boundary corner of interior angle `α` is a *boundary singularity*: the quad layout puts
`k = round(α/90°)` patch corners there and it emits `k − 1` separatrices into the domain.

| α | k | emits |
|---|---|---|
| 90° convex | 1 | 0 |
| 180° straight wall | 2 | — never reaches here; the turn is below the corner threshold |
| 270° reflex | 3 | 2 |

Two sign traps, both of which fail *silently* — the corner emits nothing and the domain just comes
out with no interior lines:

- **The interior angle is the sweep from the outgoing edge round to the incoming one**, on a CCW
  loop. Reversed, it returns the exterior angle and reads a 270° reflex corner as 90° convex.
- **The interior direction is the outgoing wall rotated by half the interior angle — not the
  bisector of the two wall vectors.** At a reflex corner those point opposite ways: `back + fwd`
  aims into the 90° wedge the domain does *not* occupy, so the probe point lands outside the mesh
  and the field sample fails.

Test the arms by whether they lie strictly inside the interior sector, not by a dot product against
the interior direction — at a 270° corner a valid arm can be 135° away from it.

### Step 6 — `repair.py` + `decomposition.py`: the network → coarse mesh

1. Snap separatrix endpoints to the boundary and to each other within a tolerance matched to the
   codebase's `geometric_key` precision.
2. Split every polyline at its intersections so the network is a proper planar graph.
3. Partition into `(boundary_polylines, other_polylines)` using the boundary test at
   [decomposition.py:211](../../src/compas_singular/algorithms/decomposition.py#L211).
4. `CoarsePseudoQuadMesh.from_polylines(...)`.
5. Fix non-quad faces — port the two cases from `solve_triangular_faces`
   ([decomposition.py:343](../../src/compas_singular/algorithms/decomposition.py#L343)):
   one boundary vertex ⇒ duplicate it; two ⇒ merge them.

**Check:** `mesh.is_quadmesh()`, `collect_strips()` succeeds, and `densification()` yields a
manifold dense mesh — i.e. the existing downstream pipeline accepts it.

**`from_polylines` cannot build a domain with no interior line.** It discards any face whose
vertices are *all* on the boundary — that is how it drops the outer face — so a square comes back
with 4 vertices and 0 faces. A square has no singularities and no reflex corners, and its correct
coarse layout is a single quad. Special-case it: build the one face directly from the boundary
arcs' endpoints.

---

## 9. Implementation status

Originally `examples/New approach/framefield/`; now `src/compas_singular/framefield/` (see
the status note at the top of §3). Runnable via `01_field.py`, `02_separatrices.py`,
`03_layout.py` in the singular312 env.

| Step | State |
|---|---|
| 1 `background.py` | **Done.** Validated on 6 domains: positive areas, loop count, Euler characteristic. |
| 2 `field.py` solve | **Done.** Sparse complex solve, hard + soft constraints. |
| 3 singularities | **Done.** Both checks pass on all 6 domains against independent predictions. |
| 4 `constraints.py` | **Done.** Boundary + guide curves; 4-fold `mode` identity asserted. |
| 5 `trace.py` | **Done.** Arm counts exact; every trace terminates legitimately. |
| 5b corner launches | **Done.** |
| 6 `repair.py` / `decomposition.py` | **Partial** — see below. |
| 7 `edges_to_curves` | **Done.** Wired through `densification`. |
| 8 comparison vs `SkeletonDecomposition` | Not started. |
| `quality.py` + `15_baseline.py` | **Done.** Two-tier gate; 64-row baseline committed and diffable. |

### Two fixes that changed the picture

**Interpolate `u`, not the angles.** `angle_in_face` originally blended the three vertex *angles*,
each defined only modulo 90°, which needs a branch choice per vertex. Near a singularity the three
genuinely differ by more than a quarter period and no branch choice is right. Walking a probe
circle around a pentagon's singularity, the unwrapped residual advanced by **2, 3 or 4 periods
depending only on the radius**, where the index says it must always be 5.

Interpolating `u` — a single-valued smooth complex field — and taking `arg(u)/4` removes the
ambiguity entirely; the branch matching disappears. This one change fixed *both* the pentagon's
resolution sensitivity and the disc, and tightened the disc's launch-angle spread from 54–200° to
105–137°. The 4th power was already doing this job; the interpolation just wasn't using it.

**Launch corner separatrices from the corner, not the probe.** The probe is inset along the
interior bisector purely to sample the field somewhere safely inside. Starting the *trace* there
displaced the whole separatrix sideways by that inset: on a T-plate the line leaving one reflex
corner ran 0.28 *below* the wall it should lie on, passed under the opposite corner instead of
terminating at it, and continued to the far edge. `from_polylines`' planar face search cannot
separate a region from a line lying on its own edge, so the plate came back as one 18-sided face.

Two reflex corners facing each other also launch the same line from both ends, so it arrives twice
reversed and has to be deduplicated.

### End-to-end results

Every number below is read from [baseline.json](baseline.json), not from memory. Regenerate and
diff it with `python 15_baseline.py`; the columns are coarse patch count, min/max face angle and
worst aspect ratio at background 0.3–0.6 with quad target 1.0.

**Rectilinear plates — exact, at every resolution.** All of these are 90.0° / 90.0° / AR 1.00,
100% coverage, no irregular interior vertex:

| Domain | Singularities | Coarse layout |
|---|---|---|
| square | none | 1 quad |
| L-plate | none (reflex corner ×1) | 3 quads — the 3-rectangle decomposition |
| T-plate | none (reflex ×2) | 4 quads |
| U-plate | none (reflex ×2) | 5 quads |
| plus-plate | none (reflex ×4) | 5 quads |
| comb-plate | none (reflex ×6) | 11 quads |
| rectangle + rectangular slot | none (reflex ×2) | 5 quads |
| square + square hole | none (8 corner launches) | 8 quads |

**Curved boundaries — all on route `field` since `arrangement.py`.** This table previously read
*"disc: falls apart at 0.4 and below"* and *"ellipse, stadium: triangular and many-sided patches at
every resolution"*. Both were the pre-arrangement picture and both are false; see §"The two
limitations that remain" for the diagnosis and `13_limits.py` LIMIT 1:

| Domain | Coarse patches @ 0.6/0.5/0.4/0.3 | min angle | worst aspect |
|---|---|---|---|
| pentagon | 5 / 5 / 5 / 5 | 65.0 – 70.7° | 1.66 – 1.71 |
| hexagon | 11 / 11 / 11 / 11 | 62.4 – 67.9° | 2.20 – 3.48 |
| disc | 5 / 5 / 9 / 9 | 51.9 – 58.7° | 1.79 – 2.76 |
| stadium | 4 / 4 / 4 / 4 | 54.6 – 60.3° | 1.65 – 2.27 |
| ellipse | 10 / 10 / 10 / 10 | 24.4 – 65.9° | 2.06 – 3.00 |

The ellipse's spread is the one outlier in the suite, and it is a single row: background 0.6, the
only resolution whose layout carries poles. See §"What the resolution sweep says" in
`15_baseline.py`.

**Still on a fallback:**

| Domain | Result |
|---|---|
| disc with a round hole | route `triangulation` — no corner anywhere and no separatrix at all, so there is no layout to repair. LIMIT 4, and a different problem from the crossings one |
| square + ring cable | route `field` at 0.6/0.5, `polygon` at 0.4/0.3, and a 180.00° face at 0.5/0.4/0.3 — the suite's one remaining hard-floor breach |

Demonstrated by [10_domains.py](10_domains.py), [11_downstream.py](11_downstream.py),
[12_cables.py](12_cables.py), [13_limits.py](13_limits.py); measured by
[15_baseline.py](15_baseline.py); see [README.md](README.md). Each of the first four opens a
compas_viewer scene via `framefield/viz.py`, following the `scene.add_group` convention the rest of
`compas_singular/examples` already uses.

### The guarantee: `quad_mesh()` always returns a mesh

**The deliverable is a mesh.** A front end that raises instead of producing one is useless however
good its field was, so `FieldDecomposition.quad_mesh()` owns the whole chain — layout, repair,
densification — and every fallback in it. Callers should not assemble `decomposition_mesh` +
`collect_strips` + `densification` by hand; each of those can fail independently.

Three routes, best first. `d.route()` says which was taken; `d.warnings()` gives the detail.

| Route | When | Element flow |
|---|---|---|
| **field** | the separatrix layout closed and densifies | follows the field — the point of all this |
| **polygon** | layout failed, domain simply connected | the domain as one n-gon, quad-split; ignores the field |
| **triangulation** | nothing else available (multiply connected, say) | follows the triangulation; ugly, high-valence, correct |

The result is validated before it is returned: all quads, manifold, and **covering 97–103% of the
domain area**. That last check is the one that matters and the one that is easy to omit — a
round-holed plate produced a manifold, all-quad mesh covering 108% of its domain because the
patches ran straight across the hole. Nothing structural was wrong with it; it simply was not a
mesh of that domain.

#### Those three checks are not a quality gate, and never were

This section used to end: *"14 domains × 4 background resolutions = 56 runs, every one a manifold
all-quad mesh at 100% coverage."* That was true and it is still true of 60 of the current 64 runs.
It is also **compatible with a mesh nobody can use**, and it is now actively misleading, because
`route()` — which used to compensate, since a bad layout normally fell back and said so — puts
every domain in `14_arrangement.py` on `'field'` since `arrangement.py` landed. Neither number
separates good from bad any more.

Measured: all-quad + manifold + 97–103% coverage certifies `rect with slot` at every resolution,
a plate whose mesh contains a **0.00° angle, a 180.00° angle, and an aspect ratio of 38**.

So quality is gated separately, in **two tiers**, because one threshold cannot do both jobs —
minimum angle on the field route spans 3.5° to 90.0° across the suite with no degeneracy among
those rows, so anything loose enough to pass the ellipse cannot notice the square slipping off 90°:

| Tier | What | Where | Absolute? |
|---|---|---|---|
| **hard floor** | no angle at or near 0° or 180°, no non-finite aspect ratio | `framefield/quality.py`, on every `quad_mesh()` | yes — permanent, and no baseline may bless a breach |
| **regression check** | every row, every metric, vs. committed numbers | `15_baseline.py` + `baseline.json` | no — relative to the recorded state |

Neither reroutes and neither raises: the "always returns a mesh" contract stands, and a failure
comes back as a named `QUALITY FAILURE` entry in `warnings()`. Rerouting on quality would make
things *worse* — the only destination left is the triangulation backstop, whose own minimum angle
is 21.5° with 802 irregular interior vertices on the round-holed disc.

`FieldDecomposition.quality()` returns the metrics as a dict so callers gate on numbers rather than
on a verdict someone typed by hand. That failure mode is not hypothetical: the report this work
replaced put the hexagon's 65.3° minimum angle into the ellipse's row.

**Pseudo-quads are the trap in every one of these metrics.** A pole is the quad `(p, a, b, p)` with
one side collapsed. Read as it is stored — a three-corner face — the ring cable's 22 pole faces
report angles from 7.4° to 117.5° and finite aspect ratios. Reconstitute the phantom fourth corner
and **all 22** report 0.00° and divide by zero. `python 15_baseline.py --poles` asserts the first
and prints the second.

Verified: **16 domains × 4 background resolutions = 64 runs**, recorded in `baseline.json`. Route
and floor split: **57 field/clear, 4 triangulation/clear, 1 field/breached, 2 polygon/breached,
0 raising.** The three remaining breaches are all `square+ring cable`.

### Reading a disc's layout, and the duplicate-trace bug

A disc's four `+1` singularities produce the canonical layout: **a central quad, ringed by four
more**. Each singularity emits `4 - 1 = 3` arms; two run to its neighbours around the ring, one runs
out to the wall. 12 traces, 8 distinct lines, 5 patches. It is the shape you would draw by hand.

Drawn on screen it looked nothing like that — the ring was doubled and visibly wobbly. **Every
sing-to-sing connection was being traced twice, once from each end, and the deduplication caught
none of the four.** Two causes:

* The duplicate test compared `pts[len(pts)//2]`. The two traces are decimated independently and
  run in opposite directions, so the same index is a different position along the curve; the four
  pairs' index-midpoints sat 0.55–1.10 apart against a 0.48 tolerance. Fixed by comparing
  **arc-length** midpoints.
* The tolerance was scaled to the background spacing. Two traces of one connection bow apart by a
  fraction of how far they *run* — both are noisiest where they leave their singularity — so a
  fixed tolerance recognised the short pairs and missed the long ones. Scaling it to the
  connection's own length fixes the rest.

Merging now **averages** the two traces along a common parameter rather than keeping one: neither
is more correct, both wobble at their launch end, and a connection between two singularities ought
to be symmetric. Disc: 12 raw traces → 8, exactly 4 ring + 4 boundary. The field route also picked
up two more cases across the test matrix (29 → 31 of 56).

### Five bugs the guarantee work exposed

Each was silent, and each was found only by checking a property rather than by watching for an
exception.

1. **Backwards faces from `from_polylines`.** It recovers faces from a planar embedding and does
   not guarantee consistent winding. A backwards face breaks strip propagation, so two edges of a
   quad that should share a strip get different densities, and `discrete_coons_patch` then indexes
   off the end of the shorter side — an `IndexError` several frames from the cause. Hexagon: 2 of
   44 faces backwards, 11 faces left mismatched.

2. **Per-face orientation fixes are wrong.** Reversing just the negative-area faces makes every
   signed area positive and makes the mesh **non-manifold** — orientation is shared across edges,
   so flipping one face and not its neighbour leaves them traversing their common edge the same
   way. Unify topologically first (a breadth-first pass over face adjacency; `Mesh.unify_cycles`
   did not return in ten minutes on a few-thousand-face mesh), then flip the whole mesh if its
   total area is negative.

3. **Quad splits chord curved walls.** The new vertex goes at the *chord* midpoint, which on a
   curved boundary lies inside the domain. Every split cuts the corner: an ellipse covered 65% of
   its own area, a disc 74%. Fixed by placing boundary midpoints on the arc — but only where the
   edge genuinely runs *along* the wall. An edge can be a **chord across a hole** with both ends on
   the rim, and snapping that one gave 114% coverage and a non-manifold mesh. Guard: snap only when
   arc length ≈ chord length.

4. **Snapped midpoints can collide.** Two distinct edges landing on the same boundary point puts
   two vertices at one position, which is non-manifold. Keep the occupied positions; fall back to
   the chord when one is taken.

5. **`edges_to_curves` must be total.** `densification` looks up *every* edge once given a mapping
   at all — there is no per-edge straight-chord branch, only a per-call one — so a partial dict
   raises `KeyError`. Fill the gaps with the chord.

### The two limitations that remain

**1. ~~Curved and multiply-connected boundaries fall back~~ — FIXED for curved boundaries; the
diagnosis below was wrong.**

This section used to read: *"A smooth boundary offers nothing to land on, so separatrices must
terminate on each other, and where they near-miss the region does not close. Confirmed not to be a
tolerance problem."* The tolerance half was right. The rest was not.

The network **does** close. `from_polylines` returns a manifold layout and `solve_non_quad_faces`
makes it all-quad. What rejected it was `repair.py`'s `area2 <= 0.0` test, and the reason traces to
`from_polylines` itself: it hands every polyline *segment* to `Mesh.from_lines`, which keys
vertices off segment endpoints and **never computes an intersection**. Two separatrices crossing in
their interiors therefore contribute no node at the crossing, the graph the face search walks is
not a planar embedding, and the faces it recovers **visit the same vertex twice**.
`topological_quad_split` turns each such face into a fan of inverted quads — 58 of 140 on the disc
at `target_length` 0.4 — and `densifiable` rejects the layout.

Measured on the raw `from_polylines` output, before any repair:

| case | interior crossings | raw faces repeating a vertex | route (before) |
|---|---|---|---|
| disc @0.6, @0.5 | 0 | 0 | field |
| disc @0.4 | 3 | 8 of 14, one 14-sided | fallback |
| hexagon @0.35 | 0 | 0 | field |
| hexagon @0.5 | 2 | 1 of 7, 19-sided | fallback |
| ellipse @0.5 | 2 | 1 of 3, 14-sided | fallback |
| stadium @0.5 | 8 | 2 of 5, one 32-sided, non-manifold | fallback |

The hexagon failing at 0.5 and 0.4 and succeeding at 0.35 is the signature: a discrete arrangement
bug comes and goes with resolution, a geometric limit does not.

The stadium had a second, independent bug. `build_network._split_loop` walks the loop between
consecutive marks by **vertex index**, and when two marks sort into the wrong relative order the
walk runs the long way round, so the arcs it emits overlap: 17 arcs covering the wall about one and
a half times, with `(0,0) → (12.93,0.15)` and `(0,0.98) → (12.93,0.15)` both present. No amount of
splitting at crossings repairs a wall that is drawn twice.

`framefield/arrangement.py` fixes all three: it partitions each loop by **arc length** so the wall
tiles exactly once, lands loose ends onto the polyline they nearly touch (per *node*, so the four
arms of a singularity move together), and splits every polyline at every crossing with one shared
node. `decomposition._build` calls it, so `decomposition_polylines`, `decomposition_mesh` and
`edges_to_curves` all see the same network. `decomposition_mesh` then asserts no recovered face
repeats a vertex and reports it through `repair_notes`.

Result: disc @0.6/0.5/0.4/0.3, hexagon @0.5/0.4/0.35, ellipse @0.5, stadium @0.5 all reach route
`field`, with zero crossings, zero dangling ends and 99–100% coverage. The rectilinear plates, the
pentagon, the square and the square-with-hole are unchanged. `14_arrangement.py` measures it.

What genuinely remains is a *quality* cost, not a failure: `solve_non_quad_faces`'s guaranteed step
is **global**, so one non-quad patch quad-splits the entire layout. The hexagon's raw layout at 0.5
is 10 quads and one triangle, and that triangle turns 11 patches into 43; at 0.35 the raw layout is
8 quads and stays 8. A local fix for the odd triangle is the next thing worth doing.

Multiply-connected domains with no corner anywhere are a different problem and are not fixed. A
disc with a round hole produces no separatrix at all, so there is no layout to repair — that is
LIMIT 4, not this one. Two radial cuts from a smooth hole were tried and did not work —
`from_polylines` could not resolve the result into an annulus decomposition, and it broke two
annulus cases that had worked. The code is kept, disabled, in `Tracer.hole_launches` with that
written down.

**2. A cable steers the field but not the mesh — the headline gap against the objective.**
Densification is a `discrete_coons_patch` between coarse edges and never consults the field. The
chain is `cable → field` ✅, `field → layout` only via singularities and corners, `layout → mesh`
field-blind. Measured on a square with a hard, wide-band cable at weight 5: the field bends onto it
to 0.0°, and the dense mesh is identical to the unguided one.

Two ways to close it, and they are not equivalent:

- *(a) make the cable change the layout.* A closed hoop forces winding the walls cannot absorb, so
  singularities appear and the layout does change. But it is a blunt instrument — it only works
  when the cable forces topology, and the patches then need repair.
- *(b) make densification field-aware* — integrate the field inside a patch instead of
  `discrete_coons_patch`. **This is the one that matches the objective**: the cable becomes an
  input to the mesh rather than a post-hoc correction. It is the larger change, since densification
  lives in `compas_singular` and knows nothing about fields.

> `guide_lines.py` is **not** the answer to (2). It snaps polyedges an existing layout already
> contains, so it can only reuse directions that layout already has — its own docstring says
> *"refinement subdivides the lines that exist, it does not invent lines in new directions"*. It is
> post-processing for final adjustment, not a way to generate a cable-aligned mesh.

**3. ~~Two rectilinear plates this section listed as robust are not~~ — FIXED, and the fix is the
interesting part.** `15_baseline.py` found both the first time it ran: `comb-plate` **raised** at
every resolution, and `rect with slot` gave **2 coarse patches for an 8-corner plate**, 100.5%
coverage, and a face carrying a 0.00° and a 180.00° angle — while reporting route `field`, all-quad
and manifold. Both are now 90.0° / 90.0° / AR 1.00 at 100% coverage, 11 and 5 patches.

One cause, and it is the **fourth instance of a single species of bug in this module**:

> a tolerance scaled to a GLOBAL quantity, applied to a domain whose LOCAL feature is smaller than
> that scale.

| # | Where | The global scale | The local feature it swallowed |
|---|---|---|---|
| 1 | `_cluster` landing tolerance | 0.40 | two hexagon landings 0.35 apart |
| 2 | duplicate test, `0.35 * span` | 5.1 on a 14-unit trace | two ellipse arms 4.7 apart |
| 3 | `boundary_corners`, every arc vertex a corner | — | (the stadium fix that introduced #4) |
| 4 | `boundary_corners` `spacing` = perimeter/25 | 2.96 | the comb's 2-unit teeth |

The two numbers that make the point sharper than the pattern does: the rect-with-slot's perimeter
is 50, so `spacing` is **exactly 2.0**, and its slot's bottom segment is **exactly 2.0** — sitting
precisely on the knife edge. And the U-plate survived only because its shortest segment (3.0)
cleared its own gap (2.24) by 1.34×. That is not a margin; that is a domain that had not fired yet.

**The fix does not tune `spacing`.** Tuning it down is the move that produced this situation: the
stadium's arc chords are 0.94 long, so a gap small enough to spare the comb re-breaks the domain
the constant was introduced for. Instead `boundary_corners` gained a second regime above
`repair.SHARP_TURN` (45°), where a turn is a corner unconditionally and is never grouped at all.
`spacing` still governs the genuinely ambiguous band below it, unchanged.

**Why an angle and not a length.** An angle has no length scale, so no small feature can be smaller
than it — the failure mode of all four instances above is structurally unavailable. And
`_densify_loop` only ever subdivides, never making a turn coarser, so a domain that passes at one
background passes at all of them. Measured across every outline in the suite the two populations do
not come close to touching — sampled curves reach 18.0° (the stadium's arc), real corners start at
60.0° (the hexagon) — putting the threshold 2.50× above the first and 1.33× below the second.

**What breaks next.** An outline sampled coarser than 45° per vertex — fewer than 8 points per full
turn — has its samples read as corners. That is a statement about *sampling resolution* rather than
domain size, which is the whole improvement, and 8 points per turn is an octagon, which is a
polygon with real corners rather than a sampled circle. Domains with genuine corners shallower than
45° fall through to the unchanged run-grouping, so that band is not newly at risk.

**The constant is tested, not latent.** `python 15_baseline.py --corners` prints the turn spectrum
per outline and fails if either population comes within 1.25× of the threshold — so the next domain
that lands near 45° announces itself instead of silently picking a regime.

### Step 7 — `edges_to_curves`

Return the traced polyline for each coarse edge so densification follows the separatrix rather than
chording it straight. Larger payoff here than for the boundary-curve case, because interior
separatrices are the alignment the field was solved for.

### Step 8 — compare against `SkeletonDecomposition`

Same input, both front ends, on the pentagon plate. Compare singularity count, worst corner angle,
and the across-guide angle measurement `guide_lines.py` already defines. **If the no-guide layout
is worse than the skeleton's, stop and fix `trace`/`repair`** — no amount of constraint work
rescues a bad tracer.
