# New approach — frame-field front end

A cross/frame-field replacement for `SkeletonDecomposition`'s medial-axis front end.

**The objective:** a quad mesh **generated from** the frame field, where cables and force vectors
are *inputs* that shape the field and the mesh comes out already aligned to them.

That is not the same thing as `src/compas_singular/guide_lines.py`, which pulls an existing mesh's
polyedges onto a curve afterwards. That is post-processing, useful for final adjustment, and it can
only ever reuse directions the layout already has — its own docstring says *"refinement subdivides
the lines that exist, it does not invent lines in new directions"*. `guide_lines.py` is not used
anywhere here. See [12_cables.py](12_cables.py) for both halves of the objective measured:
`cable -> field` in part 1, `field -> mesh` in part 2.

Nothing in `src/compas_singular` is modified: the field produces a `CoarsePseudoQuadMesh` and the
existing pipeline takes it from there.

Run everything with the `singular312` env:

```
cd "compas_singular/examples/New approach"
C:/Users/Casper/anaconda3/envs/singular312/python.exe 10_domains.py
```

Every script opens a **compas_viewer** scene, following the same convention as the rest of
`compas_singular/examples`: one `scene.add_group` per stage or domain, laid out on a grid, with the
camera framed to the whole set. Toggle groups in the scene tree to isolate a stage. Pass
`--no-view` for the numbers only.

## Documents

| File | What it is |
|---|---|
| [cross-field-quad-meshing-evaluation.md](cross-field-quad-meshing-evaluation.md) | Why a field-based front end, which representation, academic sources |
| [frame-field-front-end-module-layout.md](frame-field-front-end-module-layout.md) | Module design, the 8-step plan, status, risk register |
| [paper-diffusion-generated-crossfields.md](paper-diffusion-generated-crossfields.md) | Review of Dai/Qiao/Wang 2026 — what to take, what not to, and the measurement showing our field puts the disc's singularities in the wrong place |
| [site/](site/) | The visual overview page — what the front end can mesh today, generated from `baseline.json`. `python site/build_page.py` rebuilds it; `site/README.md` has the published URL to update rather than replace |

## Examples — start here

| File | Shows | Scene |
|---|---|---|
| **[10_domains.py](10_domains.py)** | What decomposes today, end to end | 8 domains, each with its dense mesh, patch layout and singularities |
| **[11_downstream.py](11_downstream.py)** | The output is an ordinary coarse mesh: strips, densification, grammar, JSON | one layout driving several densities |
| **[12_cables.py](12_cables.py)** | Cables as field input — the field takes them **and so does the mesh** | field with/without the cable, then the two meshes they produce |
| **[22_force_lines.py](22_force_lines.py)** | The **constraint input** explored: ten force-line patterns, the band/weight/mode knobs measured, and the load-path case. **Read before feeding your own cables in.** | 10 patterns × guides, field, layout, mesh |
| **[16_field_densify.py](16_field_densify.py)** | Field-aware densification, measured: the constant-field identity, the cable alignment, the alignment/quality trade | no scene; numbers |
| **[21_edit_coarse.py](21_edit_coarse.py)** | **Editing the coarse layout by hand** between the two halves — drag a corner, slide one along a wall, delete a patch. The field is not re-solved, and a moved edge keeps its separatrix's curvature instead of becoming a chord | 3 domains × generated + 4 edits |
| **[13_limits.py](13_limits.py)** | Every known constraint, with a reproducer. **Read before using this on your own outline.** | each failure beside the equivalent success |
| **[17_vs_skeleton.py](17_vs_skeleton.py)** | The same 16 shapes through **both** front ends — is this better than what it replaces? | field left, medial-axis skeleton right, per domain |
| **[15_baseline.py](15_baseline.py)** | The locked quality baseline — 16 domains × 4 resolutions, and the diff that catches a regression | no scene; numbers and a diff |
| **[25_cache.py](25_cache.py)** | The **solve cache** — a hit reproduces a cold solve exactly, and every way an entry can go stale, asserted one axis at a time | no scene; 44 checks |

## Checks — step-by-step verification

Tests behind the plan's per-step checks, not demonstrations.

| File | Checks |
|---|---|
| [01_field.py](01_field.py) | Triangulation validity, field solve, singularity indices vs. two independent predictions |
| [02_separatrices.py](02_separatrices.py) | Arm counts equal `4 − index`; every trace terminates legitimately |
| [03_layout.py](03_layout.py) | The network reaches `from_polylines` and the existing pipeline accepts the result |
| [19_field_accuracy.py](19_field_accuracy.py) | **The one external check** — disc singularity radius against the published 0.85. The default solver measures ~0.45 and drifts; `CrossField.solve(relax=True)` measures 0.85 |

Each also opens a viewer — the field, the raw traces, and the coarse layouts respectively.

## Reading the scenes

| Colour | Means |
|---|---|
| grey ticks | the cross field (two of the four arms; thinned for legibility) |
| black | domain walls and holes |
| orange | separatrices — the field lines that cut the domain into patches |
| blue dots | coarse patch corners |
| **red** / **blue** dot | singularity wanting valence 3 / valence 5 |
| purple dashed | a cable or force line |

## The two lines that matter

```python
d = FieldDecomposition.from_boundary(outline, guides=cables)
mesh = d.quad_mesh(target_length=0.5)      # ALWAYS returns an all-quad mesh
```

`quad_mesh` is the deliverable. Do not assemble `decomposition_mesh` +
`collect_strips` + `densification` by hand — each of those can fail on a domain whose separatrix
network did not close, and a front end that returns no mesh is useless however good its field was.
`quad_mesh` owns the whole chain and every fallback in it.

Two things to check on the result:

- **`d.route()`** — `'field'`, `'polygon'` or `'triangulation'`. Only `'field'` means the element
  flow follows the field; the other two produce a mesh of the right *shape* that says nothing about
  the cables fed in. **Not a quality signal** — every domain in `14_arrangement.py` returns
  `'field'`, including ones carrying a 180° angle.
- **`d.quality()`** — `min_angle`, `max_angle`, `aspect_max`, `share_below`, `irregular_interior`,
  `poles`, `coverage`, as a dict. **This is the one to gate on.**
- **`d.warnings()`** — everything that degraded, in full. A `QUALITY FAILURE` entry names the
  offending metric and its value.

| Route | When | What you get |
|---|---|---|
| **field** | the separatrix layout closed | element flow follows the field — the point of the front end |
| **polygon** | layout failed, domain simply connected | the domain as one n-gon, split to quads; ignores the field |
| **triangulation** | nothing else available (e.g. multiply connected) | quad-split triangulation; ugly, high-valence, correct |

Verified across 16 domains × 4 background resolutions = 64 runs, recorded row by row in
[baseline.json](baseline.json) and diffed by [15_baseline.py](15_baseline.py).

## The field travels on its own

The layout and the field come out of a decomposition separately, and each one is now a file:

```python
field = d.get_field()                       # a CrossField
coarse = d.decomposition_mesh()             # a CoarsePseudoQuadMesh

field.save_to_json('field.json')
coarse.save_to_json('coarse.json')
```

Reopened later, the field is passed to densification and steers every patch interior:

```python
field = CrossField.load_from_json('field.json')
coarse = CoarsePseudoQuadMesh.load_from_json('coarse.json')

why = field.mismatch(outline, holes, guides=cables, target_length=0.5)
if why:                                     # 'the outer boundary changed', ...
    field = CrossField.from_boundary(outline, holes, guides=cables, target_length=0.5)

coarse.collect_strips()
coarse.set_strips_density_target(0.5)
dense = coarse.densification(field=field)   # patch interiors integrated from the field
```

`save_to_json` / `load_from_json` are the same pair `Mesh` has, with the same two conveniences —
the parent folder is created, and `default=` comes back instead of raising when the file is not
there.

**Always ask `mismatch` before using a loaded field.** A field is a function of its domain, and
densifying a layout with a field solved for a *different* outline is not an error anywhere
downstream: every patch still integrates, the mesh still welds, the result still passes the
quality gate. It is simply aligned to a shape that is no longer there. Nothing else in the
pipeline can catch that, which is why the field stores what it was solved from — the outer
boundary, the holes, the guides and every solver parameter — and why `mismatch` names the first
thing that moved rather than returning a bare `False`.

Two properties worth knowing:

- **The round trip is exact, not close.** `theta` is rebuilt from `u` on load and `max |dtheta|`
  is `0.0`, which `tests/test_workflow.py` pins. Derived state is deliberately not stored: `theta`
  is recomputed, and the point locator is a bucket grid that rebuilds in milliseconds.
- **The field needs nothing else to densify.** It carries its own background triangulation and
  builds its own point locator, so a layout from anywhere — a skeleton decomposition, a hand-built
  mesh, one read back out of a Rhino document — can be densified with it. `CrossField.from_boundary`
  solves one without a decomposition around it, skipping separatrix tracing and the planar
  arrangement, which are the expensive half.

## Not re-solving what has not changed

`framefield/cache.py` stores a solved `FieldDecomposition` in memory and on disk, keyed so that it
invalidates itself. Opt in per call; `from_boundary` is untouched, so every script and harness here
still solves cold.

```python
from compas_singular.framefield import cache
d = cache.solve(outer, holes, guides=guides, target_length=0.5, build=True, verbose=True)
```

Measured on this suite: **5.11 s of solving becomes 0.039 s**, entries 98–482 KB each. The three
Rhino commands that build a decomposition — steps 3, 4 and 6 — now share one entry through
`CMD_start.get_decomposition`, so a pass over an unchanged document solves once instead of three
times, and a Rhino restart does not start over.

**What makes an entry stale.** Four things go into the key, and a miss says which one moved:

| axis | what it covers |
|---|---|
| `inputs` | the outer boundary, the holes and the guides — each **tagged with its role**, so a polyline moved from a hole to a guide is a different problem, not the same one |
| `params` | every `from_boundary` argument that reaches the solver, `relax` / `guide_weight` / `guide_band` / `field_tau` included |
| `code` | the contents of every `framefield/*.py`, so editing `repair.py` re-solves instead of serving a field the current code would not produce |
| `env` | Python, numpy, scipy and compas versions, plus `TOL.precision` |

```
field: cache miss -- framefield source changed since these inputs were last solved
field: cache miss -- solver settings changed -- target_length 0.5 -> 0.4
field: reusing a cached solve (on disk)
```

Every hit is a **freshly reconstructed object**, so two callers never share one — a
`FieldDecomposition` is mutable long after it is built (`edit_coarse` sets `_edited` and replaces
`mesh`; `quad_mesh` sets `dense`), and sharing one across commands is what made warnings report
twice and forced three files to reset `edit_notes` defensively.

Entries live under `%LOCALAPPDATA%\compas_singular\framefield-solve`; set
`COMPAS_SINGULAR_CACHE` to move them, or to an empty string for memory only. `SolveCache().clear()`
empties it, and `.entries()` lists what is there without unpickling anything.

**Turning it off.** The Rhino side reads a `"cache"` key from the document's settings blob,
resolved by `CMD_start.resolve_cache`:

| value | effect |
|---|---|
| `true` (default) | memory in front of disk — steps 3, 4 and 6 solve once between them, and a Rhino restart does not start over |
| `"memory"` | reused while Rhino stays open, **nothing written to disk** |
| `false` | solve every time |

Switching it off should be rare — the cache invalidates itself, so a changed input or a changed
library re-solves without being told. The honest reasons are timing a cold solve and suspecting the
cache rather than the code; for the latter `get_decomposition(force=True)` re-solves once without
touching the setting. Neither removes entries already written: `SolveCache().clear()` does that.

**Do not put the baseline behind it.** `15_baseline.py` and `24_symmetry.py` prove determinism by
re-solving; a cache underneath either would make that a tautology.

## State, honestly

**"A manifold all-quad mesh at 100% coverage" was never the claim it sounded like.** It was true,
it is true of all 64 runs — and it was always compatible with a mesh nobody could use.
Those three properties certified a mesh carrying a literal **180° angle** and an aspect ratio of
38. `route()` used to compensate, because a bad layout usually fell back and said so; since
`arrangement.py` put every test domain on route `'field'`, it no longer discriminates at all.

So robustness is now stated in terms of **element quality**, measured by
[15_baseline.py](15_baseline.py) against a committed per-row baseline. Two tiers, because one
threshold cannot do both jobs. Minimum angle on the field route spans **3.5° to 90.0°** across this
suite without a single degeneracy among those rows (24.4° to 90.0° if you exclude the cable
domains), so any absolute threshold loose enough to pass the ellipse is far too loose to notice the
square slipping off 90°:

| Tier | What it is | Where |
|---|---|---|
| **hard floor** | absolute, permanent: no angle at or near 0° or 180°, no non-finite aspect ratio. Degeneracies, not tastes. | `framefield/quality.py`, checked on every `quad_mesh()` |
| **regression check** | every row of the baseline, every metric, against the committed numbers | `15_baseline.py` |

`quad_mesh()` still **always returns a mesh** — a quality failure never reroutes and never raises.
It comes back as a named entry in `d.warnings()`, and `d.quality()` gives the numbers as a dict.

**Where the 64 runs actually stand:**

| | Rows |
|---|---|
| route `field`, hard floor clear | 57 |
| route `triangulation` (round-holed disc, LIMIT 4), floor clear | 4 |
| route `polygon`, **hard floor breached** | 2 |
| route `field`, **hard floor breached** | 1 |
| raises outright | 0 |

**One domain is still broken**, and the baseline records it rather than hiding it:
**`square+ring cable`** — a 180.00° face at backgrounds 0.5/0.4/0.3, and a route that degrades to
`polygon` at 0.4/0.3. It is also the suite's only source of **poles** (14 in the dense mesh at bg
0.5), which is why it is in the matrix at all: a pole is what a naive quality metric divides by zero
on. `python 15_baseline.py --poles` checks that, and prints what the naive reading would have
given — 0.00° and a division by zero on all 22 pole faces.

**Two domains the harness caught on its first run are now fixed.** `comb-plate` raised at every
resolution and `rect with slot` produced a 0.00°/180.00° face while passing every structural check
there was; both are now 90.0°/90.0°/AR 1.00 at 100% coverage, with 11 and 5 coarse patches. One
cause, and it was the fourth instance of one species: **a tolerance scaled to a global quantity,
applied to a domain whose local feature is smaller than that scale** — here
`boundary_corners`' `spacing` default of perimeter/25 swallowing the comb's 2-unit teeth. The fix
is `repair.SHARP_TURN`: above a 45° turn a vertex is a corner unconditionally and is never grouped.
An **angle has no length scale**, so no small feature can be smaller than it, which is what makes
this immune to the failure mode that produced the other three. Measured across the suite, sampled
curves reach 18.0° and real corners start at 60.0° — the threshold sits 2.50× above the first and
1.33× below the second, and `python 15_baseline.py --corners` fails if either side ever comes within
1.25×. See `boundary_corners`' docstring for what breaks next and why it is a smaller risk.

**The field route — the one that is actually worth having — covers:** rectilinear plates (square,
L, T, U, plus, comb, rectangle with a rectangular hole or slot), the pentagon, and — since
`arrangement.py` — every curved domain tested: disc at 0.6/0.5/0.4/0.3, hexagon at 0.5/0.4/0.35,
ellipse, stadium. Run [14_arrangement.py](14_arrangement.py) for the per-domain numbers.

> **Correction.** This file used to say curved boundaries fell back because "the separatrix network
> does not close", and that it was "confirmed not to be a tolerance problem". The second half was
> right and the first half was wrong. The network *does* close: `from_polylines` returns a manifold
> layout and `solve_non_quad_faces` makes it all-quad. It was rejected because `from_polylines`
> hands raw segments to `Mesh.from_lines` and **never computes an intersection**, so separatrices
> crossing in their interiors got no node at the crossing and the recovered faces repeated
> vertices — which `topological_quad_split` turned into inverted quads and `densifiable` rejected
> at `area2 <= 0.0`. The tell was there all along: the hexagon failed at `target_length` 0.5 and
> 0.4 and *succeeded* at 0.35. A geometric limit does not come and go with background resolution.
> A separate bug in `_split_loop` had the stadium's wall drawn one and a half times over.
> See [framefield/arrangement.py](framefield/arrangement.py).

**Three gaps remain, and none of them is "no mesh" any more:**

1. **Multiply-connected domains with no corner anywhere fall back**, so their meshes ignore the
   field. A disc with a round hole has no singularity and no reflex corner, hence no separatrix at
   all and no layout to repair — this is the LIMIT 4 case, not a crossing problem. Separately, the
   patch *count* is still inflated: `solve_non_quad_faces`'s guaranteed step is global, so one
   stray triangle quad-splits the whole layout (hexagon at 0.5: 11 raw patches become 43).
2. ~~**A cable steers the field but not the mesh.**~~ **Closed.**
   [framefield/densify.py](framefield/densify.py) integrates the field inside each patch instead of
   blending the patch's four edges, so a cable steers a patch interior without the layout having to
   split. On the square with a diagonal cable the mesh went from 36.9° off the cable — a plain
   axis-aligned grid, byte-identical to the unguided one — to 26.9°, and to 25.2° over the middle
   of the cable where it is clear of a wall.

   What remains is a **ceiling, not a gap**: a patch boundary is held fixed so the patches still
   weld and the strip grammar still applies, and on a plain square that boundary *is* the four
   walls. The mesh must turn to meet them however well the field is integrated. Reaching 0° would
   need the layout to change too, which is the other option in
   [12_cables.py](12_cables.py) part 3 and is a separate piece of work.

   Two deliberate exclusions, both measured rather than assumed:
   **pseudo-quad (pole) patches keep their Coons interiors**, because a pole is a collapsed side at
   a singularity and the field has no continuous branch there to integrate; and **element quality
   is spent only where a guide curve asked for it**, so the unguided domains may improve but may
   not be made worse in any respect.

3. **The constraint input is not the bottleneck; the LAYOUT is.**
   [22_force_lines.py](22_force_lines.py) puts ten force-line patterns through one plate. The field
   takes every one of them (0.0–3.8° off the curves fed in, against 36.9° unguided). The layout
   takes about half: patterns the walls can absorb are free — an orthogonal cable grid is 0.0° in
   the field *and* the mesh with minimum angle still 90° — while patterns that force a singularity
   cost minimum angle 90° → 7.2°/12.2°, and three of the ten (arch, tied arch, fan) lose the field
   route entirely.

   One mechanism accounts for the lost separatrices, and it is measured rather than inferred: a
   guide pushes singularities against the walls, and `Tracer.launch_directions` skips any probe
   circle that leaves the domain. A singularity nearer a wall than its smallest radius,
   `0.45 × target_length`, therefore launches **none** of its separatrices. Both failures in that
   file sit inside that radius (0.227 against 0.270; 0.159 against 0.225) and nothing else does.

   **On a guided domain, turn `relax=True` on.** A hard guide under the default solver builds a
   ridge in `|u|` — 1.000 on the cable, 0.604 two bands away — and since the singularities are
   where `|u|` collapses, they slide off it into the walls. `relax=True` removes the ridge, clears
   two of `square+ring cable`'s three hard-floor breaches and restores its field route at
   background 0.4, and is the **only** setting under which real principal-stress load paths beat
   the walls alone. It is still not the default, because it moves baseline rows unrelated to
   cables — that decision is open.
