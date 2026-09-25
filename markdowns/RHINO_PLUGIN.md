# compas_singular in Rhino — plan, feasibility, progress

**Started:** 2026-08-10
**Goal:** drive the whole quad-meshing workflow from Rhino — boundaries in, cables in, frame
field, coarse layout, hand edit, densities out — with the compute done by the existing
`compas_singular` + `framefield` code and only selection/conversion/display written new.

This file is the shared tracker. **§4 is the work list; §6 is the log.** Everything in §1–§2 is
measured, not assumed, and says how to re-measure it.

---

## 0. Verdict

**Feasible, and further along than the plan assumes.** The whole compute stack already runs
inside Rhino 8's own Python interpreter, unmodified — that was the question with the power to
kill the project, and the answer is no.

The remaining work is almost entirely *front end*: conversions, selection, display, and holding
state between commands. Two things in the workflow have no backend at all yet (manual frame
rotation, session state), and one existing subpackage is dead code (`compas_singular/rhino`).

---

## 1. The environment — verified, not assumed

Rhino 8 ships CPython 3.9.10, and this machine's `brg-csd` script environment already carries
the entire scientific stack the field solve needs:

| | version | needed by |
|---|---|---|
| Python | 3.9.10 | — |
| compas | 2.15.1 | everything |
| compas_rhino | 2.15.1 | conversions |
| numpy | 2.0.2 | `field.py`, `densify.py` |
| scipy | 1.13.1 | `field.py` (`factorized`), `densify.py` (`spsolve`) |

Located at `C:\Users\Casper\.rhinocode\py39-rh8\` with the env at
`site-envs\brg-csd-rjJylh8p`.

**Re-check with:**

```bash
PYTHONPATH="C:/Users/Casper/.rhinocode/py39-rh8/site-envs/brg-csd-rjJylh8p;<repo>/src" \
  C:/Users/Casper/.rhinocode/py39-rh8/python.exe -c "import framefield; print('ok')"
```

Run from `examples/New approach/`. This is also the fastest way to test Rhino compatibility
**without opening Rhino** — the interpreter is the same one, so anything that imports and
computes here does so in-app.

### The pipeline computes there

Six domains, `target_length=0.5`, run in the Rhino interpreter. Every one on route `field`:

| domain | patches | dense faces | min angle | max angle | AR | time |
|---|---|---|---|---|---|---|
| square | 1 | 144 | 90.0 | 90.0 | 1.00 | 0.2 s |
| L-plate | 3 | 108 | 90.0 | 90.0 | 1.00 | 0.1 s |
| pentagon | 5 | 157 | 64.7 | 110.4 | 1.58 | 0.5 s |
| disc | 5 | 225 | 51.3 | 145.4 | 1.73 | 0.9 s |
| square + square hole | 8 | 240 | 90.0 | 90.0 | 1.00 | 0.3 s |
| square + cable | 9 | 314 | 61.5 | 129.9 | 2.45 | 0.4 s |

These match `baseline.json` in character — numpy 2.0.2 / scipy 1.13.1 here against 2.5.1 / 1.18.0
in `singular312` makes no difference to the routes taken.

### It is fast enough to be a command, not fast enough to be a slider

A 20 × 14 m plate with one reflex corner and a cable:

| background spacing | background vertices | field + trace | layout + densify | dense faces |
|---|---|---|---|---|
| 1.0 | 258 | 0.2 s | 0.2 s | 231 |
| 0.5 | 936 | 0.4 s | 0.8 s | 853 |
| 0.3 | 2 438 | 2.1 s | 4.4 s | 2 344 |

Design consequence: **the field is solved once and reused.** Re-solving on every density change
would put a 2–6 s stall in the interaction. The field and the tracer must outlive a single
command — see §4.2.

---

## 2. The workflow, step by step

| # | Step | Backend | Rhino side | Verdict |
|---|---|---|---|---|
| 1 | define boundaries | `FieldDecomposition.from_boundary(outer, inners, target_length=)` | curve → point-list conversion, outer/inner sorting | **ready**, needs conversion |
| 2 | cables / force lines (opt.) | `guides=`, `mode=`, `guide_weight=`, `guide_band=` | selection + conversion | **ready**, needs conversion |
| 3 | generate frame field | `CrossField.solve` via `from_boundary` | no display of crosses | **ready**, needs display |
| 3b | **rotate frames manually** (opt.) | **nothing** | **nothing** | **GAP — §4.1** |
| 4 | coarse mesh — field route | `decomposition_mesh()` | draw the layout | **ready** |
| 4 | coarse mesh — skeleton route | `SkeletonDecomposition` | draw the layout | **ready**, see caveat below |
| 5 | edit coarse mesh (opt.) | `edit_coarse()` / `quad_mesh(coarse=)`, `edit.face_polylines` | bake corners, read back | **built, backend included** — `framefield/edit.py` + `FieldDecomposition.edit_coarse`, driven by `rhino_plugin/coarse_edit.py` (vertex moves). Polylines already work as an input, so topology edits need only the bake |
| 6 | set densities | `collect_strips`, `set_strips_density_target`, `set_strips_density` | per-strip selection UI | **ready**, biggest UI job |
| 6b | densities follow the field | `field_aware` + `densify()` | a checkbox | **ready**, see caveat below |

### Step 5 is not a guess — it was designed for this

`framefield/edit.py`'s own docstring says *"the generated layout goes out to Rhino, comes back
edited"*, names `compas_rhino.conversions.mesh_to_compas` as an accepted input, and
`edit_coarse` deliberately **raises instead of falling back** so a bad edit fails *"while they
are still in Rhino and can fix it"*. It accepts a mesh **or** one closed polyline per patch, and
re-welds topology from corner coordinates, so deleting or splitting a patch just works.

One rule the Rhino side must honour: **one point per corner, not one point per sample.** A patch
baked with its curved separatrix edges comes back as a 40-gon and gets fanned into 38 patches.
`face_polylines` bakes corners only; the curved separatrices go out separately as reference
geometry.

### Step 6b works across routes — and costs quality

The requirement is *"densities should follow the frame field even if the frame field was not used
to generate the coarse skeleton."* **This works today**, with no new code: solve a
`FieldDecomposition` on the same boundary, hand it the skeleton's layout, and every patch relaxes
against the field.

Measured on a square-with-square-hole plate, cable across it, skeleton layout of 12 patches:

| | min angle | max angle | AR | mean deviation from the cable |
|---|---|---|---|---|
| `field_aware=False` (Coons) | 35.5° | 145.1° | 1.99 | 32.7° |
| `field_aware=True` | 30.3° | 149.6° | 3.20 | 31.8° |

**Element quality got worse and alignment barely improved.** That is not a bug — `densify` only
spends quality where a guide asks for it, and a skeleton layout's *edges* do not follow the field,
so interior relaxation fights the patch boundaries it was handed. The field-aware toggle is
therefore a real feature on the field route and a **marginal one on the skeleton route**. Do not
default it on for skeleton layouts, and show the user both numbers.

### The skeleton route needs a densely sampled outline

`SkeletonDecomposition` on a 4-corner square raises
`KeyError: popitem(): dictionary is empty` from `Mesh.boundaries()`. Its medial axis is the
circumcentres of a Delaunay triangulation whose *only* vertices are the boundary points, so a
corners-only outline has no interior to work with. The field route has the opposite requirement
and builds its own interior — see the design doc's finding 1.1, *"the two front ends want
opposite triangulations from the same boundary."*

**Consequence for the plugin:** the outline must be sampled *before* it reaches the skeleton
route, at roughly the coarse feature size. Rhino's `rs.DivideCurve` handles this; do not send a
`rs.CurveEditPoints` outline straight through.

---

## 2b. Rhino's interpreter is long-lived — and its `sys.path` is not

**Hit on 2026-08-10, and it will be hit again.** `ff_02_field.py` died with

```
ff_common.py line 287, in get_decomposition
    from framefield.field_decomposition import FieldDecomposition
ModuleNotFoundError: No module named 'framefield'
```

while `ff_01_boundaries.py` had just run fine, and `ff_common`'s own module-level
`from compas_singular.datastructures import ...` — on a path fixed the same way — kept working.

The two halves behave differently, and that is the whole bug:

* **`sys.modules` stays warm** across script runs, so the second time a script does
  `from ff_common import ...` the module is served from cache and **its body does not re-execute**.
  `setup_paths()` therefore runs *once per Rhino session*, not once per run.
* **`sys.path` is re-initialised** by the script engine between runs.

So module-level imports survive — they were bound in run 1 and are just names now — while a
**late import inside a function**, reached for the first time in run 2, looks up a `sys.path` that
no longer has the path on it. The error names the module and says nothing about paths.

Reproduced and fixed, both verified in the Rhino interpreter:

```
run 2, before the fix:  import framefield -> No module named 'framefield'
run 2, with ensure_paths() first:  OK
```

**The fix, and the trap inside the fix.** `setup_paths()` did two things: insert the paths *and*
purge `compas_singular` from `sys.modules`. Calling it from inside `get_decomposition` — the
obvious one-liner — would purge the modules mid-session, after `ff_common` had already bound
`CoarsePseudoQuadMesh`. Measured:

```
same class object after purge: False
old instance isinstance(new) : False
```

A layout built before the call would stop being an instance of the class checked after it — a
much worse failure than the import error, and one that would surface far from its cause. So the
two jobs are now split:

* **`ensure_paths()`** — `sys.path` only, idempotent, safe anywhere. **Call it immediately before
  every late import.**
* **`setup_paths()`** — `ensure_paths()` plus the one-time purge. **Module-body use only.**

`rhino_plugin/coarse_edit.py` had the same latent bug in `_closest_on_loop` and now ensures the
path itself, derived from `__file__` rather than hard-coded.

> **§4.5 dissolved this, 2026-08-28.** `framefield` is now a subpackage of the installed
> `compas_singular`, so there is no separate path to insert for it and no framefield-specific late
> import to protect — it resolves the moment `SINGULAR_SRC` is on `sys.path`, same as any
> other `compas_singular` subpackage. `CMD00_start.ensure_paths()` still inserts `SINGULAR_SRC`
> itself (that part was never framefield-specific), and its purge in `import_compas_singular()`
> now skips `compas_singular.framefield` by name, so the in-session solve-cache tier from §4.2
> survives a purge instead of being thrown away on every command.

---

## 3. What is dead, and what has to move

### `src/compas_singular/rhino/` does not import

```
compas_singular.rhino  ->  ModuleNotFoundError: No module named 'compas_rhino.artists'
```

It is COMPAS 1.x code. `compas_rhino.artists` was replaced by `compas_rhino.scene` in COMPAS 2,
and `compas.utilities` by `compas.itertools`. Affected: `rhino/artists/*`,
`rhino/grasshopper.py`, `rhino/geometry/*`. (`rhino/constraints/*` was removed on 2026-09-16; its
surface projection now lives in `datastructures/mesh/projection.py`.)

Do **not** try to revive it wholesale. Take from it only what is still wanted — the
polyedge/strip selection helpers in `rhino/artists/` are the genuinely useful part, since
"select a strip and set its density" is step 6's core interaction and that logic is non-trivial.

The `.ghuser` set in `examples/grasshopper/user_objects/` is from the same era and is likewise
COMPAS 1.x.

### `framefield` lived in `examples/` — moved 2026-08-28

`examples/New approach/framefield/` used to be on `sys.path` only because a script run from that
folder put it there, which a plugin could not rely on. Milestone 1 held up — 64 baseline
rows said so — so the design doc's anticipated move happened: `framefield` is now
`src/compas_singular/framefield/`, imported as `compas_singular.framefield`, installed with the
rest of the package. `15_baseline.py` was re-run before and after; all 84 rows matched
byte-for-byte. 27 consumer files (this repo and the sibling `compas_topology` commands) were
repointed at the same time. See §4.5.

---

## 4. Required actions

Ordered. Each is a thing that does not exist yet.

### 4.1 Manual frame rotation — the one genuine backend gap

Nothing supports it today:

* `Constraint` is `namedtuple('Constraint', 'vkey direction weight')` — per **vertex**, `weight=None`
  meaning hard. The mechanism to pin a direction exists.
* But `FieldDecomposition.from_boundary` builds its constraint list locally at
  [decomposition.py:244](../compas_singular/examples/New%20approach/framefield/decomposition.py#L244)
  and **throws it away**. There is no `self.constraints`, no re-solve entry point, no
  constructor that takes a field.

Needed, smallest first:

1. Store the constraint list on the decomposition.
2. `add_constraints(...)` / `resolve()` — re-solve the field with extra pins, re-trace, invalidate
   the cached network. Layout and edits downstream must be explicitly discarded, not silently kept.
3. Decide what the user actually rotates: a **region** (pick a point + radius + angle, emit soft
   constraints over the vertices in it) is more usable than a single vertex and matches
   `from_curves`' existing `band` parameter.
4. Rhino side: draw the crosses (§4.3), pick, gumball or angle prompt, re-solve, redraw.

**Watch for:** a cross is 4-fold symmetric, so "rotate by 90°" is a no-op by construction. The UI
must clamp to (−45°, +45°] or the user will think it is broken.

### 4.2 Session state — **the Rhino document holds it; JSON carries the inputs**

**Decided (2026-08-10):** each step bakes its result as Rhino geometry on its own layer, and that
geometry *is* the state the next step reads. A small JSON rides along in document user text
holding the things geometry cannot express. Nothing derived is ever serialized.

Three measurements decided the shape of this:

| what | as JSON | rebuild cost |
|---|---|---|
| inputs (outline, guides, `target_length`, mode, rotations) | **222 bytes** | — |
| coarse layout, 3 patches | 304 bytes | — |
| separatrix network, 10 polylines | 3.2 KB | — |
| background + field, 936 vertices @ 0.5 spacing | **92 KB** | 0.4 s |
| same @ 0.3 spacing (2 438 vertices) | ~240 KB | 2.1 s |

**And the solve is deterministic** — two independent builds gave identical singularities and
identical coarse vertices to 9 decimals. So a rebuild reproduces the state exactly; there is no
drift to defend against.

So:

* **Never serialize the field.** 92–240 KB to store something a deterministic 0.4–2.1 s call
  reproduces exactly. Rebuild it; cache the live `FieldDecomposition` in `scriptcontext.sticky`
  keyed by a hash of the inputs JSON, so within a session it is solved once.
* **Store the 222-byte inputs blob** in document user text. It survives save/reopen, travels with
  the file, and is the whole recipe.
* **Bake everything visible**, because that is the point — and because the layout and the
  densities have to come back *edited*, which only geometry can express.

#### Layer plan

| layer | baked | read back? |
|---|---|---|
| `ff::boundary` | outline + inner boundaries | yes — inputs |
| `ff::guides` | cables / force lines | yes — inputs |
| `ff::field` | cross lines per background vertex, decimated | no — display only |
| `ff::singularities` | points, coloured by index sign | no — display only |
| `ff::separatrices` | traced polylines | no — rebuilt; bake for reference |
| `ff::coarse` | patch polylines (`edit.face_polylines`) | **yes — this is the edit seam** |
| `ff::mesh` | the dense quad mesh | no — output |

Manual frame rotations (§4.1) go in the inputs blob as `(point, radius, angle, weight)`, not as
geometry. They are constraints on a solve, not objects in a model, and baking them would make the
field's own display the thing you edit — which is the loop `edit_coarse` deliberately avoids.

#### **Densities must be keyed geometrically, never by strip index**

The trap in this design, and it is a silent one. `collect_strips` builds its list from
`self.edges()` and consumes it with `edges.pop()`, so **strip indices follow vertex insertion
order** — and a layout welded back from Rhino has whatever order the geometry came back in.

Measured: bake a 4-strip layout, shuffle the patches as a selection would, read it back —
**0 of 4 strip indices point at the same strip.** The permutation was 0→2, 1→3, 2→0, 3→1. A
density map keyed by index does not error; it silently applies the wrong numbers.

Two identities that do survive, both verified 4 of 4:

* **the set of rounded midpoints of every edge in the strip** — order-independent, use it to match
  a whole density map;
* **any one remembered edge** — store `(gkey(u), gkey(v))`, find that edge in the rebuilt mesh,
  ask which strip owns it. This is also exactly what "user picks an edge and sets its density"
  needs, so it is the natural one for the UI.

Round the same way the rest of the codebase does — 3 decimals, matching `geometric_key` and
`from_polylines`' endpoint matching.

#### Invalidation

Key the sticky cache on a hash of the inputs blob, and store the source curve GUIDs in it. A
moved boundary changes the blob, misses the cache, and re-solves — rather than silently meshing
the old outline. When the field is re-solved, an already-edited coarse layout is no longer valid
against it: say so and make the user re-bake, do not warp it silently.

### 4.3 Display

Nothing draws any of it. Needed, in the order the workflow needs them:

* **crosses** — line pairs per background vertex, decimated; prerequisite for 4.1
* **singularities** — points, coloured by index sign; the single most informative thing on screen
* **separatrices** — polylines; already available as `decomposition_polylines()`
* **coarse layout** — as a mesh *and* as `face_polylines` for baking
* **dense mesh** — `compas_rhino.conversions.mesh_to_rhino`
* **quality** — `quality()` returns a dict; colour faces by angle and put `warnings()` in the
  command line. Do not let a `QUALITY FAILURE` warning stay invisible: `quad_mesh` **always**
  returns a mesh, so a degenerate one arrives looking like a success.

### 4.4 Conversion + selection layer

The thin, unavoidable glue. `curve → list[[x,y,z]]` with sampling control; outer/inner boundary
sorting (largest area outer, or ask); guide curves; picking a strip; picking a patch.
`compas_rhino.conversions` covers the mesh directions.

### 4.5 Promote `framefield` to `src/` — **done, 2026-08-28**

Flat folder move to `src/compas_singular/framefield/`, no edits inside the package (every
intra-package import was already relative). `setup.py` switched from a hard-coded single-package
list to `find_packages(where="src")`, so the 23 subpackages — previously unshipped — are
now included too. `15_baseline.py` re-run before and after: all 84 rows byte-identical. 27
consumer files repointed at `compas_singular.framefield`, across this repo and the sibling
`compas_topology` commands (`CMD00_start.py`, `CMD_symmetry_report.py`, `coarse_edit.py`); the three
`NEW_APPROACH_PATH` `sys.path` insertions in `compas_topology` were dead code once the move
landed and were removed with it.

### 4.6 Repair or delete `compas_singular/rhino/`

Port the strip/polyedge selection helpers to COMPAS 2; delete or clearly mark the rest. Leaving a
subpackage that raises `ModuleNotFoundError` on import in the shipped package is worse than not
having it.

---

### 4.7 Carrying the workflow across two commands — what a bake loses

**Written 2026-09-02**, after `SkeletonDecomposition` and `FieldDecomposition` were
given the same shape. In a script the whole front end is four lines and every
intermediate is just a variable:

```python
decomposition = SkeletonDecomposition.from_boundary(...)   # or FieldDecomposition
coarse = decomposition.coarse_mesh()
coarse.set_strips_density_target(t=0.5)
dense  = coarse.densification(overwrite_edges_to_curves=..., field=...)
```

Lines 2 and 4 are **different Rhino commands** — step 3 and step 6 — so `coarse`,
the curves and the field have to cross a document boundary between them. This
section is what that costs and what is already paid for.

#### What a baked mesh is, and what it is not

`CMD02_coarse_mesh` bakes the layout with `rs.AddMesh`, and a Rhino mesh is
**vertices and faces**. Everything a `CoarsePseudoQuadMesh` knows lives in
`mesh.attributes`, and none of it is expressible as geometry:

| `attributes` key | what it is | survives the bake? | how step 5/6 gets it back |
|---|---|---|---|
| *(vertices, faces)* | the layout itself | yes, at **float32** — `rs.AddMesh` stores `Point3f`, ~2e-6 off on a 20-unit plate | `read_coarse` |
| `face_pole` | which corner of a pseudo-quad is collapsed | **no** | re-derived from the separate `Skeleton::Poles` point layer, via `from_vertices_and_faces_with_poles` |
| `strips` | the bands densities are defined on | **no** | re-collected — **with permuted keys**, see §4.2 |
| `strips_density` | the user's actual work | **no** | re-applied from document user text, keyed by rounded edge midpoint |
| `polyedges` | | **no** | re-collected |
| `edge_coarse_to_dense`, `vertex_coarse_to_dense` | coarse→dense correspondence | **no** | never populated by `densification` in the first place — only `from_quad_mesh` writes them |
| `quad_mesh` | the dense result | **no** | re-densified |
| — | the decomposition's branch polylines | **no** | baked to `Skeleton::Polylines`, read by `read_polylines` |
| — | the domain walls / `inputs` | **no** | re-read from the `Outer` / `Inner` layers |
| — | the cross field | **no** | re-solved, disk-cached |

So the document does hold all of it — just scattered across four layers, one
user-text blob and a disk cache, with each command re-deriving what it needs.

#### How compas_singular stored it originally: it did not have to

The library's own answer is `to_json`, and it carries **everything**. Measured
2026-09-02 on a 9-strip plate with a point feature: strips, per-strip densities
including a hand-set 9, and all 3 `face_pole` entries survived a
`to_jsonstring()` / `from_jsonstring()` round trip, in **2.2 kB**. The shipped
`examples/data/coarse_quad_mesh_british_museum_poles.json` carries the same set —
`strips`, `strips_density`, `polyedges`, `face_pole`, both coarse→dense maps.

That is worth knowing but is **not** the fix here, for the reason §4.2 already
gives: the layout has to come back **edited**, and only baked geometry can
express a user's edit. JSON would be a second source of truth for the one thing
the user changes. It is available for everything the user does *not* edit.

#### What step 6 actually still needs — one item, not three

* **The field — already solved.** `CMD00_start.get_decomposition()` re-solves and
  the disk cache (four-axis digest: inputs, params, `framefield` source, env)
  makes a hit nearly free. Measured 95 kB pickled, essentially all of it the
  background mesh, against a 0.11 s cold solve on a 12×8 plate at 0.5 (2.1 s on
  the 20×14 plate at 0.3). **Storing it in the document would be the wrong
  trade** — that was already the §4.2 verdict and it still holds. For the
  skeleton route the missing piece is smaller still: `CrossField.from_boundary`
  takes exactly what `read_boundaries` already returns.
* **The densities — already stored.** Document user text, keyed geometrically,
  and since 2026-09-02 there is **one** `apply_densities` (`density_common.py`)
  instead of a byte-identical copy in each of steps 5 and 6.
* **`edges_to_curves` — the only gap, and it is recomputation, not storage.**
  Step 6 gets it today by handing the read-back layout to
  `FieldDecomposition.edit_coarse`, which re-matches each edge to the separatrix
  it came from. That is field-only. But `coarse_edges_to_curves` re-derives the
  same mapping **from the document alone** — walls plus the baked
  `Skeleton::Polylines` — which is exactly what `CMD02_coarse_mesh` already does
  after its own bake. The work is to call it in step 6 as well.

#### The change to `CMD06_quad_mesh` — **done 2026-09-03**

Implemented as below. `CMD00_start` gained `cache_directory()` / `cache_path()` and
`read_layout()`; step 3 writes `coarse.json` (+ `field.json` on the field route)
after its bake, step 4 refreshes `coarse.json` after a commit, steps 5 and 6 read
through `read_layout()`. Step 6 no longer calls `get_decomposition()` at all —
which is what stops the skeleton route paying for a field solve it cannot use.

`read_layout()` prefers the side-car but only after checking it still describes
what is on the layer, comparing the **rounded vertex set** (a bake renumbers, so
order proves nothing; `Point3f` means exact comparison never survives). It falls
back to the document, out loud, when something moved the mesh outside these
commands.

Two things the implementation turned up:

* **`mismatch` has to be given the same resolvers step 3 solved with.** Checking
  a guided document's field against the `relax=False` default reports
  `solver settings changed -- relax True -> False` and discards a perfectly good
  field. It takes `resolve_relax(settings, guides)` and `resolve_field_symmetry(settings)`.
* **Identical quality numbers do not mean the field did nothing.** With and
  without `field=`, this domain reported the same min angle and aspect — because
  the worst face sits on a patch BOUNDARY, which `edges_to_curves` fixes either
  way. 308 of 543 vertices had moved, by up to 0.37. Compare coordinates, not
  metrics, when checking whether a field is reaching the interiors.

Round trip verified headlessly on both routes: the layout's `route` and its
hand-set density survive, the loaded field is bit-identical to the live one
(max |Δu| = 0.0), a moved outer wall is caught by `mismatch`, and the skeleton
route densifies with `chord: 0`.

#### The change as originally specified

```python
settings = get_settings()
coarse, poles = read_coarse()

wall_sampling = settings["triangulation_spacing"] * 0.25
outer_loop, inner_loops = read_boundary_loops(wall_sampling)
snap_corners_to_walls(coarse, loops=[outer_loop] + inner_loops)
curves, tally = coarse_edges_to_curves(
    coarse, loops=[outer_loop] + inner_loops,
    polylines=read_polylines("Skeleton::Polylines"))

apply_densities(coarse, settings["target_length"])

field = get_decomposition().get_field() if settings["field_aware"] else None
dense = coarse.densification(overwrite_edges_to_curves=curves, field=field)
```

Same shape as the script, and **it gives the skeleton route field densification,
which it cannot have today.**

What it trades:

* **loses `edit_coarse`'s warp branch.** An interior edge whose corner the user
  *dragged* in step 4 keeps its curvature under `edit_coarse` and falls back to a
  chord under `coarse_edges_to_curves` — the fourth branch that needs the live
  traced network. So the field route should keep `edit_coarse`; this is the
  skeleton route's path, not a replacement for both.
* **gains the wall arcs at document resolution.** `read_boundary_loops` re-reads
  the Rhino curves four times finer than the background, which is better than
  anything a decomposition can offer from its stored loops.

#### The JSON side-car — `save_to_json` / `load_from_json`

**Added 2026-09-02**, on `compas_singular.datastructures.mesh.Mesh`, so every
mesh in the library has them. `to_json` carries the whole `attributes` dict, so
this is the one mechanism that moves the layout's *knowledge* between commands
rather than only its shape:

```python
coarse.attributes['route'] = 'skeleton'                 # a plain dict; it round-trips
coarse.save_to_json(mesh_cache_path('coarse'))
coarse = CoarsePseudoQuadMesh.load_from_json(mesh_cache_path('coarse'))
```

The split is deliberate: **saving and loading are about the mesh** and live on
the class, where `load_from_json` reads as the constructor it is; **only "which
folder, for which document" is Rhino-shaped**, and that is `cache_directory()` /
`mesh_cache_path()` in `rhino/helpers/helpers.py`.

Four things they do that the bare `to_json` / `from_json` pair does not:

* **`save_to_json` creates the parent folder.** `to_json` raises
  `FileNotFoundError` on a path whose directory does not exist yet, which a
  first run always is;
* **`load_from_json` repairs the strip keys.** JSON object keys are strings, so
  `attributes['strips']`, `['strips_density']` and `['polyedges']` come back
  keyed `'0'`, `'1'`, ... The mesh still densifies — everything internal reads
  whatever keys it finds, verified identical face counts — but
  `get_strip_density(0)` raises `KeyError` and `sorted(strips())` orders `'10'`
  before `'2'`. Both silent. Digit-string keys become ints again on load;
* **`default=` instead of an exception** for a cache that has not been written
  yet, which is the normal state of a fresh document;
* **`cache_directory()` is per DOCUMENT, not per script.** A `cache/` folder next
  to the command files is shared by every model Rhino opens, so two open
  documents silently overwrite each other's layout. This is a hidden folder
  beside the `.3dm`, falling back to temp for an unsaved document.

Two things that already worked and are worth not working around: **the class
comes from the file** — compas dispatches on the stored `dtype`, so a
`CoarsePseudoQuadMesh` comes back as one however it is asked for, `Mesh` included
— and **`cls` still acts as an assertion**, because `Data.from_json` refuses
anything that is not an instance of it.

Sizes measured: a 6-patch layout **1.7 kB**, the same with densities and poles
2.3 kB, a 384-face dense mesh 31 kB, a 2400-face one 232 kB. So the LAYOUT is
small enough for document user text if travelling inside the file matters more
than being readable on disk; the dense mesh is not, and does not need to be — it
is baked.

**It is a side-car, not the source of truth.** The baked layout is what the user
edits and only geometry can express an edit, so the JSON has to be rewritten
after every bake and treated as stale otherwise. §4.2's rule stands.

#### The field is not a candidate for it — measured

`CrossField` **cannot** be JSON'd as it stands, and making it so is real work:

| | |
|---|---|
| `CrossField` is a compas `Data`? | **no** — no `to_json` |
| `BackgroundMesh`, `Symmetry`? | **no** — only `background.mesh` is |
| `field.u` | dict of **complex** — `TypeError: Object of type complex is not JSON serializable` |

So it would need `__data__`/`__from_data__` on three classes plus a real/imag
split. Against: 95 kB pickled (77 kB of it the background triangulation, 11.5 kB
the field values) for something a **0.11 s** solve reproduces exactly on a 12×8
plate at 0.5 — 2.1 s on the 20×14 plate at 0.3 — and which is already disk-cached
behind a four-axis digest (inputs, params, `framefield` source, environment) that
no hand-written key would match.

**Keep the field cached separately, as §4.2 decided.** What belongs in the JSON
is the *recipe*: the boundaries and the solve parameters. With
`CrossField.from_boundary` that recipe is now four arguments, and
`coarse.attributes['solve']` is the natural place for it.

#### What is genuinely not stored anywhere

1. **Which route built the layout.** Step 6 cannot tell a skeleton layout from a
   field one and assumes field — it calls `get_decomposition()` unconditionally,
   so choosing *Skeleton* in step 3 still pays for a field solve in step 6. One
   string in the settings blob fixes it, and it is the prerequisite for the
   change above being conditional rather than a rewrite.
2. **`polyline_features`.** `CMD02_coarse_mesh` passes the guides as
   `polyline_features` when *Skeleton* is chosen, so they cut the domain and the
   layout is built around them (see `HOW_IT_WORKS.md` §5). `CMD01_boundary_selection`
   still drops them: its `polyline_features` is a one-point stub that is never
   passed on.
3. **The domain the layout was built from.** Every step re-reads `Outer`/`Inner`
   rather than remembering, which is right — a moved curve should invalidate —
   but nothing checks that the layout still matches the walls it was built on. A
   corner off its wall shows up two steps later as a chorded boundary edge and a
   coverage below 1.

## 5. Open decisions

| # | Decision | Notes |
|---|---|---|
| D1 | **Delivery vehicle** — Rhino commands (Python scripts + `rs.Get*`), Grasshopper components, or an Eto panel | Commands suit the sequential workflow and match how the code is already shaped. GH fights step 5 (the hand edit) and step 6 (per-strip densities). A panel is the best UX and the most work. *Assumed: commands, for v1.* |
| D2 | Where the plugin code lives | Separate folder in this repo? Separate repo? Affects 4.5. |
| D3 | Which Rhino script env | `brg-csd` already has everything, but is shared with other work. A dedicated env pinned by a `# requirements:` header is cleaner. |
| D4 | Does the plugin need Rhino 7 / IronPython? | **If yes, the field route is impossible** — no numpy, no scipy. Skeleton route only. Assumed no. |
| D5 | Surfaces, or planar only? | `framefield` is planar (`background.py` returns identity tangent bases). `compas_singular` has `surface_discrete_mapping` for the skeleton route. Assumed planar for v1. |

---

## 6. Progress log

| Date | What | State |
|---|---|---|
| 2026-08-10 | Rhino 8 CPython env inventoried: compas 2.15.1, compas_rhino, numpy 2.0.2, scipy 1.13.1 | ✅ verified |
| 2026-08-10 | `framefield` + `compas_singular` import under `py39-rh8` | ✅ verified |
| 2026-08-10 | Full pipeline computes there — 6 domains, all route `field` (§1) | ✅ verified |
| 2026-08-10 | Scaling measured on a 20 × 14 m plate (§1) | ✅ verified |
| 2026-08-10 | Skeleton layout → field densification cross-route works, quality cost measured (§2) | ✅ verified |
| 2026-08-10 | `compas_singular.rhino` confirmed broken under COMPAS 2 | ❌ found |
| 2026-08-10 | Manual frame rotation confirmed absent (§4.1) | ❌ found |
| 2026-08-10 | No serialization anywhere in `framefield` (§4.2) | ❌ found |
| 2026-08-10 | Field solve is deterministic — identical singularities and coarse vertices across builds | ✅ verified |
| 2026-08-10 | State sizes measured: inputs 222 B, field 92 KB @ 0.5 spacing (§4.2) | ✅ verified |
| 2026-08-10 | **Strip indices do not survive a bake/read-back — 0 of 4.** Two geometric identities do — 4 of 4 (§4.2) | ⚠️ found |
| 2026-08-10 | Storage design decided: geometry in the document, 222-byte inputs blob in user text, field rebuilt + cached | ✅ decided |
| 2026-08-10 | Step 5 vertex-move command written — `rhino_plugin/coarse_edit.py`, edits the in-memory mesh and commits through `edit_coarse` | ✅ built |
| 2026-08-10 | Round trip verified headless on a disc: 5 patches in / 5 out, coverage 1.0, route still `field`, 49 of 49 dense boundary vertices on the wall | ✅ verified |
| 2026-08-10 | Six-step command series in `rhino_plugin/scripts/`: `ff_01_boundaries` … `ff_06_quad_mesh` + `ff_common`. Closes 4.2 (settings + densities in user text, field rebuilt and sticky-cached), 4.3 (crosses, singularities, separatrices, layout, dense mesh, quality) and 4.4 | ✅ built |
| 2026-08-10 | Full series dry-run with Rhino stubbed, disc + hole + cable: field route 36 patches → edit → per-strip density → 512 faces; skeleton route 14 patches, coverage 1.00. Density override survived the round trip, and was correctly reported LOST when the layout changed under it | ✅ verified |
| 2026-08-10 | **`edit_coarse` renumbers vertex keys — 2 of 8 unchanged on a disc.** Same weld as the bake/read-back, so §4.2's geometric keying applies after a commit too | ⚠️ found |
| 2026-08-10 | `ff_02_field.py` `ModuleNotFoundError: framefield` — Rhino keeps `sys.modules` warm but re-inits `sys.path`, so late imports break while module-level ones survive (§2b) | ✅ fixed |
| 2026-08-10 | `setup_paths` split into `ensure_paths` (safe anywhere) + `setup_paths` (body only) — the naive fix broke `isinstance`, measured (§2b) | ✅ fixed |
| 2026-08-11 | **Step 5's backend was missing** — `rhino_plugin/coarse_edit.py` and `ff_04`/`ff_06` were written against `framefield.edit` and `FieldDecomposition.edit_coarse`, neither of which existed in the tree. Written: `framefield/edit.py` (`coarse_from_skeleton`, `faces_from_geometry`, `mesh_from_faces`, `snap_to_loops`, `closest_on_loop`, `warp_polyline`, `face_polylines`) + `edit_coarse`, `quad_mesh(coarse=)`, and the warp branch in `edges_to_curves` | ✅ built |
| 2026-08-11 | Plugin's API surface exercised under Rhino 8's own CPython 3.9: the `_closest_on_loop` shim resolves, `edit_coarse` takes the in-memory mesh, every `edit_notes` key the plugin prints is present, `ff_06`'s `strict=False` polyline route densifies. Disc: route `field`, coverage 0.994, edges 7 exact / 2 boundary / 3 warped / **0 chord** | ✅ verified |
| 2026-08-11 | **The warp branch must be gated on "was this layout edited".** Left ungated it also fired on GENERATED layouts, where an unmatched edge is deliberately a chord: it moved 3 of the baseline's 84 rows and gave `square+circle @ 0.60` a 0.001° element. Gated, all 84 rows match | ⚠️ found |
| 2026-08-11 | `edit_coarse` returns a mesh with **no strips collected** — `densify()` then raises `KeyError` from `get_strip_density`. `quad_mesh` collects them itself; `ff_06` gets them from `apply_densities` | ⚠️ found |
| 2026-08-28 | `framefield` promoted to `src/compas_singular/framefield/` (§4.5) -- flat move, no internal edits; `setup.py` switched to `find_packages(where="src")`; 27 consumer files repointed; `15_baseline.py` 84/84 rows byte-identical before/after | ✅ done |
| | 4.5 promote `framefield` to `src/` — also dissolves §2b | ☑ done 2026-08-28 |
| | 4.2 implement the inputs blob + sticky cache + geometric density keying | ☑ in `ff_common` |
| | 4.4 conversion + selection layer | ☑ in `ff_common` |
| | 4.3 display | ☑ except face colouring by angle |
| | 4.1 manual frame rotation | ☐ |
| | 4.6 repair/delete `compas_singular/rhino/` | ☐ |

---

## 7. Related documents

* [frame-field front end — module layout](../examples/New%20approach/frame-field-front-end-module-layout.md) — the design, the step plan, and the implementation status of the field route
* [examples/New approach/README.md](../examples/New%20approach/README.md) — what each numbered example demonstrates
* [baseline.json](../examples/New%20approach/baseline.json) — 64 locked rows; regenerate and diff with `15_baseline.py`
* [HOW_IT_WORKS.md](HOW_IT_WORKS.md), [BASIC_COMMANDS.md](../BASIC_COMMANDS.md) — the original `compas_singular` workflow
