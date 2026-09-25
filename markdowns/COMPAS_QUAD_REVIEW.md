# `compas_quad` vs `compas_singular` — full file-by-file review

**Reviewed:** `compas_singular/compas_quad-main/compas_quad-main/` (Robin Oval, `compas_quad` 0.1.0)
against `compas_singular/src/compas_singular/` (this repo).
**Date:** 2026-08-20 · **Environment used for all runtime checks:** conda env `singular312`
(Python 3.12, **COMPAS 2.15.1**, compas_viewer 2.0.2).

**Implementation status added 2026-08-26** — see §0 below. Batches 1–2 of §7 are done; the review
text itself is left as originally written, with `> **Update …**` notes where implementing an item
proved a claim wrong.

---

## 0. Implementation status (2026-08-26)

Batches 1 and 2 of §7 are implemented and merged into `src/compas_singular/`. Everything landed
under one rule: **behavioural changes ship behind a keyword argument that defaults to today's
behaviour.** Where the old behaviour was an *exception*, the fix landed unconditionally, since no
working call site can observe the difference.

| # | Item | Status | How it landed |
|---|---|---|---|
| 1 | New `polyedge_graph` | ✅ done | `polyedge_graph(legacy=True)`; the two coloring call sites are pinned to `legacy=True` |
| 2 | `collect_polyedges` O(1) removal | ✅ done, **improved** | list-as-stack + set for membership, so polyedge/strip keys keep their order (upstream's `set()` randomises them). Also applied to `collect_strips` and the `PseudoQuadMesh` override |
| 3 | `collect_polyedge` orientation + `both_sides` | ✅ done | `collect_polyedge(u0, v0, both_sides=True, oriented=False, strict=False)` |
| 4 | `face_opposite_edge` null guard | ✅ done | unconditional, plus the `PseudoQuadMesh` override and an unpack guard in both `collect_strip`s |
| 5 | Delete the `self.times` line | ✅ done | line removed together with the now-unused `import time` |
| 6 | Polyedge accessors | ✅ done, 3 deviations | see §4.4(f) update |
| 9 | `conftest.py` for rhino/blender/ghpython | ⏭️ **deliberately skipped** | `pytest.ini` sets `testpaths = tests` with no `--doctest-modules`, so `src/compas_singular/rhino/` is never collected. Nothing to guard |
| 10 | `requirements.txt`; CI matrix + action bump | ✅ partly stale | `requirements.txt` was already correct. CI matrix → 3.9–3.12, `compas-actions.build@v1.1.1`, `checkout@v4`, `setup-python@v5`. Also deleted the dead `[tool:pytest]` block from `setup.cfg` |
| 12 | `vertex_opposite_vertex` rewrite | ✅ done | `strict=False`, threaded through `collect_polyedge` → `collect_polyedges` → `singularity_polyedge_decomposition` → `CoarseQuadMesh.from_quad_mesh` |
| 14 | `collect_strip(both_sides=False)` | ✅ done | bonus — it is in no §7 batch, but it is one line and pairs with #3 |
| §5.5 | `__all__` in the grammar modules | ✅ done, **re-diagnosed** | the review had the mechanism backwards; see the §5 update |

**Deferred, not lost:** items 7, 8 (batch 4 — `skip_singularities` threading, `quad_mesh_polyedge_subcolor`,
the four `string_generation_*` functions), 11 (batch 3 — topological welding in `densification`), and
13, 15, 16 (the ATP lizard grammar and its examples/fixtures). Everything in the §6 "do not take"
table remains rejected.

### Verification of the batch

* **Test suite: 8 → 59 passing.** New `tests/test_mesh_quad.py` (null guards, accessors, and each
  opt-in flag against the review's own annulus and 4-valent-boundary fixtures); `tests/test_import.py`
  extended with namespace assertions.
* **Default-invariance: 0 differing fields.** A harness snapshotted `attributes['polyedges']`,
  `attributes['strips']`, the `polyedge_graph` edge list, both colourings, the singularity list and
  `singularity_polyedge_decomposition` on three meshes (dense British Museum, annulus, kinked grid)
  before and after. Nothing moved at the defaults.
* **Performance:** `collect_polyedges` measured at 480 / 1540 / 3200 vertices, best-of-10 per side,
  results byte-identical each time. A quiet run gave −27 % / −26 % / −38 %; a loaded run gave
  −32 % / −35 % / −82 %. The direction is stable, the magnitude is load-sensitive — the review's
  −35 % is the right order of magnitude to quote.

### Two results that go beyond what the review could establish

1. **The merged `strict` rule is cleaner than either source.** On the annulus it yields
   `[[2,3,0,1,2], [4,0], [5,1], [6,2], [6,7,4,5,6], [7,3]]` — the correct inner and outer rings plus
   four radials, **without** the `[7,3,2,1,0,3,2,1,0]` artefact §4.4(b) notes that *both* versions
   produce. Tightening the `n == 3` test to `is_edge_on_boundary(v, nbr)` is what removes it.
   On the dense British Museum mesh `strict=True` is a **no-op** (66 polyedges either way, coarse
   V/F 20/12, 9 strips), so it only bites on meshes with holes.
2. **The legacy `polyedge_graph` does break a colouring.** §4.4(e) says "I cannot claim a colouring
   failure". On the annulus it can: `is_adjacency_two_colorable` returns
   `{0: 0, 1: -1, 2: 1, 3: 1}` — colour `-1` for polyedge 1, an invalid colouring caused by the
   self-loops. With `legacy=False` the same mesh colours cleanly.

### Known-broken, untouched

`examples/04_lizard.py`, `06_add_strips.py` and `07_lizard_grammar.py` import
`compas_singular._compat`, which no longer exists in this checkout (removed in `04e3583d`). They fail
on import, independently of this work.

---

## 1. Headline verdict

**`compas_quad` is not a newer, COMPAS-2 version of `compas_singular`. It is a *sibling fork of the same
ancestor*, and on the axis that matters most here — COMPAS 2 compatibility — `compas_singular` is
strictly further ahead.**

| | `compas_singular` (this repo) | `compas_quad-main` |
|---|---|---|
| Lines of library code | ~16 100 | ~5 400 |
| COMPAS 2 imports | complete (`_compat.py` shim layer) | partial |
| COMPAS 2 *call signatures* | complete (Mesh shim accepts both) | **not migrated** — 12 call sites still use the removed 2-arg edge API |
| `from_json` / `to_json` | works | **broken for every datastructure** |
| Scope | quad grammar + skeleton decomposition + frame field + Rhino plugin + block/dual/guide-line tooling | quad grammar only |

Importing `compas_quad.coloring` fails outright on COMPAS 2.15.1:

```
FAIL compas_quad.coloring ImportError: cannot import name 'adjacency_from_edges' from 'compas.topology'
```

So **a wholesale replacement is off the table.** What `compas_quad` *does* contain is roughly a dozen
genuine, individually valuable improvements to code that already exists here, plus one entirely new
feature (the ATP lizard grammar). Those should be cherry-picked. Section 6 has the ranked list.

---

## 2. Provenance and relationship

Both packages descend from Robin Oval's PhD work (Ecole des Ponts / ETH Zurich,
[thesis](https://pastel.archives-ouvertes.fr/tel-02917467)). `compas_quad` carries
`author_email = robin.oval@princeton.edu`, i.e. it is the later *upstream* continuation, but it was
trimmed down to the quad-grammar core. Module mapping:

| `compas_quad` | `compas_singular` equivalent |
|---|---|
| `datastructures/mesh.py` | `datastructures/mesh/mesh.py` + `_compat.mesh_substitute_vertex_in_faces` |
| `datastructures/mesh_quad.py` | `datastructures/mesh_quad/mesh_quad.py` |
| `datastructures/mesh_quad_coarse.py` | `datastructures/mesh_quad_coarse/mesh_quad_coarse.py` |
| `datastructures/mesh_quad_pseudo.py` | `datastructures/mesh_quad_pseudo/mesh_quad_pseudo.py` |
| `datastructures/mesh_quad_pseudo_coarse.py` | `datastructures/mesh_quad_pseudo_coarse/mesh_quad_pseudo_coarse.py` |
| `grammar/addition.py` | `datastructures/mesh_quad/grammar/add_strip.py` |
| `grammar/deletion.py` | `datastructures/mesh_quad/grammar/delete_strip.py` (+ parts of `grammar_pattern.py`) |
| `grammar/lizard.py` | `datastructures/lizard/lizard.py` |
| `grammar/addition2.py` | **nothing — entirely new** |
| `coloring/coloring.py` | `topology/coloring.py` |
| `coloring/coloring_mesh.py` | `datastructures/mesh/coloring.py` |
| `coloring/coloring_quadmesh.py` | `datastructures/mesh_quad/coloring.py` + `mesh_quad_coarse/coloring.py` |
| `coloring/twocolorer.py` | `algorithms/twocoloring.py` |
| `utilities.py` | `utilities/lists.py` |

**In `compas_singular` and absent from `compas_quad` entirely:** `algorithms/` (skeleton
decomposition, isomorphism, layout, mapping, propagation, triangulation), `geometry/`, `editing/coarse_layout.py`,
`blocks.py` and `guide_lines.py` (both since removed), `rhino/` (the whole plugin: `coarse_curves.py`, `dual_mesh.py`,
`mesh_ui.py`, `grasshopper.py`, constraints, artists), `datastructures/skeleton/`, `datastructures/network/`,
`utilities/pareto.py`, `_compat.py`, `_viewer_patches.py`, and the `framefield/` example package.
341 top-level names exist only on this side.

---

## 3. Environment facts established by testing

Everything below was measured, not inferred.

```
compas 2.15.1
compas.topology.adjacency_from_edges          -> ABSENT (renamed vertex_adjacency_from_edges)
compas.utilities.pairwise                     -> ABSENT (moved to compas.itertools)
compas.datastructures.Network                 -> present (alias of Graph); .nodes() yes, .vertices() NO
compas.datastructures.Mesh.join / .weld       -> present
compas.itertools.linspace(0, 1, n)            -> yields n values; raises for n < 2
Mesh objects are always truthy (no __len__/__bool__)
```

Import smoke test of `compas_quad` in that environment:

```
OK    compas_quad
OK    compas_quad.utilities
OK    compas_quad.datastructures
OK    compas_quad.grammar
FAIL  compas_quad.coloring   ImportError: adjacency_from_edges
```

`datastructures` imports only because `setuptools` still shims `distutils` on 3.12 (used by
`mesh_quad_pseudo.py`); that shim is scheduled to go away.

---

## 4. File-by-file analysis

### 4.1 `src/compas_quad/__init__.py`, `__main__.py`, `setup.py`

Boilerplate, identical in shape to this repo's. `setup.py` declares `packages=["compas_quad"]`
(no subpackages listed — relies on `include_package_data`), Python ≥3.6.
`requirements.txt` = `compas` (this repo's is empty — arguably worth copying that one line).

**Verdict: nothing to take.**

---

### 4.2 `utilities.py` ⟷ `utilities/lists.py`

Same three functions. Differences:

* Parameter renamed `l` → `in_list` (kills flake8 E741). Cosmetic but nice.
* **`__all__` has a missing comma**, so it actually evaluates to
  `['are_items_in_listlist_split', 'sublist_from_to_items_in_closed_list']` — verified at runtime.
  A bug in `compas_quad`, not here.
* Drops `common_items` and `remove_isomorphism_in_integer_list`, which this repo's
  `algorithms/isomorphism.py` needs.

**Verdict: do not take.** Optionally apply the `l` → `in_list` rename locally.

---

### 4.3 `datastructures/mesh.py` ⟷ `datastructures/mesh/mesh.py`

> **Update 2026-08-26:** every reference to `_compat.py` in this review is stale. The module was
> deleted in `04e3583d` ("Direct integration of compas v2 functions"); the COMPAS 2 shims now live on
> the `Mesh` subclass at [mesh.py:28-48](../src/compas_singular/datastructures/mesh/mesh.py#L28-L48),
> and `mesh_substitute_vertex_in_faces` is imported straight from `compas.datastructures`.

`compas_quad`'s `Mesh` is minimal: a `move(vector)` method plus the free function
`mesh_substitute_vertex_in_faces` (identical to the one in [_compat.py:152](../src/compas_singular/_compat.py#L152)).

This repo's `Mesh` is much richer: `boundaries()`, `is_boundary_vertex_kink()`, `boundary_kinks()`,
`vertex_centroid()`, `to_vertices_and_faces(keep_keys=True)` **and — critically — the COMPAS 2
compatibility shims** at [mesh.py:26-48](../src/compas_singular/datastructures/mesh/mesh.py#L26-L48)
that let `edge_midpoint(u, v)`, `edge_length(u, v)`, `edge_faces(u, v)`, `is_edge_on_boundary(u, v)`
keep their old two-argument form.

That shim is exactly what `compas_quad` lacks. Proof:

```python
>>> QuadMesh(compas_quad).vertex_opposite_vertex(0, 4)
TypeError: Mesh.is_edge_on_boundary() takes 2 positional arguments but 3 were given
```

12 call sites in `compas_quad` still use the removed signature (`mesh_quad.py` ×7,
`mesh_quad_pseudo.py` ×2, `coloring_mesh.py` ×2, `mesh_quad_coarse.py` ×1).

**Verdict: keep ours.** Only `move()` is worth adding as a convenience (this repo has the free
function `mesh_move_by` in `datastructures/mesh/operations.py`, which is equivalent).

---

### 4.4 `datastructures/mesh_quad.py` ⟷ `datastructures/mesh_quad/mesh_quad.py`

This is the most productive file in the whole comparison. Eight substantive changes:

#### (a) `face_opposite_edge` — null guard ✅ TAKE

```python
fkey = self.halfedge[u][v]
if fkey is None:          # ← new
    return None
```
Ours raises `TypeError` on `face_vertex_descendant(None, v)` when handed a boundary halfedge.
One-line fix, no downside.

#### (b) `vertex_opposite_vertex` — rewritten ⚠️ TAKE HALF

Ours:
```python
if self.is_vertex_singular(v): return None
elif self.is_vertex_on_boundary(v):
    if not self.is_vertex_on_boundary(u): return None
    else: return [nbr for nbr in self.vertex_neighbors(v) if nbr != u and self.is_vertex_on_boundary(nbr)][0]
else: ...nbrs[nbrs.index(u) - 2]
```
Theirs:
```python
nbrs = self.vertex_neighbors(v, ordered=True); n = len(nbrs)
if n == 4: return nbrs[nbrs.index(u) - 2]
if self.is_edge_on_boundary(u, v) and n == 3:
    for nbr in nbrs:
        if nbr != u and self.is_vertex_on_boundary(nbr): return nbr
return None
```

Two independent changes hide in there, and they pull in opposite directions.

**The `n == 3` clause is a real bug fix — take it.** Ours only checks that *vertex* `u` is on a
boundary, not that the *edge* `(u, v)` is. On an annulus (4 quads round a square hole), walking from
outer vertex 0 into inner vertex 4 across an interior edge, ours turns onto the hole boundary:

```
singular  opp(0, 4) = 5      collect_polyedge(0,4) = [0, 4, 5, 1, 0]
quad-rule opp(0, 4) = None   collect_polyedge(0,4) = [4, 0]
```

and the resulting polyedge sets are correspondingly garbage vs. clean:

```
singular : [[6,7,3,2,1,0,4,5,1], [5,6,2,1,0,4,5], [4,7,3,2,1,0,4], [0,3,2,1,0]]
quad-rule: [[6,7,4,5,6], [7,3], [6,2], [5,1], [4,0], [7,3,2,1,0,3,2,1,0]]
```

Ours mixes inner ring (4-7) and outer ring (0-3) vertices in a single polyedge. Theirs recovers the
correct inner loop plus four radial polyedges.

**The `n == 4` clause is a regression — do not take it as written.** It drops the
`is_vertex_singular` test, so a *boundary* vertex of degree 4 (which is singular: regular boundary
valency is 3) is now walked straight through. Measured on a 2×2 grid with one extra quad hanging off
the mid-bottom vertex (vertex 1, ordered nbrs `[0, 4, 2, 9]`, boundary edges `(1,0)` and `(1,9)`):

```
singular : opp(0,1) = None   opp(4,1) = None
quad-rule: opp(0,1) = 2      opp(4,1) = 9     ← boundary polyedge diverted into the interior,
                                                and an interior polyedge turned onto the boundary
```

**Recommended merge** — take their structure, guard the `n == 4` branch, and fix the residual bug
they left in the `n == 3` branch (it tests `is_vertex_on_boundary(nbr)`, which cannot distinguish
boundary *loops*; it should test `is_edge_on_boundary(v, nbr)`):

```python
def vertex_opposite_vertex(self, u, v):
    nbrs = self.vertex_neighbors(v, ordered=True)
    n = len(nbrs)
    if n == 4 and not self.is_vertex_on_boundary(v):
        return nbrs[nbrs.index(u) - 2]
    if n == 3 and self.is_edge_on_boundary(u, v):
        for nbr in nbrs:
            if nbr != u and self.is_edge_on_boundary(v, nbr):
                return nbr
    return None
```
That last change fixes the remaining `[7,3,2,1,0,3,2,1,0]` artefact above, which *both* versions produce.

> ⚠️ This function is the engine under `collect_polyedge` → `collect_polyedges` →
> `singularity_polyedge_decomposition` → `CoarseQuadMesh.from_quad_mesh`. Changing it changes coarse
> layouts. Re-run the 84-row frame-field baseline before adopting.

#### (c) `collect_polyedge` — orientation preservation + `both_sides` ✅ TAKE

Adds a `flipped` flag and reverses back at the end, so the returned polyedge runs in the direction of
the seed halfedge. Ours silently returns the reverse whenever the walk hits an extremity first:

```
seed (2, 1):  singular -> [0, 1, 2, 6, 5, 1, 0]     (contains "1,2" — seed direction lost)
              quad     -> [0, 1, 5, 6, 2, 1, 0]     (contains "2,1" — seed direction kept)
```

Also adds `both_sides=False` for a one-sided walk (`[2, 1, 0]`), which is what a directional
tracer wants. Pure addition; nothing here depends on the old flipping behaviour.

#### (d) `collect_polyedges` — `set` instead of `list` ✅ TAKE

`edges = set(self.edges())` instead of `list(...)`, so removal is O(1) rather than O(n).
Measured on a 1369-vertex / 1296-face dense mesh: **0.241 s → 0.156 s (−35 %)**.
Caveat: `set.pop()` is arbitrary-order, so polyedge *keys* will be assigned in a different order.
Nothing here indexes polyedges by key across runs, but worth a grep before merging.

> **Update 2026-08-26 — taken, but not as written.** The caveat is avoidable: keeping the
> deterministic `list` as the seed order and adding a `set` purely for membership gets the same
> speedup with **byte-identical output**, keys included. That removes the need for the grep and made
> the whole batch verifiable by straight comparison. Applied to `collect_polyedges`,
> `collect_strips`, and the `PseudoQuadMesh.collect_strips` override.

#### (e) `polyedge_graph` — rewritten ✅ TAKE

Ours is O(n_polyedges²) and structurally wrong: for every vertex it appends `(key, key_2)` for the
*first* polyedge containing that vertex — which is usually the polyedge itself. Measured on the same
1369-vertex mesh:

```
singular : 0.031 s, 2730 edges, of which 1365 are self-loops (u, u)
quad     : 0.005 s, 1365 edges, 0 self-loops, 4 vertices skipped (the true singularities)
```

**6× faster and half the edges it emits are meaningless.** In my test the downstream
`quad_mesh_polyedge_2_coloring` still produced the same colouring on both (the self-loops happen to
be absorbed by `is_adjacency_two_colorable`), so I cannot claim a colouring failure — but the graph
is wrong and the new one is right.

> **Update 2026-08-26 — the colouring failure is claimable.** On the annulus fixture from §4.4(b),
> the legacy graph makes `is_adjacency_two_colorable` return `{0: 0, 1: -1, 2: 1, 3: 1}` — colour
> `-1`, i.e. not a colouring at all. The self-loops are not always absorbed. Ported as
> `polyedge_graph(legacy=True)`; on the dense British Museum mesh `legacy=False` gives 472 edges and
> 0 self-loops against 950 and 473, with the same colouring result.

Take it, but **replace the `print('exception in polyedge graph', pkeys)` with silence or a debug
log** — it fires once per singularity on every call.

#### (f) New accessors ✅ TAKE

`number_of_polyedges()`, `polyedge_vertices(pkey)`, `polyedge_edges(pkey)`,
`polyedge_midpoint(pkey)`, `polyline(pkey)`, `polyedge_length(polyedge)`.
All trivial, all useful. Note the naming inconsistency in upstream: every other `polyedge_*` takes a
`pkey`, but `polyedge_length` takes a vertex list. Fix that on the way in.

> **Update 2026-08-26 — ported, with three deviations.** "All trivial" was too generous:
>
> 1. `polyedge_midpoint` uses `Polyline(...).point(0.5)` upstream. **`Polyline.point` does not exist
>    in COMPAS 2.15.1** — verified, it raises `AttributeError`; the method is `point_at`. Ported as
>    `point_at(0.5)`.
> 2. `polyedge_edges` returns `pairwise(...)` upstream, and `compas.itertools.pairwise` returns a
>    `zip` — a one-shot iterator, which is a trap in an accessor. Ported returning a `list`.
> 3. `polyedge_length` takes a `pkey` here, as recommended above.

#### (g) `collect_strip(both_sides=False)` ⚠️ TAKE WITH FIX

```python
if self.halfedge[u0][v0] is None:
    if not both_sides:
            return (u0, v0)      # ← returns a TUPLE, not a list of edges
    u0, v0 = v0, u0
```
Every other return path yields `[(u, v), ...]`. Should be `return [(u0, v0)]`. (The stray extra
indentation is legal but sloppy.)

#### (h) `singularity_polyedge_decomposition` — `boundaries()` → `vertices_on_boundaries()`

Behaviourally equivalent here; ours has its own `boundaries()` on the `Mesh` subclass.
No reason to change. **Skip.**

---

### 4.5 `datastructures/mesh_quad_coarse.py` ⟷ `mesh_quad_coarse/mesh_quad_coarse.py`

#### (a) The getter/setter → callable-accessor rename ❌ DO NOT TAKE

```
get_quad_mesh()/set_quad_mesh()             ->  dense_mesh(mesh=None)
get_polygonal_mesh()/set_polygonal_mesh()   ->  pattern(pattern=None)
get_strip_density()/set_strip_density()     ->  strip_density(skey, d=None)
get_strip_densities()/set_strips_density()  ->  strips_density(d=None)
attributes['quad_mesh'] / ['polygonal_mesh'] ->  attributes['dense_mesh'] / ['pattern']
```

Cost, counted across `compas_singular` + `compas_topology` (excluding the vendored copy):

| symbol | call sites |
|---|---:|
| `set_strips_density` | 106 |
| `get_quad_mesh` | 79 |
| `set_strip_density` | 30 |
| `get_strip_density` | 16 |
| `set_quad_mesh` | 7 |
| `get_strip_densities` / `set_polygonal_mesh` / `get_polygonal_mesh` | 5 |
| **total** | **~243** |

That includes the Rhino `CMD_*` and `ff_*` command scripts in the sibling `compas_topology` repo,
which are not covered by this repo's tests.

And the new form carries a **live bug**: the setters test `if d:`, so density `0` is silently ignored.
Verified:

```
strip_density(skey, 0)  -> still 4        (the write was dropped)
strips_density(0)       -> attributes['strips_density'] left EMPTY -> KeyError on next read
```

`set_strips_density(d, skeys=None)` also loses its `skeys` argument in the rename.

**Verdict: reject the rename.** If the shorter names are wanted for ergonomics, add them as thin
aliases that delegate to the existing methods and use `if d is not None:`.

#### (b) The new `densification()` — topological welding ✅ TAKE (for `CoarseQuadMesh`)

This is the single best algorithmic idea in `compas_quad`. Ours builds one Coons patch per coarse
face and then calls `meshes_join_and_weld`, which welds by **geometric key** — i.e. any two vertices
that land on the same coordinates get fused, whether or not they are topologically the same vertex.
Theirs builds an explicit `parent2child` index map per coarse halfedge and welds by **topology**.

Measured on a coarse mesh containing a deliberate slit (two topologically distinct vertices at
identical xyz):

```
singular : V/F 17/8   manifold = False     ← the slit was welded shut
quad     : V/F 18/8   manifold = True      ← slit preserved
```

On a clean 2×2 coarse mesh at density 4 both give V/F 81/64, manifold, in ~3 ms. So: same output when
geometry is unambiguous, strictly better when it is not.

**But it is not a drop-in.** Three things must be merged in:

1. **`edges_to_curves` is missing from `CoarseQuadMesh.densification` upstream.** That parameter is
   ours (the curved-boundary work; see [mesh_quad_coarse.py:252](../src/compas_singular/datastructures/mesh_quad_coarse/mesh_quad_coarse.py#L252))
   and dropping it would faceted-ise every curved boundary in the Rhino pipeline.
2. **The Coons argument order changes** from `discrete_coons_patch(ab, bc, cd[::-1], da[::-1])` to
   `discrete_coons_patch(da[::-1], cd[::-1], bc, ab)`, so that `vertices` comes out row-major and the
   `parent2child` index arithmetic (`n`, `n-1+k*n`, `len(v)-n …`) is valid. The curve sampling has to
   be re-mapped to match.
3. `parent2child` is computed and then **thrown away**. It is precisely the coarse-edge → dense-polyedge
   map that `attributes['edge_coarse_to_dense']`, `edges_to_curves`, and the frame-field densifier all
   want. Store it while porting.

Upstream keeps the old implementation as `densification_old()`. Do the same here during migration.

#### (c) `from_quad_mesh` — `adjacency_from_edges` → `vertex_adjacency_from_edges`

Ours already imports the COMPAS 2 name via `_compat`. Otherwise identical. **Skip.**

---

### 4.6 `datastructures/mesh_quad_pseudo.py` ⟷ `mesh_quad_pseudo/mesh_quad_pseudo.py`

**Ours is strictly better. Take nothing.**

`compas_quad` carries a hand-rolled COMPAS 1.x `data` property setter (~85 lines) that branches on
`LooseVersion(compas.__version__)` and reconstructs vertices/faces/edgedata by hand. On COMPAS 2 that
whole path is dead code — and it imports `from distutils.version import LooseVersion`, which is
removed in Python 3.12 (it only resolves here because `setuptools` still shims `distutils`).

Ours replaces all of it with 6 lines of `__from_data__` override
([mesh_quad_pseudo.py:23-30](../src/compas_singular/datastructures/mesh_quad_pseudo/mesh_quad_pseudo.py#L23))
that restores the integer face keys of `face_pole` after JSON round-trip.

Verified round-trip:

```
singular  face_pole before {0:4, 1:4, 2:4, 3:4}  ->  after {0:4, 1:4, 2:4, 3:4}  (int keys)
quad      face_pole before {0:4, 1:4, 2:4, 3:4}  ->  FAILED
          TypeError: PseudoQuadMesh.__init__() got an unexpected keyword argument 'default_vertex_attributes'
```

That last failure is the `__init__(self)` signature — upstream dropped `*args, **kwargs` from every
constructor (`QuadMesh`, `CoarseQuadMesh`, `PseudoQuadMesh`, `CoarsePseudoQuadMesh`), which breaks
COMPAS 2 deserialisation for **all** of them, not just the pseudo mesh:

```
singular CoarseQuadMesh json round-trip OK
quad     CoarseQuadMesh json round-trip FAILED (same TypeError)
```

---

### 4.7 `datastructures/mesh_quad_pseudo_coarse.py` ⟷ `mesh_quad_pseudo_coarse/mesh_quad_pseudo_coarse.py`

**The new pole densification is broken. Do not take it.**

Same square-with-a-centre-pole input, density 4, both packages:

```
singular : V/F 65/64  manifold=True   face_pole entries=16  tri faces=16   6.5 ms
quad     : V/F 74/64  manifold=False  face_pole entries=0   tri faces=4    3.5 ms
```

Three defects:

1. **Non-manifold output** — the pole vertices are never merged (74 vertices instead of 65).
2. **`face_pole` is lost entirely.** The old implementation built `face_pole_map` from geometric keys
   and wrote `dense_mesh.attributes['face_pole'] = face_pole`; the new one deleted that block and never
   replaced it. Everything downstream that asks "is this face a pseudo-quad?" — `vertex_index`,
   `strip_faces`, `is_strip_closed`, `has_strip_poles`, the whole pole/block workflow — goes blind.
3. **It re-introduces the `linspace` off-by-one this repo already fixed.** Upstream still has
   `[curve.point_at(t) for t in linspace(0, 1, d)]`; ours is `linspace(0, 1, d + 1)`. Confirmed
   semantics: `linspace(0,1,3) -> [0.0, 0.5, 1.0]` (n points, not n subdivisions) and
   `linspace(0,1,1)` **raises** `ValueError`, i.e. upstream crashes at density 1.

Upstream's `densification_old` is equivalent to ours modulo that `linspace` bug, so it adds nothing.

*If* the topological weld is wanted for pseudo meshes too, it has to be written here: keep upstream's
`weld_map` skeleton, add pole handling (their `u == v` branch is incomplete), and re-attach `face_pole`.
That is a real piece of work, not a port.

---

### 4.8 `grammar/addition.py` ⟷ `mesh_quad/grammar/add_strip.py`

The core `add_strip` is **byte-for-byte identical** apart from line wrapping. Differences:

* Adds `add_strips(mesh, polyedges, callback=None, callback_args=None)` — a plain loop. Ours
  (`grammar_pattern.add_strips`) is **more capable**: it updates the remaining polyedges after each
  insertion via `strip_polyedge_update`, which upstream leaves as a `# update polyedges` TODO.
  **Keep ours.** The `callback` hook is a reasonable addition to ours.
* Keeps `add_strip_old` (the pre-rewrite version). Ours lives in `grammar_pattern.add_strip` — same
  code plus a `boundary_kinks(pi/12)` call. Already covered.
* Populates `__all__` (ours is an empty list in `grammar/add_strip.py`, so `from .add_strip import *`
  re-exports everything by accident). **Worth fixing here.**

> Note: this repo currently exports `add_strip` from `grammar/add_strip.py` (returns
> `(n, old_to_new)`) but `add_strips` from `grammar_pattern.py` (which internally calls the *old*
> `add_strip` returning `(skey, left, right)`). Verified via `inspect.getmodule`. That inconsistency
> is ours, not upstream's, and is worth resolving independently of this review.

---

### 4.9 `grammar/addition2.py` — **entirely new** ⚠️ TAKE WITH A FIX

438 lines with no counterpart here. Implements strip addition as a **turtle walk** rather than as an
operation on a pre-computed polyedge:

* `lizard_atp(mesh, lizard, movements)` — interprets an `a`/`t`/`p` string, where `a` toggles
  "recording" mode and the recorded `t`/`p` substring is played back as a strip insertion.
* `add_strip_lizard` (`'0'`/`'1'` alphabet) and `add_strip_lizard_2` (`'p'`/`'t'` alphabet) — otherwise
  identical twins; the duplication is upstream sloppiness.
* Inserts triangles as it advances and retro-converts the previous triangle to a quad, so the mesh is
  quad-valid at each completed step, and closes the strip on return to the boundary.

This is a genuinely different generative mechanism from `add_strip(mesh, polyedge)` and is the basis
of upstream's newest examples (`99_lizard3.py`, `99_lizard4.py`). It is the one thing in `compas_quad`
that is *new capability* rather than *a better version of what we have*.

**Blocking issue:** its `pivot`/`turn` helpers use `nbrs[i + 1 - len(nbrs)]`, the COMPAS 0.x/1.x
neighbour winding. Under COMPAS 2.15.1 that winding is reversed — see §4.10 for the measurement.
Every one of these helpers needs `nbrs[i - 1]` before the module will behave.

Also: the `if __name__ == '__main__'` block imports `compas.numerical.fd_numpy` and `compas_view2`,
both gone in COMPAS 2; and `t1, t2` are computed and never used.

---

### 4.10 `grammar/lizard.py` ⟷ `datastructures/lizard/lizard.py`

The `Lizard` class is identical **except** for exactly the thing this repo already fixed:

```python
# compas_quad
new_head = nbrs[i + 1 - len(nbrs)]
# compas_singular  (with the comment explaining why)
new_head = nbrs[i - 1]   # COMPAS 2.x orders vertex neighbors with the opposite winding
```

Measured, 2×2 coarse grid, `initiate(tail=0, head=1)`:

```
singular (i-1)   'atta'    -> V/F 12/6  quad=True  manifold=True
singular (i-1)   'attta'   -> V/F 13/7  quad=True  manifold=True
singular (i-1)   'atttta'  -> V/F 14/8  quad=True  manifold=True
quad     (i+1)   'atta'    -> FAILED IndexError: list index out of range
quad     (i+1)   'attta'   -> V/F 13/7  quad=True  manifold=True
quad     (i+1)   'atttta'  -> V/F 13/8  quad=True  manifold=False
```

**Do not take the class.** ❌

**Do take the four string generators** ✅ — they are pure functions with no COMPAS coupling and no
counterpart here:

* `string_generation_brute(characters, length)` — exhaustive `itertools.product`.
* `string_generation_random(characters, number, length, ratios=None)` — weighted random.
* `string_generation_structured(characters, number, length)` — biased so that `a…a` brackets are
  balanced and enclose a polyedge of length ≥ 2, i.e. it emits far fewer invalid strings.
* `string_generation_evolution(characters, number, length)` — single-character mutation chain, for
  hill-climbing over the grammar.

---

### 4.11 `grammar/deletion.py` ⟷ `mesh_quad/grammar/delete_strip.py`

* `network_disconnected_nodes(network)` → `[nwk.vertices() for nwk in network.exploded()]`.
  **This is broken twice over on COMPAS 2:** `Graph`/`Network` exposes `.nodes()`, not `.vertices()`
  (verified: `Network.vertices -> False`), *and* `.nodes()` returns a generator, which the loop below
  consumes twice (once for `centroid_points`, once for the substitution loop) — the second pass would
  see an empty iterator. Ours calls `connected_components(network.adjacency)` via `_compat` and
  returns lists. **Keep ours.**
* Adds `collateral_strip_deletions` and `total_boundary_deletions` — but **we already have both**, in
  `grammar_pattern.py` ([:565](../src/compas_singular/datastructures/mesh_quad/grammar_pattern.py#L565)
  and [:573](../src/compas_singular/datastructures/mesh_quad/grammar_pattern.py#L573)), and
  `algorithms/twocoloring.py` already imports them. Identical code.
* Populates `__all__` (ours is empty — same issue as `add_strip.py`).

**Verdict: take nothing but the `__all__` fix.**

---

### 4.12 `coloring/coloring.py` ⟷ `topology/coloring.py`

Identical modulo whitespace. **Nothing to take.**

---

### 4.13 `coloring/coloring_mesh.py` ⟷ `datastructures/mesh/coloring.py`

Identical except upstream imports `compas.topology.adjacency_from_edges` (removed in COMPAS 2 — ours
uses `vertex_adjacency_from_edges`) and calls `is_edge_on_boundary(u, v)` without a shim.
**Ours is ahead. Nothing to take.**

---

### 4.14 `coloring/coloring_quadmesh.py` ⟷ `mesh_quad/coloring.py` + `mesh_quad_coarse/coloring.py`

Upstream merges our two modules into one and adds:

* **`skip_singularities` parameter** on `quad_mesh_polyedge_2_coloring`, threaded through to
  `polyedge_graph`. ✅ Take — it pairs with §4.4(e) and lets you ask for the colouring *with*
  singularity connections counted, which we currently cannot express.
* **`quad_mesh_polyedge_subcolor(quad_mesh, color=0)`** — new. Re-colours the subset of polyedges
  sharing one colour, by their adjacency through faces. ✅ Take; self-contained.
* Broken imports: `compas.topology.adjacency_from_edges` **and** `compas.utilities.pairwise`, both
  gone in COMPAS 2 (verified). Fix on the way in.

Whether to merge our two coloring modules into one is a taste call; the split (`mesh_quad/` vs
`mesh_quad_coarse/`) mirrors our package layout and I would leave it.

---

### 4.15 `coloring/twocolorer.py` ⟷ `algorithms/twocoloring.py`

Renamed `TwoColourableProjection` → `TwoColorer`. Body otherwise identical apart from reflow, with one
real change:

**`projection_0` ends with `self.times = (t1 - t0, t2 - t0, t3 - t0)` — and `t0`…`t3` are never
assigned anywhere in the file.** Guaranteed `NameError` at the end of every `projection_0()` call.
Verified: [twocoloring.py:450](../src/compas_singular/algorithms/twocoloring.py#L450), with `import time`
at line 5 and no assignment in between. Upstream deleted the line (and the unused import).

✅ **Take that deletion.** It is a real bug fix in our code.
❌ Skip the class rename — `TwoColourableProjection` is referenced in our docs and `algorithms/__init__.py`.

> **Update 2026-08-26 — done.** Both the `self.times` assignment and the unused `import time` are
> deleted, so the line numbers quoted above no longer resolve. `self.times = None` in `__init__` is
> kept. Verified: `projection_0()` on the British Museum coarse mesh now returns 9 result entries
> instead of raising. Note for anyone with a script around this call — it used to raise `NameError`
> at the very end of every *successful* run, so a `try/except` wrapper would have swallowed the
> completion; that wrapper is now dead code and the script will fall through instead.
> The class rename was skipped, as recommended.

---

### 4.16 `examples/`

* `00_densification.py` — uses COMPAS 2's `compas_viewer.viewer.Viewer`. Roughly equivalent to our
  `examples/00_densification.py`.
* `01_redensification.py`, `02_densification_poles.py` — still on `compas_view2.app.App` (COMPAS 1).
* `99_lizard.py` … `99_lizard4.py` — the interesting ones. `99_lizard3/4.py` are the most modern code
  in the whole package: `compas_fd.solvers.fd_numpy` (correct for COMPAS 2), `compas_viewer.Viewer`,
  a clean form-finding `postprocessing()` that maps the boundary to a circle and relaxes with FDM.
  `99_lizard4.py` demonstrates the full `lizard_atp` → pole extraction → `CoarsePseudoQuadMesh` →
  densify → form-find loop.
  ⚠️ `99_lizard.py` and `99_lizard2.py` still use `compas.numerical.fd_numpy` (removed in COMPAS 2)
  and `mesh.key_index()` / `index_key()` (renamed to `vertex_index()` / `index_vertex()`).
  `99_lizard4.py` also has a hardcoded `to_json('/Users/roval/Desktop/mesh4.json')`.

✅ Worth adapting `99_lizard3.py`/`99_lizard4.py` as our `examples/` entry point for the ATP grammar,
*after* §4.9's winding fix.

**`examples/jsons/`** — three fixtures (British Museum coarse quad mesh, the same with poles, and the
dense mesh). We have none. ⚠️ They are **COMPAS 0.x/1.x raw-dict JSON** (no `dtype`/`data` wrapper):

```
CoarseQuadMesh.from_json(...) -> TypeError: The data in the file is not a CoarseQuadMesh
```

They would need a one-off conversion pass before they are useful as test data — but a real-world
coarse/dense pair with poles is exactly what `tests/test_pipeline.py` lacks, so it is worth doing.

---

### 4.17 Project infrastructure

| item | assessment |
|---|---|
| `tasks.py` | Upstream adds `lint`, `format`, `testdocs`, `linkcheck` invoke tasks and moves docs from `docsource/` → `docs/` + `dist/docs`. ✅ Modest improvement; adopting the docs move would churn our `docsource/` tree, so probably just take `lint`/`format`. |
| `.github/workflows/build.yml` | Python matrix `3.6–3.10` (ours: `3.6–3.8`); `compas-actions.build@v1.1.1` (ours: `v1.0.0`). ✅ Both worth taking — though on a 3.12 target the matrix should really be `3.9–3.12`. |
| `.github/workflows/pr-checks.yml` | New: enforces a `CHANGELOG.md` entry per PR. Optional. |
| `conftest.py` | Adds doctest namespace fixtures and skips `rhino`/`blender`/`ghpython` paths during collection. ✅ **Genuinely useful here** — we have `rhino/` modules that cannot import outside Rhino. We have no `conftest.py` at all. |
| `setup.cfg` | `testpaths` includes `src/` (doctest collection). Ours has a `[black]` section theirs lacks. Minor. |
| `requirements.txt` | Theirs: `compas`. Ours: **empty**. ✅ Trivially worth fixing. |
| `docs/` | A single `automodule:: compas_quad` stub. Ours (`docsource/`, 7 chapters) is far richer. Skip. |

> **Update 2026-08-26.** Three rows of this table are stale or were reversed on implementation:
>
> * **`requirements.txt` is not empty** in this checkout — it pins `compas>=2.0,<3` and
>   `networkx>=2.0`, and `pyproject.toml` reads it via `[tool.setuptools.dynamic]`. No change made.
> * **`conftest.py` was deliberately skipped.** `pytest.ini` sets `testpaths = tests` and does not
>   enable `--doctest-modules`, so `src/compas_singular/rhino/` is never collected in the first place.
>   Adding the upstream conftest would guard against nothing.
> * **`setup.cfg`** turned out to carry a `[tool:pytest]` block that `pytest.ini` silently overrides —
>   dead config containing `--strict`, which pytest 8 removed. The block was deleted so it cannot
>   spring back on whoever eventually drops `pytest.ini`.
>
> `.github/workflows/build.yml` was updated as recommended: matrix → 3.9–3.12,
> `compas-actions.build@v1.1.1`, `checkout@v4`, `setup-python@v5`.

---

## 5. Bugs found *in this repo* while doing the comparison

Independent of whether anything is ported:

1. **`algorithms/twocoloring.py:450`** — `self.times = (t1 - t0, t2 - t0, t3 - t0)` with all four names
   undefined. `projection_0()` raises `NameError` on every successful completion.
2. **`polyedge_graph`** emits a self-loop `(pkey, pkey)` for roughly every non-singular vertex —
   1365 of 2730 edges on a 1369-vertex test mesh.
3. **`vertex_opposite_vertex`** turns interior polyedges onto boundaries when `u` merely *is* a
   boundary vertex without the edge `(u, v)` being a boundary edge (annulus test, §4.4b).
4. **`face_opposite_edge`** raises instead of returning `None` for a boundary halfedge.
5. **`grammar/add_strip.py` and `grammar/delete_strip.py` both have `__all__ = []`**, so
   `from .add_strip import *` re-exports every imported symbol including `pairwise` and
   `breadth_first_paths`.
6. **`add_strip` vs `add_strips` resolve to different implementations** — `add_strip` from
   `grammar/add_strip.py`, `add_strips` from `grammar_pattern.py` (which calls its own older
   `add_strip`). Verified with `inspect.getmodule`.
7. **`requirements.txt` is empty.**

> **Update 2026-08-26 — status of these seven, and a re-diagnosis of #5.**
>
> Fixed: **1** (line and the unused `import time` deleted — `projection_0()` now completes instead of
> raising at the end of every successful run), **3** (behind `strict=True`), **4** (unconditional).
> **7** is stale: `requirements.txt` already pins `compas>=2.0,<3` and `networkx>=2.0`.
> **2** is addressed behind `polyedge_graph(legacy=False)`.
>
> **#5 has the mechanism backwards, and #6 is a symptom of it.** `__all__ = []` re-exports
> *nothing* — when `__all__` is defined, `import *` imports exactly those names. The actual leak is
> the `__all__ = [name for name in dir() if not name.startswith('_')]` idiom in the package
> `__init__`s: `from .add_strip import *` binds `add_strip` as a **submodule attribute**, and that
> module object then gets re-exported and shadows the *function* of the same name. Measured before
> the fix:
>
> ```
> type(compas_singular.datastructures.add_strip)    -> module   (not callable)
> type(compas_singular.datastructures.delete_strip) -> module   (not callable)
> import compas_singular.datastructures.mesh_quad.grammar
>         -> ImportError: cannot import name 'grammar' from '...mesh_quad.mesh_quad'
> ```
>
> That was live breakage, not cosmetics: [algorithms/mapping.py:7](../src/compas_singular/algorithms/mapping.py#L7)
> imports `delete_strip` from the package namespace and calls it at lines 54, 172, 179 and 188 —
> `TypeError: 'module' object is not callable` on every one.
>
> **Populating the two `__all__`s as written above would have caused a regression.** The two
> implementations share names but not signatures, and `mesh_quad/__init__.py` imports `.grammar`
> *last*:
>
> | name | `grammar/…` | `grammar_pattern.py` |
> |---|---|---|
> | `add_strip` | `(mesh, polyedge) -> (n, old_to_new)` | `(mesh, polyedge) -> (skey, left, right)` |
> | `add_strips` | polyedges **not** updated afterwards | updates the remaining polyedges |
> | `delete_strip` | `(mesh, skey, update_data=True)` | `(mesh, skey, preserve_boundaries=False)`, also fixes `face_pole` |
> | `delete_strips` | `(mesh, skeys, callback=None, …)` | `(mesh, skeys, preserve_boundaries=False)` |
>
> Exporting them would have rebound `datastructures.delete_strips` to the callback version and broken
> [twocoloring.py:423](../src/compas_singular/algorithms/twocoloring.py#L423)
> (`delete_strips(copy_mesh, combination, preserve_boundaries=True)`).
>
> **What was done instead:** leave `__all__ = []` in the two grammar modules untouched and fix the
> root cause — exclude module objects from the computed `__all__`, in all 19 package `__init__.py`
> files that use the idiom. Result: all four names resolve to grammar_pattern's richer functions,
> `add_strips`/`delete_strips` are **unchanged**, `mapping.py` starts working, and
> `import compas_singular.datastructures.mesh_quad.grammar` succeeds. Regression tests for all of it
> are in `tests/test_import.py`.
>
> #6 remains open as a *design* question — two different `add_strip` contracts still exist under one
> name. That inconsistency is ours, not upstream's, and is worth resolving independently.

---

## 6. Ranked recommendation

### Take (high value, low risk)

| # | Change | Where | Evidence | Status |
|---|---|---|---|---|
| 1 | New `polyedge_graph` | `mesh_quad/mesh_quad.py` | 6× faster, removes 1365 spurious self-loops | ✅ as `legacy=False` |
| 2 | `collect_polyedges` uses a `set` | `mesh_quad/mesh_quad.py` | 0.241 s → 0.156 s | ✅ list+set, keys preserved |
| 3 | `collect_polyedge` orientation fix + `both_sides` | `mesh_quad/mesh_quad.py` | seed direction preserved | ✅ as `oriented=True` |
| 4 | `face_opposite_edge` null guard | `mesh_quad/mesh_quad.py` | prevents `TypeError` | ✅ unconditional |
| 5 | Delete the `self.times` line | `algorithms/twocoloring.py` | removes a guaranteed `NameError` | ✅ done |
| 6 | Polyedge accessors (`polyline`, `polyedge_vertices`, `polyedge_edges`, `polyedge_midpoint`, `polyedge_length`, `number_of_polyedges`) | `mesh_quad/mesh_quad.py` | pure additions | ✅ 3 deviations, see §4.4(f) |
| 7 | `quad_mesh_polyedge_subcolor` + `skip_singularities` flag | `mesh_quad/coloring.py` | new capability | ⏳ deferred (batch 4) |
| 8 | The four `string_generation_*` functions | new `datastructures/lizard/generation.py` | new capability, zero coupling | ⏳ deferred (batch 4) |
| 9 | `conftest.py` skipping `rhino`/`blender`/`ghpython` | repo root | we have none | ⏭️ skipped — not needed, see §0 |
| 10 | `requirements.txt` = `compas`; CI matrix + action bump | repo root | trivial | ✅ CI done; requirements already correct |

### Take with modification (high value, needs work)

| # | Change | Required fix | Status |
|---|---|---|---|
| 11 | **Topological welding in `CoarseQuadMesh.densification`** | Merge in `edges_to_curves`; re-map the Coons argument order; **store** `parent2child`; keep the old one as `densification_old` | ⏳ deferred (batch 3) |
| 12 | `vertex_opposite_vertex` rewrite | Guard `n == 4` with `not is_vertex_on_boundary(v)`; change the `n == 3` boundary test from `is_vertex_on_boundary(nbr)` to `is_edge_on_boundary(v, nbr)` | ✅ as `strict=True`, both fixes applied |
| 13 | **`grammar/addition2.py` (ATP lizard)** | Change every `nbrs[i + 1 - len(nbrs)]` to `nbrs[i - 1]`; de-duplicate `add_strip_lizard` / `add_strip_lizard_2`; drop the `__main__` block | ⏳ deferred (batch 4) |
| 14 | `collect_strip(both_sides=False)` | `return [(u0, v0)]`, not `(u0, v0)` | ✅ done |
| 15 | `examples/99_lizard3.py` / `99_lizard4.py` | Adapt after #13; they already use `compas_fd` + `compas_viewer` correctly | ⏳ deferred (batch 4) |
| 16 | `examples/jsons/*` fixtures | Convert from COMPAS 1.x raw-dict JSON before use | ⏳ deferred — and partly stale: `examples/data/coarse_quad_mesh_british_museum*.json` already exist here and are used by `tests/test_pipeline.py` |

### Do not take

| # | Change | Why |
|---|---|---|
| 17 | The `get_*`/`set_*` → callable-accessor rename | ~243 call sites incl. the `compas_topology` Rhino commands; and `if d:` silently drops density 0 |
| 18 | `CoarsePseudoQuadMesh.densification` (new) | Non-manifold output, `face_pole` lost, `linspace` off-by-one re-introduced |
| 19 | `PseudoQuadMesh.data` property setter | COMPAS 1.x deserialiser + `distutils`; our `__from_data__` is correct |
| 20 | `__init__(self)` without `*args, **kwargs` | Breaks `from_json` for every datastructure |
| 21 | `Lizard.turn` / `.pivot` winding | Verified `IndexError` / non-manifold under COMPAS 2.15.1 |
| 22 | `network.exploded()` in `delete_strip` | `Graph` has no `.vertices()`; generator consumed twice |
| 23 | `adjacency_from_edges` / `compas.utilities.pairwise` imports | Removed in COMPAS 2 |
| 24 | `collateral_strip_deletions`, `total_boundary_deletions` | Already in `grammar_pattern.py` |
| 25 | `TwoColourableProjection` → `TwoColorer` rename | Referenced in our docs and `algorithms/__init__.py` |
| 26 | `add_strips` (upstream version) | Ours updates the remaining polyedges; theirs has it as a TODO |
| 27 | `utilities.py` | Loses `common_items` / `remove_isomorphism_in_integer_list`; has an `__all__` typo |

---

## 7. Suggested sequencing

1. **Zero-risk batch first** — items 4, 5, 6, 9, 10 and the `__all__` fixes from §5. No behavioural
   change to any layout; run the existing tests and move on. ✅ **done 2026-08-26**
2. **Polyedge batch** — items 1, 2, 3, 12. These *do* change coarse layouts. Run the 84-row frame-field
   baseline before and after and diff it; expect movement, and decide row by row.
   ✅ **done 2026-08-26** — but not the way this line assumes, on two counts. There is no
   `framefield/` package in this checkout, so the 84-row baseline does not exist; a purpose-built
   invariance harness over three meshes was used instead. And because everything landed behind
   defaulted keywords, **no coarse layout moved at all** — the diff is empty by construction, and the
   row-by-row decision is deferred to whoever flips `strict` / `legacy`. §0 records what those flags
   do when switched on.
3. **Densification batch** — item 11. Highest value, highest care. Gate it behind
   `densification_old()` and diff dense meshes across the example suite. Verify the curved-boundary
   coverage number (0.996) has not regressed. ⏳ **not started**
4. **New capability** — items 7, 8, 13, 15, 16. Additive; can land any time after batch 1.
   ⏳ **not started**
5. **Never** — everything in the "do not take" table. Worth recording so this comparison does not get
   redone in six months. ✅ still rejected; nothing from it was taken.

### Suggested next step

Batches 1–2 delivered the fixes but left them dormant, which is the deliberate consequence of the
opt-in policy. The natural follow-up is a **decision pass** on the three flags now that the machinery
and the A/B harness exist: `polyedge_graph(legacy=False)` looks safe to flip (it is strictly more
correct and the colourings agree on real meshes), `collect_polyedge(oriented=True)` needs a check of
the two Grasshopper call sites, and `strict=True` needs the coarse-layout diff on real project meshes
— it is a no-op on the British Museum fixture but changes everything on meshes with holes.

---

## 8. Appendix — reproduction

All measurements above came from scripts run with
`C:\Users\Casper\anaconda3\envs\singular312\python.exe`, with `compas_quad` on the path via

```python
sys.path.insert(0, r'...\compas_singular\compas_quad-main\compas_quad-main\src')
```

Scratch scripts (session-local): `t1_densify.py` (densification parity + density-0 bug),
`t2_poles.py` (pole densification), `t3_weld.py` (geometric vs topological weld),
`t5_opp2.py` (`vertex_opposite_vertex`), `t6_lizard.py` (winding),
`t7_perf.py` (polyedge benchmarks), `t9_json.py` / `t11_misc.py` (serialisation, `__all__`),
`t10_orient.py` (polyedge orientation).

### Appendix B — reproducing the batch 1–2 verification (2026-08-26)

The evidence for §0 is checked in, not session-local. Note that `compas_singular` may resolve to a
different clone on this machine, so put this tree first:

```powershell
$env:PYTHONPATH = "<repo>\src"
python -m pytest tests -q          # 59 passing
```

* `tests/test_mesh_quad.py` — null guards, the six accessors, and each opt-in flag, asserted against
  the annulus and 4-valent-boundary fixtures from §4.4(b). Both the fixed *and* the historical
  behaviour are pinned, so flipping a default will fail a test rather than pass silently.
* `tests/test_import.py` — the namespace assertions from the §5 update: every subpackage imports, no
  package re-exports a submodule, all four grammar names are callable, and `delete_strip(s)` keeps
  `preserve_boundaries`.

The default-invariance and A/B harnesses were session-local: `invariance.py` (snapshot/diff of
polyedges, strips, graph edges, colourings, singularities and decomposition), `ab.py` (flag
characterisation), `bench.py` (`collect_polyedges` old vs new at three mesh sizes). Their conclusions
are recorded in §0; re-create them from that description if the numbers need re-checking.
