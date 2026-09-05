# Cross-Field Quad Meshing as an Alternative to compas_singular

**Date:** 2026-08-03
**Question:** Should the new workflow use cross fields to generate quad meshes inside an arbitrary
boundary, with mesh edges constrained perpendicular to guide curves?
Is there a better approach?

**Input assumption (important):** guide curves are **not** always force lines. They may also be a
**cable layout** — a construction-driven set of curves with no guarantee of mutual orthogonality.
The whole evaluation below is structured around that.

---

## 1. Verdict

**Direction-field-guided meshing is correct.** It is the only method family that admits interior
directional constraints at all — paving, indirect (Blossom-Quad / Catmull–Clark) and quadtree
overlay have nowhere to put them.

But **"cross field" is the wrong commitment.** A cross field is 4-fold rotationally symmetric,
which *hard-codes orthogonality of the two edge families*. That is exact for force lines and
over-constrained for a general cable layout.

The right commitment is one level up: **frame fields** (Panozzo et al. 2014), of which cross fields
are the orthogonal, unit-length special case. This costs almost nothing — see §3.

### Why force lines are the easy case

In plane stress the two families of principal directions are perpendicular everywhere except at
isotropic points, and those isotropic points are exactly index-±½ singularities. A principal stress
field **is** a cross field. You are not approximating the constraint — you are transcribing it, and
the field can be constrained directly from the stress tensor rather than solved for smoothness and
then hoped to stay aligned.

### Why a cable layout may not be

A cable layout carries no such guarantee. Whether a cross field still suffices depends entirely on
how many cable families there are and how they meet — see the case analysis next.

---

## 2. Input taxonomy — this determines everything downstream

| Case | Input | Is a cross field sufficient? |
|---|---|---|
| **A** | Two orthogonal families (principal stress trajectories) | **Yes, exactly.** Both cross directions are pinned and mutually consistent. |
| **B** | **One** family (a single cable set: parallel-ish, radial, fanned) | **Yes, exactly.** Only one of the two directions per point is constrained; the perpendicular family is *derived* by the smoothness solve. |
| **C** | Two or more families meeting at non-90° angles (e.g. a cable net crossing at 60°) | **No.** Over-constrained. Needs a frame field. |
| **D** | Families that cross at *varying* / locally contradictory angles | **No, and neither does a frame field** at the crossings. Needs constraint prioritisation or domain cutting (§5.2). |

**Case B is the one most likely to be misclassified as hard.** "Mesh edges perpendicular to the
cables" with a *single* cable family is precisely the situation a cross field handles best: you pin
one direction, and perpendicularity of the other is not an extra constraint but the definition of
the object. If the cable layout turns out to be single-family, the original cross-field plan stands
unmodified.

### A semantic point that only bites in cases C/D

Under 4-fold symmetry, *"edges perpendicular to the guide curve"* and *"edges tangent to the guide
curve"* are **the same constraint**, differing by a 90° phase. So the perpendicular/tangent
distinction is free for cases A and B.

It stops being free once the frame is non-orthogonal, and then it becomes a real design question
that must be answered per input type:

- **Force lines** → the mesh edges should *follow* them; the perpendicular family comes along
  automatically.
- **Cable layout** → in a real cable net the cables *are* mesh edge chains, bounding the quads,
  rather than being perpendicular to them.

Both of these are **alignment**, not perpendicularity. Worth confirming that "perpendicular" in the
original brief was not shorthand for "aligned" — the answer changes what gets pinned.

---

## 3. Recommended representation: frame field, not cross field

Panozzo et al. give the clean route, and it is why the generalisation is cheap:

1. Solve for a **frame field** satisfying a sparse set of constraints (non-orthogonal,
   non-unit-length, smoothly varying linear transformations on tangent spaces).
2. Compute a **surface deformation that warps the frame field into a cross field**.
3. Run the ordinary cross-field quad-meshing pipeline **in the deformed domain**.
4. Map the resulting mesh back.

So the entire pipeline in §4 survives; it gains a warp step in front and an unwarp step behind.
Cases A and B simply produce an identity-ish deformation and reduce to the original plan.

Bonus: the same machinery controls **element density and anisotropy**, not just alignment — likely
useful for grading mesh density along cables.

Alternative representation if the deformation route proves awkward: **N-PolyVector fields**
(Diamanti et al. 2014), which drop orthogonality and rotational symmetry and are still computable
by a **sparse linear solve with no integer variables**.

---

## 4. Recommended pipeline

A path named in neither branch of the earlier discussion, and the one that fits this codebase:

```
guide curves (force lines and/or cables)
  → frame field  [ → warp to cross field, if non-orthogonal ]
    → singularity graph
      → separatrix tracing
        → coarse quad layout (four-sided patches)
          → compas_singular densification / strip machinery
            [ → unwarp ]
```

The mixed-integer parametrization is **skipped entirely**. Separatrices traced from each
singularity partition the domain into four-sided patches; each patch is meshed by structured /
transfinite mapping.

This is Viertel & Osting's pipeline (with provable bounds on singularity count and location for
planar domains), and it is the topological half of Gmsh's `quadqs`.

### Why this fits compas_singular specifically

compas_singular's data model **is already** a coarse quad mesh plus strip densification:
`add_strip`, strip-based densification, `edges_to_curves` for boundary recovery.

Only the **front end** is replaced — the medial-axis skeleton decomposition that currently produces
the coarse layout gives way to a field-derived layout. Everything downstream survives, including
the strip-editing design space that is the main reason to use compas_singular at all. Oval's own
later work moves the same direction (the 2023 vector-encoding paper decouples singularity topology
from geometry).

### Note on difficulty — a correction to the earlier framing

The earlier analysis treated the field as the interesting stage and the parametrization as merely
"technically delicate." In engineering-risk terms that is backwards:

| Stage | Real difficulty |
|---|---|
| Boundary-aligned cross/frame field on a planar domain | **Low.** Sparse linear solve in the power-4 representation, Ginzburg–Landau PDE formulation, or PolyVector complex-polynomial formulation. |
| Mixed-integer parametrization + quantization | **High.** Seam placement, integer rounding, injectivity, degenerate zero-width strips. |
| Mesh extraction from the integer grid | Moderate; mostly bookkeeping once the parametrization is valid. |

Bommes' *Mixed-Integer Quadrangulation* and *Integer-Grid Maps* exist precisely because the middle
row is hard; Lyon et al.'s quantization work exists because IGM still was not reliable enough.

> **Do not reimplement MIQ.** If the full parametrization route is ever taken, drive **libigl** or
> **Gmsh** instead of writing it from scratch. The separatrix route above avoids the question.

---

## 5. Two decisions to make before writing code

### 5.1 Which input case is being targeted — A, B, C or D?

Decide this first; it fixes the field representation (§2, §3). If both force lines and cable
layouts must be supported by one code path, implement the **frame field** and let cases A/B fall
out as the degenerate orthogonal instance. Do not implement a cross field and retrofit
non-orthogonality later — the retrofit is the deformation step, and it is easier to have it present
from the start than to thread it in.

### 5.2 Must the guide curves be *exactly* mesh polylines, or approximately?

Field alignment is **soft**. Even with hard field constraints, extraction can shift an isoline off
the guide curve by a fraction of an element.

If exact incidence is required:

- **Option A — feature-curve constraints in the parametrization.** MIQ's feature alignment, or
  Campen & Kobbelt's aligned parameterization. Correct but expensive.
- **Option B — pre-cut the domain along the guide curves.** *(recommended, and now more strongly
  than before)*

  Interior guide curves become internal boundaries. The standard boundary-tangency constraint then
  gives exact incidence **and** perpendicularity for free, and each subdomain meshes independently.

  **This matters more given non-orthogonal input.** Cutting absorbs the non-orthogonality at the
  patch *corners*, where it is a local topology problem, instead of fighting it globally in the
  field. For case D (contradictory crossings) it may be the only tractable option.

  compas_singular already has curve-feature handling in its skeleton decomposition, so this is also
  the shortest path to something working.

  **Caveat:** the headless feature-curve path is fragile and leaves interior triangles. Budget time
  to fix that before relying on it.

---

## 6. Method families considered and rejected

| Method | Why not |
|---|---|
| **Instant Meshes / QuadriFlow** | Do support orientation-field constraints and are extremely robust, but produce **no coarse layout** and do **not conform exactly** to boundaries. Wrong output type for structural design. |
| **Direct streamline tracing** (guide curves + orthogonal trajectories; Alliez et al.'s anisotropic remeshing) | Tempting and simple, but degenerates near isotropic points with no clean way to close the mesh there. Also assumes orthogonality outright. Fine for visualization, not for an analyzable mesh. |
| **Paving / advancing front** (Blacker & Stephenson) | Boundary-perpendicular rows for free, but edge flow is dictated entirely by boundary offsets. No place to inject interior guide curves. |
| **Indirect** (triangulate → Blossom-Quad / midpoint subdivision) | Robust and simple, but edge flow is whatever the triangulation gave. Will not follow cables or force lines. |
| **Grid overlay / quadtree** | Fast, poor boundary elements, no alignment control. |

---

## 7. Sources

### Non-orthogonal fields — now the primary reference

- [Frame Fields: Anisotropic and Non-Orthogonal Cross Fields (Panozzo, Puppo, Tarini, Sorkine-Hornung, SIGGRAPH 2014)](https://igl.ethz.ch/projects/frame-fields/)
  — the warp-to-cross-field construction; also controls density and anisotropy
- [Designing N-PolyVector Fields with Complex Polynomials (Diamanti et al. 2014)](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.12426)
  — sparse linear solve, no integer variables
- [Designing 2D and 3D Non-Orthogonal Frame Fields (Corman et al.)](https://members.loria.fr/ECorman/Papers/ff_nonortho.pdf)

### Field-guided pipeline — foundations

- [An Approach to Quad Meshing Based on Harmonic Cross-Valued Maps and the Ginzburg–Landau Theory (Viertel & Osting)](https://arxiv.org/pdf/1708.02316)
  — the separatrix-partition route, planar domains, with guarantees
- [Quasi-structured quadrilateral meshing in Gmsh (Reberol, Georgiadis, Remacle 2021)](https://arxiv.org/abs/2103.04652)
  — production-grade version of the recommended pipeline
- [Computing two-dimensional cross fields — a PDE approach based on the Ginzburg–Landau theory (Beaufort, Remacle et al.)](https://www.researchgate.net/publication/317356724_Computing_two_dimensional_cross_fields_-_A_PDE_approach_based_on_the_Ginzburg-Landau_theory)
- [Ginzburg–Landau energy and placement of singularities in generated cross fields](https://arxiv.org/pdf/2010.16381)
- [Mixed-Integer Quadrangulation / Integer-Grid Maps (Bommes et al.)](https://www.researchgate.net/publication/262352155_Integer-Grid_Maps_for_Reliable_Quad_Meshing)
  — the reference formulation for alignment-constrained quad meshing
- [Parametrization Quantization with Free Boundaries for Trimmed Quad Meshing (Lyon, Bommes, Kobbelt 2019)](https://cgg.unibe.ch/media/papers/1269/a51-lyon.pdf)
  — why quantization is the hard part

### Structural / force-line alignment

- [Aligning principal stress and curvature directions (Pellis & Pottmann, AAG 2018)](https://www.geometrie.tuwien.ac.at/geom/ig/publications/principalstress/principalstress.pdf)
- [Optimized Quad Gridshell from Stress Field and Curvature Field](https://www.researchgate.net/publication/326540880_Optimized_Quad_Gridshell_from_Stress_Field_and_Curvature_Field)
- [Stress Line Generation for Structurally Performative Architectural Design (Tam & Mueller, ACADIA 2015)](https://papers.cumincad.org/data/works/att/acadia15_095.pdf)

### compas_singular lineage — read before replacing the front end

- [A vector encoding for topology finding of structured quad-based patterns (Oval, Mesnil, Van Mele, Block, Baverel 2023)](https://journals.sagepub.com/doi/10.1177/09560599231207650)
- [Similarity-driven topology finding of surface patterns for structural design (Oval et al. 2024)](https://block.arch.ethz.ch/brg/files/OVAL_2024_CAD_similarity-driven-topology-finding-of-surface-patterns-for-structural-design_1720030018.pdf)
- [Two-Colour Topology Finding of Quad-Mesh Patterns (Oval et al. 2021)](https://www.sciencedirect.com/science/article/pii/S0010448521000415)
- [Feature-based Topology Finding of Patterns for Shell Structures (Oval et al. 2019)](https://www.researchgate.net/publication/331064073_Feature-based_Topology_Finding_of_Patterns_for_Shell_Structures)

### General

- [Quad meshing survey / paper compilation (continuously updated)](https://github.com/Bigger-and-Stronger/quad-meshing-survey)
