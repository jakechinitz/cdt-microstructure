# The admissibility-closure coupling (`v6_closure_run.py`) — model record

This documents the second theory coupling: the paper's **own** cell weighting
(Sections 5–6, Appendix B of the Entropic Scalar EFT draft) attached to the
v6 CDT substrate. It exists because the paper explicitly states its
weighting is *not* a spin-foam vertex amplitude — "it solves a finite local
boundary-counting problem and weights configurations by admissibility
closure" (§26.8) — so the EPRL run, even fixed, tests a neighboring theory.
This coupling tests the stated one.

## 1. What the paper fixes (no freedom here)

Verbatim from Appendix B, implemented in `build_energy_table()` /
`ensemble_report()` and **verified numerically in this repo**:

- Cell = tetrahedron; each face carries a label m ∈ {−3,…,3} (seven-state
  face sector, j_eff = 3).
- Injectivity: the four face labels are distinct. With parity doubling the
  isolated-cell state count is Ω_tet = 2·P(7,4) = 1680.
- Closure invariant: K²(m) = 48 − (S² − Σ²)/3, S = Σmᵢ, Σ² = Σmᵢ².
- Admissibility weight: p_η(b) ∝ exp(−η K²(b)) at the paper-fixed precision
  η\* = 0.0298668443935 (stationary normalized closure evidence).
- Checks reproduced exactly by this code: ⟨K²⟩_η\* = 50.223 = 3/(2η\*), and
  the admissibility-weighted entropy g_share,eff = **7.4198** (paper's
  value; raw ceiling ln 1680 = 7.4265). Every run prints this fidelity line
  at startup.

## 2. The identification (the one modeling step)

The theory's volumetric cells are identified with the **spatial tetrahedra
of the CDT time slices**; their triangular faces are the shared faces
through which neighboring cells correlate. This is the natural reading —
the paper's cells are 3D volumetric boundary cells, and the CDT slices are
exactly the spatial 3-manifolds of the geometry — and it is what makes the
coupling *geometric*: labels live on slice triangles, each shared by
exactly two cells, so the joint admissibility statistics depend on the
slice's gluing pattern. That dependence is the entire coupling channel.

**Precondition (engine fix, same audit cycle):** each triangle shared by
exactly two spatial tets requires the slices to be closed simplicial
3-manifolds. The pre-fix v6 move set violated this (see §5 below); closure
runs enforce `--causal-slices` (default) and the run refuses to be read on
defective slices.

## 3. The coupled model

    pi(g, m) ∝ exp( −S_Regge(g) − eps·(N41−target)²
                    − β · Σ_tets [ E(m_faces(tet)) − mu ] )

    E(m₁..m₄) = η·K²(m) + λ_inj · (# equal-label pairs)

with uniform base measure 1/7 per triangle label (per-face normalized label
sum). Sampling reuses the exact-sampler pattern from the fixed EPRL run:
label heat-bath at the same exponent β as the geometry action; uniform
label births at geometry moves (exactly detailed-balanced under the
normalized base measure); mu = per-tet label free energy by thermodynamic
integration, stored with a calibration context in the checkpoint.

The geometry marginal is w(g) ∝ E_{m~unif}[e^{−βΣE}]^(centered) — the
slice-dependent **closure free energy**. Geometries whose slice gluing
makes admissibility-closed labelings entropically easier are favored: the
theory's "weights configurations by admissibility closure", exactly.

**The theory-designated point is β = 1, η = η\*, λ_inj → hard.** β = 0 is
the bare control; the β sweep is a diagnostic ramp, not the claim.

### Declared forks (modeling choices the paper does not pin down)

- **F1 — soft injectivity.** The paper's injectivity is a hard constraint;
  the run uses a penalty λ_inj per colliding pair (default 3.0) because a
  hard 7-color constraint at max conflict degree 6 does not guarantee an
  ergodic single-face heat-bath (q = Δ+1 is the marginal case). The hard
  model is the λ_inj → ∞ limit; runs report the collision fraction (≈9% at
  β=1, λ=3 on small tests) so the softness is visible. Sensitivity in
  λ_inj should be checked before strong claims.
- **F2 — face orientation.** The closure sum is over *oriented* faces; on a
  shared face the two cells see opposite normals, suggesting the neighbor
  reads −m where the cell reads m. Implementing that consistently requires
  oriented-slice bookkeeping through all Pachner moves; the current model
  uses the unoriented reading (both cells see m). K² is not invariant under
  flipping one m, so F2 is a real fork — flagged for follow-up, not silently
  resolved.
- **F3 — parity.** Parity doubling contributes exactly ln 2 per cell —
  extensive in the pinned N41 — so it drops out of the coupled dynamics and
  is omitted.

## 4. Why the EPRL audit findings don't recur (by construction)

| EPRL audit finding | Closure model |
|---|---|
| 1. Centering not volume-neutral | mu = label free energy by TI, and both N_tets (= N41/2) and N_triangles (= 2·N_tets) are **pinned by the N41 penalty** — the extensive part cannot drive volume at all; mu mainly keeps the (4,1)/(3,2) mix unbiased. |
| 2. Heat-bath ignored β | β-consistent from the first commit (`_heatbath_pass(..., beta, ...)`). |
| 3. Label births outside Hastings | Uniform births exactly detailed-balanced under the normalized base measure (same theorem as the EPRL fix). |
| 4. Slot convention noise | E is permutation-symmetric in the four faces by construction — there is no slot. |
| 5. Sign problem | The weight is positive. Nothing is discarded; the sampled model IS the stated model (up to F1–F2). |

Placebo arm (`--placebo`): the 210 sorted-label-orbit energies are shuffled
— same value pool and permutation symmetry, closure structure destroyed.
Only a real-vs-placebo difference at matched β is attributable to the
closure weighting.

## 5. Engine precondition: the causal-slice fix (found during this build)

While validating the cell identification, a **bare-engine defect** was
found: the v6 move set preserves the 4D manifold conditions but not the CDT
foliation. A (2,4) move on a (2,2)-interface tetrahedron with both apexes
at the same time creates a spatial tet whose two pentachora point the same
way in time (a slice "fold"); (3,3) has analogous cases. Measured on
pre-fix runs: ~5% of slice triangles off the closed-3-manifold condition
(incidence 1/3/4 instead of 2). This means pre-fix runs — including the
bare-baseline volume scans — sampled a *generalized* ensemble, not the
standard AJL CDT ensemble the protocol cites.

Fix: `causal_slice_ok()` (v6_cdt_moves.py) rejects any move that would
break the foliation — for every spatial triangle of a new pentachoron,
require exactly two containing spatial tets, each with one future and one
past apex. The check is complete (an affected triangle's link circle
necessarily passes through a new pentachoron) and preserves detailed
balance exactly as `strictness_ok` does. Verified: clean over 42k
proposals with undo cycles, all seven moves still fire, growth to target
works, and defective old checkpoints **heal monotonically** on resume
(101 → 81 defects in 6 sweeps in the test) since defects can be destroyed
but never re-created. Every run now prints a `slice` column and a final
foliation verdict, and the flag `--no-causal-slices` reproduces the old
ensemble for comparison.

**Consequence for existing results:** the bare d_H ≈ 3.37 @ 20k baseline
and all prior sweeps were measured on the generalized ensemble. The bare
volume scan should be redone with the fixed engine before the next round
of theory claims (expect quantitative shifts; direction unknown).

## 6. Protocol

```bash
# arms resume from a CLEAN bare checkpoint grown with the fixed engine:
./run_closure_sweep.sh "0.0 0.3 1.0" 20000       # + placebo at beta=1 automatically
```

Gates, in order: (G1) validity + foliation CLEAN on every arm; (G2)
thermalized; (G3′) matched N4 across arms; (G4) β=0 reproduces bare;
(G5) chain not frozen. Then read d_H / blob / cos³ vs the β=0 control and
the β=1 placebo. Success band mirrors the pre-registered EPRL reading:
flat-or-better for the real arm **with** real-vs-placebo separation if any
movement is claimed.

## 7. Machinery acceptance test (small volume — NOT a physics result)

Three arms, K=16, N41=2000, eps=0.01, 400 sweeps each, all resumed from one
CLEAN bare base grown with the fixed engine (single seed, base still slowly
growing — machinery check only):

| arm            | final N4 | d_H  | blob | ⟨E⟩/tet | tets w/ collision |
|----------------|----------|------|------|---------|-------------------|
| β=0 control    | 3640     | 2.48 | 2.10 | 3.92    | 0.631             |
| β=1 real       | 3420     | 2.47 | 2.62 | 1.73    | **0.078**         |
| β=1 placebo    | 3374     | 2.30 | 3.32 | 1.97    | 0.733             |

Checks: foliation CLEAN and audits at ~1e-12 drift on all arms; N4 matched
within a few % (no volume leak); the collision fraction behaves exactly as
designed — 0.631 at β=0 (uniform prediction 1 − 840/2401 = 0.650), driven
to 0.078 by the real closure weight at β=1, while the placebo (scrambled
energies, no injectivity structure) stays high (0.733). The d_H/blob spread
between arms at this size is within single-seed noise of a non-equilibrated
base — read nothing into it; the production sweep with gates is the
measurement.
