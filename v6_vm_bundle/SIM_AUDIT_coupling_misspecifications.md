# Sim audit — where the EPRL coupling is misspecified (2026-07)

Audit of the v6 theory run (`v6_theory_run.py` + `v6_run_lib.py`) prompted by
the question: *is the matched-volume β-sweep testing what we think it is?*
Answer: **no — as shipped, the swept signal is dominated by β-dependent
sampler artifacts that the β=0 control is structurally blind to.** Five
concrete misspecifications, each quantified below on this repo's own tensor
and engine. Findings 1–3 are fixable bugs; findings 4–5 bound what any
frozen-j run can claim even after the fixes.

The control/regge/ours design, the G1–G5 gates, and the prereg discipline are
sound. The problem is that every artifact below **vanishes identically at
β=0 and grows monotonically with β** — i.e. each one masquerades exactly as
the signal the sweep is designed to detect. A clean β=0 control (G4) cannot
catch any of them.

---

## 1. The centering is not volume-neutral (the big one)

`Centering.calibrate()` fixes mu **once, on the first action evaluation** —
which on a fresh run or a bare-checkpoint resume happens with
**uniform-random labels**. Measured: uniform-label mean cost = **30.53**
(the prereg's mu ≈ 30.54 confirms the production sweep calibrated exactly
there). After the heat-bath equilibrates the labels, the mean cost drops to
**≈ 28.0**. So from sweep 1 onward the "centered, volume-neutral" term is
actually an **extensive offset ≈ −2.5 per pentachoron**, i.e. an
uncontrolled effective cosmological-constant shift

    k4_eff ≈ k4 − 2.5·β      (bare k4 = 0.9)

| β    | k4 shift (as shipped) | as % of bare k4 |
|------|----------------------|-----------------|
| 0.05 | −0.13                | 14%             |
| 0.1  | −0.25                | 28%             |
| 0.3  | −0.76                | 84%             |

The eps-penalty pins **N41 only**; the N32 sector is unpinned, so the shift
drives the N32/N41 mix — a size→shape confound in exactly the observables
(d_H, blob) the sweep reads. This is very likely why the earlier runs failed
the G3 matched-volume gate.

Reproduce: `python centering_beta_audit.py`.

**Fix:** calibrate mu per β from label-equilibrated configurations (run the
heat-bath to equilibrium at fixed geometry *before* calibrating), and/or gate
on matched **N4** (not just N41) across the sweep. Do not center on a running
mean — that breaks detailed balance.

## 2. The heat-bath ignores β

The geometry moves feel `β·S_EPRL`, but `make_heatbath()` resamples labels
with weight `Π|A|^1` at **every** β (`lw += np.log(amp)` — no β). The joint
(geometry, labels) chain therefore has no consistent stationary measure: the
labels always sit at the maximally-concentrated β=1 conditional. Measured
consequence (β-consistent heat-bath on a fixed geometry, `centering_beta_audit.py`):

| β    | correct offset·β | shipped offset·β | inflation |
|------|------------------|------------------|-----------|
| 0.05 | −0.005           | −0.128           | ~25×      |
| 0.1  | −0.032           | −0.256           | ~8×       |
| 0.3  | −0.247           | −0.769           | ~3×       |

**Fix:** one line — `logw[c] = beta * lw` (using the same β as the geometry
action; in `eprl_only` mode β=1 is correct). Then recalibrate mu per β
(finding 1), because the correct offset is itself β-dependent.

## 3. Label births are outside the Hastings ratio

Tetrahedra created by a geometry move get a fresh **uniform** label, and
neighbours' shared-face labels are refreshed, with neither the proposal
density (D^±Δtets) nor the fresh-label cost (uniform labels cost ≈ +2.5 above
equilibrated) folded into the acceptance ratio. `hastings_log()` corrects the
geometry proposal only. The docstring triages this as "adequate for the
steer-direction comparison" — but the bias is **β-scaled and acts precisely
on the volume-changing moves**, so it contaminates the swept signal rather
than a nuisance direction. Combined with finding 1 it also creates a ratchet
(a region grows at near-neutral EPRL cost with fresh labels, the heat-bath
then relaxes them, and shrinking the region afterwards costs +2.5·β per
pentachoron).

**Fix:** draw new labels from their heat-bath conditional and include the
proposal ratio, or gate volume moves through a label-marginalized weight.

## 4. The slot assignment is convention-noise (measured, large)

The vertex tensor is **strongly non-symmetric** under permutations of its
five intertwiner slots: adjacent-swap relative asymmetry ≈ **1.20**, worst
permutation ≈ **1.35** (order-unity — the permuted tensor is as different
from the original as a random other tensor). But `Intertwiners.faces()`
assigns slot k to "neighbour across the face opposite vertex k" in the
engine's internal vertex ordering — which has no relation to the sl2cfoam
recoupling/pairing basis, and the same shared tetrahedron enters its two
pentachora at unrelated slot positions. The quantity being coupled is
therefore **not the EPRL amplitude of the triangulation** (which would need
the convention-checked {15j} gluing flagged in `VERTEX_PROVENANCE.md`); it
is the tensor sampled at essentially arbitrary slot assignments. Also
unresolved at the same level: whether sl2cfoam-next's gluing convention
requires (2i+1) edge-dimension weights on the intertwiner sum, and the
missing face-amplitude factor ((2j+1)^{N2} at frozen j — extensive in
(N0, N4), i.e. another β-dependent tilt of the effective bare couplings if
the theory term is meant to be the full EPRL measure rather than the vertex
part alone).

## 5. |A| discards an order-unity interference structure

49.6% of tensor entries are negative (49.5% of the |amplitude| mass). For
the true gluing, the shared-tet intertwiner sum Σ_i A·A has measured
per-edge coherence |Σ A·A| / Σ|A||A| with **median 0.44**, and **19% of
edges are >90% sign-cancelled**. The positivized Π|A|^β model is therefore
not a mild approximation to the EPRL weighting — it is a different theory,
and the difference compounds over thousands of internal tets. (Positivization
is unavoidable for Metropolis; the claim just has to be about the positivized
model, and peaked-j makes this *worse*, not better — semiclassical vertex
amplitudes oscillate like e^{iS_Regge}.)

---

## The placebo control: shuffled tensor reproduces the "steer"

Decisive check (new: `make_shuffled_control.py`): rerun the sweep with an
**entry-shuffled** copy of the tensor — identical value distribution
(same mean/std of −log|A|, same "teeth" σ), all EPRL correlation structure
destroyed. Every artifact above depends only on entry statistics, so:

* shuffled reproduces the real-tensor movement at the same β ⇒ the movement
  is machinery, not EPRL;
* only a real-vs-shuffled **difference** at matched β is attributable to the
  amplitude.

Small-volume result from this audit (K=16, N41=2000, eps=0.01, k4=0.9,
400 sweeps, matched N41 ≈ 2000):

| run                  | final N4 | final blob | blob (late-run) | cos3err (late) |
|----------------------|----------|------------|-----------------|----------------|
| β=0 control          | 2962     | 2.70       | ~2.0–2.7        | ~0.7–0.8       |
| β=0.3, real tensor   | 2714     | 3.21       | ~3.2–4.0        | ~0.5–0.6       |
| β=0.3, **shuffled**  | 2660     | 3.20       | ~3.1–3.8        | ~0.6–0.7       |

The shuffled placebo reproduces the real-tensor signature (N4 suppressed at
matched N41, blob elevated) — at this size the entire visible β-effect is
coupling machinery. **The production sweep should not be interpreted until a
matched-volume shuffled-control arm runs at 20k and comes back different from
the real-tensor arm.**

---

## What this means for the frozen-j vs peaked-j question

Frozen-j is not just "possibly too weak" — it is structurally unable to
express the main channel by which a spin-foam amplitude could steer CDT
geometry: with every triangle at the same spin, all areas are equal, so the
amplitude cannot see hinge deficit angles; only intertwiner correlations
remain, and `eprl_term_diagnostic.py` already measured those to be type-blind
(|split| ≈ 0.05 ≪ σ ≈ 0.87). So even a perfectly implemented frozen-j
coupling is expected to be flat *for reasons that say nothing about the
theory*.

Peaked-j is testable in principle — sl2cfoam-next computes mixed-spin
vertices, so per-triangle spins j_t ∈ {2,3,4} could be sampled jointly with
intertwiners under a Gaussian peaking weight, with vertex tensors cached on
demand (≤3^10 spin combos per vertex, far fewer in practice). But it inherits
findings 2–5 and sharpens the sign problem. Fix the machinery first, and
demand the shuffled placebo separates from the real tensor, before spending
on peaked-j.

Also worth stating: per the paper (§26.8, App. B), the theory's own
weighting is **admissibility-closure counting on the 1680-state boundary
ensemble** — explicitly *not* a spin-foam vertex amplitude — and its "7" is
the seven **face labels** (m = −3…3 of j_eff = 3), not the seven-dimensional
**intertwiner space** of the EPRL j=3 vertex that the sim sums over. The two
sevens are different objects. The theory's actual per-cell weight is
positive, local, and cheap — it has no sign problem, no sl2cfoam convention
dependence, and no frozen-j truncation issue. Coupling *that* weight to the
CDT substrate is both closer to the paper's claim and far easier to get
right than EPRL.

---

*Files added with this audit: `centering_beta_audit.py` (findings 1–2),
`make_shuffled_control.py` (placebo arm). Both run on the shipped bundle
with no extra dependencies.*
