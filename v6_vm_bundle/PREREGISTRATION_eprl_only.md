# Pre-registration — `eprl_only` (no-Regge) run

Written **before** the result. The `eprl_only` mode removes the Einstein-Regge
action entirely and lets the (faithful, sl2cfoam-derived) EPRL amplitude be the
only thing selecting geometry — subject to the volume penalty and CDT's move-set.
This bears more directly on "does the substrate generate spacetime?" than the
regge_plus_eprl sweep, because no GR is inserted. Committed for timestamp.

## Established facts going in (do not re-litigate)

- Tensor is FAITHFUL: validated=1, j=3, γ=1.2, shells=5, rho_gj; stats match all
  prior diagnostics. The geometry-blindness is a property of the CORRECT
  amplitude (type-split 0.054 ≪ σ=0.87), not a convention artifact. No escape
  hatch.
- Volume IS still pinned in eprl_only: run() always adds eps·(N41−target)². So
  this is "EPRL action + volume penalty + CDT moves," NOT "geometry with nothing
  holding it." Predicted failure is therefore NOT volume runaway / expander
  collapse (the condensate-model prediction does not transfer — CDT's manifold
  constraint forbids arbitrary expander graphs).
- CENTERING is the experiment's fork:
  * `--no-center-eprl`: dominated by the amplitude's μ·N4 piece (an
    artefact-grade effective cosmological constant) → crushes volume to minimum,
    freezes. This re-tests a known artifact; NOT the clean experiment.
  * `--center-eprl` (USE THIS): removes the trivial μ·N4 scale, eps holds volume,
    and the non-trivial EPRL FLUCTUATIONS select shape. This is the clean test of
    "does the amplitude's structure pick 4D geometry."

## Pre-registered interpretation rule (fix BEFORE the read)

Run: centered eprl_only, matched volume (eps holds N41), thermalized, valid
manifold throughout. Read the d_s(σ) flow and the cos³ blob, not d_H at sweep 20.

| Outcome | Reading |
|---|---|
| **de Sitter**: cos³ blob forms, d_s(σ) flows up toward ~4 with scale, valid manifold | **POSITIVE — the substrate's amplitude generates 4D geometry.** Striking, publishable, supports the generation claim. (Low prior given the blindness diagnostics, but it is exactly the result that would mean what we want.) |
| **entropic collapse**: blob→1 (flat/stalk), low d_H, d_s flat/low, possibly crumpled | **NEGATIVE — the substrate does NOT generate 4D on its own at frozen j=3.** The 4D in regge_plus_eprl came from the Einstein-Regge action, not the amplitude. This is a REAL finding, NOT a neutral null, NOT rescued by "Regge works / the tetrahedra are fine." |
| **non-manifold / freezes / N41 won't hold** | setup failure (or acceptance crushed) — not a physics result; diagnose, don't interpret. |

The negative is the likely outcome and it is the hard-to-want one. Pre-committing:
an entropic collapse under centered eprl_only is the **direct negative answer** to
"can the substrate's dynamics generate spacetime," and will be written as such.

## Bound on the claim either way (state in the writeup)

Even eprl_only is NOT "substrate dynamics from scratch." CDT's MOVE-SET still
proposes geometry changes and the manifold constraint defines legal geometries —
that kinematics is AJL's. The amplitude only SELECTS among CDT-proposed
geometries. So:
  * a POSITIVE is "the EPRL amplitude selects 4D de Sitter among CDT-proposed
    geometries with no Regge action" — strong, but not "geometry from nothing."
  * a NEGATIVE is "the amplitude fails to select 4D among CDT-proposed
    geometries; volume-fixed entropy + frozen-j amplitude relax to the
    branched-polymer/crumpled phase."
Do not state the generation claim more broadly than "amplitude selects among
CDT-proposed geometries."

## Frozen-j caveat (unchanged, applies here too)

eprl_only freezes j=3 and evolves the INTERTWINER labels (i∈{0..6}) by heat-bath.
Whether "the substrate's dynamics" in the paper's sense is faithfully captured by
"freeze j, evolve intertwiners under the EPRL vertex, CDT moves for geometry" is a
real question the paper does not settle. A negative here is a negative about the
FROZEN-j realization, and could in principle be lifted by peaked-j (which does not
exist yet) — same standing caveat as the sweep.

## Timing asymmetry + transient trap (pre-register BEFORE reading)

Collapse and de Sitter have very different detection times, and the inconvenient
one is fast:
- **Collapse is FAST & decisive.** A geometry-blind amplitude that can't hold 4D
  falls toward the entropic phase quickly (uncentered smoke reacted ~20 sweeps).
  "Did it fail to hold a manifold / stable shape" — answerable in hours, decisive
  if NO.
- **De Sitter is SLOW & deceptive to confirm.** d_s approaches 4 FROM BELOW with
  volume (β=0 scan: d_s_peak 2.36→2.65 over N4 2k→18k, peak σ marching out). So a
  genuine de Sitter result reads LOW and still-climbing at small N4/short therm —
  indistinguishable early from "not de Sitter." Confirmation needs volume +
  equilibration, overnight+, not hours.

START FRESH (no Regge-equilibrated resume). Reason: resuming from a Regge-shaped
checkpoint means the first hours are a TRANSIENT decay of the old GR-built de
Sitter state, not the EPRL equilibrium. An early "looks de Sitter" from a Regge
seed is the LEAST trustworthy reading — a false positive that fades. A de Sitter
result from a fresh (non-Regge) start is the one with no asterisk, because it
can't be the starting config talking.

Pre-committed: **the early "still looks like de Sitter" reading is the transient,
NOT the verdict.** Do not bank it. The verdict is the thermalized, large-N d_s
flow rising toward 4 with the peak marching outward — same bar as the β=0 scaling
test. Early reads can only decisively show the NEGATIVE (collapse).

## What to actually run / read

1. Launch centered eprl_only at matched volume, long thermalization. Watch that it
   stays a valid manifold (links/valid ok) and N41 holds — if not, it's a setup
   failure, not a result.
2. Read the d_s(σ) FLOW (ds_scaling_check-style) + cos³ blob, after thermalization.
   d_H alone is insufficient (and is anyway nearly blind to this term).
3. Apply the table above. Do not soften an entropic collapse into a null.
