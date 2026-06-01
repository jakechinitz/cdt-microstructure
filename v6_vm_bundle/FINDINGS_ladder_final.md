# Findings — eprl_only ladder (final read)

Written after the ladder runs hit their natural limits. This is the honest
accounting of what the frozen-j=3 EPRL simulation can and cannot establish. None
of it reaches "the substrate generates 4D spacetime"; the readable results are a
viability fact, a debunked artifact, and three instrument-limit findings. The
paper does not depend on any of this (it stands on the static weak-field chain);
this was a falsification side-probe of one auxiliary dynamical claim.

## What IS established (solid)

1. **Viability.** Under pure-amplitude weighting (eprl_only, centered, no Regge),
   the geometry holds a valid 4-manifold and keeps N41 pinned at target for
   thousands of sweeps — it does not poof, freeze-to-minimum, or go non-manifold
   on its own. Basic geometric viability of the frozen-j=3 amplitude as a CDT
   weighting. (Independent of all matching/equilibration gates.)

2. **The "sharper de Sitter" blob is an ARTIFACT, not a prediction.** EPRL's
   blob_score ~5.0 looked like a sharper de Sitter profile than Regge (~2.8). It
   is not: blob = max(V3)/mean(V3) is a peakedness ratio, and EPRL's cos3_relerr
   is HIGH (~0.84, i.e. the profile does NOT fit cos^3) while the d_s(sigma) flow
   at matched N4~84k is IDENTICAL to Regge to +/-0.02 across all scales (rise
   +1.61 vs +1.59). Two independent channels (cos^3 fit quality; direct d_s flow)
   agree the high blob is a non-de-Sitter profile feature, not sharper de Sitter.
   Do not publish blob ~5.0 as a de Sitter result.

## Instrument-limit findings (real, but about the apparatus, not the ontology)

3. **Floor and ceiling cannot be volume-and-equilibration matched (Gate 0 fail).**
   N41-pinning lets each action settle at its own N4 (no_action -> N4/N41 ~12+ via
   N32 entropy; Regge -> ~4 via its k4 term). Their equilibrated N4 ranges do not
   overlap (floor ~230-257k; Regge ~67-95k), so the crisp three-rung d_s
   comparison at matched volume is not constructible on this machine with this
   volume-fixing scheme.

4. **eprl_only does not equilibrate volume at feasible runtime.** EPRL's N4 grew
   monotonically with NO deceleration over 2575 sweeps (rate steady ~500/25sw to
   the end; ratio climbed 2.9->4.16->4.69->5.01 straight through Regge's value with
   no ceiling). The loud frozen-j term suppresses the volume-reducing moves enough
   that the chain cannot reach a stationary N4 in available time. So the treatment
   rung could not be brought to a readable equilibrium; "EPRL ~ Regge in the ratio"
   was a transient waypoint, NOT an equilibrium, and is explicitly NOT banked.

5. **Regge ceiling lost manifold validity at large N4.** The Regge run verified OK
   for ~8650 sweeps then went BAD (S^3 vertex-link failures; gluing still passes)
   around N4~93.6k — the same large-N pseudo-manifold drift seen earlier in the
   10k bare run. Ceiling d_s data before sweep 8650 is usable; after is not.

## Bottom line

The frozen-j=3 EPRL simulation, on this hardware with N41-pinned CDT moves, cannot
deliver the matched-volume equilibrated d_s comparison the compatibility test
requires: the rungs don't co-equilibrate (3), the treatment doesn't equilibrate at
all (4), and the ceiling degrades at large N (5). What survives is viability (1)
and the debunking of the blob artifact (2). The compatibility question ("is the
j=3 structure CDT-like or expander-like") is therefore left UNTESTED by a clean
comparison — not answered. The generation question is further out still and needs
the spin-sum (a cluster-scale computation that does not exist here). These are
honest negatives about what is measurable; they are not evidence against the
ontology, and the paper's static-chain results are untouched.

## Why not just run longer

(4) is not a "wait longer" — the growth rate showed zero deceleration over
thousands of sweeps, so there is no evidence more runtime reaches equilibrium; the
honest statement is non-equilibration at feasible runtime. (3) is structural (the
actions' divergent volume preferences) and longer runs make it worse, not better.
(5) would need a smaller step / link-repair, a separate engineering effort. The
decisive next step for the actual physics is the spin-sum amplitude, not more of
this run.
