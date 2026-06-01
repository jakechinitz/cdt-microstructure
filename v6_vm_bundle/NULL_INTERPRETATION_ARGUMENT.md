# The r_nn argument — does a flat blob mean the GOOD null or the EMPTY null?

Written **before** the matched-volume blob delta is read. This is the
load-bearing reasoning the scripts cannot perform: `blob_compare.py` will *remind*
you to make this argument, but making it is intellectual work, and it may come out
against us. Committed for an immutable timestamp so the conclusion can't be
back-fit to whichever null flatters the result.

## The two nulls (only one is publishable as a physics finding)

If the blob is flat at β=0.3 vs β=0 at matched volume:

- **GOOD null** — *"frozen-j=3 EPRL intertwiner fluctuations are largely
  orthogonal to the de Sitter geometry: the spin labels and the macroscopic shape
  decouple."* This is a finding about quantum gravity ONLY IF the term had a real
  chance to move the geometry and didn't.
- **EMPTY null** — *"we coupled to a quantity nearly orthogonal to the observable,
  so nothing moved."* This is a finding about our variable choice, not physics.

The diagnostic places us uncomfortably near the empty null: the centered
amplitude's variance is σ≈0.87 but the (4,1)/(3,2) type-split is only 0.054
(split/σ≈0.06). What rescues us toward the GOOD null — IF it holds — is the
dual-graph nearest-neighbour correlation **r_nn≈0.14 (~4σ vs shuffled null)**: a
real, geometrically-structured channel that is not type-aligned.

## The argument we OWE (and its falsifiable threshold)

Claim to defend: *the r_nn≈0.14 channel, amplified by β over ~N₄ pentachora,
constituted a real opportunity for the geometry to move — so a flat blob is
SURPRISING (informative), not predetermined.*

To make this non-circular, fix the threshold BEFORE the read:

1. **Channel strength is real, not noise.** r_nn=0.14 is ~4σ above the shuffled
   null (null std 0.035). PASS — established.

2. **The channel is dynamically non-negligible, not just statistically nonzero.**
   This is the hard part and the one that can fail. A 14% nearest-neighbour
   correlation means the amplitude field has coherent patches on the dual graph;
   coupled at β=0.3 with per-pentachoron action scale β·σ≈0.26 (O(1) on the
   Metropolis scale), over ~10⁵ pentachora, the *collective* bias on the
   geometry is not obviously sub-threshold. The honest test: does β·σ·(coherence
   length set by r_nn) reach the scale of a single accepted Pachner move's action
   change? If yes, the term could have steered and a flat result is informative.

3. **Pre-registered quantitative bar.** We call the GOOD null only if BOTH hold:
   (a) the β=0.3 per-move action contribution from the EPRL term, measured on the
   thermalized config, is ≥ O(1) for a non-negligible fraction of attempted moves
   (i.e. the term actually entered accept/reject decisions, not washed out); AND
   (b) the blob is flat within the gated tolerance anyway.
   If (a) FAILS — the term's per-move contribution was sub-threshold on essentially
   all moves — then we are in the EMPTY null: the result is "our coupling was too
   weak/orthogonal to matter," a statement about setup, and the finding is smaller.

## UPDATE (smoke-test surprise): the term is LOUD, not weak — reframes the empty null

A quick `per_move_eprl_scale.py` run (small N4=672, thin — qualitative only)
found |β·dS_EPRL| median ≈ 5 vs |dS_Regge| ≈ 1.8 at β=0.3. So the term is NOT
washed out per-move — it is *louder* than gravity in accept/reject. This splits
the analysis cleanly and changes criterion 3:

- The empty null is NOT "coupling too weak." The term is loud.
- The live distinction is now **loud-and-steering (GOOD)** vs **loud-but-
  orthogonal (EMPTY)**: a term can dominate the action yet, with type-split/σ≈0.06,
  push almost entirely along label directions the blob can't see. Loud noise is
  still noise. So criterion 3a ("did the term enter decisions?") is necessary but
  NOT sufficient — orthogonality (criterion i) is the real gate.
- NEW concern (prereg gate G5): a loud term *suppresses acceptance*
  (~e^{-|β·dS|}). Because Pachner moves are local (fixed footprint), this need not
  ease at large N4. β=0.3 may be a partially FROZEN chain — which is not a result.
  Before trusting any flat β=0.3 read, confirm from the live log that N4/blob are
  still evolving (not stuck). If β=0.3 is frozen, the usable comparison may be
  β=0 vs a *lower* β (0.05/0.1) where the chain still moves.

## This may come out against us — and that's fine, stated honestly

It is entirely possible that working through (2)/(3) we conclude the channel was
too weak: that β=0.3 was not enough to give r_nn=0.14 real dynamical leverage. In
that case the honest writeup is the EMPTY null, explicitly: *"at frozen-j=3 and
β≤0.3, the EPRL term's coupling to the de Sitter-relevant geometry was below the
threshold to steer it; we cannot distinguish decoupling-in-principle from
too-weak-in-this-setup."* That is a smaller, but still honest, result — and it
points to the fix (higher β until per-move contribution is O(1), or a
geometrically-aligned observable), not to a claim about quantum gravity.

The discipline: which null we are in is decided by the **measured per-move EPRL
action contribution at β=0.3** (criterion 3a), NOT by how much we'd like the
result to be a physics finding. Measure 3a when the run is thermalized; it is a
diagnostic the current logs do not yet contain.

## What to actually do at the read

1. Run `blob_compare.py` — gates first, then the blob delta.
2. If gates PASS and blob is FLAT: measure criterion 3a (per-move EPRL action
   contribution distribution at β=0.3 on the thermalized config). This decides
   GOOD vs EMPTY null. It is the real intellectual work; the blob delta was the
   easy part.
3. If gates PASS and blob MOVED: relaunch β=0.05, 0.1 to confirm smooth/monotonic
   response before believing it (a lone 0→0.3 jump is not a signal).
4. If gates say NOT READABLE: more thermalization / matched volume. No workaround.

## Hierarchy reminder (do not relead with the satisfying number)

bare 40k climbing = instrument validation (known AJL/Kommu phenomenology), NOT
the headline. β=0.3-vs-β=0 blob = the experiment. When the bare climb looks great
and the sweep is flat, the pull will be to lead with the climb. Resist it: the
climb confirms the engine works; the flat sweep (read through the null argument
above) is the part that is actually ours.
