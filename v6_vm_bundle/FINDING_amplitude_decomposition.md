# Finding: the trustworthy / geometry-driving split in frozen-j=3 EPRL

This is a structural result about whether the frozen-j=3 EPRL amplitude can be the
spacetime-generation mechanism the theory needs. It was derived by analysis (not
by any overnight run), and it is reported here as a finding in its own right —
arguably the clearest thing this investigation has produced — rather than buried
as a caveat about flag settings. Credit: surfaced by adversarial review pushing
back on a "let the comparison decide" framing that wrongly treated two unlike
confounds as symmetric.

## The decomposition

Write the per-pentachoron EPRL contribution as `-log|A_v| = mu + delta_p`, where
`mu` is the configuration-mean (~30.5 for vertex_j3.npz) and `delta_p` is the
fluctuation. The EPRL action over a triangulation then splits into two channels:

1. **mu*N4 — the geometry-DRIVING channel.** Depends on the volume; with the
   per-pentachoron scale ~30 it is by far the larger piece and is what can exert
   real force on geometry. BUT this absolute scale `mu` is **convention-dependent
   and artefact-grade at frozen j** (see VERTEX_PROVENANCE.md): the overall
   normalization of the j=3 vertex is not a trustworthy physical quantity until
   the amplitude is summed over spins properly.

2. **Sum of delta_p — the TRUSTWORTHY, convention-independent channel.** The
   fluctuations about the mean are robust to normalization. BUT they were
   measured to be **nearly geometry-blind**: type-split 0.054 << sigma 0.87,
   dual-graph r_nn ~ 0.14 (weak). This is the part `--center-eprl` keeps.

## Why this matters more than any single run

The generation claim needs *some* part of the amplitude to drive geometry toward
4D in a trustworthy way. The decomposition says the two requirements live in
**different parts of the amplitude**:

  * the part with enough teeth to drive geometry (mu*N4) is the part you can't
    trust at frozen j;
  * the part you can trust (the fluctuations) is nearly geometry-blind.

If that split is correct, then **neither channel can support the generation claim
at frozen j=3**: the driving channel is untrustworthy, the trustworthy channel
can't drive. That is a precise statement of why frozen-j=3 EPRL is structurally
unable to be the generation mechanism — independent of centered-vs-uncentered,
independent of the overnight runs, which can only corroborate it.

## The honest consequence: the clean test is not currently runnable

`--center-eprl` and `--no-center-eprl` at frozen j are both PROXIES that fail in
complementary ways:
  * centered  = trustworthy channel only -> tests whether geometry-blind
    fluctuations can build shape (epistemically weak: a null is near-predetermined
    *in principle*, though see the empirical caveat below);
  * uncentered = driving channel included -> any 4D it produces is driven by the
    untrustworthy absolute scale, so a positive is not a trustworthy positive.

The genuinely clean generation test is **EPRL summed over spins j** (where the
absolute scale becomes physical), which is a different and much larger computation
that does not currently exist (the multi-j / peaked-j tensors were never built).
So the accurate status is: *the clean frozen-j-free generation test is not
runnable today; both frozen-j runs are complementary proxies, neither decisive.*

## Empirical head-to-head (small scale, both modes, same target N41=2500)

```
centered   eprl_only:  N41 held ~2550;  d_H 1.74->2.31;  blob ~2.2;  valid throughout
uncentered eprl_only:  N41 held ~1020;  d_H 1.36->1.97;  blob ~2.8->3.55; valid throughout
```

Three results:

1. **Volume-pin prediction CONFIRMED.** Uncentered's mu*N4 term relocates the
   volume pin DOWN by ~mu/(2*eps) ~ 1400 simplices: target 2500 - 1400 = ~1100,
   observed ~1020. So the driving channel acts as a cosmological constant that
   *shifts where volume sits*, NOT a runaway. (Reviewer's "volume wanders" refuted;
   the eps penalty regulates in both modes, just balanced at different N4.)

2. **Neither mode collapsed.** The strong "geometry-blind => entropic collapse by
   construction" prediction FAILED for BOTH centered and uncentered: both held
   valid manifolds with persistent blobs (2.2 and 3.5) and rising d_H. At this
   scale neither relaxes to a structureless stalk. So "centered is degenerate by
   construction" is not borne out, and the decomposition must NOT be read as
   predicting collapse — it is about trust+drive, not about the dynamics dying.

3. **The head-to-head is itself volume-confounded** (the irony): because the two
   modes pinned at different N4 (~2550 vs ~1020), uncentered's higher blob /
   lower d_H are partly just small-volume effects, not a clean shape difference.
   So this comparison CANNOT be read as "uncentered builds a better blob." To
   compare shape cleanly you must match N4 across modes (raise uncentered's target
   by ~1400, or retune k4), which has not been done.

Caveat on all of the above: N4 ~ 1000-3800 is tiny and finite-size dominated;
read trends, not absolute d_H. None of this overturns the analytic decomposition;
it sharpens it: the driving channel demonstrably moves volume (confirming it has
teeth on N4), the trustworthy channel demonstrably does not kill the geometry
(both hold blobs), and a clean shape test remains blocked by volume-matching plus
the frozen-j trust problem.

## What to do (the meta-point)

Stop launching variants; read what is in flight. The decomposition is the result
and it needed no overnight run to establish — the runs corroborate, they cannot
produce anything more decisive than the analysis already has. Specifically:
  * keep the centered fresh run going (informative: the convention-independent
    channel, and empirically not collapsing);
  * do NOT launch uncentered as "the fix" — it is not one;
  * when the uncentered head-to-head lands, read it as CONFIRMATION of the
    decomposition (collapse, or 4D-driven-by-untrustworthy-scale), not as a
    tiebreaker between two viable theories;
  * put this decomposition in the writeup as a primary finding.
