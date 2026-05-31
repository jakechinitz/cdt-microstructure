# Pre-registration — matched-volume EPRL β-sweep (frozen-j=3)

Written **before** reading the results of the `eps=0.01` matched-volume sweep, so
the interpretation is fixed in advance and cannot be reverse-engineered to fit
whatever comes out. Committed to git for an immutable timestamp.

## What is running

Four runs, each resumed from the same thermalized bare 20k checkpoint
(`old/scan_20000.json`), identical except for the EPRL coupling β:

| run | β (EPRL coupling) | role |
|---|---|---|
| `thy_b0.0_20k`  | 0.0  | **control** (pure Regge through the theory harness) |
| `thy_b0.05_20k` | 0.05 | coupling on, weak |
| `thy_b0.1_20k`  | 0.1  | coupling on, medium |
| `thy_b0.3_20k`  | 0.3  | coupling on, strong |

Fixed for all: centering ON (mu≈30.54), k4=0.9, eps=0.01, K=80, N41 target=20000.
Bare reference at this volume: **d_H ≈ 3.37**, valid manifold, de Sitter blob.

## Gates that must pass BEFORE any result is read

If any gate fails, the output is `NOT READABLE` — fix the setup and rerun, do not
interpret d_H.

- **G1 — Validity.** `valid = ok` (and links `ok`/`gen`, never `BAD`) on all four.
- **G2 — Thermalization.** Action S / N4 has stopped trending on each run
  (flat, not still drifting) over the measured window.
- **G3 — Matched volume.** N41 pinned near 20,000 AND N4 comparable across all
  four β (within a few %). This is the gate the two prior runs failed. Without
  it, d_H differences are the size→dimension artifact, not the theory.
- **G4 — Control faithful.** β=0.0 reproduces the bare value, d_H ≈ 3.37
  (± noise). If the control does not land near bare, the harness itself perturbs
  geometry and NO β>0 point is trustworthy.
- **G5 — Chain not frozen.** Acceptance non-trivial; N4/d_H actually move
  sweep-to-sweep (β=0.3 most at risk). A frozen chain is not a result.

## Observables, in priority order

1. **PRIMARY — de Sitter volume profile.** `blob` score, `active` slices,
   `cos3err`. The cos³ blob is the canonical CDT phase-C signature and the most
   trustworthy shape observable here.
2. **SECONDARY — d_H** (Hausdorff, shell-growth). Trusted, but a single number;
   read alongside the profile, not alone.
3. **CROSS-CHECK ONLY — d_s** (`ds_crosscheck.py`). Known to under-read on CDT
   dual graphs (see `VERTEX_PROVENANCE.md` §5 / `ds_crosscheck.py` header).
   Used only to confirm direction-of-motion agrees with d_H/blob; never a verdict.

## Pre-registered interpretation (success = FLAT OR BETTER)

The premise: the EPRL amplitude should be **compatible with, or support**,
physical 4D de Sitter geometry. So at matched volume, thermalized, with the
control passing (G4), as β increases:

| Outcome across β | Reading |
|---|---|
| d_H / blob **hold near bare** (d_H≈3.37, blob persists) | **SUCCESS — theory survives.** Amplitude is compatible with 4D; coupling it on does not destroy the de Sitter geometry. (The expected, hoped-for result.) |
| d_H / blob **climb** (d_H↑ toward 4, blob sharpens) | **STRONG SUCCESS.** Amplitude actively supports 4D. Would be a surprise at frozen-j. |
| d_H / blob **degrade** (d_H↓ toward ~2, blob dissolves), monotonic in β | **CONCERNING — amplitude fights 4D.** Investigate; do NOT pre-excuse it. |

### Strength of claim — bounded honestly

- A **flat/better** result is a real, defensible statement: *"the frozen-j=3 EPRL
  amplitude, coupled at matched volume, does not destabilize the emergent 4D de
  Sitter geometry of the bare CDT substrate."* It is a **compatibility /
  survival** result, NOT a proof the theory predicts 4D.
- A **degrading** result is provisional: because this is the frozen-j=3
  truncation (not peaked-j), a clean negative would motivate building the
  peaked-j amplitude before any claim of falsification. Frozen-j may not
  represent the full amplitude. (This caveat applies to a NEGATIVE only — it is
  not a license to expect or dismiss a negative in advance.)
- Either way this run **cannot, by itself, falsify or confirm the theory.** It is
  a controlled compatibility check on the frozen-j placeholder. The headline
  physics result remains the bare-Regge 4D emergence (the substrate test).

## Already-settled supporting facts (logged, not pending)

- EPRL term is **not inert**: centered per-pentachoron σ≈0.87 ("HAS TEETH"),
  variance largely orthogonal to the (4,1)/(3,2) geometric DOF (type-blind).
  → a flat sweep is a genuine finding, not a dead-term artifact.
  (`eprl_term_diagnostic.py`, reproduced on the VM.)
- Amplitude provenance and its validation scope are recorded in
  `VERTEX_PROVENANCE.md` (small-spin gate passed; asymptotic gate trusted-not-run;
  unpinned sl2cfoam-next version is a known reproducibility gap).
