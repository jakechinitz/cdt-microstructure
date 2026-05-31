# Pre-registration — matched-volume EPRL β-sweep (frozen-j=3)

Written **before** reading the results of the `eps=0.01` matched-volume sweep, so
the interpretation is fixed in advance and cannot be reverse-engineered to fit
whatever comes out. Committed to git for an immutable timestamp.

## Evidence hierarchy (do not relabel the substrate as the headline)

The bare-Regge run (β=0, and the separate bare 40k continuation) reproduces known
AJL/Kommu CDT phenomenology: d_H climbing toward ~3.5–4 on a pure-Regge scaffold.
That is a **validation checkpoint that the CDT instrument works** — necessary,
expected, NOT novel, and it contains **no theory**. It cannot be evidence for
"the theory generates 4D from the graph," because at β=0 nothing of the theory is
in the run.

The novel claim — *your* EPRL microstructure doing the generating/steering — lives
ONLY in the β>0 sweep, which the diagnostics show is weakly geometry-coupled and
whose most probable honest outcome is a decoupling null. So the hierarchy is:

  * Bare 40k climb / β=0  → "the substrate/instrument works." Figure 1. Bank as
    validation, never as the headline.
  * β=0.3 blob vs β=0 blob at matched volume → THE experiment. The one number the
    sweep exists to produce; most likely a decoupling null, which is a real
    finding stated as such — and stops being one if dressed as "graph makes 4D."

## Mechanical note on the bare 40k continuation (read first rows with care)

The 40k bare run was resumed from `old/scan_40000.json`, which holds N4≈190,498 —
an EARLIER/smaller checkpoint than the N4≈223k at which the d_H≈3.51 reading
occurred. So this is NOT a seamless continuation of that climb: with
`extra_state absent` (thermalization not carried over) and a smaller N4, the
first logged rows will read LOWER than 3.51 and must re-equilibrate/re-grow back
through ~223k before d_H is comparable. Do not read the opening dip as a result.
(A nearer-continuous seed exists in the contaminated-era `results/scan_40000.json`
at N4≈223k — valid geometry for a bare continuation despite its double-writer
history — if a true continuation is wanted later.)

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

## AMENDMENT (post-diagnostic, pre-results): a flat d_H is WEAK evidence

Before any sweep numbers, the term-teeth diagnostic was extended (see
`eprl_term_diagnostic.py`) and found: the centered EPRL variance is real
(σ≈0.87) but **largely orthogonal to the geometry** — type-split 0.054
(split/σ≈0.06) and dual-graph nearest-neighbour correlation r_nn≈+0.14 (weak,
~4σ). This was verified NOT to be a centering artifact (raw and centered
type-splits identical to 4e-15). Consequence: **a flat d_H across β is now
expected almost regardless of the theory**, so flat d_H carries little
information and must NOT be banked as "theory tolerates 4D." Two corrections:
- The decisive observable is the **cos³ BLOB / profile**, not d_H — does the
  weak-but-real geometric channel move the de Sitter shape as β→0.3?
- d_s is NOT a usable arbiter: a measured β=0 d_s(σ) flow peaks at ~2.4 on a
  known-good bare config (the documented CDT-dual-graph under-read), so d_s would
  falsely condemn the substrate. Do not gate on d_s.

## Pre-registered interpretation (success = FLAT OR BETTER, read via BLOB)

The premise: the EPRL amplitude should be **compatible with, or support**,
physical 4D de Sitter geometry. At matched volume, thermalized, control passing
(G4), as β increases — judged PRIMARILY on the cos³ blob/profile, with d_H
secondary and read in light of the weak-coupling amendment above:

| Outcome across β | Reading |
|---|---|
| **blob persists** (and d_H holds ~3.37) | **SURVIVAL, weakly informative.** Consistent with "amplitude doesn't destroy 4D," BUT given the measured orthogonality this is close to the expected-anyway outcome; state it as a decoupling result, not a tolerance proof. |
| blob **sharpens** / d_H climbs | **STRONG SUCCESS.** The geometric channel actively supports 4D — would be notable given how weak the coupling is. |
| blob **dissolves** / d_H↓ toward ~2, monotonic in β | **CONCERNING — amplitude fights 4D.** Investigate; do NOT pre-excuse. |
| nothing moves at all, even β=0.3 | **DECOUPLING NULL (the likely result).** "Frozen-j=3 EPRL intertwiner fluctuations are orthogonal to the de Sitter geometry; spin labels and macroscopic shape decouple in this discretization." A real, publishable null — a statement about variable choice + frozen-j, not about QG endorsing de Sitter. |

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
