# Pre-registration — the three-run ladder (floor / treatment / ceiling)

Written **before** the matched-volume read. This supersedes the scattered runs
(regge sweep, bare 40k, tiny eprl_only pilots) with the minimal sufficient set to
extract what frozen-j=3 can say. Committed for timestamp.

## The reframe that makes j=3 answerable

The paper does NOT fix the absolute amplitude scale — it specifies the *structure*
(seven states, the counting, closure). So asking j=3 to settle a *generation*
claim (which needs the scale) was the wrong question. The right, scale-free
question j=3 CAN answer: **is the seven-state structure geometrically COMPATIBLE
with 4D — does it permit, favor, or fight de Sitter, relative to no structure?**

That is read as the treatment's position on a floor→ceiling ladder.

## Floor is NOT zero (measured, refutes the "nothing emerges" intuition)

No-action (action=0, CDT moves + eps volume-fix only), grown to N4≈11k–17k:
```
d_H ≈ 3.19,  blob ≈ 2.0–3.9,  d_s_peak ≈ 2.81 @ σ≈22,  valid manifold,  N41 held
```
The bare scaffold produces a valid manifold with a blob and climbing d_H/d_s with
NO action at all — the v4 near-Ramanujan-expander resurfacing. So "eprl_only
emerged a manifold/blob" is WORTHLESS without this floor to subtract: the scaffold
emerges that on its own. This is the control the analysis kept demanding.

## CORRECTION to a prior claim (must not stand uncorrected)

`ds_scaling_check.py` calls propose_and_apply with NO accept/reject — it is
accept-all, i.e. it IS the no-action floor. So the earlier "substrate validated:
d_s_peak climbs 2.36→2.65, 4D-trending" was measuring the EXPANDER FLOOR, not
Regge-driven de Sitter. That d_s rise is the scaffold's, not evidence the
substrate's *dynamics* produce 4D. The substrate-validation claim is hereby
downgraded: what was shown is that the floor (uniform measure over CDT-legal
volume-fixed triangulations) has d_s climbing with N4 — which is exactly why the
floor is a non-trivial baseline, not a vindication of any action.

## The three runs — identical except the action

Same moves, same eps=1e-2, same K, same target N41, same seed-handling. Only the
action varies. All grown to the SAME volume (~20k N41-scale, i.e. N4 ~ 90–100k) —
matched volume is non-negotiable (it killed the centered/uncentered comparison)
and d_s needs N4 ≳ 15–20k to separate 4D from small-N ambiguity.

1. **FLOOR — no action.** geometry_action = 0. Uniform measure baseline.
2. **TREATMENT — centered eprl_only.** Centered: for a STRUCTURAL test you want the
   untrustworthy absolute scale removed, leaving the convention-independent
   structure the paper actually specifies.
3. **CEILING — Regge.** Known to build 4D de Sitter. "What success looks like."
   Reuse a matched-volume β=0 / bare checkpoint if available rather than rerun.

## Pre-registered read (decided before seeing it)

Observable: d_s flow at matched volume; d_H/blob secondary. Treatment's position
on the floor→ceiling ladder is the answer:

| Treatment vs floor/ceiling | Reading |
|---|---|
| eprl_only ≈ floor (no-action) | j=3 structure is **INERT** — geometry is the scaffold's doing, amplitude rode along. Modest, honest. |
| eprl_only between floor and ceiling, toward Regge | structure does **GR-like geometric work** — structural compatibility with (and mild favoring of) 4D. The interesting positive. |
| eprl_only **worse** than floor | structure **FIGHTS** 4D — a real negative for the ontology. |

## What this can and cannot claim (writeup discipline)

- CAN claim: **structural compatibility** of the j=3 seven-state structure with 4D
  geometry, measured as treatment-vs-floor at matched volume.
- CANNOT claim: **generation**. That needs the absolute scale the paper rightly
  doesn't fix → needs the spin-sum → needs a cluster. Defer explicitly as future
  work. A positive on the ladder is "the structure is compatible with / favors 4D,"
  NOT "generates spacetime."
- STANDING LIMIT: even outcome-2 is "EPRL structure SELECTS among CDT-PROPOSED
  geometries" — the move-set is still CDT's. Not "substrate dynamics from scratch."
  State this so the claim isn't broader than the run.

## Cost / timing

Three runs grown to N4~90–100k = the overnight window, not hours. Collapse (if
treatment falls below floor) shows fast; compatibility (treatment ≥ floor, d_s
climbing) needs the volume + thermalization, same bar as the floor's own d_s scan.
