# Runbook: the remaining production runs

Audited 2026-07-04 (code paths checked end-to-end; analyzers verified on
synthetic data with known answers). All commands run from `v6_vm_bundle/`.
`BASE` below = the clean fixed-engine bare checkpoint
(`scan_20000_causal.json`); closure-thermalized checkpoints where noted.

**Timing calibration:** every engine log prints a `wall_s` column. After
~20 minutes of any run, `sweeps/hour = (rows so far) x measure_every /
(wall_s/3600)` extrapolates the finish time exactly. Estimates below are
anchored to the earlier PC session's pace at N41=20000.

---

## 1. Stage-4 v1.2 (T3/T4 gates) — ~6–12 h — SUPERSEDES the v1.1 run

The v1.1 run's T4 "screening" verdict (xi ~ 1.6) measured TWO artifacts,
both now fixed (see the addendum in CAPACITY_CONSERVATION.md):

- **The rectifier (the big one).** Vacuum n_coll is 92% zero / 8% one, so
  `max(0, n_coll - nbar)` removes only ~8% of vacuum absorption: the v1.1
  "excess" mode had intrinsic xi = sqrt(D/(kappa*0.071)) ~ 7.5 BY
  CONSTRUCTION. Fix: absorption is driven by the PERSISTENT failure
  average (`--persist`, EWMA, default 200 sweeps; rectified vacuum rate
  falls ~1/sqrt(W), pushing intrinsic xi to ~23) — theory-faithful, since
  the paper defines mass as persistent closure failure and the vacuum's
  flicker turnover is already the anchor budget, not commitment.
- **Deficit reference.** Deficits are now measured against the same-slice
  far field (new `far_f` CSV column), removing the closed-slice zero-mode
  surplus that steepened the v1.1 fit (and explains measured 1.6 < 7.5);
  the analyzer prints the level-0 dressing profile with sign (v1.1's
  d=1 elevation to 7.76 is a real under-absorption surplus, expected to
  collapse under v1.2). Verified on synthetic data with both
  contaminations planted: true 1/d recovered exactly, gate PASS.

```
python v6_capacity_run.py --resume clo_b1.0_20k_K40.json \
    --absorb excess --persist 200 --pins-per-level 25 --sweeps 4000 \
    --out cap_v12 > logs/cap_v12.log 2>&1
python v6_capacity_run.py --analyze cap_v12
```

Run it on a **K=40 base** (thicker slices, usable radius ~7-8 vs ~5): grow
`scan_20000_causal_K40.json` bare, thermalize a beta=1 closure arm on it,
resume from that. With intrinsic xi ~ 23 >> radius ~ 8, a screened verdict
at xi ~ radius is now REAL physics against Route B, and a flat
ln(deficit*d) out to the radius is an honest T4 PASS. v1.1's M1/M2 remain
claim-grade (near-field sink quantities, insensitive to the far-field
artifacts), though M2's absolute value will shift under v1.2.

Gates: T1 conservation exact (run aborts otherwise); T2 <f> at 7.4198;
M1 dressing-subtracted proportionality; M2 junction constancy; T3 ~ d^-1;
T4 massless over the measurable range.

## 2. Stage-3 replication seeds — ~4–8 h per seed-pair

Per seed S in {1, 2, 3}: two arms, SAME resume, out/checkpoint names carry
seed+arm (the checkpoint default now derives from --out, so concurrent
arms can no longer clobber each other):

```
python v6_defect_run.py --resume clo_b1.0_20k.json --pin-mode fail \
    --seed S --out def_s{S}_real   > logs/def_s{S}_real.log 2>&1 &
python v6_defect_run.py --resume clo_b1.0_20k.json --pin-mode closed \
    --seed S --out def_s{S}_placebo > logs/def_s{S}_placebo.log 2>&1 &
```

Both arms of a seed share pin locations (pin RNG = seed+31 on the same
geometry): the real-minus-placebo difference is anchored identically.
Analyze pooled: `python v6_defect_run.py --analyze def_s*_real.csv
def_s*_placebo.csv`. Hard gate per arm: "pins alive N/N" (an arm printing
PIN LOSS is void). Error bars in analyze ignore autocorrelation — for the
paper, quote across-seed scatter, not within-run SEM.

## 3. Bridge scaling family — demo volume first (~2–4 h), production later

```
BASE=<clean demo checkpoint, N41=2000> EPSLIST="0.005 0.01 0.02" K=16 \
    ./run_bridge_sweep.sh "0.25 0.5 1.0" 2000
# when arms converge:
python bridge_predict.py --analyze logs/bridge_b*_2k.log \
    --pred results/bridge_pred_2k.csv
```

Predictions are written to the CSV before any arm launches. The analyzer
prints two prediction columns: `pre-reg` (independent TI, the
pre-registration) and `exact` (-beta*mu_arm/(4*eps) from the mu the
centered arm actually subtracted, read from its log — the pair difference
isolates exactly that constant, so `ratio` uses it and predictor TI noise
cannot masquerade as a family failure). Gates: ratio ~ 1 everywhere;
concave beta-curve (S1); 1/eps collapse (S2); placebo on its own mu_plc
line, 28% from the real line (S3). Production repeat at 20000 with
EPSLIST="5e-4 1e-3 2e-3" once the shape holds (~1–2 days, 20 arms).

## 4. Autotuned phase grid — ~6–12 h at --jobs 4 (may already be running)

```
python phase_grid.py --tune-k4 3 --target-n41 10000 --K 80 \
    --sweeps 400 --jobs 4 --out grid10k_tuned
python phase_grid.py --plot grid10k_tuned.csv
```

CSV is append-safe (resume by re-running; finished points are skipped).
Read `pin` (residual) and the recorded post-tune `k4` per point; a point
with large |pin| after tuning is still not a fair test. Success = a region
with blob >> 1, active >= 3, d_H rising toward 4, foliation CLEAN.

## 5. No-action floor row — ~2–4 h

```
python v6_closure_run.py --local-closure --beta-closure 0.0 \
    --resume scan_20000_causal.json --checkpoint results/floor_b0_20k.json \
    --target-n41 20000 --K 80 --eps 1e-3 --max-sweeps 4000 \
    > logs/floor_b0_20k.log 2>&1
```

The beta=0 fixed-engine control row for the comparison tables (exact-bare
by construction: uniform births are detailed-balanced at beta=0).

## 6. Rung-2 curvature displacement — script READY; launch after the
##    volume-pair verdict confirms the boundary region

```
python rung2_boundary.py --predict            # c0 predictions on record
python rung2_boundary.py --k0-list "2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4" \
    --Delta 0.4 --seeds 3 --sweeps 800 --jobs 4 --out rung2.csv
python rung2_boundary.py --analyze rung2.csv
```

The test: in bare-k0 coordinates the coupled ensemble's boundary must
shift by +c0*r (c0_real = +0.0187/slice-edge, r = N1/N0 measured from the
runs, ~ +0.12 at r~6.5) while the PLACEBO arm's boundary must not move
(c0_plc = +0.0003, seventy times smaller despite a LARGER volume-sector
coefficient — the shuffle destroys structure, and only structure prices
curvature). Structure-vs-magnitude, where the bridge family tested
magnitude. Per (k0, seed): k4 tuned on bare bursts, then bare/real/placebo
arms share the tuned checkpoint and k4; centering removes the volume tilt,
so only the curvature sector can move the crossing.

Verified on synthetic scans: planted shifts +0.12/+0.002 recovered at
+0.107/+0.007 with 3 seeds (resolution ~0.015 vs the 0.12 signal, ~8:1).
Cost: 8 k0 x 3 seeds x (tune + 3x800 sweeps) at 10k volume ~ 2-3x the
tuned grid, 1-2 days at --jobs 4. Gates: real shift at prediction, placebo
consistent with zero, foliation CLEAN, shared tuned k4 per point.
CSV is append-safe (rerun to resume).

---

Suggested order: check #4's CSV first (it may be done); then #1 and the
first seed of #2 together; #3 demo family alongside (cheap); remaining
seeds + #5 as slots free; #6 when the boundary is known; #3 production
last. Total ~3–5 days of PC time with 3–4 concurrent processes; the
paper-critical results (#1 gates, first #2 seed, #3 family shape) land in
the first ~24 h.
