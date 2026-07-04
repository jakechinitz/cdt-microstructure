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

## 1. Stage-4 v1.1 (T3/T4 gates) — ~6–12 h

```
python v6_capacity_run.py --resume clo_b1.0_20k.json \
    --absorb excess --pins-per-level 25 --sweeps 4000 \
    --out cap_v11 > logs/cap_v11.log 2>&1
python v6_capacity_run.py --analyze cap_v11
```

Notes from the audit (both fixed in code, defaults are now safe):
- In `--absorb excess` mode absorption is now OFF during therm (field held
  exactly at the vacuum anchor) and measurements start only after a field
  burn-in (`--post-therm`, default = therm). Without this, the screened
  therm-transient contaminated the profile and biased T4 toward a FALSE
  screening verdict.
- The analyzer now computes T3 (power fit; massless Green ~ d^-1) and T4
  (screening slope; gate = consistent with zero) itself, verified on
  synthetic Coulomb (recovers -1.00, PASS) and Yukawa xi=2 (recovers 2.0,
  FAIL).
- Usable shell range is capped by the slice radius (~125 cells/slice at
  20k/K=80); if T3/T4 report "insufficient range", raise N41 or lower K
  before reading anything into it.

Gates: T1 conservation exact (run aborts otherwise); T2 <f> at 7.4198;
M1 dressing-subtracted proportionality; M2 junction constancy; T3 ~ d^-1;
T4 massless.

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

## 6. Rung-2 curvature displacement — BLOCKED on #4

Needs the phase boundary located by the grid first. Then: bare vs
beta=1-coupled arms straddling the boundary; prediction c0 ~ +0.019
(closure) / ~0.08 per cell (seven-channel field) as a boundary shift.
Pre-registered in induced_couplings.py + paper J.4; no script yet — write
it against the measured boundary geometry when #4 lands.

---

Suggested order: check #4's CSV first (it may be done); then #1 and the
first seed of #2 together; #3 demo family alongside (cheap); remaining
seeds + #5 as slots free; #6 when the boundary is known; #3 production
last. Total ~3–5 days of PC time with 3–4 concurrent processes; the
paper-critical results (#1 gates, first #2 seed, #3 family shape) land in
the first ~24 h.
