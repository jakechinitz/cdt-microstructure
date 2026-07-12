# The 100k scaling program — one launch

`run_100k_program.py` runs the full pre-registered battery for the
production-volume (N41=100k) graph as a single resumable launch: it grows
the bases, thermalizes the coupled arms, runs the experiments, applies the
gates between stages, and writes `REPORT.md` with autocorrelation-corrected
numbers. Re-running the same command resumes; finished jobs never repeat.

## What changed relative to the v5 runs (and why)

1. **Tuned phase point.** All arms run at (k0=2.0, Delta=0.4) with k4
   auto-tuned per volume (`v6_verify_run.py --tune-k4`), not the borrowed
   literature point (2.2, 0.6, k4=0.9). The autotuned phase grid showed
   fixed k4=0.9 fights the volume pin across the grid; the candidate
   extended region peaks at (2.0, 0.4) with tuned k4~0.7. The paper
   (App. J.1) pre-registers the confirmation as a volume pair at that
   peak — stage A of this program IS that confirmation run.
2. **Autocorrelation-corrected errors** (`autocorr.py`, Sokal windowing;
   self-test included). Every analyzer now quotes IAT-corrected within-run
   SEMs and across-seed/across-pin scatter. This removes the standing
   "errors not thinned for autocorrelation" caveat from the defect and
   capacity results. On the existing 20k stage-3 data the correction is
   material for the geometry observables (tau ~ 100 measurement steps:
   coordination/volume errors were ~10x underestimated; the shell-1
   capacity separation survives at ~30 sigma, the geometric separations
   at ~4-5 sigma).
3. **Replication is default.** Two defect seed pairs (real+placebo,
   shared pin anchors) run automatically; the analyzer quotes across-seed
   scatter when it has more than one seed.
4. **Capacity instrument gets the radius it needs.** The K=40 sibling
   base at 100k carries ~1250 cells/slice (vs 250 at 20k/K40), pushing the
   usable graph-Coulomb radius to ~12 steps: the T3 profile exponent
   (-0.65±0.14 at 40 slices, ideal -1) and the T4 massless bound both
   sharpen. `--rmax-cap 8` is the default measurement range.

## Launch

```bash
cd v6_vm_bundle
# sanity: plan only
python run_100k_program.py --dir runs/p100k --dry-run
# tiny end-to-end pipeline check (minutes)
python run_100k_program.py --dir runs/smoke --smoke --jobs 4
# the real thing (needs ~6 cores free; several days wall time)
nohup python -u run_100k_program.py --dir runs/p100k --jobs 6 \
      > runs/p100k.launch.log 2>&1 &
```

Watch: `tail -f runs/p100k.launch.log` and the per-job logs in
`runs/p100k/logs/`. State: `runs/p100k/state.json`.

## Stages, dependencies, gates

| stage | job | what | gate |
|---|---|---|---|
| A | base_K80 | bare 100k, K=80, tuned point | manifold PASS + foliation CLEAN |
| A | base_K40 | bare 100k, K=40 (capacity substrate) | same |
| A | base_ladder | bare 40k, K=80, same point (volume pair) | same |
| B | clo_b0 / clo_b1 / clo_plc | closure arms on base_K80 (control / theory point / orbit-shuffled placebo) | manifold + foliation + no audit drift |
| B | clo_b1_K40 | theory-point arm on base_K40 | same |
| C | def_s{1,2}_{real,placebo} | pinned-defect seed pairs on clo_b1 | pins alive N/N |
| C | cap_v12 | conservation instrument (absorb=excess, persist=200) on clo_b1_K40 | exact conservation (abort otherwise) + completion |
| D | report | analyzers + gate table -> REPORT.md | — |

A failed gate blocks the failed arm's dependents only; siblings continue.
`--retry-failed` clears failed jobs so a relaunch reruns them.

## Reading REPORT.md

* **Volume pair**: blob and d_H must RISE from 40k to 100k at the tuned
  point — the pre-registered acceptance for the candidate extended region.
* **Closure arms**: collision fraction must order (b0 ~0.65 -> b1 ordered,
  placebo high); K80 arms' N4 spread < 2% (matched volume); geometric
  observables statistically unchanged.
* **Defect**: real-vs-placebo separation table with z-scores
  (autocorrelation-corrected; across-seed where seeds > 1).
* **Capacity**: M1 dressing-subtracted linearity, M2 junction constancy,
  T3 profile exponent, T4 screening bound — the same pre-registered gates
  as CAPACITY_CONSERVATION.md, at 2-3x the usable radius.

## Deliberately NOT in the one-launch

* **Rung-2 curvature displacement** (`rung2_boundary.py`): pre-registered
  to launch only after the volume-pair verdict confirms the boundary
  region. Once REPORT.md's ladder verdict is positive:
  `python rung2_boundary.py --predict` then the k0 scan (see
  RUNBOOK_REMAINING.md #6).
* **Bridge production family** (`run_bridge_sweep.sh`): its discipline is
  demo-volume-first with predictions written before launch; run it as its
  own pass (RUNBOOK_REMAINING.md #3).
* **Capacity-through-moves (v2 instrument)** and the **oriented-closure
  fork**: physics-code extensions, to be built and validated on small
  volumes with planted-answer tests before they touch production data —
  not bolted on in the same pass as the scaling run.

## Cost anchor

40k/K80 ran ~9 s/sweep at N4~92k (2026-07 logs). Per-sweep cost scales
~linearly with N4, so 100k (N4~235k) is ~20-25 s/sweep: stage A ~1-1.2
days/arm, stage B similar, stage C defects ~1 day/arm and capacity is
cheaper (frozen geometry). With `--jobs 6` the whole program is ~3-4 days
wall; RAM per 100k process is the main sizing check (a 20k checkpoint JSON
is ~50 MB; in-memory ~a few GB — verify early with `top`).
