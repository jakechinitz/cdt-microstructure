# cdt-microstructure

Lattice program for the entropic scalar EFT paper: the tetrahedral
admissibility ensemble (Part II of the paper) coupled to a dynamical
causal-dynamical-triangulations host, with the pre-registered
compatibility, defect, bridge, and capacity-conservation experiments
(paper Section 23 / Appendix J).

## Layout

| path | what |
|---|---|
| `paper/action_workspace_v5.tex` | the manuscript (v5 action workspace) |
| `v6_vm_bundle/` | the simulation: engine, runs, analyzers, docs. **All commands run from here.** |
| `results/2026-07-03/` | stage 1–4 production runs at N41=20k (bare rescan, closure sweep, defect, capacity v1.1, K40 pair, phase grid) |
| `results/2026-07-05/` | replication + follow-ups (defect seed 1, bridge family at 2k, tuned grid, capacity v1.1 rerun) |
| `archive/` | superseded/duplicate uploads, kept for provenance (nothing is deleted) |

## Quick start — the 100k scaling program (one launch)

```bash
cd v6_vm_bundle
pip install numpy scipy
python run_100k_program.py --dir runs/p100k --dry-run        # inspect the plan
python run_100k_program.py --dir runs/smoke --smoke --jobs 4 # minutes-long e2e check
nohup python -u run_100k_program.py --dir runs/p100k --jobs 6 \
      > runs/p100k.launch.log 2>&1 &                         # the real thing (~days)
```

The program grows the 100k bases at the tuned phase point, thermalizes the
closure arms, runs the defect and capacity experiments with gates between
stages, and writes `runs/p100k/REPORT.md`. It is resumable: re-run the same
command after any interruption. See `v6_vm_bundle/PROGRAM_100K.md`.

Key docs inside `v6_vm_bundle/`: `RUN_ON_VM.md` (engine + bare scans),
`CLOSURE_MODEL.md` (the coupled measure), `CAPACITY_CONSERVATION.md`
(the transport instrument), `RUNBOOK_REMAINING.md` (remaining passes:
bridge production family, rung-2 after the ladder verdict).
