#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_eprl_sweep.sh -- matched-volume EPRL coupling sweep for the v6 CDT engine
#
# Runs v6_theory_run.py at several beta_eprl values, each resuming from the SAME
# thermalized bare checkpoint and pinned to the SAME target volume, so the only
# thing that changes across runs is the EPRL coupling strength.
#
# POST-AUDIT PROTOCOL (see SIM_AUDIT_coupling_misspecifications.md):
#   * Each arm auto-calibrates its own volume-neutral mu by thermodynamic
#     integration on the base checkpoint (a few minutes, printed as "# [mu TI]"
#     lines at startup) and stores it in its checkpoint for resumes.
#   * The label heat-bath runs at the arm's own beta (beta-consistent sampler).
#   * The vertex tensor is slot-symmetrized by default (ordering-invariant
#     amplitude). SYMM=0 restores the raw tensor.
#   * A PLACEBO ARM runs automatically at the largest beta with an
#     entry-shuffled tensor (same statistics, EPRL structure destroyed).
#
# READING THE RESULT (gates first, then the comparison):
#   G3': matched volume now means N4 (BOTH simplex types), not just the pinned
#        N41, comparable across all arms. If N4 drifts with beta, the coupling
#        is not volume-neutral -- fix before reading shape observables.
#   PLACEBO: only a real-vs-shuffled DIFFERENCE at the same beta is
#        attributable to the EPRL amplitude's structure. If the placebo
#        reproduces the real arm, the movement is machinery/entry-statistics.
#   Then: d_H / blob / cos3err vs beta against the bare run at the same volume.
#   (beta=0.0 remains the harness-faithfulness control; note it cannot catch
#   beta-dependent artifacts, which is what the placebo arm is for.)
#
# USAGE (from the v6_vm_bundle directory):
#   ./run_eprl_sweep.sh                       # defaults: betas + N41=20000
#   ./run_eprl_sweep.sh "0.0 0.05 0.1 0.3"    # custom beta list
#   ./run_eprl_sweep.sh "0.0 0.1" 40000       # custom betas + volume
#
# Knobs via environment variables (sensible defaults shown):
#   BASE=old/scan_20000.json   bare checkpoint to resume each run from
#   K=80  K0=2.2  DELTA=0.6  K4=0.9  EPS=1e-3   action / volume-clamp settings
#   MODE=regge_plus_eprl       (or eprl_only)
#   AUDIT=25                   incremental-action self-audit cadence (0=off)
#   WALL=                      per-run wall-hours cap (e.g. WALL=12); empty=none
#   MAXSWEEPS=1000000          hard sweep cap
#   LAUNCH=auto                auto|systemd|nohup process supervisor
#   CENTER=1                   1=centered (volume-neutral); 0=raw (--no-center-eprl)
#   SYMM=1                     1=slot-symmetrized tensor; 0=raw tensor
#   PLACEBO=1                  1=also run the shuffled-tensor arm at max beta
#
# NOTE ON CORES: each beta is one single-threaded process (the placebo arm is
# one more). Launch at most as many concurrent arms as you have free vCPUs.
# ---------------------------------------------------------------------------
set -euo pipefail

BETAS="${1:-0.0 0.02 0.05 0.1 0.3}"
N41="${2:-20000}"
NK="$((N41/1000))k"

BASE="${BASE:-old/scan_${N41}.json}"
K="${K:-80}"; K0="${K0:-2.2}"; DELTA="${DELTA:-0.6}"; K4="${K4:-0.9}"; EPS="${EPS:-1e-3}"
MODE="${MODE:-regge_plus_eprl}"; AUDIT="${AUDIT:-25}"; MAXSWEEPS="${MAXSWEEPS:-1000000}"
WALL="${WALL:-}"; LAUNCH="${LAUNCH:-auto}"; CENTER="${CENTER:-1}"
SYMM="${SYMM:-1}"; PLACEBO="${PLACEBO:-1}"

PY="$PWD/.venv/bin/python3"; [ -x "$PY" ] || PY="$(command -v python3)"
mkdir -p logs results

if [ ! -f "$BASE" ]; then
  echo "ERROR: bare checkpoint '$BASE' not found." >&2
  echo "       Point BASE= at a thermalized bare run at N41=$N41 (see RUN_ON_VM.md step 2)." >&2
  exit 1
fi

# Pick a process supervisor.
if [ "$LAUNCH" = auto ]; then
  if command -v systemd-run >/dev/null 2>&1; then LAUNCH=systemd; else LAUNCH=nohup; fi
fi

center_flag="--center-eprl"; [ "$CENTER" = 0 ] && center_flag="--no-center-eprl"
symm_flag="--symmetrize-vertex"; [ "$SYMM" = 0 ] && symm_flag="--no-symmetrize-vertex"
wall_flag=""; [ -n "$WALL" ] && wall_flag="--wall-hours $WALL"

echo "# EPRL sweep: betas=[$BETAS]  N41=$N41  base=$BASE"
echo "# mode=$MODE  K=$K  k0=$K0 Delta=$DELTA k4=$K4 eps=$EPS  centering=$([ "$CENTER" = 0 ] && echo off || echo 'on (mu via TI per arm)')  symmetrize=$([ "$SYMM" = 0 ] && echo off || echo on)  placebo=$([ "$PLACEBO" = 0 ] && echo off || echo on)  launch=$LAUNCH"
[ "$CENTER" = 0 ] && echo "# !! centering OFF: you must retune --k4 per beta or the volume will drift."

launch_one () {  # $1=tag  $2=extra python args
  local tag="$1"; shift
  local ckpt="results/thy_${tag}.json"
  local log="logs/thy_${tag}.log"
  local cmd="$PY -u v6_theory_run.py --mode $MODE --local-eprl \
--audit-every $AUDIT $center_flag $symm_flag --resume $BASE --checkpoint $ckpt \
--target-n41 $N41 --K $K --k0 $K0 --Delta $DELTA --k4 $K4 --eps $EPS \
$wall_flag --max-sweeps $MAXSWEEPS $*"
  echo ">> $tag  -> $log"
  if [ "$LAUNCH" = systemd ]; then
    systemd-run --unit="cdt-eprl-${tag}" --collect \
      /bin/bash -c "cd '$PWD' && $cmd > '$log' 2>&1"
  else
    nohup bash -c "cd '$PWD' && $cmd > '$log' 2>&1" >/dev/null 2>&1 &
  fi
}

BMAX=""
for B in $BETAS; do
  launch_one "b${B}_${NK}" "--beta-eprl $B"
  BMAX="$B"   # betas are listed ascending; last = largest
done

# --- placebo arm: entry-shuffled tensor at the largest beta -----------------
if [ "$PLACEBO" != 0 ] && [ -n "$BMAX" ] && [ "$BMAX" != 0.0 ] && [ "$BMAX" != 0 ]; then
  SHUF="vertex_j3_shuffled.npz"
  [ -f "$SHUF" ] || "$PY" make_shuffled_control.py --out "$SHUF"
  launch_one "b${BMAX}shuf_${NK}" "--beta-eprl $BMAX --vertex $SHUF"
fi

echo
echo "# launched. watch with:"
echo "#   tail -n 4 logs/thy_b*_${NK}.log"
echo "#   grep -hE 'mu TI|centering|audit ok|Traceback|BAD' logs/thy_b*_${NK}.log | tail -n 30"
echo "# gates: (1) every arm's N4 comparable at matched N41; (2) real-vs-shuffled"
echo "#        difference at beta=$BMAX is the EPRL signal -- no difference = no signal."
