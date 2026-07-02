#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_closure_sweep.sh -- matched-volume ADMISSIBILITY-CLOSURE sweep
#
# Couples the paper's own cell weighting (K^2 admissibility at eta*, soft
# injectivity; see CLOSURE_MODEL.md) to the v6 CDT substrate. Arms resume
# from the SAME thermalized bare checkpoint at the same pinned volume.
# beta=1 is the THEORY-DESIGNATED point; beta=0 the control; a placebo arm
# (orbit-shuffled energy table) runs at the largest beta automatically.
#
# REQUIRES a bare checkpoint grown with the CAUSAL-SLICE-FIXED engine
# (foliation CLEAN) -- the closure cells live on the slices. Old defective
# checkpoints heal under the fixed engine but start clean if you can.
#
# USAGE (from the v6_vm_bundle directory):
#   ./run_closure_sweep.sh                     # betas "0.0 0.3 1.0", N41=20000
#   ./run_closure_sweep.sh "0.0 1.0" 20000
#
# Env knobs (defaults shown):
#   BASE=scan_${N41}_causal.json  bare checkpoint (fixed engine!)
#   K=80 K0=2.2 DELTA=0.6 K4=0.9 EPS=1e-3
#   ETA=            (empty = the paper's eta*)   LAMBDA=3.0 (soft injectivity)
#   AUDIT=25  MAXSWEEPS=1000000  WALL=  LAUNCH=auto  PLACEBO=1
# ---------------------------------------------------------------------------
set -euo pipefail

BETAS="${1:-0.0 0.3 1.0}"
N41="${2:-20000}"
NK="$((N41/1000))k"

BASE="${BASE:-scan_${N41}_causal.json}"
K="${K:-80}"; K0="${K0:-2.2}"; DELTA="${DELTA:-0.6}"; K4="${K4:-0.9}"; EPS="${EPS:-1e-3}"
ETA="${ETA:-}"; LAMBDA="${LAMBDA:-3.0}"
AUDIT="${AUDIT:-25}"; MAXSWEEPS="${MAXSWEEPS:-1000000}"
WALL="${WALL:-}"; LAUNCH="${LAUNCH:-auto}"; PLACEBO="${PLACEBO:-1}"

PY="$PWD/.venv/bin/python3"; [ -x "$PY" ] || PY="$(command -v python3)"
mkdir -p logs results

if [ ! -f "$BASE" ]; then
  echo "ERROR: bare checkpoint '$BASE' not found." >&2
  echo "       Grow one with the FIXED engine: python v6_verify_run.py" >&2
  echo "       --target-n41 $N41 --K $K --checkpoint $BASE  (foliation must be CLEAN)" >&2
  exit 1
fi

if [ "$LAUNCH" = auto ]; then
  if command -v systemd-run >/dev/null 2>&1; then LAUNCH=systemd; else LAUNCH=nohup; fi
fi

eta_flag=""; [ -n "$ETA" ] && eta_flag="--eta $ETA"
wall_flag=""; [ -n "$WALL" ] && wall_flag="--wall-hours $WALL"

echo "# CLOSURE sweep: betas=[$BETAS]  N41=$N41  base=$BASE"
echo "# eta=${ETA:-eta*}  lambda_inj=$LAMBDA  K=$K k0=$K0 Delta=$DELTA k4=$K4 eps=$EPS  placebo=$([ "$PLACEBO" = 0 ] && echo off || echo on)  launch=$LAUNCH"

launch_one () {  # $1=tag  $2...=extra python args
  local tag="$1"; shift
  local ckpt="results/clo_${tag}.json"
  local log="logs/clo_${tag}.log"
  local cmd="$PY -u v6_closure_run.py --local-closure --audit-every $AUDIT \
$eta_flag --lambda-inj $LAMBDA --resume $BASE --checkpoint $ckpt \
--target-n41 $N41 --K $K --k0 $K0 --Delta $DELTA --k4 $K4 --eps $EPS \
$wall_flag --max-sweeps $MAXSWEEPS $*"
  echo ">> $tag  -> $log"
  if [ "$LAUNCH" = systemd ]; then
    systemd-run --unit="cdt-clo-${tag}" --collect \
      /bin/bash -c "cd '$PWD' && $cmd > '$log' 2>&1"
  else
    nohup bash -c "cd '$PWD' && $cmd > '$log' 2>&1" >/dev/null 2>&1 &
  fi
}

BMAX=""
for B in $BETAS; do
  launch_one "b${B}_${NK}" "--beta-closure $B"
  BMAX="$B"
done

if [ "$PLACEBO" != 0 ] && [ -n "$BMAX" ] && [ "$BMAX" != 0.0 ] && [ "$BMAX" != 0 ]; then
  launch_one "b${BMAX}plc_${NK}" "--beta-closure $BMAX --placebo"
fi

echo
echo "# launched. watch:  tail -n 4 logs/clo_*_${NK}.log"
echo "#   grep -hE 'mu TI|centering|audit ok|foliation|Traceback|BAD' logs/clo_*_${NK}.log | tail -n 30"
echo "# gates: foliation CLEAN on every arm; matched N4; then real-vs-placebo at beta=$BMAX."
