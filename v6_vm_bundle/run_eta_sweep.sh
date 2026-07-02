#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_eta_sweep.sh -- is the paper's eta* dynamically special?
#
# Sweeps the admissibility precision eta AROUND the theory-fixed value
# eta* = 0.0298668443935 at full coupling (beta=1), all arms resuming from
# the same clean bare checkpoint at the same pinned volume, plus a beta=0
# control and a placebo arm at eta*. If the geometry response (d_H / blob /
# cos3err) is extremal or otherwise distinguished AT eta*, that is a
# dynamical echo of the closure-evidence stationarity argument (App. B) --
# the sharpest single figure this program can produce. If eta* is nothing
# special dynamically, that is also worth knowing and does not by itself
# contradict the paper (eta* is fixed by closure evidence, not by geometry).
#
# USAGE (from the v6_vm_bundle directory):
#   ./run_eta_sweep.sh                                   # default grid, N41=20000
#   ./run_eta_sweep.sh "0.0075 0.015 0.0299 0.06 0.12" 20000
#
# Env knobs: BASE (clean fixed-engine checkpoint), K K0 DELTA K4 EPS LAMBDA
#            AUDIT MAXSWEEPS WALL LAUNCH PLACEBO  (as in run_closure_sweep.sh)
#
# Afterwards:
#   python plot_sweep.py logs/clo_eta*_${NK}.log --x eta --logx \
#          --control logs/clo_eta_ctrl_${NK}.log --mark-x 0.0298668443935 \
#          --out eta_sweep_${NK}
# ---------------------------------------------------------------------------
set -euo pipefail

ETAS="${1:-0.0075 0.015 0.0298668443935 0.06 0.12}"
N41="${2:-20000}"
NK="$((N41/1000))k"
ETA_STAR="0.0298668443935"

BASE="${BASE:-scan_${N41}_causal.json}"
K="${K:-80}"; K0="${K0:-2.2}"; DELTA="${DELTA:-0.6}"; K4="${K4:-0.9}"; EPS="${EPS:-1e-3}"
LAMBDA="${LAMBDA:-3.0}"
AUDIT="${AUDIT:-25}"; MAXSWEEPS="${MAXSWEEPS:-1000000}"
WALL="${WALL:-}"; LAUNCH="${LAUNCH:-auto}"; PLACEBO="${PLACEBO:-1}"

PY="$PWD/.venv/bin/python3"; [ -x "$PY" ] || PY="$(command -v python3)"
mkdir -p logs results

if [ ! -f "$BASE" ]; then
  echo "ERROR: bare checkpoint '$BASE' not found (grow with the FIXED engine)." >&2
  exit 1
fi
if [ "$LAUNCH" = auto ]; then
  if command -v systemd-run >/dev/null 2>&1; then LAUNCH=systemd; else LAUNCH=nohup; fi
fi
wall_flag=""; [ -n "$WALL" ] && wall_flag="--wall-hours $WALL"

echo "# ETA sweep at beta=1: etas=[$ETAS] (eta*=$ETA_STAR)  N41=$N41  base=$BASE"
echo "# lambda_inj=$LAMBDA  K=$K k0=$K0 Delta=$DELTA k4=$K4 eps=$EPS  placebo@eta*=$([ "$PLACEBO" = 0 ] && echo off || echo on)"

launch_one () {
  local tag="$1"; shift
  local ckpt="results/clo_${tag}.json"
  local log="logs/clo_${tag}.log"
  local cmd="$PY -u v6_closure_run.py --local-closure --audit-every $AUDIT \
--lambda-inj $LAMBDA --resume $BASE --checkpoint $ckpt \
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

launch_one "eta_ctrl_${NK}" "--beta-closure 0.0"
for E in $ETAS; do
  launch_one "eta${E}_${NK}" "--beta-closure 1.0 --eta $E"
done
[ "$PLACEBO" != 0 ] && launch_one "eta${ETA_STAR}plc_${NK}" \
  "--beta-closure 1.0 --eta $ETA_STAR --placebo"

echo
echo "# launched. plot when finished:"
echo "#   python plot_sweep.py logs/clo_eta0*_${NK}.log --x eta --logx \\"
echo "#          --control logs/clo_eta_ctrl_${NK}.log --mark-x $ETA_STAR --out eta_sweep_${NK}"
