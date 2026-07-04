#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_bridge_sweep.sh -- Regge-bridge SCALING-FAMILY sweep
#
# Upgrades the single -62 measurement to a family test of
#     Delta N41(beta, eps) = -beta*mu(beta) / (4*eps)
# (see bridge_predict.py / REGGE_BRIDGE.md / paper App. J.6). For every
# (beta, eps) point it launches a matched centered/uncentered PAIR resuming
# from the same clean base checkpoint; the pair difference isolates the
# mu-term. A placebo pair (orbit-shuffled table, with its OWN predicted
# mu_plc) runs at the largest beta and first eps.
#
# The analytic predictions are computed FIRST (bridge_predict.py, one TI
# pass per table on the base geometry) and written to results/, so every
# number is on record before the arms start.
#
# USAGE (from the v6_vm_bundle directory):
#   ./run_bridge_sweep.sh                       # betas "0.25 0.5 1.0", N41=20000
#   ./run_bridge_sweep.sh "0.5 1.0" 20000
#
# Env knobs (defaults shown):
#   BASE=scan_${N41}_causal.json   bare checkpoint (FIXED engine, foliation CLEAN)
#   EPSLIST="1e-3"                 add e.g. "5e-4 1e-3 2e-3" for the S2 eps-law
#   K=80 K0=2.2 DELTA=0.6 K4=0.9
#   ETA= (empty = eta*)  LAMBDA=3.0  AUDIT=25  MAXSWEEPS=1000000
#   WALL=  LAUNCH=auto  PLACEBO=1
#
# Afterwards:
#   python bridge_predict.py --analyze logs/bridge_*_${NK}.log \
#          --pred results/bridge_pred_${NK}.csv
# ---------------------------------------------------------------------------
set -euo pipefail

BETAS="${1:-0.25 0.5 1.0}"
N41="${2:-20000}"
NK="$((N41/1000))k"

BASE="${BASE:-scan_${N41}_causal.json}"
EPSLIST="${EPSLIST:-1e-3}"
K="${K:-80}"; K0="${K0:-2.2}"; DELTA="${DELTA:-0.6}"; K4="${K4:-0.9}"
ETA="${ETA:-}"; LAMBDA="${LAMBDA:-3.0}"
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
eta_flag=""; [ -n "$ETA" ] && eta_flag="--eta $ETA"
wall_flag=""; [ -n "$WALL" ] && wall_flag="--wall-hours $WALL"

echo "# BRIDGE sweep: betas=[$BETAS]  eps=[$EPSLIST]  N41=$N41  base=$BASE"
echo "# eta=${ETA:-eta*}  lambda_inj=$LAMBDA  K=$K k0=$K0 Delta=$DELTA k4=$K4  placebo=$([ "$PLACEBO" = 0 ] && echo off || echo on)  launch=$LAUNCH"

# --- predictions on record BEFORE the arms launch ---------------------------
PRED="results/bridge_pred_${NK}.csv"
echo "# computing analytic predictions -> $PRED"
$PY -u bridge_predict.py --resume "$BASE" --betas "$BETAS" \
  --eps-list "$EPSLIST" $eta_flag --lambda-inj "$LAMBDA" \
  --out "$PRED" | tee "logs/bridge_pred_${NK}.log"

launch_one () {  # $1=tag  $2=eps  $3...=extra python args
  local tag="$1" eps="$2"; shift 2
  local ckpt="results/${tag}.json"
  local log="logs/${tag}.log"
  local cmd="$PY -u v6_closure_run.py --local-closure --audit-every $AUDIT \
$eta_flag --lambda-inj $LAMBDA --resume $BASE --checkpoint $ckpt \
--target-n41 $N41 --K $K --k0 $K0 --Delta $DELTA --k4 $K4 --eps $eps \
$wall_flag --max-sweeps $MAXSWEEPS $*"
  echo ">> $tag  -> $log"
  if [ "$LAUNCH" = systemd ]; then
    systemd-run --unit="cdt-${tag}" --collect \
      /bin/bash -c "cd '$PWD' && $cmd > '$log' 2>&1"
  else
    nohup bash -c "cd '$PWD' && $cmd > '$log' 2>&1" >/dev/null 2>&1 &
  fi
}

BMAX=""; E1=""
for B in $BETAS; do
  case "$B" in 0|0.0) echo "# skip beta=0 (centered == uncentered)"; continue;; esac
  for E in $EPSLIST; do
    [ -z "$E1" ] && E1="$E"
    launch_one "bridge_b${B}_e${E}_cen_${NK}" "$E" "--beta-closure $B"
    launch_one "bridge_b${B}_e${E}_unc_${NK}" "$E" "--beta-closure $B --no-center-closure"
  done
  BMAX="$B"
done

if [ "$PLACEBO" != 0 ] && [ -n "$BMAX" ]; then
  launch_one "bridge_b${BMAX}_e${E1}_cenplc_${NK}" "$E1" "--beta-closure $BMAX --placebo"
  launch_one "bridge_b${BMAX}_e${E1}_uncplc_${NK}" "$E1" "--beta-closure $BMAX --placebo --no-center-closure"
fi

echo
echo "# launched. analyze when the arms have converged:"
echo "#   python bridge_predict.py --analyze logs/bridge_b*_${NK}.log --pred $PRED"
echo "# gates: ratio~1 across the family; concave beta-curve (S1);"
echo "#        1/eps collapse (S2); placebo tracks mu_plc not mu_real (S3);"
echo "#        foliation CLEAN on every arm."
