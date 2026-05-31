#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_eprl_sweep.sh -- matched-volume EPRL coupling sweep for the v6 CDT engine
#
# Runs v6_theory_run.py at several beta_eprl values, each resuming from the SAME
# thermalized bare checkpoint and pinned to the SAME target volume, so the only
# thing that changes across runs is the EPRL coupling strength. With EPRL
# centering on (the default in v6_theory_run.py) the amplitude term is
# volume-neutral, so kappa_4 stays at its bare value for EVERY beta -- no
# per-beta retuning. The result you read is d_H / blob / cos3err vs beta,
# compared against the bare run at the same volume:
#   * a smooth, monotonic move with beta  => the EPRL amplitude steers geometry
#   * flat (== bare) across all beta       => EPRL leaves the 4D geometry alone
# Both are real, defensible results. (beta=0.0 is the control: it must reproduce
# the bare run -- if it does, the theory harness is faithful.)
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
#
# NOTE ON CORES: each beta is one single-threaded process. Launch at most as many
# concurrent betas as you have free vCPUs (stop the bare scans first if needed).
# ---------------------------------------------------------------------------
set -euo pipefail

BETAS="${1:-0.0 0.02 0.05 0.1 0.3}"
N41="${2:-20000}"
NK="$((N41/1000))k"

BASE="${BASE:-old/scan_${N41}.json}"
K="${K:-80}"; K0="${K0:-2.2}"; DELTA="${DELTA:-0.6}"; K4="${K4:-0.9}"; EPS="${EPS:-1e-3}"
MODE="${MODE:-regge_plus_eprl}"; AUDIT="${AUDIT:-25}"; MAXSWEEPS="${MAXSWEEPS:-1000000}"
WALL="${WALL:-}"; LAUNCH="${LAUNCH:-auto}"; CENTER="${CENTER:-1}"

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
wall_flag=""; [ -n "$WALL" ] && wall_flag="--wall-hours $WALL"

echo "# EPRL sweep: betas=[$BETAS]  N41=$N41  base=$BASE"
echo "# mode=$MODE  K=$K  k0=$K0 Delta=$DELTA k4=$K4 eps=$EPS  centering=$([ "$CENTER" = 0 ] && echo off || echo on)  launch=$LAUNCH"
[ "$CENTER" = 0 ] && echo "# !! centering OFF: you must retune --k4 per beta or the volume will drift."

for B in $BETAS; do
  tag="b${B}_${NK}"
  ckpt="results/thy_${tag}.json"
  log="logs/thy_${tag}.log"
  cmd="$PY -u v6_theory_run.py --mode $MODE --beta-eprl $B --local-eprl \
--audit-every $AUDIT $center_flag --resume $BASE --checkpoint $ckpt \
--target-n41 $N41 --K $K --k0 $K0 --Delta $DELTA --k4 $K4 --eps $EPS \
$wall_flag --max-sweeps $MAXSWEEPS"
  echo ">> beta=$B  -> $log"
  if [ "$LAUNCH" = systemd ]; then
    systemd-run --unit="cdt-eprl-${tag}" --collect \
      /bin/bash -c "cd '$PWD' && $cmd > '$log' 2>&1"
  else
    nohup bash -c "cd '$PWD' && $cmd > '$log' 2>&1" >/dev/null 2>&1 &
  fi
done

echo
echo "# launched. watch with:"
echo "#   tail -n 4 logs/thy_b*_${NK}.log"
echo "#   grep -hE 'centering|audit ok|Traceback|BAD' logs/thy_b*_${NK}.log | tail"
