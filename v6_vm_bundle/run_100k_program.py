#!/usr/bin/env python3
"""ONE-LAUNCH 100k PROGRAM -- grow the production graph, thermalize the
coupled arms, and run the full pre-registered test battery, with gates
enforced between stages. Designed so that once launched, the main outputs
of the scaling program arrive without further babysitting:

  stage A  BASES (parallel)
    A1 base_K80    bare growth+thermalization at the TUNED phase point
                   (k0=2.0, Delta=0.4, k4 auto-tuned per volume -- the
                   phase-grid finding: fixed k4=0.9 fights the pin)
    A2 base_K40    thick-slice sibling (double the cells per slice =
                   ~2x usable far-field radius for the capacity gates)
    A3 base_ladder lower-volume rung at the SAME point (default N41/2.5):
                   the pre-registered volume-pair confirmation -- blob and
                   d_H must RISE with volume for the candidate region
  stage B  CLOSURE ARMS (need A; matched-volume compatibility test)
    B1 clo_b0      beta=0 control on base_K80 (exact-bare by construction)
    B2 clo_b1      beta=1 theory point on base_K80
    B3 clo_plc     beta=1 orbit-shuffled placebo on base_K80
    B4 clo_b1_K40  beta=1 on base_K40 (the capacity instrument's substrate)
  stage C  EXPERIMENTS (need B)
    C*  defect seed pairs (real=fail + placebo=closed per seed, shared pin
        locations) resuming clo_b1  [replication: >= 2 seeds default]
    C+  capacity v1.2 (absorb=excess, persistent-failure window) resuming
        clo_b1_K40  [T1-T4 + M1/M2 gates]
  stage D  ANALYSIS -> REPORT.md
    autocorrelation-corrected analyzers (autocorr.py), gate table, ladder
    verdict, matched-volume check, collision-fraction ordering.

GATES between stages (an arm that fails its gate blocks its dependents,
never the siblings): manifold PASS + foliation CLEAN on every geometry
arm; pin survival N/N on defect arms; conservation (the capacity run
aborts on violation).

USAGE
  python run_100k_program.py --dir runs/p100k --jobs 6            # launch
  python run_100k_program.py --dir runs/p100k --jobs 6            # resume
  python run_100k_program.py --dir runs/p100k --dry-run           # plan
  python run_100k_program.py --dir runs/p100k --smoke --jobs 4    # tiny e2e
  python run_100k_program.py --dir runs/p100k --report-only

State lives in <dir>/state.json (job -> status); finished jobs are never
re-run, so the script can be re-launched after a crash/reboot and it
continues where it stopped. --retry-failed clears failed jobs first.
"""
from __future__ import annotations
import argparse
import json
import os
import re
import subprocess
import sys
import time

BUNDLE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable or "python3"


# ---------------------------------------------------------------------------
# job table
# ---------------------------------------------------------------------------

def build_jobs(a):
    """Returns {name: job} with job = dict(cmd, log, deps, gate, stage).
    Paths are absolute inside a.dir. Commands that need the tuned k4 of
    their base carry the placeholder '@K4:<base-job>' resolved at launch
    time from the base's log."""
    d = a.dir
    j = {}

    def base_job(name, K, n41, stage="A"):
        ck = os.path.join(d, f"{name}.json")
        j[name] = dict(stage=stage, deps=[], ckpt=ck, cmd=[
            PY, "-u", os.path.join(BUNDLE, "v6_verify_run.py"),
            "--target-n41", str(n41), "--K", str(K),
            "--k0", str(a.k0), "--Delta", str(a.Delta), "--k4", str(a.k4),
            "--eps", str(a.eps_grow), "--seed", str(a.seed),
            "--tune-k4", str(a.tune_k4), "--tune-burst", str(a.tune_burst),
            "--max-sweeps", str(a.sweeps_base),
            "--measure-every", str(a.measure_every_base),
            "--checkpoint", ck],
            log=os.path.join(d, "logs", f"{name}.log"),
            gate=dict(require=[r"manifold check .*: PASS",
                               r"FINAL foliation: .*CLEAN"],
                      forbid=[r"Traceback", r"verify FAILED"]))

    base_job("base_K80", a.K, a.n41)
    base_job("base_K40", a.K_cap, a.n41)
    if a.ladder:
        base_job("base_ladder", a.K, a.n41_ladder)

    def closure_job(name, base, beta, K, placebo=False):
        ck = os.path.join(d, f"{name}.json")
        cmd = [PY, "-u", os.path.join(BUNDLE, "v6_closure_run.py"),
               "--local-closure", "--audit-every", str(a.audit_every),
               "--lambda-inj", str(a.lambda_inj),
               "--beta-closure", str(beta),
               "--resume", j[base]["ckpt"], "--checkpoint", ck,
               "--target-n41", str(a.n41), "--K", str(K),
               "--k0", str(a.k0), "--Delta", str(a.Delta),
               "--k4", f"@K4:{base}", "--eps", str(a.eps_couple),
               "--seed", str(a.seed),
               "--max-sweeps", str(a.sweeps_closure),
               "--measure-every", str(a.measure_every_closure)]
        if placebo:
            cmd += ["--placebo"]
        j[name] = dict(stage="B", deps=[base], ckpt=ck, cmd=cmd,
                       log=os.path.join(d, "logs", f"{name}.log"),
                       gate=dict(require=[r"manifold check: PASS",
                                          r"foliation: CLEAN"],
                                 forbid=[r"Traceback", r"BAD"]))

    closure_job("clo_b0", "base_K80", 0.0, a.K)
    closure_job("clo_b1", "base_K80", 1.0, a.K)
    closure_job("clo_plc", "base_K80", 1.0, a.K, placebo=True)
    closure_job("clo_b1_K40", "base_K40", 1.0, a.K_cap)

    for s in range(1, a.seeds + 1):
        for mode, tag in (("fail", "real"), ("closed", "placebo")):
            name = f"def_s{s}_{tag}"
            out = os.path.join(d, name)
            j[name] = dict(stage="C", deps=["clo_b1"], csv=out + ".csv", cmd=[
                PY, "-u", os.path.join(BUNDLE, "v6_defect_run.py"),
                "--resume", j["clo_b1"]["ckpt"], "--pin-mode", mode,
                "--pins", str(a.pins), "--pin-sep", str(a.pin_sep),
                "--rmax", str(a.rmax_defect), "--seed", str(s),
                "--sweeps", str(a.sweeps_defect),
                "--therm", str(a.therm_defect),
                "--target-n41", str(a.n41), "--K", str(a.K),
                "--k0", str(a.k0), "--Delta", str(a.Delta),
                "--k4", "@K4:base_K80", "--eps", str(a.eps_couple),
                "--out", out],
                log=os.path.join(d, "logs", f"{name}.log"),
                gate=dict(require=[r"pins alive"],
                          forbid=[r"PIN LOSS", r"Traceback"]))

    cap_out = os.path.join(d, "cap_v12")
    j["cap_v12"] = dict(stage="C", deps=["clo_b1_K40"], csv=cap_out + ".csv",
                        cmd=[
        PY, "-u", os.path.join(BUNDLE, "v6_capacity_run.py"),
        "--resume", j["clo_b1_K40"]["ckpt"],
        "--absorb", "excess", "--persist", str(a.persist),
        "--pins-per-level", str(a.pins_per_level),
        "--pin-sep", str(a.pin_sep_cap), "--rmax", str(a.rmax_cap),
        "--sweeps", str(a.sweeps_cap), "--therm", str(a.therm_cap),
        "--measure-every", str(a.measure_every_cap),
        "--seed", str(a.seed), "--out", cap_out],
        log=os.path.join(d, "logs", "cap_v12.log"),
        gate=dict(require=[r"STAGE 4 RUN DONE"],
                  forbid=[r"conservation violated", r"Traceback"]))
    return j


# ---------------------------------------------------------------------------
# scheduler
# ---------------------------------------------------------------------------

def load_state(path):
    return json.load(open(path)) if os.path.exists(path) else {}


def save_state(path, st):
    tmp = path + ".tmp"
    json.dump(st, open(tmp, "w"), indent=1)
    os.replace(tmp, path)


def tuned_k4(log, default):
    """Last '# [k4-tune] final k4 = X' in a base log, else the default."""
    if not os.path.exists(log):
        return default
    val = default
    for line in open(log, errors="replace"):
        m = re.search(r"\[k4-tune\] final k4 = ([-\d.]+)", line)
        if m:
            val = m.group(1)
    return val


def gate_check(job):
    """True if the job's log passes its require/forbid patterns."""
    log = job["log"]
    if not os.path.exists(log):
        return False, "log missing"
    txt = open(log, errors="replace").read()
    for pat in job["gate"].get("forbid", []):
        if re.search(pat, txt):
            return False, f"forbidden pattern '{pat}'"
    for pat in job["gate"].get("require", []):
        if not re.search(pat, txt):
            return False, f"missing required pattern '{pat}'"
    return True, "ok"


def run_jobs(jobs, args):
    state_path = os.path.join(args.dir, "state.json")
    st = load_state(state_path)
    if args.retry_failed:
        for k, v in list(st.items()):
            if v.get("status") == "failed":
                del st[k]
        save_state(state_path, st)
    os.makedirs(os.path.join(args.dir, "logs"), exist_ok=True)
    running = {}                            # name -> (Popen, t0, logf)

    def ready(name):
        job = jobs[name]
        if name in running or st.get(name, {}).get("status") in ("done",
                                                                 "failed"):
            return False
        return all(st.get(dep, {}).get("status") == "done"
                   for dep in job["deps"])

    def blocked(name):
        """Permanently blocked: some dependency failed."""
        return any(st.get(dep, {}).get("status") == "failed"
                   for dep in jobs[name]["deps"])

    def resolve(cmd):
        out = []
        for tok in cmd:
            if tok.startswith("@K4:"):
                base = tok[4:]
                out.append(str(tuned_k4(jobs[base]["log"], args.k4)))
            else:
                out.append(tok)
        return out

    pending = [n for n in jobs if st.get(n, {}).get("status") != "done"]
    print(f"# scheduler: {len(pending)} job(s) to run, "
          f"{len(jobs) - len(pending)} already done, jobs={args.jobs}")
    while True:
        # launch what we can
        for name in jobs:
            if len(running) >= args.jobs:
                break
            if ready(name) and not blocked(name):
                cmd = resolve(jobs[name]["cmd"])
                logf = open(jobs[name]["log"], "a")
                logf.write(f"\n# === launch {time.strftime('%F %T')}: "
                           f"{' '.join(cmd)}\n")
                logf.flush()
                p = subprocess.Popen(cmd, stdout=logf, stderr=subprocess.STDOUT,
                                     cwd=BUNDLE)
                running[name] = (p, time.time(), logf)
                print(f"[{time.strftime('%H:%M:%S')}] LAUNCH {name} "
                      f"(pid {p.pid}) -> {jobs[name]['log']}", flush=True)
        # mark permanently-blocked jobs
        for name in jobs:
            if (name not in running and name not in st and blocked(name)):
                st[name] = {"status": "failed",
                            "note": "dependency failed"}
                print(f"[{time.strftime('%H:%M:%S')}] BLOCKED {name} "
                      f"(dependency failed)", flush=True)
                save_state(state_path, st)
        # reap
        done_now = []
        for name, (p, t0, logf) in running.items():
            rc = p.poll()
            if rc is None:
                continue
            logf.close()
            ok, why = (gate_check(jobs[name]) if rc == 0
                       else (False, f"exit code {rc}"))
            st[name] = {"status": "done" if ok else "failed",
                        "rc": rc, "gate": why,
                        "wall_s": int(time.time() - t0)}
            save_state(state_path, st)
            print(f"[{time.strftime('%H:%M:%S')}] "
                  f"{'DONE ' if ok else 'FAIL '}{name}  rc={rc}  "
                  f"gate={why}  wall={st[name]['wall_s']}s", flush=True)
            done_now.append(name)
        for name in done_now:
            del running[name]
        # finished?
        unfinished = [n for n in jobs
                      if st.get(n, {}).get("status") not in ("done", "failed")]
        if not unfinished and not running:
            break
        if not running and not any(ready(n) and not blocked(n)
                                   for n in unfinished):
            # deadlock guard (shouldn't happen)
            for n in unfinished:
                st[n] = {"status": "failed", "note": "unschedulable"}
            save_state(state_path, st)
            break
        time.sleep(args.poll)
    return st


# ---------------------------------------------------------------------------
# analysis + report
# ---------------------------------------------------------------------------

def _grep1(log, pat, group=1, default=None):
    if not os.path.exists(log):
        return default
    out = default
    for line in open(log, errors="replace"):
        m = re.search(pat, line)
        if m:
            out = m.group(group)
    return out


def _run_capture(cmd):
    try:
        r = subprocess.run(cmd, cwd=BUNDLE, capture_output=True, text=True,
                           timeout=3600)
        return r.stdout + (("\n[stderr]\n" + r.stderr) if r.returncode else "")
    except Exception as e:                                # noqa: BLE001
        return f"(analyzer failed: {e})"


def write_report(jobs, st, args):
    d = args.dir
    L = []
    L.append(f"# 100k program report -- {time.strftime('%F %T')}")
    L.append(f"\nrun dir: `{d}`   N41={args.n41} (ladder {args.n41_ladder})"
             f"  K={args.K}/{args.K_cap}  point: k0={args.k0} "
             f"Delta={args.Delta} (k4 auto-tuned)\n")
    # --- job/gate table ------------------------------------------------------
    L.append("## Jobs and gates\n")
    L.append("| job | stage | status | gate | wall |")
    L.append("|---|---|---|---|---|")
    for name, job in jobs.items():
        s = st.get(name, {})
        L.append(f"| {name} | {job['stage']} | {s.get('status', 'pending')} "
                 f"| {s.get('gate', s.get('note', ''))} "
                 f"| {s.get('wall_s', '')} |")
    # --- tuned k4 ------------------------------------------------------------
    L.append("\n## Tuned couplings\n")
    for b in ("base_K80", "base_K40", "base_ladder"):
        if b in jobs:
            L.append(f"- {b}: k4 = {tuned_k4(jobs[b]['log'], args.k4)}")
    # --- ladder verdict ------------------------------------------------------
    if "base_ladder" in jobs:
        L.append("\n## Volume-pair (phase-region confirmation)\n")
        L.append("Pre-registered acceptance: blob score and d_H must RISE "
                 "with volume at the tuned point.\n")
        L.append("| base | N41 | d_H | blob | active |")
        L.append("|---|---|---|---|---|")
        for b, n in (("base_ladder", args.n41_ladder), ("base_K80", args.n41)):
            log = jobs[b]["log"]
            dh = _grep1(log, r"d_H = ([\d.]+)")
            blob = _grep1(log, r"blob score\s*= ([\d.]+)")
            act = _grep1(log, r"active slices = (\d+)")
            L.append(f"| {b} | {n} | {dh} | {blob} | {act} |")
        try:
            d1 = float(_grep1(jobs["base_ladder"]["log"], r"d_H = ([\d.]+)"))
            d2 = float(_grep1(jobs["base_K80"]["log"], r"d_H = ([\d.]+)"))
            b1 = float(_grep1(jobs["base_ladder"]["log"],
                              r"blob score\s*= ([\d.]+)"))
            b2 = float(_grep1(jobs["base_K80"]["log"],
                              r"blob score\s*= ([\d.]+)"))
            L.append(f"\nVERDICT: d_H {'RISES' if d2 > d1 else 'does NOT rise'}"
                     f" ({d1} -> {d2}); blob "
                     f"{'RISES' if b2 > b1 else 'does NOT rise'}"
                     f" ({b1} -> {b2}).")
        except (TypeError, ValueError):
            L.append("\nVERDICT: incomplete (missing numbers in logs).")
    # --- closure compatibility -----------------------------------------------
    L.append("\n## Closure arms (matched-volume compatibility)\n")
    L.append("| arm | final N4 | N41 | d_H | collision fraction | mu |")
    L.append("|---|---|---|---|---|---|")
    n4s = {}
    for name in ("clo_b0", "clo_b1", "clo_plc", "clo_b1_K40"):
        if name not in jobs:
            continue
        log = jobs[name]["log"]
        n4 = _grep1(log, r"final N4=(\d+)")
        n41 = _grep1(log, r"final N4=\d+\s+N41=(\d+)")
        dh = _grep1(log, r"d_H = ([\d.]+)")
        coll = _grep1(log, r"tets with label collision: ([\d.]+)")
        mu = _grep1(log, r"mu = ([\d.]+) per cell") or \
            _grep1(log, r"reusing mu = ([\d.]+)")
        if n4 and name != "clo_b1_K40":
            n4s[name] = int(n4)
        L.append(f"| {name} | {n4} | {n41} | {dh} | {coll} | {mu} |")
    if len(n4s) >= 2:
        spread = (max(n4s.values()) - min(n4s.values())) / max(n4s.values())
        L.append(f"\nmatched-volume spread across K80 arms: {100 * spread:.1f}%"
                 f" ({'PASS (<2%)' if spread < 0.02 else 'CHECK'})")
    L.append("\nExpected ordering (paper 23.5): collision fraction "
             "b0 ~ 0.65 >> b1 (ordered), placebo stays high ~ 0.73; "
             "geometry statistically unchanged across arms.")
    # --- defect analysis -------------------------------------------------
    reals = [jobs[n]["csv"] for n in jobs if n.startswith("def_s")
             and n.endswith("_real") and st.get(n, {}).get("status") == "done"]
    plcs = [jobs[n]["csv"] for n in jobs if n.startswith("def_s")
            and n.endswith("_placebo")
            and st.get(n, {}).get("status") == "done"]
    L.append("\n## Defect experiment (pooled, autocorrelation-corrected)\n")
    if reals and plcs:
        out = _run_capture([PY, os.path.join(BUNDLE, "v6_defect_run.py"),
                            "--analyze"] + reals + plcs)
        L.append("```\n" + out.strip() + "\n```")
    else:
        L.append("(defect arms incomplete -- no analysis)")
    # --- capacity analysis -----------------------------------------------
    L.append("\n## Capacity conservation instrument (v1.2 gates)\n")
    if st.get("cap_v12", {}).get("status") == "done":
        out = _run_capture([PY, os.path.join(BUNDLE, "v6_capacity_run.py"),
                            "--analyze", os.path.join(d, "cap_v12")])
        L.append("```\n" + out.strip() + "\n```")
    else:
        L.append("(capacity run incomplete -- no analysis)")
    L.append("\n---\nnext steps: if the volume pair confirms the region, "
             "launch rung-2 (`rung2_boundary.py --predict` then the k0 scan) "
             "and the production bridge family (`run_bridge_sweep.sh`).")
    path = os.path.join(d, "REPORT.md")
    open(path, "w").write("\n".join(L) + "\n")
    print(f"\n# report -> {path}")
    return path


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="runs/p100k",
                    help="run directory (state, checkpoints, logs, report)")
    ap.add_argument("--jobs", type=int, default=4,
                    help="max concurrent processes (one core each)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--report-only", action="store_true")
    ap.add_argument("--retry-failed", action="store_true")
    ap.add_argument("--poll", type=float, default=20.0)
    # volumes / geometry
    ap.add_argument("--n41", type=int, default=100000)
    ap.add_argument("--n41-ladder", type=int, default=40000)
    ap.add_argument("--ladder", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="run the lower-volume rung for the volume-pair check")
    ap.add_argument("--K", type=int, default=80)
    ap.add_argument("--K-cap", type=int, default=40,
                    help="slice count for the capacity sibling (thicker "
                         "slices = longer usable far-field radius)")
    # tuned phase point (phase-grid peak; k4 is a starting guess, auto-tuned)
    ap.add_argument("--k0", type=float, default=2.0)
    ap.add_argument("--Delta", type=float, default=0.4)
    ap.add_argument("--k4", type=float, default=0.72)
    ap.add_argument("--tune-k4", type=int, default=3)
    ap.add_argument("--tune-burst", type=int, default=60)
    ap.add_argument("--eps-grow", type=float, default=1e-4)
    ap.add_argument("--eps-couple", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=0)
    # sweeps
    ap.add_argument("--sweeps-base", type=int, default=4000)
    ap.add_argument("--sweeps-closure", type=int, default=4000)
    ap.add_argument("--sweeps-defect", type=int, default=3000)
    ap.add_argument("--sweeps-cap", type=int, default=4000)
    ap.add_argument("--measure-every-base", type=int, default=50)
    ap.add_argument("--measure-every-closure", type=int, default=20)
    ap.add_argument("--measure-every-cap", type=int, default=5)
    # closure / defect / capacity knobs
    ap.add_argument("--lambda-inj", type=float, default=3.0)
    ap.add_argument("--audit-every", type=int, default=25)
    ap.add_argument("--seeds", type=int, default=2,
                    help="defect replication seeds (real+placebo per seed)")
    ap.add_argument("--pins", type=int, default=100)
    ap.add_argument("--pin-sep", type=int, default=5)
    ap.add_argument("--rmax-defect", type=int, default=3)
    ap.add_argument("--therm-defect", type=int, default=200)
    ap.add_argument("--pins-per-level", type=int, default=25)
    ap.add_argument("--pin-sep-cap", type=int, default=6)
    ap.add_argument("--rmax-cap", type=int, default=8,
                    help="capacity shell range (K40 at 100k reaches ~12)")
    ap.add_argument("--therm-cap", type=int, default=400)
    ap.add_argument("--persist", type=int, default=200)
    ap.add_argument("--smoke", action="store_true",
                    help="tiny end-to-end pipeline test (minutes, not days)")
    args = ap.parse_args()

    if args.smoke:
        args.n41, args.n41_ladder = 300, 150
        args.K, args.K_cap = 8, 6
        args.sweeps_base, args.sweeps_closure = 30, 30
        args.sweeps_defect, args.sweeps_cap = 40, 120
        args.measure_every_base = 10
        args.measure_every_closure = 10
        args.measure_every_cap = 2
        args.tune_k4, args.tune_burst = 1, 10
        args.seeds = 1
        args.pins, args.pin_sep, args.rmax_defect = 4, 4, 2
        args.therm_defect = 5
        args.pins_per_level, args.pin_sep_cap, args.rmax_cap = 2, 4, 3
        args.therm_cap, args.persist = 20, 10
        args.poll = 2.0

    args.dir = os.path.abspath(args.dir)
    os.makedirs(args.dir, exist_ok=True)
    jobs = build_jobs(args)

    if args.dry_run:
        print(f"# plan: {len(jobs)} jobs -> {args.dir}")
        for name, job in jobs.items():
            deps = ",".join(job["deps"]) or "-"
            print(f"\n[{job['stage']}] {name}  (deps: {deps})")
            print("  " + " ".join(job["cmd"]))
        return

    if not args.report_only:
        st = run_jobs(jobs, args)
    else:
        st = load_state(os.path.join(args.dir, "state.json"))
    write_report(jobs, st, args)
    bad = [n for n, s in st.items() if s.get("status") != "done"]
    if bad:
        print(f"# {len(bad)} job(s) not done: {', '.join(bad)}")
        print("# rerun the same command to resume; --retry-failed to clear "
              "failures first")
    else:
        print("# ALL JOBS DONE -- see REPORT.md")


if __name__ == "__main__":
    main()
