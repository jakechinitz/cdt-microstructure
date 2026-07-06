#!/usr/bin/env python3
"""RUNG 2 -- the curvature-sector test: does the substrate displace the
host's phase boundary by the precomputed amount?

THE TEST. The ring transfer matrix (induced_couplings.py, App. J.4) prices
the closure weighting's induced curvature-sector coupling at
c0 = +0.0187 per slice edge at the theory point. Expressed against the CDT
action's -(k0+6*Delta)*N0 term via the measured slice edge-per-vertex ratio
r = N1/N0, the coupled ensemble at bare k0 behaves like the bare ensemble
at k0 - c0*r: the phase boundary, measured in bare-k0 coordinates, must
shift UP by

    delta_k0(real) = + c0_real * r   ~ +0.12   (r ~ 6.5)

The orbit-shuffled placebo table carries c0_plc = +0.0003 -- SEVENTY times
smaller -- despite a LARGER volume-sector coupling (c1 = 1.20 vs 0.53).
The shuffle preserves the energy pool and destroys closure structure, and
only structure prices curvature. So the pre-registered pair is:

    delta_k0(placebo) = + c0_plc * r ~ +0.002  (consistent with zero)

A real-arm shift at the predicted size with a placebo shift at zero is the
first measurement in this program of the substrate steering the HOST's
phase structure -- the curvature sector, where geometry's stiffness lives
-- and it dissociates closure STRUCTURE from energy MAGNITUDE (the bridge
family tested magnitude; this tests structure).

DESIGN. Fine k0 scan across the boundary the tuned grid located (coherent
to k0 ~ 2.0-2.5 at Delta = 0.4, dead by 3.0-3.5). Per (k0, seed): k4 is
auto-tuned on the BARE system (short bursts, pin residual as feedback),
then three arms run from the SAME tuned checkpoint with the same k4 --
bare (beta=0), real (beta=1, centered), placebo (beta=1, centered,
shuffled table). Centering removes the volume-sector tilt (J.2/J.6), so
what remains to move the boundary is the curvature sector. The order
parameter is the slice-volume coherence (slicecorr) that mapped the
boundary; the crossing point k0* is interpolated per arm at the midpoint
of its plateaus, and the verdict is k0*(arm) - k0*(bare) against the
predictions above, with r measured from the runs' own configurations.

USAGE:
  python rung2_boundary.py --predict                 # write predictions
  python rung2_boundary.py --k0-list "2.0 2.2 ... 3.4" --jobs 4 \
         --out rung2.csv                             # launch scan (resumable)
  python rung2_boundary.py --analyze rung2.csv       # crossings + verdict
"""
from __future__ import annotations
import argparse
import csv
import os
import sys
import time
from collections import defaultdict
from itertools import combinations

import numpy as np

ETA_STAR = 0.0298668443935
PLACEBO_SEED = 12345
ARMS = ("bare", "real", "plc")

FIELDS = ["k0", "Delta", "arm", "seed", "k4_tuned", "beta", "mu",
          "N4", "N41", "slicecorr", "d_H", "blob", "active", "n1n0",
          "pin", "foliation", "wall_s"]


# ---------------------------------------------------------------------------
# Predictions (pure math -- on record before any lattice runs)
# ---------------------------------------------------------------------------

def predict(eta, lam, out):
    from induced_couplings import ring_analysis
    import induced_couplings as ic
    import v6_closure_run as vc
    real = ring_analysis(1.0, eta, lam, verbose=False)
    orig = ic.build_energy_table
    ic.build_energy_table = lambda e, l: vc.build_energy_table(e, l, PLACEBO_SEED)
    plc = ring_analysis(1.0, eta, lam, verbose=False)
    ic.build_energy_table = orig
    print(f"# rung-2 predictions (ring transfer matrix, beta=1, "
          f"eta={eta:.6g}, lambda={lam})")
    print(f"#   c0_real = {real['c0']:+.4f} per slice edge  "
          f"(c1 = {real['c1']:+.4f})")
    print(f"#   c0_plc  = {plc['c0']:+.4f} per slice edge  "
          f"(c1 = {plc['c1']:+.4f})   ratio {plc['c0']/real['c0']:.3f}")
    print(f"#   delta_k0(arm) = + c0_arm * r,  r = N1/N0 measured from the "
          f"runs themselves (r ~ 6.5 -> real +{real['c0']*6.5:.3f}, "
          f"placebo +{plc['c0']*6.5:.3f})")
    if out:
        with open(out, "w") as f:
            f.write("table,c0,c1\n")
            f.write(f"real,{real['c0']:.6f},{real['c1']:.6f}\n")
            f.write(f"placebo,{plc['c0']:.6f},{plc['c1']:.6f}\n")
        print(f"# wrote {out}")
    return real["c0"], plc["c0"]


# ---------------------------------------------------------------------------
# One (k0, seed) point: tune k4 bare, then three arms from the tuned state
# ---------------------------------------------------------------------------

def slice_ratio(T):
    """Measured N1/N0 on the spatial slices (edges per vertex)."""
    from v6_closure_run import tet_of
    verts, edges = set(), set()
    for vs in T.pent.values():
        tet = tet_of(T, vs)
        if tet is None:
            continue
        tt = sorted(tet)
        verts.update(tt)
        edges.update(combinations(tt, 2))
    return len(edges) / max(1, len(verts))


def run_point(task):
    (k0, Delta, k4, eps, target, K, sweeps, tune_rounds, seed,
     eta, lam, beta) = task
    import tempfile
    from v6_run_lib import (run, load_checkpoint, dual_adjacency,
                            hausdorff_dim, volume_profile, profile_metrics)
    from phase_grid import slice_corr
    from v6_closure_run import (build_energy_table, FaceLabels,
                                IncrementalClosure, make_heatbath,
                                calibrate_mu_ti)
    from v6_theory_run import Centering

    t0 = time.time()
    tmp = []

    def mktmp():
        fd, p = tempfile.mkstemp(suffix=".json")
        os.close(fd)
        tmp.append(p)
        return p

    # --- bare k4 tuning, shared by all three arms (phase_grid protocol) ----
    tuned = mktmp()
    k4t = k4
    resume = None
    if tune_rounds > 0:
        burst = max(20, sweeps // (2 * tune_rounds))
        for _ in range(int(tune_rounds)):
            T = run(f"tune k0={k0}", k0=k0, Delta=Delta, k4=k4t,
                    target_N41=target, K=K, eps=eps, seed=seed,
                    max_sweeps=burst, measure_every=burst, checkpoint=tuned,
                    resume=resume, verbose=False)
            resume = tuned
            k4t = k4t + 2 * eps * (T.type_counts()[0] - target)
    else:
        T = run(f"grow k0={k0}", k0=k0, Delta=Delta, k4=k4t,
                target_N41=target, K=K, eps=eps, seed=seed,
                max_sweeps=max(20, sweeps // 4),
                measure_every=sweeps, checkpoint=tuned, verbose=False)

    rows = []
    for arm in ARMS:
        ck = mktmp()
        mu = None
        if arm == "bare":
            T = run(f"r2 {arm} k0={k0}", k0=k0, Delta=Delta, k4=k4t,
                    target_N41=target, K=K, eps=eps, seed=seed + 7,
                    max_sweeps=sweeps, measure_every=max(1, sweeps // 4),
                    checkpoint=ck, resume=tuned, verbose=False)
        else:
            Etab = build_energy_table(
                eta, lam, PLACEBO_SEED if arm == "plc" else None)
            rng = np.random.default_rng(seed + 7)
            labels = FaceLabels(rng)
            T0, _, _ = load_checkpoint(tuned)
            mu = calibrate_mu_ti(T0, Etab, beta, seed=seed + 777,
                                 verbose=False)
            del T0
            centering = Centering(enabled=True, mu_fixed=mu, label="closure")
            labels.mu_saved = mu
            inc = IncrementalClosure(labels, Etab, k0, Delta, k4t, beta,
                                     centering=centering)
            T = run(f"r2 {arm} k0={k0}", k0=k0, Delta=Delta, k4=k4t,
                    target_N41=target, K=K, eps=eps, seed=seed + 7,
                    max_sweeps=sweeps, measure_every=max(1, sweeps // 4),
                    checkpoint=ck, resume=tuned, extra_state=labels,
                    delta_action=inc, extra_hook=make_heatbath(inc),
                    audit_every=max(25, sweeps // 8), verbose=False)
        ids, adj = dual_adjacency(T)
        prof = volume_profile(T)
        pm = profile_metrics(prof)
        n41 = T.type_counts()[0]
        sl_bad, _, sl_ss = getattr(T, "_final_slices", (None, {}, None))
        rows.append({
            "k0": k0, "Delta": Delta, "arm": arm, "seed": seed,
            "k4_tuned": round(k4t, 4),
            "beta": 0.0 if arm == "bare" else beta,
            "mu": (round(mu, 4) if mu is not None else ""),
            "N4": T.n_pent(), "N41": n41,
            "slicecorr": round(slice_corr(prof), 4),
            "d_H": round(hausdorff_dim(adj), 3),
            "blob": round(pm["blob_score"], 3),
            "active": pm["active_slices"],
            "n1n0": round(slice_ratio(T), 4),
            "pin": round(2 * eps * (n41 - target), 4),
            "foliation": "CLEAN" if (sl_bad == 0 and sl_ss == 0) else "DEFECTIVE",
            "wall_s": int(time.time() - t0),
        })
    for p in tmp:
        try:
            os.unlink(p)
        except OSError:
            pass
    return rows


# ---------------------------------------------------------------------------
# Analysis: boundary crossing per arm, displacement vs prediction
# ---------------------------------------------------------------------------

def crossing(k0s, vals):
    """k0 where vals crosses the midpoint of its plateaus, descending."""
    k0s, vals = np.asarray(k0s, float), np.asarray(vals, float)
    order = np.argsort(k0s)
    k0s, vals = k0s[order], vals[order]
    thr = (vals.max() + vals.min()) / 2
    for i in range(len(k0s) - 1):
        a, b = vals[i], vals[i + 1]
        if (a >= thr) and (b < thr):
            return float(k0s[i] + (a - thr) / (a - b) * (k0s[i + 1] - k0s[i]))
    return None


def analyze(path, pred_csv=None):
    rows = list(csv.DictReader(open(path)))
    if not rows:
        sys.exit("no rows")
    preds = {}
    if pred_csv and os.path.exists(pred_csv):
        for r in csv.DictReader(open(pred_csv)):
            preds[r["table"]] = float(r["c0"])
    by = defaultdict(list)              # (arm, k0) -> slicecorr values
    seeds = defaultdict(set)
    rbar = []
    for r in rows:
        by[(r["arm"], float(r["k0"]))].append(float(r["slicecorr"]))
        seeds[r["arm"]].add(r["seed"])
        rbar.append(float(r["n1n0"]))
    rbar = float(np.mean(rbar))
    print(f"# measured slice edge/vertex ratio r = N1/N0 = {rbar:.3f}")
    print(f"{'arm':>6} {'k0*':>8}   (order parameter: slicecorr, "
          f"midpoint crossing)")
    stars = {}
    for arm in ARMS:
        ks = sorted({k for (a, k) in by if a == arm})
        if len(ks) < 3:
            print(f"{arm:>6} {'--':>8}   (needs >=3 k0 points)")
            continue
        vals = [float(np.mean(by[(arm, k)])) for k in ks]
        st = crossing(ks, vals)
        stars[arm] = st
        prof = "  ".join(f"{k:g}:{v:.3f}" for k, v in zip(ks, vals))
        print(f"{arm:>6} {(f'{st:.3f}' if st else 'none'):>8}   [{prof}]")
    if "bare" in stars and stars["bare"] is not None:
        for arm, tbl in (("real", "real"), ("plc", "placebo")):
            if stars.get(arm) is None:
                continue
            dk = stars[arm] - stars["bare"]
            p = preds.get(tbl)
            ptxt = f"   predicted {p*rbar:+.3f}" if p is not None else ""
            print(f"delta_k0({arm}) = {dk:+.3f}{ptxt}")
        print("# gates: real shift at the predicted size AND placebo shift "
              "consistent with zero -- structure, not magnitude, moves the "
              "boundary. Errors: rerun with more --seeds and jackknife; "
              "crossing resolution is set by the k0 step and run length.")


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predict", action="store_true",
                    help="compute + write the c0 predictions and exit")
    ap.add_argument("--analyze", default=None, metavar="CSV")
    ap.add_argument("--pred", default="rung2_pred.csv",
                    help="prediction CSV (written by --predict, read by "
                         "--analyze)")
    ap.add_argument("--k0-list", default="2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4")
    ap.add_argument("--Delta", type=float, default=0.4)
    ap.add_argument("--k4", type=float, default=0.72,
                    help="k4 starting guess (tuned per point; the tuned "
                         "grid found 0.67-0.75 in the coherent region)")
    ap.add_argument("--eps", type=float, default=1e-3)
    ap.add_argument("--target-n41", type=int, default=10000)
    ap.add_argument("--K", type=int, default=80)
    ap.add_argument("--sweeps", type=int, default=800,
                    help="per arm, after tuning (the grid mapped at 400; "
                         "the crossing needs better statistics)")
    ap.add_argument("--tune-rounds", type=int, default=3)
    ap.add_argument("--seeds", type=int, default=1,
                    help="independent seeds per k0 (>=3 for jackknife errors)")
    ap.add_argument("--eta", type=float, default=ETA_STAR)
    ap.add_argument("--lambda-inj", type=float, default=3.0)
    ap.add_argument("--beta", type=float, default=1.0)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--out", default="rung2.csv")
    args = ap.parse_args()

    if args.analyze:
        analyze(args.analyze, args.pred)
        return
    if args.predict:
        predict(args.eta, args.lambda_inj, args.pred)
        return

    predict(args.eta, args.lambda_inj, args.pred)   # always on record first

    done = set()
    if os.path.exists(args.out):
        for r in csv.DictReader(open(args.out)):
            done.add((float(r["k0"]), r["arm"], int(r["seed"])))
        print(f"# resuming: {len(done)} rows already in {args.out}")
    tasks = []
    for k0 in [float(x) for x in args.k0_list.split()]:
        for s in range(args.seeds):
            if all((k0, a, s) in done for a in ARMS):
                continue
            tasks.append((k0, args.Delta, args.k4, args.eps, args.target_n41,
                          args.K, args.sweeps, args.tune_rounds, s,
                          args.eta, args.lambda_inj, args.beta))
    print(f"# rung-2 scan: {len(tasks)} (k0, seed) points x 3 arms, "
          f"jobs={args.jobs}, out={args.out}")

    new = not os.path.exists(args.out)
    f = open(args.out, "a", newline="")
    w = csv.DictWriter(f, fieldnames=FIELDS)
    if new:
        w.writeheader()

    def emit(rows):
        for r in rows:
            w.writerow(r)
        f.flush()
        r0 = rows[0]
        print(f"# done k0={r0['k0']} seed={r0['seed']}  "
              + "  ".join(f"{r['arm']}:{r['slicecorr']}" for r in rows)
              + f"  (k4={r0['k4_tuned']}, {rows[-1]['wall_s']}s)", flush=True)

    if args.jobs > 1:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = [ex.submit(run_point, t) for t in tasks]
            for fut in as_completed(futs):
                emit(fut.result())
    else:
        for t in tasks:
            emit(run_point(t))
    f.close()
    print(f"# scan complete. analyze:  python rung2_boundary.py "
          f"--analyze {args.out}")


if __name__ == "__main__":
    main()
