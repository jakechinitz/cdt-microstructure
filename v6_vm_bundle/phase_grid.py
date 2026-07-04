#!/usr/bin/env python3
"""Locate phase C in THIS engine's conventions: a coarse (k0, Delta) grid of
short bare runs at fixed pinned volume, with per-point phase diagnostics.

WHY: the borrowed literature point (k0=2.2, Delta=0.6) produced a blob score
that FALLS with volume on the fixed-foliation engine (10k/20k/40k rescan) --
the noise-peak signature -- while d_H and cos3err trended correctly. Phase
coordinates do not transfer between CDT codes; the de Sitter phase C must be
located in this engine's own conventions.

PER-POINT DIAGNOSTICS (all from the final thermalized stretch):
  blob        max/mean of V3(tau)          -- C large, A ~1+noise, B ~K
  active      slices clearly above stalk   -- C: several, B: ~1
  slicecorr   Pearson corr of adjacent V3(tau),V3(tau+1) -- A ~0, C >> 0
  maxfrac     largest slice / total        -- B-collapse detector
  d_H         dual-graph Hausdorff
  n0n4        N0/N4 (the standard k0 order parameter)
  pin         eps-pin residual  2*eps*<N41-target>  -- large |pin| means the
              point is far from pseudo-critical k4 and is being held by the
              pin, not by the action; treat its diagnostics with suspicion.

PHASE-C CANDIDATE = high slicecorr + blob well above 1 + active >= 3 +
maxfrac small + d_H high. The grid does NOT decide by itself: a single-volume
blob cannot distinguish signal from a noise peak. ACCEPTANCE for any
candidate point is the volume-pair check -- blob must RISE with volume:

    python v6_verify_run.py --k0 <K0> --Delta <D> --target-n41 10000 ... ckptA
    python v6_verify_run.py --k0 <K0> --Delta <D> --target-n41 20000 ... ckptB
    python plot_sweep.py logA logB --x n41 --out confirm_<K0>_<D>

USAGE:
  python phase_grid.py --k0 "1.0 1.5 2.0 2.5 3.0 3.5" --Delta "0.2 0.4 0.6 0.8" \
         --target-n41 10000 --K 80 --sweeps 400 --jobs 4 --out grid10k
  python phase_grid.py --plot grid10k.csv          # heatmaps (needs matplotlib)

Results append to <out>.csv as each point finishes; rerunning skips points
already present (crash/interrupt-safe; delete the row to redo a point).
"""
from __future__ import annotations
import argparse
import csv
import os
import sys
import time

import numpy as np


# ---------------------------------------------------------------------------
# One grid point = one short bare run + diagnostics
# ---------------------------------------------------------------------------

def slice_corr(prof):
    """Pearson correlation of adjacent slice volumes (cyclic). Phase A
    (uncorrelated slices) ~ 0; phase C (coherent blob) strongly positive."""
    a = np.asarray(prof, float)
    if len(a) < 4 or a.std() == 0:
        return 0.0
    b = np.roll(a, 1)
    return float(np.corrcoef(a, b)[0, 1])


def run_point(args_tuple):
    (k0, Delta, k4, eps, target, K, sweeps, seed, tune_rounds) = args_tuple
    # import inside the worker (multiprocessing-safe)
    from v6_run_lib import run, dual_adjacency, hausdorff_dim, volume_profile, \
        profile_metrics
    import os
    import tempfile
    t0 = time.time()
    # --- optional k4 auto-tuning to pseudo-criticality ----------------------
    # Real CDT tunes k4 per (k0, Delta) point; a fixed k4 makes the eps pin
    # fight real volume pressure in parts of the grid (visible as a large
    # |pin| residual) and those cells are not fair tests. The pin force IS
    # the mistuning signal, so use it as the feedback error:
    #     k4  <-  k4 + 2*eps*(N41 - target)
    # over a few short bursts before the measurement run.
    ckpt = None
    if tune_rounds > 0:
        fd, ckpt = tempfile.mkstemp(suffix=".json")
        os.close(fd)
        burst = max(20, sweeps // (2 * tune_rounds))
        resume = None
        for _ in range(int(tune_rounds)):
            T = run(f"tune k0={k0} D={Delta}", k0=k0, Delta=Delta, k4=k4,
                    target_N41=target, K=K, eps=eps, seed=seed,
                    max_sweeps=burst, measure_every=burst, checkpoint=ckpt,
                    resume=resume, verbose=False)
            resume = ckpt
            k4 = k4 + 2 * eps * (T.type_counts()[0] - target)
    T = run(f"grid k0={k0} D={Delta}", k0=k0, Delta=Delta, k4=k4,
            target_N41=target, K=K, eps=eps, seed=seed, max_sweeps=sweeps,
            measure_every=max(1, sweeps // 4), checkpoint=None,
            resume=ckpt, verbose=False)
    if ckpt:
        try:
            os.unlink(ckpt)
        except OSError:
            pass
    ids, adj = dual_adjacency(T)
    dH = hausdorff_dim(adj)
    prof = volume_profile(T)
    pm = profile_metrics(prof)
    a = np.asarray(prof, float)
    n41 = T.type_counts()[0]
    okf, repf = getattr(T, "_final_verify", (None, {}))
    sl_bad, _, sl_ss = getattr(T, "_final_slices", (None, {}, None))
    return {
        "k0": k0, "Delta": Delta, "k4": round(k4, 4), "K": K,
        "target_n41": target,
        "sweeps": sweeps, "seed": seed,
        "N4": T.n_pent(), "N41": n41,
        "n0n4": round(len(T.vinc) / max(1, T.n_pent()), 4),
        "d_H": round(dH, 3),
        "blob": round(pm["blob_score"], 3),
        "active": pm["active_slices"],
        "slicecorr": round(slice_corr(prof), 3),
        "maxfrac": round(float(a.max() / max(1.0, a.sum())), 4),
        "cos3err": (round(pm["cos3_relerr"], 3)
                    if pm["cos3_relerr"] is not None else ""),
        "pin": round(2 * eps * (n41 - target), 4),
        "valid": bool(okf),
        "foliation": "CLEAN" if (sl_bad == 0 and sl_ss == 0) else "DEFECTIVE",
        "wall_s": int(time.time() - t0),
    }


FIELDS = ["k0", "Delta", "k4", "K", "target_n41", "sweeps", "seed", "N4",
          "N41", "n0n4", "d_H", "blob", "active", "slicecorr", "maxfrac",
          "cos3err", "pin", "valid", "foliation", "wall_s"]


# ---------------------------------------------------------------------------
# Heatmap rendering (matplotlib optional)
# ---------------------------------------------------------------------------

def plot_grid(csv_path):
    rows = list(csv.DictReader(open(csv_path)))
    if not rows:
        sys.exit(f"{csv_path} is empty")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("matplotlib not installed (pip install matplotlib)")
    k0s = sorted({float(r["k0"]) for r in rows})
    Ds = sorted({float(r["Delta"]) for r in rows})
    panels = [("slicecorr", "adjacent-slice correlation", "higher = phase C"),
              ("blob", "blob score (max/mean)", "single-volume; confirm by pair"),
              ("d_H", "d_H (dual graph)", ""),
              ("maxfrac", "max slice fraction", "high = phase-B collapse")]
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.6), dpi=150)
    fig.patch.set_facecolor("#fcfcfb")
    for ax, (key, label, note) in zip(axes.ravel(), panels):
        M = np.full((len(Ds), len(k0s)), np.nan)
        for r in rows:
            try:
                M[Ds.index(float(r["Delta"])), k0s.index(float(r["k0"]))] = \
                    float(r[key])
            except (ValueError, KeyError):
                pass
        im = ax.imshow(M, origin="lower", aspect="auto", cmap="Blues")
        ax.set_xticks(range(len(k0s)), [f"{v:g}" for v in k0s], fontsize=8)
        ax.set_yticks(range(len(Ds)), [f"{v:g}" for v in Ds], fontsize=8)
        ax.set_xlabel("k0", fontsize=9)
        ax.set_ylabel("Delta", fontsize=9)
        ttl = label + (f"\n({note})" if note else "")
        ax.set_title(ttl, fontsize=9.5, color="#0b0b0b")
        for i in range(len(Ds)):
            for j in range(len(k0s)):
                if not np.isnan(M[i, j]):
                    ax.text(j, i, f"{M[i, j]:g}", ha="center", va="center",
                            fontsize=7.5, color="#0b0b0b")
        fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle(f"bare-engine phase grid — {os.path.basename(csv_path)}",
                 fontsize=11, color="#0b0b0b")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.splitext(csv_path)[0] + ".png"
    fig.savefig(out, facecolor=fig.get_facecolor())
    print(f"wrote {out}")


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--k0", default="1.0 1.5 2.0 2.5 3.0 3.5")
    ap.add_argument("--Delta", default="0.2 0.4 0.6 0.8")
    ap.add_argument("--k4", type=float, default=0.9)
    ap.add_argument("--eps", type=float, default=1e-3)
    ap.add_argument("--target-n41", type=int, default=10000)
    ap.add_argument("--K", type=int, default=80)
    ap.add_argument("--sweeps", type=int, default=400,
                    help="sweeps per point (growth + a short thermal tail)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=1,
                    help="grid points run in parallel (one core each)")
    ap.add_argument("--tune-k4", type=int, default=0, metavar="ROUNDS",
                    help="auto-tune k4 to pseudo-criticality per point over "
                         "N short bursts before measuring (recommended: 3). "
                         "Uses the eps-pin residual as the feedback error; "
                         "the tuned k4 is recorded in the CSV.")
    ap.add_argument("--out", default="phase_grid")
    ap.add_argument("--plot", default=None, metavar="CSV",
                    help="render heatmaps from an existing results CSV and exit")
    args = ap.parse_args()

    if args.plot:
        plot_grid(args.plot)
        return

    csv_path = args.out + ".csv"
    done = set()
    if os.path.exists(csv_path):
        for r in csv.DictReader(open(csv_path)):
            done.add((float(r["k0"]), float(r["Delta"])))
        print(f"# resuming: {len(done)} points already in {csv_path}")

    points = [(float(a), float(b))
              for a in args.k0.split() for b in args.Delta.split()
              if (float(a), float(b)) not in done]
    print(f"# phase grid: {len(points)} points to run "
          f"(N41={args.target_n41}, K={args.K}, {args.sweeps} sweeps, "
          f"jobs={args.jobs})")
    if not points:
        print("# nothing to do -- all points present; use --plot to render")
        return

    tasks = [(k0, D, args.k4, args.eps, args.target_n41, args.K,
              args.sweeps, args.seed, args.tune_k4) for k0, D in points]
    new_file = not os.path.exists(csv_path)
    f = open(csv_path, "a", newline="")
    w = csv.DictWriter(f, fieldnames=FIELDS)
    if new_file:
        w.writeheader()
        f.flush()

    def record(res):
        w.writerow(res)
        f.flush()
        print(f"  k0={res['k0']:g} D={res['Delta']:g}: blob={res['blob']} "
              f"active={res['active']} slicecorr={res['slicecorr']} "
              f"d_H={res['d_H']} maxfrac={res['maxfrac']} pin={res['pin']} "
              f"foliation={res['foliation']} ({res['wall_s']}s)", flush=True)

    if args.jobs > 1:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = [ex.submit(run_point, t) for t in tasks]
            for fut in as_completed(futs):
                record(fut.result())
    else:
        for t in tasks:
            record(run_point(t))
    f.close()
    print(f"# done -> {csv_path}   (render: python phase_grid.py --plot {csv_path})")
    print("# REMEMBER: a good-looking point is only a CANDIDATE until the "
          "volume-pair check passes (blob must RISE with volume).")


if __name__ == "__main__":
    main()
