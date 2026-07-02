#!/usr/bin/env python3
"""Summary figure for a sweep of v6 runs (closure eta/beta sweeps, EPRL beta
sweeps, or bare volume scans).

Parses the final RESULT block of each run log, writes a CSV always, and -- if
matplotlib is installed (optional; the VM protocol needs only numpy/scipy) --
renders small-multiple panels of the reliable observables (d_H, blob score,
cos^3 fit error) against the swept variable. Placebo arms are drawn as their
own open markers, the beta=0 / bare control as a dashed reference band, and
the theory point (eta* for eta sweeps) as a vertical marker.

Design rules (kept deliberately boring): one panel per observable -- the three
observables have different scales and NEVER share a y-axis; one hue for the
real arms, one reserved hue for placebo arms; neutral dashed lines for
references; direct labels instead of legend boxes where they fit.

Usage:
  python plot_sweep.py logs/clo_*_20k.log --x eta --control logs/clo_b0.0_20k.log \
         --mark-x 0.0298668443935 --out eta_sweep_20k
  python plot_sweep.py logs/thy_b*_20k.log --x beta --out eprl_sweep_20k
  python plot_sweep.py scan_10k.log scan_20k.log scan_40k.log --x n41 --out volume_scan
"""
from __future__ import annotations
import argparse
import csv
import re
import sys

ETA_STAR = 0.0298668443935

# text/marks palette (validated reference palette; see repo audit docs)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e5e4e0"
SERIES = "#2a78d6"          # real arms
PLACEBO = "#e34948"         # placebo arms (reserved, never reused)


def parse_log(path):
    """Pull the swept parameters + final observables out of one run log."""
    txt = open(path, errors="replace").read()
    row = {"log": path}

    def grab(pat, cast=float, group=1):
        m = None
        for m in re.finditer(pat, txt):
            pass                                  # keep the LAST match (final block)
        return cast(m.group(group)) if m else None

    row["beta"] = grab(r"\bbeta(?:_eprl)?=([0-9.eE+-]+)")
    row["eta"] = grab(r"(?<![A-Za-z_])eta=([0-9.eE+-]+)")
    row["placebo"] = bool(re.search(r"placebo=True|PLACEBO", txt))
    row["d_H"] = grab(r"d_H = ([0-9.]+)")
    row["blob"] = grab(r"blob score\s*=\s*([0-9.]+)")
    row["cos3err"] = grab(r"rel\.?\s*RMS err\s*=?\s*([0-9.]+)")
    row["N4"] = grab(r"final N4=(\d+)", int)
    row["N41"] = grab(r"final N4=\d+\s+N41=(\d+)", int)
    row["rail2"] = grab(r"2t=([0-9.]+)")
    row["rail4"] = grab(r"4t=([0-9.]+)")
    m = re.search(r"(?:=>|foliation:)\s*(CLEAN|DEFECTIVE)", txt)
    row["foliation"] = m.group(1) if m else None
    row["mean_E"] = grab(r"<E>/tet=([0-9.]+)")
    row["coll_frac"] = grab(r"collision:\s*([0-9.]+)")
    return row


OBS = [("d_H", "d_H (Hausdorff, dual graph)"),
       ("blob", "blob score (max/mean of V3)"),
       ("cos3err", "cos³ fit rel. RMS error")]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("logs", nargs="+")
    ap.add_argument("--x", choices=["eta", "beta", "n41"], default="eta")
    ap.add_argument("--control", default=None,
                    help="log of the bare / beta=0 control run, drawn as a "
                         "dashed reference line in each panel")
    ap.add_argument("--mark-x", type=float, default=None,
                    help="vertical marker (e.g. eta* for eta sweeps)")
    ap.add_argument("--mark-label", default="η*")
    ap.add_argument("--logx", action="store_true",
                    help="log-scale x (natural for eta grids)")
    ap.add_argument("--out", default="sweep")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    xkey = {"eta": "eta", "beta": "beta", "n41": "N41"}[args.x]
    rows = [parse_log(p) for p in args.logs]
    rows = [r for r in rows if r.get(xkey) is not None and r.get("d_H") is not None]
    if not rows:
        sys.exit("no parseable finished runs among the given logs "
                 "(need the final RESULT block)")
    rows.sort(key=lambda r: (r[xkey], r["placebo"]))
    ctrl = parse_log(args.control) if args.control else None

    csv_path = args.out + ".csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
        if ctrl:
            w.writerow(ctrl)
    print(f"wrote {csv_path}  ({len(rows)} arms"
          f"{' + control' if ctrl else ''})")
    for r in rows + ([ctrl] if ctrl else []):
        flag = " PLACEBO" if r["placebo"] else ""
        print(f"  {xkey}={r[xkey]}{flag}: d_H={r['d_H']} blob={r['blob']} "
              f"cos3err={r['cos3err']} N4={r['N4']} foliation={r['foliation']}")
        if r["foliation"] == "DEFECTIVE":
            print("    !! DEFECTIVE foliation -- this arm is not readable")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed -- CSV only "
              "(pip install matplotlib for the figure)")
        return

    real = [r for r in rows if not r["placebo"]]
    plac = [r for r in rows if r["placebo"]]
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(7.0, 8.4), dpi=150)
    fig.patch.set_facecolor("#fcfcfb")
    xlabel = {"eta": "admissibility precision η",
              "beta": "coupling β",
              "n41": "target N41"}[args.x]

    for ax, (key, label) in zip(axes, OBS):
        ax.set_facecolor("#fcfcfb")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.grid(axis="y", color=GRID, linewidth=0.6)
        ax.tick_params(colors=INK2, labelsize=9)
        ax.set_ylabel(label, fontsize=9.5, color=INK)

        xs = [r[xkey] for r in real]
        ys = [r[key] for r in real]
        ax.plot(xs, ys, "-o", color=SERIES, linewidth=2, markersize=6.5,
                zorder=3)
        if plac:
            ax.plot([r[xkey] for r in plac], [r[key] for r in plac], "s",
                    markerfacecolor="none", markeredgecolor=PLACEBO,
                    markeredgewidth=2, markersize=8, zorder=3)
        if ctrl and ctrl.get(key) is not None:
            ax.axhline(ctrl[key], color=INK2, linestyle="--", linewidth=1.2,
                       zorder=2)
        if args.mark_x is not None:
            ax.axvline(args.mark_x, color=INK2, linestyle=":", linewidth=1.2,
                       zorder=1)
        # torus rails on the d_H panel only (context, very recessive)
        if key == "d_H" and real and real[0].get("rail2") and real[0].get("rail4"):
            for rk, rl in (("rail2", "2-torus rail"), ("rail4", "4-torus rail")):
                v = real[0][rk]
                ax.axhline(v, color=GRID, linewidth=1.0, zorder=1)
                ax.annotate(rl, xy=(1.0, v), xycoords=("axes fraction", "data"),
                            xytext=(-4, 3), textcoords="offset points",
                            ha="right", fontsize=8, color=INK2)
        if args.logx:
            ax.set_xscale("log")

    # direct labels once, on the top panel (identity never color-alone)
    top = axes[0]
    if real:
        top.annotate("real arms", xy=(real[0][xkey], real[0]["d_H"]),
                     xytext=(0, 10), textcoords="offset points", ha="left",
                     fontsize=9, color=SERIES)
    if plac:
        top.annotate("placebo", xy=(plac[-1][xkey], plac[-1]["d_H"]),
                     xytext=(8, -4), textcoords="offset points",
                     fontsize=9, color=PLACEBO)
    if ctrl and ctrl.get("d_H") is not None:
        top.annotate("β=0 control", xy=(0.02, ctrl["d_H"]),
                     xycoords=("axes fraction", "data"), xytext=(0, -12),
                     textcoords="offset points", ha="left",
                     fontsize=9, color=INK2)
    if args.mark_x is not None:
        top.annotate(args.mark_label, xy=(args.mark_x, 1.0),
                     xycoords=("data", "axes fraction"), xytext=(4, -2),
                     textcoords="offset points", fontsize=10, color=INK2)

    axes[-1].set_xlabel(xlabel, fontsize=10, color=INK)
    ttl = args.title or f"{xlabel} sweep — matched-volume observables"
    import textwrap
    fig.suptitle("\n".join(textwrap.wrap(ttl, 68)), fontsize=11, color=INK,
                 y=0.998, va="top")
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    png_path = args.out + ".png"
    fig.savefig(png_path, facecolor=fig.get_facecolor())
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
