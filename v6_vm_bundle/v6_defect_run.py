#!/usr/bin/env python3
"""STAGE 3 -- defect response: does geometry (and the surrounding capacity
state) respond to persistent closure failure, the theory's definition of mass?

THE EXPERIMENT. Pin N well-separated spatial tetrahedra as DEFECTS: freeze
their four face labels into a maximally-failing configuration (all faces the
same label: 6 collisions, K^2 = 48) and forbid geometry moves that would
destroy the carrier tet. The pinned labels are shared with neighbor cells, so
the failure couples outward; everything else stays fully dynamical. Then
measure, sweep by sweep, shell-averaged observables around each defect:

    mean_E     mean closure energy of cells at graph distance d (capacity
               strain -- the theory's gravity-analog field)
    coll       fraction of cells at distance d with a label collision
    mean_q     mean slice-edge coordination in the shell (= local discrete
               CURVATURE: deficit angle ~ 2*pi - q*theta)
    shell_n    number of cells in the shell (local volume growth)

THREE-ARM DESIGN (all measured identically, in the same run where possible):
  * real pins    -- failing frozen labels (the defect).
  * placebo pins -- frozen labels in a BEST-CLOSED injective configuration
                    (K^2 = 40.67, zero collisions). Identical anchoring and
                    move-blocking footprint, no failure content. Anything the
                    real pins do beyond the placebo is attributable to the
                    CLOSURE FAILURE itself, not to pinning.
  * virtual refs -- freshly sampled unpinned cells far from all pins, every
                    measurement. The vacuum baseline.

Run the real and placebo arms as separate processes from the same thermalized
checkpoint (they cannot share a lattice); virtual references are measured
inside both.

MECHANICS. Pin protection is folded into the incremental action: a move whose
reconciliation removes a pinned carrier returns dS = +1e30, so Metropolis
rejects it and the standard undo/reconcile path restores state (audit-safe).
The heat-bath skips pinned triangles. mu is resolved exactly as in the
closure run (pins are ~1% of cells; their fixed contribution is a constant).

READING THE RESULT (write-up discipline):
  * Capacity response = mean_E / coll shell profiles of REAL separated from
    PLACEBO and decaying toward VIRTUAL with distance.
  * Geometric response = mean_q / shell_n separation (the stronger claim:
    local geometry deforms around mass).
  * Expect SHORT RANGE: the measured label correlation length is under one
    lattice step (induced_couplings.py), so shells 1-2 carry the signal.
    This is the microscopic seed of the paper's source->strain mechanism,
    NOT Newton's law; no 1/r claim is available at these scales.

USAGE:
  # from a thermalized closure checkpoint (beta=1 arm of stage 2):
  python v6_defect_run.py --resume clo_b1.0_20k.json --pin-mode fail \
         --pins 100 --target-n41 20000 --K 80 --sweeps 3000 --out def_real
  python v6_defect_run.py --resume clo_b1.0_20k.json --pin-mode closed \
         --pins 100 --target-n41 20000 --K 80 --sweeps 3000 --out def_placebo
  # analyze both:
  python v6_defect_run.py --analyze def_real.csv def_placebo.csv
"""
from __future__ import annotations
import argparse
import csv
import sys
from collections import defaultdict

import numpy as np

from v6_run_lib import run, load_checkpoint
from v6_theory_run import Centering
from v6_closure_run import (ETA_STAR, build_energy_table, tet_of, tris_of,
                            FaceLabels, IncrementalClosure, calibrate_mu_ti,
                            _heatbath_pass, _tri_map)

# frozen label configurations (labels 0..6 <-> m = l-3)
FAIL_LABELS = (3, 3, 3, 3)        # all m=0: 6 collisions, K^2 = 48
CLOSED_LABELS = (3, 4, 5, 6)      # m = {0,1,2,3}: injective, K^2 = 40.67 (min)


class PinnedClosure(IncrementalClosure):
    """IncrementalClosure + carrier protection: any move that removes a
    pinned defect tet is made unacceptable (dS = +1e30 -> Metropolis rejects
    -> normal undo/reconcile restores the cache)."""
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.pins = []            # list of tet frozensets

    def delta(self, T, added, removed):
        dS = super().delta(T, added, removed)
        for tet in self.pins:
            if tet not in self.tet_pids:
                return 1e30       # carrier destroyed -> forbid
        return dS


def select_anchors(inc, n, min_dist, rng, exclude=(), exclude_dist=0):
    """Greedy selection of up to n mutually well-separated tets (BFS distance
    in the slice tet-adjacency graph >= min_dist), avoiding everything within
    exclude_dist of the tets in `exclude`."""
    adj = _adjacency(inc)
    forbidden = set()
    for e in exclude:
        forbidden |= _ball(adj, e, exclude_dist)
    tets = [t for t in inc.tet_pids if t not in forbidden]
    order = rng.permutation(len(tets))
    chosen = []
    for i in order:
        if len(chosen) >= n:
            break
        cand = tets[i]
        near = _ball(adj, cand, max(0, min_dist - 1))
        if any(c in near for c in chosen):
            continue
        chosen.append(cand)
    return chosen


def _adjacency(inc):
    adj = defaultdict(set)
    for tri, tets in _tri_map(inc.tet_pids).items():
        for a in tets:
            for b in tets:
                if a != b:
                    adj[a].add(b)
    return adj


def _ball(adj, start, r):
    seen = {start}
    frontier = [start]
    for _ in range(r):
        nxt = []
        for u in frontier:
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    nxt.append(w)
        frontier = nxt
    return seen


def _bfs_dist(adj, a, b, cap):
    if a == b:
        return 0
    seen = {a}
    frontier = [a]
    for d in range(1, cap + 1):
        nxt = []
        for u in frontier:
            for w in adj[u]:
                if w == b:
                    return d
                if w not in seen:
                    seen.add(w)
                    nxt.append(w)
        frontier = nxt
    return cap + 1


def apply_pins(inc, pins, labels4):
    for tet in pins:
        for tri, l in zip(tris_of(tet), labels4):
            inc.labels.lab[tri] = int(l)
            inc.labels.pinned.add(tri)
    inc.pins = list(pins)


def measure_shells(inc, anchors, kind, sweep, writer, rmax=3,
                   exclude_cells=frozenset()):
    """Shell-averaged observables around each anchor; one CSV row per
    (anchor-average, shell). `exclude_cells` (the pins themselves) are
    traversed for distance but never counted -- they are sources, not
    medium; without this, pins separated by exactly `shell` contaminate
    each other's outer-shell statistics."""
    adj = _adjacency(inc)
    # slice-edge coordination map (edge -> #tets containing it), whole slice
    ecount = defaultdict(int)
    from itertools import combinations
    for tet in inc.tet_pids:
        for e in combinations(sorted(tet), 2):
            ecount[e] += 1
    mu = inc.centering.value()
    acc = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0, 0])  # shell -> sums
    for a in anchors:
        if a not in inc.tet_pids:
            continue                       # (virtual ref destroyed)
        seen = {a}
        frontier = [a]
        for d in range(1, rmax + 1):
            nxt = []
            for u in frontier:
                for w in adj[u]:
                    if w not in seen:
                        seen.add(w)
                        nxt.append(w)
            s = acc[d]
            for c in nxt:
                if c in exclude_cells:
                    continue
                l = [inc.labels.lab[t2] for t2 in tris_of(c)]
                s[0] += float(inc.E[l[0], l[1], l[2], l[3]])
                s[1] += 1.0 if len(set(l)) < 4 else 0.0
                qs = [ecount[e] for e in combinations(sorted(c), 2)]
                s[2] += float(np.mean(qs))
                s[3] += 1.0
            s[4] += 1                      # anchors contributing
            frontier = nxt
    for d, s in sorted(acc.items()):
        if s[3] > 0:
            writer.writerow({
                "sweep": sweep, "kind": kind, "shell": d,
                "n_cells": s[3], "n_anchors": s[4],
                "mean_E": round(s[0] / s[3], 5),
                "coll": round(s[1] / s[3], 5),
                "mean_q": round(s[2] / s[3], 5),
                "shell_n": round(s[3] / max(1, s[4]), 4),
            })


CSV_FIELDS = ["sweep", "kind", "shell", "n_cells", "n_anchors", "mean_E",
              "coll", "mean_q", "shell_n"]


def analyze(paths):
    rows = []
    for p in paths:
        rows += list(csv.DictReader(open(p)))
    if not rows:
        sys.exit("no rows")
    agg = defaultdict(list)
    for r in rows:
        for key in ("mean_E", "coll", "mean_q", "shell_n"):
            agg[(r["kind"], int(r["shell"]), key)].append(float(r[key]))
    kinds = sorted({r["kind"] for r in rows})
    shells = sorted({int(r["shell"]) for r in rows})
    print(f"{'kind':>10} {'shell':>5} " + "".join(
        f"{k:>20}" for k in ("mean_E", "coll", "mean_q", "shell_n")))
    for kind in kinds:
        for sh in shells:
            cells = []
            for key in ("mean_E", "coll", "mean_q", "shell_n"):
                v = agg.get((kind, sh, key))
                if v:
                    v = np.array(v)
                    cells.append(f"{v.mean():>12.4f} ±{v.std()/max(1,len(v))**0.5:.4f}")
                else:
                    cells.append(f"{'--':>20}")
            print(f"{kind:>10} {sh:>5} " + "".join(f"{c:>20}" for c in cells))
    print("\nREAD: real vs placebo at each shell = the closure-failure effect"
          " (anchoring subtracted);\n      real vs virtual = total effect;"
          " decay with shell = the strain profile.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--analyze", nargs="+", default=None, metavar="CSV",
                    help="aggregate + compare finished runs' CSVs and exit")
    ap.add_argument("--resume", default=None,
                    help="thermalized closure checkpoint (beta=1 arm)")
    ap.add_argument("--pins", type=int, default=100)
    ap.add_argument("--pin-mode", choices=["fail", "closed", "none"],
                    default="fail",
                    help="fail = defect arm; closed = placebo arm; "
                         "none = pure control (virtual refs only)")
    ap.add_argument("--pin-sep", type=int, default=5,
                    help="min BFS separation between pins (slice graph); "
                         "keep >= rmax+2 so shells stay pin-free")
    ap.add_argument("--rmax", type=int, default=3)
    ap.add_argument("--beta-closure", type=float, default=1.0)
    ap.add_argument("--eta", type=float, default=ETA_STAR)
    ap.add_argument("--lambda-inj", type=float, default=3.0)
    ap.add_argument("--sweeps", type=int, default=3000)
    ap.add_argument("--measure-every", type=int, default=2,
                    help="shell-measurement cadence (sweeps)")
    ap.add_argument("--therm", type=int, default=200,
                    help="sweeps after pinning before measurements start")
    ap.add_argument("--target-n41", type=int, default=20000)
    ap.add_argument("--K", type=int, default=80)
    ap.add_argument("--k0", type=float, default=2.2)
    ap.add_argument("--Delta", type=float, default=0.6)
    ap.add_argument("--k4", type=float, default=0.9)
    ap.add_argument("--eps", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--checkpoint", default=None,
                    help="engine checkpoint path (default: <out>_ckpt.json, "
                         "so concurrently running arms never clobber each "
                         "other's checkpoints)")
    ap.add_argument("--out", default="defect")
    ap.add_argument("--wall-hours", type=float, default=None)
    args = ap.parse_args()

    if args.analyze:
        analyze(args.analyze)
        return
    if not args.resume:
        sys.exit("--resume with a thermalized closure checkpoint is required "
                 "(run stage 2 first)")
    if args.checkpoint is None:
        args.checkpoint = args.out + "_ckpt.json"

    Etab = build_energy_table(args.eta, args.lambda_inj)
    rng = np.random.default_rng(args.seed)
    labels = FaceLabels(rng)

    # mu: reuse from checkpoint when context matches, else TI (same protocol
    # as v6_closure_run; pins are ~1% of cells -- constant offset, harmless)
    import hashlib
    mu_ctx = {"beta_eff": float(args.beta_closure),
              "etab_sha": hashlib.sha256(Etab.tobytes()).hexdigest()[:12]}
    T0, _, extra0 = load_checkpoint(args.resume)
    mu = None
    if isinstance(extra0, dict) and extra0.get("kind") == "closure" \
            and extra0.get("mu") is not None and extra0.get("mu_ctx") == mu_ctx:
        mu = float(extra0["mu"])
        print(f"# [defect] reusing mu = {mu:.4f} from checkpoint", flush=True)
    if mu is None:
        mu = calibrate_mu_ti(T0, Etab, args.beta_closure, seed=args.seed + 777)
    centering = Centering(enabled=True, mu_fixed=mu, label="closure")
    labels.mu_saved = mu
    labels.mu_ctx = mu_ctx

    inc = PinnedClosure(labels, Etab, args.k0, args.Delta, args.k4,
                        args.beta_closure, centering=centering)
    inc.full(T0)                                  # prime on the start state
    if isinstance(extra0, dict) and extra0.get("kind") == "closure":
        labels.load(extra0)                       # thermalized labels
        inc.full(T0)

    # --- pin selection + application ---------------------------------------
    if args.pin_sep < args.rmax + 2:
        print(f"# !! pin-sep {args.pin_sep} < rmax+2 = {args.rmax + 2}: outer "
              f"shells of different pins will overlap; increase --pin-sep or "
              f"lower --rmax", flush=True)
    pin_rng = np.random.default_rng(args.seed + 31)
    pins = []
    if args.pin_mode != "none" and args.pins > 0:
        pins = select_anchors(inc, args.pins, args.pin_sep, pin_rng)
        cfg = FAIL_LABELS if args.pin_mode == "fail" else CLOSED_LABELS
        apply_pins(inc, pins, cfg)
        inc.full(T0)
        print(f"# [defect] pinned {len(pins)} tets in mode '{args.pin_mode}' "
              f"(requested {args.pins}, separation >= {args.pin_sep})",
              flush=True)
    del T0                                        # run() reloads the ckpt

    csv_f = open(args.out + ".csv", "w", newline="")
    writer = csv.DictWriter(csv_f, fieldnames=CSV_FIELDS)
    writer.writeheader()

    pin_set = frozenset(pins)
    state = {"sweep": 0}

    def hook(T, hb_rng):
        state["sweep"] += 1
        _heatbath_pass(inc.labels, inc.tet_pids, inc.E, inc.beta, hb_rng)
        sw = state["sweep"]
        if sw > args.therm and sw % args.measure_every == 0:
            if pins:
                measure_shells(inc, pins, args.pin_mode, sw, writer,
                               rmax=args.rmax, exclude_cells=pin_set)
            # fresh virtual references far from pins, every measurement
            vrng = np.random.default_rng(args.seed + 1000 + sw)
            refs = select_anchors(inc, max(20, args.pins // 2), 2, vrng,
                                  exclude=pins, exclude_dist=args.rmax + 2)
            measure_shells(inc, refs, "virtual", sw, writer, rmax=args.rmax)
            csv_f.flush()

    protect = set()
    for tet in pins:
        protect |= tet
    T = run(f"DEFECT [{args.pin_mode}, pins={len(pins)}]",
            k0=args.k0, Delta=args.Delta, k4=args.k4,
            target_N41=args.target_n41, K=args.K, eps=args.eps,
            seed=args.seed, max_sweeps=args.sweeps, measure_every=50,
            checkpoint=args.checkpoint, resume=args.resume,
            extra_state=labels, delta_action=inc, extra_hook=hook,
            audit_every=50, causal=True, protect_vertices=protect,
            wall_budget_s=(args.wall_hours * 3600 if args.wall_hours else None))
    csv_f.close()

    # pin survival is a HARD gate: every pin must survive or the arm is void
    alive = sum(1 for p in pins if p in inc.tet_pids)
    print(f"\n# DEFECT RUN DONE  mode={args.pin_mode}  pins alive: "
          f"{alive}/{len(pins)}  measurements -> {args.out}.csv")
    if pins and alive != len(pins):
        print("# !! PIN LOSS -- this arm is NOT READABLE (carrier protection "
              "failed); investigate before using the CSV")
    print(f"# analyze:  python v6_defect_run.py --analyze {args.out}.csv "
          f"<other arms' csvs>")


if __name__ == "__main__":
    main()
