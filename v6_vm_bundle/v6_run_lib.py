#!/usr/bin/env python3
"""Shared run harness for the v6 4D CDT engine -- used by BOTH the bare-Regge
verification run (v6_verify_run.py) and the with-theory EPRL run
(v6_theory_run.py), so they share one protocol, one set of (reliable)
measurements, and one checkpoint format.

RELIABLE OBSERVABLES (the d_s diffusion estimator is NOT trustworthy on CDT dual
graphs -- a known extended 3D CDT read d_s~2 -- so we use these instead):
  * Hausdorff dimension d_H from dual-graph shell growth V(r) ~ r^d_H, compared
    against matched-size 2/3/4-torus rails measured the SAME way.
  * Spatial-volume profile V3(tau) = N_41 per time slice -- the de Sitter
    signature is the cos^3 "blob" (vs a flat "stalk" in the collapsed regime).

PROTOCOL (matches the CDT literature, AJL / Gorlich thesis):
  fix N_41 via a quadratic penalty, kappa_4 ~ pseudo-critical, grow to a target
  volume, thermalize, then measure periodically and checkpoint to disk so a long
  VM run can resume.
"""
from __future__ import annotations
import json
import time
import numpy as np
from collections import defaultdict

from v6_cdt import build_s1xs3, verify
from v6_cdt_run import propose_and_apply, hastings_log


# ---------------------------------------------------------------------------
# Reliable geometry observables
# ---------------------------------------------------------------------------

def dual_adjacency(T):
    """Neighbour lists of the spacetime dual graph (pentachora = nodes, shared
    tetrahedral face = edge). Returns (id_list, adj) with contiguous indices."""
    ids = list(T.pent.keys())
    idx = {p: i for i, p in enumerate(ids)}
    adj = [[] for _ in ids]
    for p in ids:
        for q in T.nbr[p]:
            if q in idx:
                adj[idx[p]].append(idx[q])
    return ids, adj


def hausdorff_dim(adj, n_centers=24, maxr=64, seed=0):
    """d_H from cumulative dual-graph ball growth V(r) ~ r^d_H on the rising
    part (before finite-size saturation)."""
    nt = len(adj)
    if nt < 50:
        return float("nan")
    rng = np.random.default_rng(seed)
    shells = np.zeros(maxr + 1)
    centers = rng.choice(nt, size=min(n_centers, nt), replace=False)
    for c in centers:
        seen = np.zeros(nt, bool); seen[c] = True
        frontier = [c]
        for r in range(1, maxr + 1):
            nxt = []
            for u in frontier:
                for w in adj[u]:
                    if not seen[w]:
                        seen[w] = True; nxt.append(w)
            if not nxt:
                break
            shells[r] += len(nxt); frontier = nxt
    shells /= len(centers)
    V = np.cumsum(shells); r = np.arange(maxr + 1)
    peak = int(np.argmax(shells[1:])) + 1 if shells[1:].any() else 3
    lo, hi = 2, max(4, peak)
    m = (r >= lo) & (r <= hi) & (V > 0)
    if m.sum() < 3:
        return float("nan")
    return float(np.polyfit(np.log(r[m]), np.log(V[m]), 1)[0])


def torus_rails(n_target):
    """d_H of matched-size 2/3/4-tori measured with the SAME estimator, as a
    calibration-free dimension ruler."""
    import step3_linkA_harness as linkA
    out = {}
    for dim in (2, 3, 4):
        s = max(3, round(n_target ** (1.0 / dim)))
        G = linkA.torus(s, dim).tolil()
        adj = [list(map(int, G.rows[i])) for i in range(s ** dim)]
        out[dim] = (s ** dim, hausdorff_dim(adj))
    return out


def volume_profile(T):
    """V3(tau): number of (4,1)/(1,4) simplices whose 4-vertex spatial face sits
    in slice tau. The de Sitter blob is a cos^3 bump; a stalk is ~flat-minimal."""
    K = T.K
    prof = [0] * K
    for vs in T.pent.values():
        ts = [T.vtime[v] for v in vs]
        vals, counts = np.unique(ts, return_counts=True)
        if len(vals) == 2 and set(counts.tolist()) == {1, 4}:
            t4 = int(vals[np.argmax(counts)])   # slice holding the 4-vertex side
            prof[t4] += 1
    return prof


def blob_score(prof):
    """A crude 'is there a localized blob' score: max/mean of the slice profile.
    ~1 => flat (stalk/uniform); >>1 => a localized de Sitter-like blob."""
    a = np.array(prof, float)
    return float(a.max() / a.mean()) if a.mean() > 0 else 0.0


def profile_metrics(prof):
    """Richer de Sitter diagnostics on the spatial-volume profile V3(tau).
    The de Sitter universe is a localized blob V3 ~ A*cos^3((t-t0)/B) sitting on
    a minimal 'stalk'. Returns a dict:
      blob_score   : max/mean (>1 => localized)
      max_slice    : peak slice volume
      mean_stalk   : mean of the bottom-half ('stalk') slices
      active_slices: # slices clearly above the stalk (blob extent)
      cos3_width   : fitted blob half-width B (slices); None if no fit
      cos3_relerr  : relative RMS error of the A*cos^3 fit over the blob (lower
                     = more de-Sitter-like); None if no fit
      centered     : profile rolled so the peak is centered (for plotting)
    """
    a = np.array(prof, float)
    K = len(a)
    out = {"blob_score": blob_score(prof), "max_slice": float(a.max()) if K else 0.0,
           "mean_stalk": 0.0, "active_slices": 0, "cos3_width": None,
           "cos3_relerr": None, "centered": []}
    if K == 0 or a.sum() == 0:
        return out
    stalk = float(np.median(a))
    out["mean_stalk"] = float(a[a <= stalk].mean()) if (a <= stalk).any() else stalk
    margin = stalk + 0.25 * (a.max() - stalk)
    out["active_slices"] = int((a > margin).sum())
    # center on the peak
    shift = (K // 2) - int(np.argmax(a))
    centered = np.roll(a, shift)
    out["centered"] = [round(float(x), 1) for x in centered]
    # cos^3 fit of the (stalk-subtracted) blob, centered at K//2
    try:
        import warnings
        from scipy.optimize import curve_fit, OptimizeWarning
        y = np.clip(centered - out["mean_stalk"], 0, None)
        x = np.arange(K) - K // 2

        def cos3(x, A, B):
            u = np.clip(x / B, -np.pi / 2, np.pi / 2)
            return A * np.cos(u) ** 3

        A0 = max(y.max(), 1.0); B0 = max(out["active_slices"] / 2.0, 1.0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            popt, _ = curve_fit(cos3, x, y, p0=[A0, B0], maxfev=5000)
        resid = y - cos3(x, *popt)
        denom = np.sqrt((y ** 2).mean()) or 1.0
        out["cos3_width"] = float(abs(popt[1]))
        out["cos3_relerr"] = float(np.sqrt((resid ** 2).mean()) / denom)
    except Exception:
        pass
    return out


# ---------------------------------------------------------------------------
# Checkpointing (so a multi-day VM run can resume)
# ---------------------------------------------------------------------------

def save_checkpoint(T, path, meta=None, extra=None):
    """Checkpoint the triangulation (+ optional `extra` payload, e.g. the
    serialized intertwiner labels, so a theory run is a FAITHFUL resume)."""
    data = {
        "K": T.K,
        "vtime": {str(k): v for k, v in T.vtime.items()},
        "pent": {str(k): list(v) for k, v in T.pent.items()},
        "nbr": {str(k): v for k, v in T.nbr.items()},
        "next_pid": T._next_pid, "next_vid": T._next_vid,
        "meta": meta or {},
        "extra": extra,
    }
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f)
    import os
    os.replace(tmp, path)   # atomic: a crash mid-write can't corrupt the ckpt


def load_checkpoint(path):
    """Returns (Triangulation, meta, extra)."""
    from v6_cdt import Triangulation
    d = json.load(open(path))
    T = Triangulation(K=d["K"])
    T.vtime = {int(k): v for k, v in d["vtime"].items()}
    T._next_vid = d["next_vid"]; T._next_pid = d["next_pid"]
    for pid_s, vs in d["pent"].items():
        pid = int(pid_s)
        vs = tuple(vs)
        T.pent[pid] = vs
        T.nbr[pid] = list(d["nbr"][pid_s])
        for v in vs:
            s = T.vinc.get(v)
            if s is None:
                T.vinc[v] = {pid}; T.vert_bag.add(v)
            else:
                s.add(pid)
        T.pent_bag.add(pid)
        setattr(T, T._classify(vs), getattr(T, T._classify(vs)) + 1)
    return T, d.get("meta", {}), d.get("extra")


# ---------------------------------------------------------------------------
# The run protocol (geometry-only; action is pluggable via the driver)
# ---------------------------------------------------------------------------

def run(action_label, k0, Delta, k4, target_N41, K=21, eps=1e-4, beta=1.0,
        seed=0, max_sweeps=20000, measure_every=50, checkpoint=None,
        strictness=2, causal=True, protect_vertices=None, wall_budget_s=None,
        verbose=True, resume=None, extra_hook=None, geometry_action=None,
        extra_state=None, delta_action=None, audit_every=0,
        link_check_every=10):
    """Run CDT to `target_N41`, then thermalize, logging d_H + the volume
    profile and checkpointing.

    `geometry_action(T) -> float` (optional) REPLACES the bare Regge action --
    the with-theory run passes the EPRL log-amplitude here. The N_41 volume
    penalty is always added on top. `extra_hook(T, rng)` (optional) runs once per
    sweep -- the EPRL run uses it for an intertwiner heat-bath pass.

    `delta_action` (optional) is a LOCAL incremental-action object that supplies
    the FULL geometry action (Regge + theory) as O(footprint) deltas instead of
    O(N) recomputes. When given it overrides `geometry_action`. It must expose:
      full(T)->float, delta(T, added, removed)->dS, accept(), reject(T, a, r),
      refresh(T) [after a heat-bath], and audit(T)->(running, recomputed).
    `audit_every` (>0) periodically cross-checks the running delta sum against a
    from-scratch recompute and aborts on drift -- the safety net for the fast path.

    Returns the final Triangulation."""
    from v6_cdt import regge_action
    if geometry_action is None and delta_action is None:
        def geometry_action(T):
            return regge_action(T, k0, Delta, k4)

    rng = np.random.default_rng(seed)
    if resume:
        T, _, extra = load_checkpoint(resume)
        if extra_state is not None and extra is not None:
            extra_state.load(extra)        # faithful resume of intertwiner labels
        if verbose:
            print(f"[resume] N4={T.n_pent()} N41={T.type_counts()[0]} "
                  f"(extra_state {'restored' if extra is not None else 'absent'})")
    else:
        T = build_s1xs3(K=K)

    def action():
        n41 = T.type_counts()[0]
        return geometry_action(T) + eps * (n41 - target_N41) ** 2

    if delta_action is not None:
        delta_action.full(T)          # prime the incremental cache
        S = None                      # delta path is self-contained
    else:
        S = action()
    t0 = time.time()
    if verbose:
        print(f"# v6 run [{action_label}]  target N41={target_N41}  K={K}  "
              f"(k0={k0}, D={Delta}, k4={k4}, eps={eps})")
        print(f"{'sweep':>6} {'N4':>7} {'N41':>7} {'d_H':>6} {'blob':>6} "
              f"{'active':>6} {'cos3err':>7} {'valid':>6} {'links':>5} "
              f"{'slice':>5} {'wall_s':>8}")

    for sw in range(1, max_sweeps + 1):
        for _ in range(max(1, T.n_pent())):
            N4_b, N0_b = T.n_pent(), len(T.vinc)
            if delta_action is not None:
                n41_b = T.type_counts()[0]
                T.begin_record()
                mt, ok, undo = propose_and_apply(T, rng, strictness, causal,
                                                 protect_vertices)
                if not ok:
                    # NOT a no-op: a filter-triggered undo (strictness/causal)
                    # re-creates pentachora with FRESH pids. Absorb the
                    # apply+undo log so pid-keyed incremental state stays
                    # honest (pre-fix this silently ghosted cache entries).
                    ua, ur = T.take_record()
                    if ua or ur:
                        delta_action.reject(T, ua, ur)
                    continue
                added, removed = T.take_record()
                dS_geom = delta_action.delta(T, added, removed)
                n41_a = T.type_counts()[0]
                dS_pen = eps * ((n41_a - target_N41) ** 2 - (n41_b - target_N41) ** 2)
                H = hastings_log(mt, N4_b, N0_b)
                if rng.random() < np.exp(min(50.0, -beta * (dS_geom + dS_pen) + H)):
                    delta_action.accept()
                else:
                    T.begin_record()
                    undo()
                    ua, ur = T.take_record()     # the undo's footprint
                    delta_action.reject(T, ua, ur)
            else:
                mt, ok, undo = propose_and_apply(T, rng, strictness, causal,
                                                 protect_vertices)
                if not ok:
                    continue
                S_new = action()
                H = hastings_log(mt, N4_b, N0_b)
                if rng.random() < np.exp(min(50.0, -beta * (S_new - S) + H)):
                    S = S_new
                else:
                    undo()
        # Audit the ACCUMULATED incremental state (before any full refresh), so
        # this genuinely checks the move-loop deltas, not a fresh rebuild.
        if delta_action is not None and audit_every and sw % audit_every == 0:
            run_s, recomp_s = delta_action.audit(T)
            drift = abs(run_s - recomp_s)
            if drift > 1e-6 * max(1.0, abs(recomp_s)):
                raise RuntimeError(
                    f"[delta-action audit] running S={run_s:.6f} != "
                    f"recompute {recomp_s:.6f} (drift {drift:.2e}) at sweep {sw}")
            elif verbose:
                print(f"# [audit ok] sweep {sw}: S={run_s:.4f} matches recompute "
                      f"(drift {drift:.1e})", flush=True)
        if extra_hook is not None:
            extra_hook(T, rng)
            if delta_action is not None:
                delta_action.refresh(T)          # heat-bath changed labels
            elif S is not None:
                S = action()   # global path: cached S is stale after the hook
        if sw % measure_every == 0:
            ids, adj = dual_adjacency(T)
            dH = hausdorff_dim(adj)
            prof = volume_profile(T)
            pm = profile_metrics(prof)
            meas_i = sw // measure_every
            # full S^3-link verification periodically (it's the expensive check);
            # the gluing-based necessary conditions run every measurement.
            do_links = bool(link_check_every) and (meas_i % link_check_every == 0)
            okv, repv = verify(T, check_links=do_links)
            lf = repv["link_failures"]
            link_s = ("-" if lf is None else "ok" if lf == 0 else
                      "gen" if lf == "unreliable(generalized)" else "BAD")
            from v6_cdt_moves import causal_slice_report
            sl_bad, _, sl_ss = causal_slice_report(T)
            slice_s = "ok" if (sl_bad == 0 and sl_ss == 0) else f"{sl_bad + sl_ss}"
            if verbose:
                cerr = pm["cos3_relerr"]
                cerr_s = f"{cerr:>7.3f}" if cerr is not None else f"{'--':>7}"
                print(f"{sw:>6} {T.n_pent():>7} {T.type_counts()[0]:>7} "
                      f"{dH:>6.2f} {pm['blob_score']:>6.2f} "
                      f"{pm['active_slices']:>6} {cerr_s} "
                      f"{'ok' if okv else 'BAD':>6} {link_s:>5} {slice_s:>5} "
                      f"{time.time()-t0:>8.0f}", flush=True)
            if not okv:
                print(f"# !! verify FAILED at sweep {sw}: gluing={repv['gluing']}",
                      flush=True)
            if checkpoint:
                ex = extra_state.serialize() if extra_state is not None else None
                save_checkpoint(T, checkpoint,
                                meta={"sweep": sw, "d_H": dH, "profile": prof,
                                      "profile_metrics": pm,
                                      "verify_ok": bool(okv), "links": link_s,
                                      "slice_defects": int(sl_bad + sl_ss)},
                                extra=ex)
        if wall_budget_s and time.time() - t0 > wall_budget_s:
            if verbose:
                print(f"# wall budget reached at sweep {sw}")
            break
    # FINAL full verification (incl. S^3 vertex links) -- the gate must pass this
    okf, repf = verify(T, check_links=True)
    from v6_cdt_moves import causal_slice_report
    sl_bad, sl_hist, sl_ss = causal_slice_report(T)
    T._final_verify = (okf, repf)            # entry points surface this
    T._final_slices = (sl_bad, sl_hist, sl_ss)
    if verbose:
        ga = repf["gluing"]
        print(f"# FINAL verify: ok={okf}  gluing_ok={repf['ok_gluing_only']}  "
              f"links={repf['link_failures']}  simplicial={ga['is_simplicial']}  "
              f"chi={repf['euler_char']}  n_tets={ga['n_tets']} "
              f"(expect {5*T.n_pent()//2})  connected={ga['connected']}", flush=True)
        print(f"# FINAL foliation: slice-triangle incidence {sl_hist} "
              f"(all-2 = closed 3-manifold slices)  same-side tets={sl_ss}  "
              f"=> {'CLEAN (standard CDT ensemble)' if sl_bad == 0 and sl_ss == 0 else 'DEFECTIVE (generalized ensemble -- see causal_slice_ok)'}",
              flush=True)
    return T
