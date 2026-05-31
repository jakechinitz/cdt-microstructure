#!/usr/bin/env python3
"""WITH-THEORY run: v6 4D CDT geometry coupled to the theory's EPRL / closure
spin-foam amplitude.

This is the ACTUAL theory test. Each tetrahedron (a shared tetrahedral face
between two pentachora, in the v6 gluing representation) carries an intertwiner
label; each pentachoron (4-simplex) contributes the rank-5 vertex tensor
T[i_0,...,i_4] on its five tetrahedral faces; the EPRL action is
  S_EPRL = -sum_pentachora log|T[face intertwiners]|.
Intertwiners are updated by a heat-bath between geometry sweeps.

TWO COUPLING MODES (the key fix -- replacing Regge is the riskier test):
  * --mode regge_plus_eprl  (DEFAULT, the safer first real test):
        S_geom = regge_action(T,k0,Delta,k4) + beta_eprl * S_EPRL
    The known-good Einstein-Regge dynamics still drives the geometry; the EPRL
    amplitude is a *correction* whose strength is beta_eprl. Turn beta_eprl up
    from ~0 and watch whether the amplitude steers d_H / the volume profile.
  * --mode eprl_only:
        S_geom = S_EPRL
    Pure amplitude geometry (no Regge). Only interpretable once the amplitude
    fidelity is validated.

>>> READ BEFORE TRUSTING ANY NUMBER FROM THIS FILE <<<
  1. RUN v6_verify_run.py FIRST. If the bare engine does NOT reproduce 4D CDT
     (d_H -> ~4, de Sitter blob), a result here is uninterpretable -- you cannot
     tell tool artefact from theory. Bare verification is the gate.
  2. AMPLITUDE FIDELITY: this loads vertex_j3.npz (the frozen-j=3 vertex tensor).
     Until the full five-tetrahedron {15j} contraction is convention-checked,
     ABSOLUTE d_H / d_s numbers are ARTEFACT-grade. What IS robust
     (finite-size-cancelling) is the COMPARISON at matched volume: bare vs
     regge_plus_eprl (sweep beta_eprl) -- does the amplitude move the geometry?
  3. DETAILED BALANCE: tetrahedra created by a geometry move get a fresh
     intertwiner (uniform), not folded into the Hastings ratio. Adequate for the
     steer-direction comparison; not for precision sampling. The heat-bath
     re-equilibrates labels each sweep.
  4. kappa_4 = 0.9 is a STARTING GUESS in our conventions, not gospel. If N_41
     runs away or freezes, retune kappa_4 (and/or eps) -- the volume penalty is
     the real anchor; kappa_4 only sets the pseudo-critical baseline.

USAGE (after verification passes):
  # additive correction (safe first test), global O(N) action:
  python v6_theory_run.py --mode regge_plus_eprl --beta-eprl 0.1 \
         --target-n41 4000 --K 24 --checkpoint thy.json
  # same, but fast local incremental action with periodic self-audit:
  python v6_theory_run.py --mode regge_plus_eprl --beta-eprl 0.1 --local-eprl \
         --audit-every 20 --target-n41 8000 --K 40 --checkpoint thy.json
"""
from __future__ import annotations
import argparse
import numpy as np

from v6_run_lib import (run, dual_adjacency, hausdorff_dim, volume_profile,
                        profile_metrics, torus_rails)
from v6_cdt import regge_action


# ---------------------------------------------------------------------------
# Intertwiner labels (per shared tetrahedron) -- with checkpoint (de)serialize
# ---------------------------------------------------------------------------

class Intertwiners:
    """Per-tetrahedron intertwiner labels. A tetrahedron is the shared face
    between two pentachora P and Q (= P.nbr[i]); keyed by frozenset({P,Q})
    (with strictness, two pentachora share at most one face, so this is unique).
    """
    def __init__(self, T, D, rng):
        self.D = D
        self.rng = rng
        self.lab = {}
        if T is not None:
            self.sync(T)

    def ensure(self, T):
        """Add labels for any new tetrahedra but NEVER drop existing ones. Used
        inside action evaluation so a REJECTED geometry proposal cannot mutate
        the label state of surviving tetrahedra (only adds harmless stale entries
        for the proposed-then-undone tets, which sync() prunes later)."""
        for p in T.pent:
            for q in T.nbr[p]:
                if q in T.pent:
                    k = frozenset((p, q))
                    if k not in self.lab:
                        self.lab[k] = int(self.rng.integers(self.D))

    def sync(self, T):
        """Assign labels to new tetrahedra (uniform) AND drop stale ones. Only
        safe to call on an ACCEPTED state (the heat-bath / full rebuild), never
        mid-proposal."""
        cur = set()
        for p in T.pent:
            for q in T.nbr[p]:
                if q in T.pent:
                    cur.add(frozenset((p, q)))
        for k in cur:
            if k not in self.lab:
                self.lab[k] = int(self.rng.integers(self.D))
        for k in list(self.lab):
            if k not in cur:
                del self.lab[k]

    def faces(self, T, p):
        return tuple(self.lab[frozenset((p, T.nbr[p][i]))] for i in range(5))

    # --- checkpoint (de)serialization (Fix: faithful theory-run resume) ------
    def serialize(self):
        """JSON-able snapshot of the labels (keyed by sorted pid pairs)."""
        out = []
        for k, c in self.lab.items():
            p, q = sorted(k)
            out.append([p, q, int(c)])
        return {"D": int(self.D), "lab": out}

    def load(self, d):
        """Restore labels from a serialize() payload. Pids are preserved by the
        geometry checkpoint, so the frozenset keys remain valid."""
        self.D = int(d["D"])
        self.lab = {frozenset((int(a), int(b))): int(c) for a, b, c in d["lab"]}


# ---------------------------------------------------------------------------
# EPRL centering: remove the trivial cosmological-constant (volume) shift
# ---------------------------------------------------------------------------

class Centering:
    """The raw EPRL action S_EPRL = sum_p (-log|amp_p|) is EXTENSIVE: its bulk is
    just (mean per-pentachoron cost) * N4, i.e. a renormalization of kappa_4 (the
    cosmological constant). That trivial piece collapses/inflates the universe and
    is normalization-dependent (artefact-grade for the frozen-j3 tensor) -- it is
    NOT the physics. Subtracting a fixed reference mu per pentachoron,

        S_EPRL -> sum_p (-log|amp_p| - mu),

    leaves only the FLUCTUATIONS of the amplitude across configurations -- the
    part that can actually steer geometry SHAPE -- and makes the term
    volume-neutral, so kappa_4 stays at its bare value across the whole beta
    sweep (no per-beta retuning). mu is auto-calibrated once from the starting
    configuration (or fixed via --eprl-mu)."""
    def __init__(self, enabled=True, mu_fixed=None):
        self.enabled = enabled
        self.mu_fixed = mu_fixed
        self.mu = None                      # set once by calibrate()

    def value(self):
        return (self.mu or 0.0) if self.enabled else 0.0

    def calibrate(self, raw_costs):
        """Fix mu on first call only (so it never drifts mid-run)."""
        if not self.enabled:
            self.mu = 0.0
            return
        if self.mu is not None:
            return
        if self.mu_fixed is not None:
            self.mu = float(self.mu_fixed)
        else:
            costs = list(raw_costs)
            self.mu = (sum(costs) / len(costs)) if costs else 0.0
        print(f"# [EPRL centering] mu = {self.mu:.4f} per pentachoron "
              f"({'fixed' if self.mu_fixed is not None else 'auto-calibrated'}); "
              f"EPRL term now measures fluctuations about the mean "
              f"(volume-neutral -- keep k4 at the bare value)", flush=True)


# ---------------------------------------------------------------------------
# Global (O(N) per move) EPRL action + heat-bath  -- the reference path
# ---------------------------------------------------------------------------

def make_eprl_action(intw, Ttensor, centering):
    """-sum_p (log|T[faces of p]| + mu) recomputed over ALL pentachora (mu is the
    centering constant; mu=0 => raw extensive action). Correct but O(N) per move;
    use IncrementalEPRL for large runs. Uses ensure() (add-only) so a rejected
    proposal never drops a surviving label."""
    def eprl(T):
        intw.ensure(T)
        raw = []
        for p in T.pent:
            amp = abs(Ttensor[intw.faces(T, p)])
            raw.append(-np.log(amp if amp > 1e-300 else 1e-300))
        centering.calibrate(raw)            # fixes mu on the first evaluation
        mu = centering.value()
        return sum(c - mu for c in raw)
    return eprl


def make_geometry_action(mode, intw, Ttensor, k0, Delta, k4, beta_eprl, centering):
    """Build the (global) geometry action for the chosen coupling mode."""
    eprl = make_eprl_action(intw, Ttensor, centering)
    if mode == "eprl_only":
        return eprl
    if mode == "regge_plus_eprl":
        def geom(T):
            return regge_action(T, k0, Delta, k4) + beta_eprl * eprl(T)
        return geom
    raise ValueError(f"unknown mode {mode!r}")


def make_heatbath(intw, Ttensor):
    """One heat-bath pass: resample each tet's intertwiner conditional on the two
    pentachora that share it. (Equilibrates labels between geometry sweeps.)"""
    D = intw.D

    def hook(T, rng):
        intw.sync(T)
        for key in list(intw.lab.keys()):
            ps = [p for p in key if p in T.pent]
            if len(ps) != 2:
                continue
            logw = np.zeros(D)
            old = intw.lab[key]
            for c in range(D):
                intw.lab[key] = c
                lw = 0.0
                for p in ps:
                    amp = abs(Ttensor[intw.faces(T, p)])
                    lw += np.log(amp if amp > 1e-300 else 1e-300)
                logw[c] = lw
            logw -= logw.max()
            w = np.exp(logw); tot = w.sum()
            intw.lab[key] = old if tot < 1e-300 else int(rng.choice(D, p=w / tot))
    return hook


# ---------------------------------------------------------------------------
# Local incremental EPRL (Fix 3): full geometry action as O(footprint) deltas
# ---------------------------------------------------------------------------

class IncrementalEPRL:
    """Maintains the FULL geometry action (Regge + beta_eprl * S_EPRL) as a
    running value, updating only the pentachora touched by each move via the
    engine's change-log. The Regge part is O(1) from running counts; the EPRL
    part is recomputed only for added pentachora and their immediate neighbours
    (whose shared-face label changed). A periodic audit() recomputes from scratch
    and the driver aborts on any drift -- the safety net for this fast path.

    The cache is reconciled FORWARD on both accept and reject (undo re-creates
    pentachora with fresh pids, so we never 'roll back' -- we just absorb the
    undo's own change-log), so it always matches the current triangulation.
    """
    def __init__(self, intw, Ttensor, mode, k0, Delta, k4, beta_eprl, rng,
                 centering=None):
        self.intw = intw
        self.T = Ttensor
        self.D = intw.D
        self.rng = rng
        self.use_regge = (mode == "regge_plus_eprl")
        self.beta = beta_eprl if self.use_regge else 1.0
        self.k0, self.Delta, self.k4 = k0, Delta, k4
        self.centering = centering if centering is not None else Centering(enabled=False)
        self.contrib = {}      # pid -> (-log|amp| - mu)  [centered]
        self.S_eprl = 0.0
        self.regge_prev = 0.0
        self._pending_regge = 0.0

    # --- low level ---------------------------------------------------------
    def _regge(self, T):
        return regge_action(T, self.k0, self.Delta, self.k4) if self.use_regge else 0.0

    def _raw_contrib(self, T, p):
        faces = tuple(self.intw.lab[frozenset((p, T.nbr[p][i]))] for i in range(5))
        amp = abs(self.T[faces])
        return -np.log(amp if amp > 1e-300 else 1e-300)

    def _contrib(self, T, p):
        """Centered per-pentachoron contribution (raw -log|amp| minus mu)."""
        return self._raw_contrib(T, p) - self.centering.value()

    def _ensure_labels(self, T, pids):
        for p in pids:
            for q in T.nbr[p]:
                if q in T.pent:
                    k = frozenset((p, q))
                    if k not in self.intw.lab:
                        self.intw.lab[k] = int(self.rng.integers(self.D))

    def _reconcile(self, T, added, removed):
        """Bring contrib/S_eprl into agreement with the current T after a set of
        add/remove ops captured by the change-log."""
        for rec in removed:
            pid = rec[0]
            c = self.contrib.pop(pid, None)
            if c is not None:
                self.S_eprl -= c
        add_present = [p for p in added if p in T.pent]
        affected = set(add_present)
        for p in add_present:
            for q in T.nbr[p]:
                if q in T.pent:
                    affected.add(q)         # neighbour's shared-face label changed
        self._ensure_labels(T, affected)
        for p in affected:
            new = self._contrib(T, p)
            old = self.contrib.get(p, 0.0)
            self.S_eprl += new - old
            self.contrib[p] = new

    # --- driver interface --------------------------------------------------
    def full(self, T):
        """Recompute everything from scratch and prime the cache. On the very
        first call this also calibrates the centering constant mu from the
        starting configuration."""
        self.intw.sync(T)
        raw = {p: self._raw_contrib(T, p) for p in T.pent}
        self.centering.calibrate(raw.values())     # fixes mu once
        mu = self.centering.value()
        self.contrib = {p: c - mu for p, c in raw.items()}
        self.S_eprl = sum(self.contrib.values())
        self.regge_prev = self._regge(T)
        return self.regge_prev + self.beta * self.S_eprl

    def delta(self, T, added, removed):
        """dS for the just-applied move (move already mutated T)."""
        pre = self.regge_prev + self.beta * self.S_eprl       # pre-move value
        self._reconcile(T, added, removed)                    # -> post-move cache
        regge_now = self._regge(T)
        post = regge_now + self.beta * self.S_eprl
        self._pending_regge = regge_now
        return post - pre

    def accept(self):
        self.regge_prev = self._pending_regge

    def reject(self, T, undo_added, undo_removed):
        """Move rejected: undo() already ran. Absorb the undo's change-log so the
        cache matches the restored (re-pid'd) triangulation."""
        self._reconcile(T, undo_added, undo_removed)
        self.regge_prev = self._regge(T)

    def refresh(self, T):
        """After a heat-bath pass mutates labels, rebuild the cache."""
        self.full(T)

    def audit(self, T):
        """(running_S, from_scratch_S) for the driver's drift check."""
        running = self.regge_prev + self.beta * self.S_eprl
        s = 0.0
        for p in T.pent:
            s += self._contrib(T, p)
        recomputed = self._regge(T) + self.beta * s
        return running, recomputed


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", choices=["regge_plus_eprl", "eprl_only"],
                   default="regge_plus_eprl",
                   help="regge_plus_eprl (default, safer): Regge + beta_eprl*S_EPRL; "
                        "eprl_only: pure amplitude geometry")
    p.add_argument("--beta-eprl", type=float, default=0.1,
                   help="EPRL coupling strength in regge_plus_eprl mode "
                        "(sweep from ~0 upward to see the steer)")
    p.add_argument("--local-eprl", action="store_true",
                   help="use the O(footprint) incremental action (fast); else "
                        "the global O(N) recompute (simple/reference)")
    p.add_argument("--audit-every", type=int, default=25,
                   help="with --local-eprl, recompute-and-check the action every "
                        "N sweeps (0 disables; aborts on drift)")
    p.add_argument("--center-eprl", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="subtract the mean per-pentachoron EPRL cost so the "
                        "amplitude term is volume-neutral (removes the trivial "
                        "cosmological-constant shift). Keep --k4 at the bare value "
                        "across the whole beta sweep. Use --no-center-eprl for the "
                        "raw extensive action (then retune --k4 per beta).")
    p.add_argument("--eprl-mu", type=float, default=None,
                   help="fixed centering constant mu per pentachoron; default "
                        "auto-calibrates from the starting configuration")
    p.add_argument("--target-n41", type=int, default=4000)
    p.add_argument("--K", type=int, default=24)
    p.add_argument("--k0", type=float, default=2.2)
    p.add_argument("--Delta", type=float, default=0.6)
    p.add_argument("--k4", type=float, default=0.9,
                   help="STARTING GUESS in our conventions; retune if N_41 "
                        "runs away/freezes")
    p.add_argument("--eps", type=float, default=1e-3)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-sweeps", type=int, default=20000)
    p.add_argument("--measure-every", type=int, default=20)
    p.add_argument("--vertex", type=str, default="vertex_j3.npz")
    p.add_argument("--checkpoint", type=str, default="v6_theory_ckpt.json")
    p.add_argument("--resume", type=str, default=None)
    p.add_argument("--wall-hours", type=float, default=None)
    args = p.parse_args()

    print("!! WITH-THEORY (EPRL) run. Verify the BARE engine first "
          "(v6_verify_run.py). Amplitude is frozen-j3 (artefact-grade absolutes);"
          " the robust signal is the bare-vs-EPRL comparison at matched volume.")
    print(f"   mode={args.mode}  beta_eprl={args.beta_eprl}  "
          f"local_eprl={args.local_eprl}  center_eprl={args.center_eprl}\n")

    from vertex_tensor import FaithfulVertex
    Ttensor = FaithfulVertex.load(args.vertex).dense_tensor().astype(np.float64)
    D = Ttensor.shape[0]
    rng = np.random.default_rng(args.seed)
    intw = Intertwiners(None, D, rng)             # labels filled lazily on first T
    centering = Centering(enabled=args.center_eprl, mu_fixed=args.eprl_mu)

    common = dict(k0=args.k0, Delta=args.Delta, k4=args.k4,
                  target_N41=args.target_n41, K=args.K, eps=args.eps,
                  seed=args.seed, max_sweeps=args.max_sweeps,
                  measure_every=args.measure_every, checkpoint=args.checkpoint,
                  resume=args.resume, extra_state=intw,
                  wall_budget_s=(args.wall_hours * 3600 if args.wall_hours else None))

    if args.local_eprl:
        inc = IncrementalEPRL(intw, Ttensor, args.mode, args.k0, args.Delta,
                              args.k4, args.beta_eprl, rng, centering=centering)
        heatbath = make_heatbath(intw, Ttensor)
        T = run(f"EPRL [{args.mode}, local]", delta_action=inc,
                extra_hook=heatbath, audit_every=args.audit_every, **common)
    else:
        geom = make_geometry_action(args.mode, intw, Ttensor, args.k0, args.Delta,
                                    args.k4, args.beta_eprl, centering)
        heatbath = make_heatbath(intw, Ttensor)
        T = run(f"EPRL [{args.mode}, global]", geometry_action=geom,
                extra_hook=heatbath, **common)

    ids, adj = dual_adjacency(T)
    dH = hausdorff_dim(adj); prof = volume_profile(T); pm = profile_metrics(prof)
    rails = torus_rails(T.n_pent())
    print("\n" + "=" * 64)
    print("  WITH-THEORY (EPRL) RESULT  -- compare to bare verification run")
    print("=" * 64)
    okf, repf = getattr(T, "_final_verify", (None, {}))
    cen = ("off" if not args.center_eprl else
           f"mu={centering.mu:.4f}" if centering.mu is not None else "on")
    print(f"  mode={args.mode}  beta_eprl={args.beta_eprl}  centering={cen}")
    print(f"  manifold check (gluing-based + S^3 links): "
          f"{'PASS' if okf else 'FAIL'}  [links={repf.get('link_failures')}, "
          f"simplicial={repf.get('gluing', {}).get('is_simplicial')}]")
    print(f"  final N4={T.n_pent()}  N41={T.type_counts()[0]}")
    print(f"  d_H = {dH:.2f}   rails: 2t={rails[2][1]:.2f} 3t={rails[3][1]:.2f} "
          f"4t={rails[4][1]:.2f}")
    print(f"  blob score = {pm['blob_score']:.2f}   active slices = "
          f"{pm['active_slices']}/{T.K}   max slice = {pm['max_slice']:.0f}")
    if pm['cos3_relerr'] is not None:
        print(f"  cos^3 fit: width={pm['cos3_width']:.1f}  "
              f"rel.RMS err={pm['cos3_relerr']:.3f}")
    print("  (compare d_H / blob / cos^3-err to the bare-Regge run at the SAME N"
          " to see if the EPRL amplitude steers geometry toward 4D.)")


if __name__ == "__main__":
    main()
