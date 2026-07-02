#!/usr/bin/env python3
"""WITH-THEORY run #2: v6 4D CDT geometry coupled to the paper's OWN
admissibility-closure weighting (Section 5-6 / Appendix B) -- NOT a spin-foam
amplitude. See CLOSURE_MODEL.md for the derivation and the fidelity notes.

THE CELL AND ITS WEIGHT (verbatim from the paper, verified numerically here):
the UV cell is a tetrahedron whose four faces carry labels m in {-3..3} (the
seven-state face sector, j_eff = 3), injectively (all four distinct), with the
closure invariant

    K^2(m) = 48 - (S^2 - Sigma^2)/3 ,   S = sum m_i,  Sigma^2 = sum m_i^2,

and admissibility weight exp(-eta* K^2) at the paper-fixed precision
eta* = 0.0298668443935 (stationary closure evidence <K^2> = 3/(2 eta*)).
This module reproduces the paper's numbers exactly from that spec:
<K^2>_eta* = 50.223 and g_share,eff = 7.4198.

THE IDENTIFICATION (the one modeling step): the theory's volumetric cells are
the SPATIAL TETRAHEDRA of the CDT time slices; their triangular faces are the
shared faces through which neighboring cells correlate -- the "sharing" of the
sharing entropy. Labels live on slice triangles; each triangle is shared by
exactly two spatial tetrahedra (this REQUIRES the causal-slice engine fix --
slices must be closed 3-manifolds; the run enforces it by default and refuses
to trust results on defective slices).

THE COUPLED MODEL (same exact-sampler pattern as the fixed EPRL run):

    pi(g, m)  propto  exp( -S_Regge(g) - eps*(N41-target)^2
                           - beta * sum_tets [ E(m_faces(tet)) - mu ] )
    E(m1..m4) = eta * K^2(m)  +  lambda_inj * (# equal-label pairs)

with uniform base measure 1/7 per triangle label. Injectivity is enforced as
a soft penalty lambda_inj (the paper's hard constraint is the
lambda_inj -> infinity limit; soft keeps the label chain ergodic and makes
beta a genuine coupling dial). beta = 1 with eta = eta* is the
theory-designated point; beta = 0 recovers bare CDT exactly.

WHY THE EPRL AUDIT FINDINGS DO NOT RECUR HERE (by construction):
  * positive weight -- no sign problem, nothing discarded (finding 5);
  * E is permutation-symmetric in the four faces -- no slot convention
    (finding 4);
  * heat-bath at the same exponent beta as the geometry action from day one
    (finding 2), uniform label births exactly detailed-balanced under the
    normalized base measure (finding 3);
  * mu is the label free energy by thermodynamic integration (finding 1) --
    and both N_tets (= N41/2) and N_triangles (= 2 N_tets) are pinned by the
    N41 volume penalty, so the extensive part of the closure term cannot
    drive volume at all; mu only keeps the action numerics clean and the
    (4,1)/(3,2) mix unbiased.
  * placebo arm: --placebo shuffles the 210 label-orbit energies (same value
    pool, closure structure destroyed) -- only a real-vs-placebo difference
    at matched beta is attributable to the closure weighting.

USAGE:
  python v6_closure_run.py --beta-closure 1.0 --local-closure --audit-every 25 \
         --resume scan_20k.json --target-n41 20000 --K 80 --checkpoint clo.json
"""
from __future__ import annotations
import argparse
from itertools import combinations, product
import numpy as np

from v6_run_lib import (run, dual_adjacency, hausdorff_dim, volume_profile,
                        profile_metrics, torus_rails)
from v6_cdt import regge_action
from v6_theory_run import Centering

ETA_STAR = 0.0298668443935          # paper App. B: stationary closure evidence


# ---------------------------------------------------------------------------
# The per-cell energy table (7^4, symmetric): eta*K^2 + lambda_inj*collisions
# ---------------------------------------------------------------------------

def build_energy_table(eta=ETA_STAR, lambda_inj=3.0, placebo_seed=None):
    """E[l1,l2,l3,l4] over labels l=0..6 (m=l-3). Permutation-symmetric by
    construction. placebo_seed != None shuffles the energies among the 210
    sorted-label orbits: same value pool and symmetry, closure structure
    destroyed -- the placebo arm."""
    E = np.empty((7, 7, 7, 7), dtype=np.float64)
    for idx in product(range(7), repeat=4):
        m = [i - 3 for i in idx]
        S = float(sum(m)); Sig2 = float(sum(x * x for x in m))
        K2 = 48.0 - (S * S - Sig2) / 3.0
        coll = sum(1 for a in range(4) for b in range(a + 1, 4)
                   if m[a] == m[b])
        E[idx] = eta * K2 + lambda_inj * coll
    if placebo_seed is not None:
        orbits = {}
        for idx in product(range(7), repeat=4):
            orbits.setdefault(tuple(sorted(idx)), []).append(idx)
        keys = sorted(orbits)
        vals = np.array([E[orbits[k][0]] for k in keys])
        perm = np.random.default_rng(placebo_seed).permutation(len(keys))
        for k, pi in zip(keys, perm):
            for idx in orbits[k]:
                E[idx] = vals[pi]
    return E


def ensemble_report(eta=ETA_STAR):
    """Reproduce the paper's App. B numbers from the injective ensemble --
    the fidelity check that this module implements the stated weighting."""
    from itertools import permutations
    K2 = []
    for m in permutations(range(-3, 4), 4):
        S = sum(m); Sig2 = sum(x * x for x in m)
        K2.append(48.0 - (S * S - Sig2) / 3.0)
    K2 = np.array(K2)
    w = np.exp(-eta * K2); Z = w.sum()
    mean_K2 = float((K2 * w).sum() / Z)
    g_eff = float(np.log(2 * Z) + eta * mean_K2)   # entropy of p_eta, 1680 states
    return mean_K2, 3.0 / (2 * eta), g_eff


# ---------------------------------------------------------------------------
# Spatial-slice geometry helpers
# ---------------------------------------------------------------------------

def tet_of(T, vs):
    """The spatial tetrahedron of a pentachoron vertex-tuple (frozenset of the
    4 same-time vertices) or None for (3,2)-types. Only valid for PRESENT
    pentachora ((8,2) deletes the removed vertex's vtime, so removed records
    cannot be re-classified -- IncrementalClosure keeps a pid->tet map
    instead)."""
    ts = [T.vtime[v] for v in vs]
    lo = min(ts)
    n_lo = sum(1 for t in ts if t == lo)
    if n_lo == 4:
        return frozenset(v for v, t in zip(vs, ts) if t == lo)
    if n_lo == 1:
        hi = max(ts)
        if sum(1 for t in ts if t == hi) == 4:
            return frozenset(v for v, t in zip(vs, ts) if t == hi)
    return None


def tris_of(tet):
    return [frozenset(c) for c in combinations(sorted(tet), 3)]


# ---------------------------------------------------------------------------
# Face labels (per slice triangle) -- with checkpoint (de)serialize
# ---------------------------------------------------------------------------

class FaceLabels:
    """m-labels (stored 0..6, m = l-3) on the spatial triangles of the slices.
    Keys are frozensets of 3 vertex ids -- stable across moves (vertex ids
    persist; a re-created spatial tet maps onto the same triangle keys)."""
    def __init__(self, rng):
        self.rng = rng
        self.lab = {}
        self.pinned = set()       # triangles whose labels are FROZEN (defects)
        self.mu_saved = None
        self.mu_ctx = None

    def ensure(self, tris):
        for tri in tris:
            if tri not in self.lab:
                self.lab[tri] = int(self.rng.integers(7))

    def prune(self, live_tris):
        for k in list(self.lab):
            if k not in live_tris:
                del self.lab[k]

    def serialize(self):
        out = [[*sorted(k), int(c)] for k, c in self.lab.items()]
        d = {"kind": "closure", "lab": out}
        if self.pinned:
            d["pinned"] = [sorted(k) for k in self.pinned]
        if self.mu_saved is not None:
            d["mu"] = float(self.mu_saved)
            if self.mu_ctx is not None:
                d["mu_ctx"] = dict(self.mu_ctx)
        return d

    def load(self, d):
        """Tolerant: ignores foreign payloads (e.g. an EPRL intertwiner
        payload from a v6_theory_run checkpoint). Does not touch mu_saved --
        the driver resolves mu before the run starts."""
        if not isinstance(d, dict) or d.get("kind") != "closure":
            print("# [closure labels] resume payload is not a closure payload "
                  "-- starting with fresh labels", flush=True)
            return
        pin_labels = {k: self.lab[k] for k in self.pinned if k in self.lab}
        self.lab = {frozenset((int(a), int(b), int(c))): int(l)
                    for a, b, c, l in d["lab"]}
        # pins applied BEFORE a resume survive the load: frozen labels win
        # over the payload, and pin sets merge (stage 3 pins a fresh copy of
        # a thermalized stage-2 checkpoint).
        self.lab.update(pin_labels)
        self.pinned |= {frozenset(map(int, k)) for k in d.get("pinned", [])}


# ---------------------------------------------------------------------------
# Incremental closure action (delta_action interface for the run harness)
# ---------------------------------------------------------------------------

class IncrementalClosure:
    """Full geometry action  Regge + beta * sum_tets (E(tet) - mu)  maintained
    incrementally from the engine change-log. Tets are keyed by their vertex
    frozenset, so undo's re-pid'd pentachora reconcile onto the same keys.
    Geometry moves never change EXISTING labels (labels live on triangles),
    so only created/destroyed tets need reconciling -- simpler than EPRL."""
    def __init__(self, labels, Etab, k0, Delta, k4, beta, centering=None):
        self.labels = labels
        self.E = Etab
        self.k0, self.Delta, self.k4 = k0, Delta, k4
        self.beta = beta
        self.centering = centering if centering is not None else Centering(enabled=False)
        self.tet_pids = {}      # tet -> set of pids containing it
        self.pid_tet = {}       # pid -> tet (only {4,1}-type pids appear)
        self.contrib = {}       # tet -> centered energy
        self.S_clo = 0.0
        self.regge_prev = 0.0
        self._pending_regge = 0.0

    # --- low level ---------------------------------------------------------
    def _regge(self, T):
        return regge_action(T, self.k0, self.Delta, self.k4)

    def _energy(self, tet):
        l = [self.labels.lab[tri] for tri in tris_of(tet)]
        return float(self.E[l[0], l[1], l[2], l[3]])

    def _contrib(self, tet):
        return self._energy(tet) - self.centering.value()

    # --- driver interface --------------------------------------------------
    def full(self, T, prune=False):
        """Recompute from scratch. prune=False (default) is add-only on the
        label state so a REJECTED proposal evaluated through here (global
        path) cannot mutate surviving labels; refresh() prunes on accepted
        states only."""
        self.tet_pids = {}
        self.pid_tet = {}
        for p, vs in T.pent.items():
            tet = tet_of(T, vs)
            if tet is not None:
                self.tet_pids.setdefault(tet, set()).add(p)
                self.pid_tet[p] = tet
        live = set()
        for tet in self.tet_pids:
            ts = tris_of(tet)
            live.update(ts)
            self.labels.ensure(ts)
        if prune:
            self.labels.prune(live)
        raw = {tet: self._energy(tet) for tet in self.tet_pids}
        self.centering.calibrate(raw.values())
        mu = self.centering.value()
        self.contrib = {tet: e - mu for tet, e in raw.items()}
        self.S_clo = sum(self.contrib.values())
        self.regge_prev = self._regge(T)
        return self.regge_prev + self.beta * self.S_clo

    def _reconcile(self, T, added, removed):
        for rec in removed:
            pid = rec[0]
            tet = self.pid_tet.pop(pid, None)
            if tet is None:
                continue
            s = self.tet_pids.get(tet)
            if s is None:
                continue
            s.discard(pid)
            if not s:
                del self.tet_pids[tet]
                self.S_clo -= self.contrib.pop(tet)
        for pid in added:
            if pid not in T.pent:
                continue
            tet = tet_of(T, T.pent[pid])
            if tet is None:
                continue
            self.pid_tet[pid] = tet
            s = self.tet_pids.setdefault(tet, set())
            s.add(pid)
            if tet not in self.contrib:
                self.labels.ensure(tris_of(tet))
                c = self._contrib(tet)
                self.contrib[tet] = c
                self.S_clo += c

    def delta(self, T, added, removed):
        pre = self.regge_prev + self.beta * self.S_clo
        self._reconcile(T, added, removed)
        regge_now = self._regge(T)
        post = regge_now + self.beta * self.S_clo
        self._pending_regge = regge_now
        return post - pre

    def accept(self):
        self.regge_prev = self._pending_regge

    def reject(self, T, undo_added, undo_removed):
        self._reconcile(T, undo_added, undo_removed)
        self.regge_prev = self._regge(T)

    def refresh(self, T):
        """After a heat-bath pass changes labels, rebuild energies (and prune
        stale labels -- called on accepted states only)."""
        self.full(T, prune=True)

    def audit(self, T):
        running = self.regge_prev + self.beta * self.S_clo
        s = 0.0
        seen = set()
        for p, vs in T.pent.items():
            tet = tet_of(T, vs)
            if tet is not None and tet not in seen:
                seen.add(tet)
                s += self._contrib(tet)
        recomputed = self._regge(T) + self.beta * s
        return running, recomputed

    # --- closure observables (for logging) -----------------------------
    def stats(self):
        if not self.contrib:
            return {}
        mu = self.centering.value()
        e = np.array([c + mu for c in self.contrib.values()])
        colls = 0
        for tet in self.tet_pids:
            l = [self.labels.lab[tri] for tri in tris_of(tet)]
            colls += (4 - len(set(l))) > 0
        return {"n_tets": len(self.contrib), "mean_E": float(e.mean()),
                "std_E": float(e.std()),
                "frac_tets_with_collision": colls / len(self.contrib)}


# ---------------------------------------------------------------------------
# Label heat-bath (beta-consistent from day one) + mu thermodynamic integration
# ---------------------------------------------------------------------------

def _tri_map(tet_pids):
    m = {}
    for tet in tet_pids:
        for tri in tris_of(tet):
            m.setdefault(tri, []).append(tet)
    return m

def _heatbath_pass(labels, tet_pids, Etab, beta, rng):
    """Resample every triangle label from its exact conditional
    prod_{tets containing it} exp(-beta E). PINNED triangles (defect
    carriers, stage 3) are never resampled."""
    tmap = _tri_map(tet_pids)
    for tri, tets in tmap.items():
        if tri in labels.pinned:
            continue
        logw = np.zeros(7)
        old = labels.lab[tri]
        for c in range(7):
            labels.lab[tri] = c
            s = 0.0
            for tet in tets:
                l = [labels.lab[t2] for t2 in tris_of(tet)]
                s += Etab[l[0], l[1], l[2], l[3]]
            logw[c] = -beta * s
        logw -= logw.max()
        w = np.exp(logw); tot = w.sum()
        labels.lab[tri] = old if tot < 1e-300 else int(rng.choice(7, p=w / tot))


def make_heatbath(inc):
    """Per-sweep hook: one heat-bath pass over all triangle labels at the
    run's coupling beta. (The harness calls inc.refresh(T) right after.)"""
    def hook(T, rng):
        _heatbath_pass(inc.labels, inc.tet_pids, inc.E, inc.beta, rng)
    return hook


def calibrate_mu_ti(T, Etab, beta, seed, n_points=6, equil=3, measure=2,
                    verbose=True):
    """Volume-neutral centering constant: the per-tet label FREE ENERGY
    mu(beta) = (1/beta) int_0^beta <E_bar>_s ds by annealed thermodynamic
    integration on the starting geometry (same construction as the EPRL fix;
    see v6_theory_run.calibrate_mu_ti). At beta<=0: uniform-label mean."""
    rng = np.random.default_rng(seed)
    labels = FaceLabels(rng)
    tet_pids = {}
    for p, vs in T.pent.items():
        tet = tet_of(T, vs)
        if tet is not None:
            tet_pids.setdefault(tet, set()).add(p)
    for tet in tet_pids:
        labels.ensure(tris_of(tet))

    def mean_E():
        s = 0.0
        for tet in tet_pids:
            l = [labels.lab[tri] for tri in tris_of(tet)]
            s += Etab[l[0], l[1], l[2], l[3]]
        return s / max(1, len(tet_pids))

    if beta <= 0:
        return mean_E()
    ss = np.linspace(0.0, float(beta), int(n_points) + 1)
    fs = []
    if verbose:
        print(f"# [mu TI] closure: {len(tet_pids)} spatial tets, "
              f"{n_points + 1} grid points x {equil}+{measure} passes", flush=True)
    for s in ss:
        for _ in range(int(equil)):
            _heatbath_pass(labels, tet_pids, Etab, s, rng)
        acc = []
        for _ in range(int(measure)):
            _heatbath_pass(labels, tet_pids, Etab, s, rng)
            acc.append(mean_E())
        fs.append(sum(acc) / len(acc))
        if verbose:
            print(f"# [mu TI]   s={s:.4f}  <E_bar>_s = {fs[-1]:.4f}", flush=True)
    trapz = getattr(np, "trapezoid", None) or np.trapz
    mu = float(trapz(fs, ss) / beta)
    if verbose:
        print(f"# [mu TI] mu({beta}) = {mu:.4f}  (uniform {fs[0]:.4f} -> "
              f"equilibrated {fs[-1]:.4f}; mu is the free energy, between)",
              flush=True)
    return mu


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--beta-closure", type=float, default=1.0,
                   help="coupling strength; beta=1 with eta=eta* is the "
                        "theory-designated point, beta=0 the bare control")
    p.add_argument("--eta", type=float, default=ETA_STAR,
                   help="admissibility precision (default: the paper's eta*)")
    p.add_argument("--lambda-inj", type=float, default=3.0,
                   help="soft injectivity penalty per equal-label pair "
                        "(paper's hard constraint = infinity limit)")
    p.add_argument("--placebo", action="store_true",
                   help="PLACEBO arm: shuffle the 210 label-orbit energies "
                        "(closure structure destroyed, value pool kept)")
    p.add_argument("--placebo-seed", type=int, default=12345)
    p.add_argument("--local-closure", action="store_true",
                   help="use the O(footprint) incremental action (recommended); "
                        "else a global O(N) recompute cross-check path")
    p.add_argument("--audit-every", type=int, default=25)
    p.add_argument("--center-closure", action=argparse.BooleanOptionalAction,
                   default=True)
    p.add_argument("--closure-mu", type=float, default=None,
                   help="fixed centering constant (else TI-calibrated / "
                        "reused from a matching checkpoint)")
    p.add_argument("--recalibrate-mu", action="store_true")
    p.add_argument("--mu-ti-points", type=int, default=6)
    p.add_argument("--mu-ti-equil", type=int, default=3)
    p.add_argument("--mu-ti-measure", type=int, default=2)
    p.add_argument("--causal-slices", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="enforce the CDT foliation. The closure model lives ON "
                        "the slices -- do not turn this off for closure runs.")
    p.add_argument("--target-n41", type=int, default=4000)
    p.add_argument("--K", type=int, default=24)
    p.add_argument("--k0", type=float, default=2.2)
    p.add_argument("--Delta", type=float, default=0.6)
    p.add_argument("--k4", type=float, default=0.9)
    p.add_argument("--eps", type=float, default=1e-3)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-sweeps", type=int, default=20000)
    p.add_argument("--measure-every", type=int, default=20)
    p.add_argument("--checkpoint", type=str, default="v6_closure_ckpt.json")
    p.add_argument("--resume", type=str, default=None)
    p.add_argument("--wall-hours", type=float, default=None)
    args = p.parse_args()

    mk2, mk2_pred, g_eff = ensemble_report(args.eta)
    print("!! WITH-THEORY (ADMISSIBILITY-CLOSURE) run -- the paper's own cell "
          "weighting, not a spin-foam amplitude.")
    print(f"   beta={args.beta_closure}  eta={args.eta:.6g}  "
          f"lambda_inj={args.lambda_inj}  placebo={args.placebo}")
    print(f"   [fidelity] injective ensemble at this eta: <K^2>={mk2:.3f} "
          f"(closure relation 3/(2 eta)={mk2_pred:.3f})  entropy={g_eff:.4f} "
          f"(paper: 7.4198 at eta*)")
    if not args.causal_slices:
        print("   !! --no-causal-slices: closure labels assume closed-3-manifold"
              " slices; results on the generalized ensemble are NOT readable.")

    Etab = build_energy_table(args.eta, args.lambda_inj,
                              args.placebo_seed if args.placebo else None)
    if args.placebo:
        print(f"# [placebo] orbit energies shuffled (seed {args.placebo_seed}) "
              f"-- this arm measures machinery + value statistics, not closure.")

    rng = np.random.default_rng(args.seed)
    labels = FaceLabels(rng)
    beta_eff = args.beta_closure

    # --- resolve mu: CLI > matching checkpoint > TI calibration -------------
    import hashlib
    mu_ctx = {"beta_eff": float(beta_eff),
              "etab_sha": hashlib.sha256(Etab.tobytes()).hexdigest()[:12]}
    mu = args.closure_mu
    T0 = None
    if args.center_closure and mu is None and args.resume:
        from v6_run_lib import load_checkpoint
        T0, _, extra0 = load_checkpoint(args.resume)
        if isinstance(extra0, dict) and extra0.get("kind") == "closure" \
                and extra0.get("mu") is not None and not args.recalibrate_mu:
            if extra0.get("mu_ctx") == mu_ctx:
                mu = float(extra0["mu"])
                print(f"# [closure centering] reusing mu = {mu:.4f} from the "
                      f"resume checkpoint", flush=True)
            else:
                print(f"# [closure centering] checkpoint mu context mismatch "
                      f"-- recalibrating", flush=True)
    if args.center_closure and mu is None:
        if T0 is None:
            if args.resume:
                from v6_run_lib import load_checkpoint
                T0, _, _ = load_checkpoint(args.resume)
            else:
                from v6_cdt import build_s1xs3
                T0 = build_s1xs3(K=args.K)
                print("# [mu TI] fresh run: calibrating on the initial (thin) "
                      "geometry; production should resume from a thermalized "
                      "checkpoint.", flush=True)
        mu = calibrate_mu_ti(T0, Etab, beta_eff, seed=args.seed + 777,
                             n_points=args.mu_ti_points,
                             equil=args.mu_ti_equil,
                             measure=args.mu_ti_measure)
    del T0
    centering = Centering(enabled=args.center_closure, mu_fixed=mu,
                          label="closure")
    labels.mu_saved = mu
    labels.mu_ctx = mu_ctx

    inc = IncrementalClosure(labels, Etab, args.k0, args.Delta, args.k4,
                             beta_eff, centering=centering)
    heatbath = make_heatbath(inc)
    common = dict(k0=args.k0, Delta=args.Delta, k4=args.k4,
                  target_N41=args.target_n41, K=args.K, eps=args.eps,
                  seed=args.seed, max_sweeps=args.max_sweeps,
                  measure_every=args.measure_every, checkpoint=args.checkpoint,
                  resume=args.resume, extra_state=labels,
                  causal=args.causal_slices,
                  wall_budget_s=(args.wall_hours * 3600 if args.wall_hours else None))

    if args.local_closure:
        T = run(f"CLOSURE [beta={beta_eff}{', PLACEBO' if args.placebo else ''}]",
                delta_action=inc, extra_hook=heatbath,
                audit_every=args.audit_every, **common)
    else:
        def geom(Tg):
            return inc.full(Tg)
        T = run(f"CLOSURE [beta={beta_eff}, global"
                f"{', PLACEBO' if args.placebo else ''}]",
                geometry_action=geom, extra_hook=heatbath, **common)
        inc.full(T)

    ids, adj = dual_adjacency(T)
    dH = hausdorff_dim(adj); prof = volume_profile(T); pm = profile_metrics(prof)
    rails = torus_rails(T.n_pent())
    st = inc.stats()
    okf, repf = getattr(T, "_final_verify", (None, {}))
    sl_bad, sl_hist, sl_ss = getattr(T, "_final_slices", (None, {}, None))
    print("\n" + "=" * 64)
    print("  WITH-THEORY (CLOSURE) RESULT -- compare to bare + placebo arms")
    print("=" * 64)
    print(f"  beta={beta_eff}  eta={args.eta:.6g}  lambda_inj={args.lambda_inj}"
          f"  placebo={args.placebo}  mu={centering.mu if centering.mu is not None else 'off'}")
    print(f"  manifold check: {'PASS' if okf else 'FAIL'}  "
          f"[links={repf.get('link_failures')}]   foliation: "
          f"{'CLEAN' if sl_bad == 0 and sl_ss == 0 else f'DEFECTIVE ({sl_bad}+{sl_ss})'}")
    print(f"  final N4={T.n_pent()}  N41={T.type_counts()[0]}")
    print(f"  d_H = {dH:.2f}   rails: 2t={rails[2][1]:.2f} 3t={rails[3][1]:.2f} "
          f"4t={rails[4][1]:.2f}")
    print(f"  blob score = {pm['blob_score']:.2f}   active slices = "
          f"{pm['active_slices']}/{T.K}   max slice = {pm['max_slice']:.0f}")
    if pm['cos3_relerr'] is not None:
        print(f"  cos^3 fit: width={pm['cos3_width']:.1f}  "
              f"rel.RMS err={pm['cos3_relerr']:.3f}")
    if st:
        print(f"  closure term: <E>/tet={st['mean_E']:.4f}  std={st['std_E']:.4f}"
              f"  tets with label collision: {st['frac_tets_with_collision']:.3f}"
              f"  (n_tets={st['n_tets']})")
    print("  (read against the beta=0 control and the --placebo arm at the "
          "same beta; only a real-vs-placebo difference is closure physics.)")


if __name__ == "__main__":
    main()
