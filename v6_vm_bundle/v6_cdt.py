#!/usr/bin/env python3
"""v6 -- faithful 4D Causal Dynamical Triangulations engine (the substrate).

WHY v6 EXISTS: v5 represented a triangulation as a set of vertex-TUPLES with
collision checks. That can only encode strict simplicial complexes, which is the
WRONG ensemble for CDT -- it is dominated by branched polymers (we proved v5
sits at d_s~2 at every coupling/scale, with branched spatial slices). Genuine
CDT samples GENERALIZED triangulations, where the same vertex set can label
distinct simplices glued differently, and degenerate identifications are
controlled by `strictness`. The validated reference (JorenB) does exactly this
with explicit face-gluing pointers + strictness checks.

v6 ARCHITECTURE (mirrors the reference, in Python so the EPRL amplitude layers
on natively):
  - A pentachoron (4-simplex) is an OBJECT identified by an integer id, NOT by
    its vertex set. Two distinct pentachora may share a vertex set.
  - Each pentachoron stores vs[0..4] (vertex ids, sorted as labels only) and
    nbr[0..4] (neighbor pentachoron ids); nbr[i] is glued across the tetrahedral
    face OPPOSITE vs[i] (i.e. the 4 vertices vs minus vs[i]).
  - Vertices carry integer time labels; time is periodic (S^1) with period K.
  - CDT simplex types (4,1)/(3,2) by the time-split of the 5 vertices.

This module: core structure + S^1xS^3 init (bootstrapped from a valid
vertex-tuple triangulation, then converted to gluings) + INDEPENDENT
verification + pluggable action. Moves are in v6_cdt_moves.py.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from itertools import combinations
from typing import Optional


def _faces(vs):
    """The 5 tetrahedral faces of a pentachoron, face i = vs without vs[i]
    (opposite vs[i]). Returns list of 5 sorted 4-tuples."""
    return [tuple(vs[j] for j in range(5) if j != i) for i in range(5)]


class Bag:
    """Indexable set: O(1) add / discard / uniform random pick / len."""
    __slots__ = ("items", "pos")

    def __init__(self):
        self.items = []
        self.pos = {}

    def add(self, x):
        if x in self.pos:
            return
        self.pos[x] = len(self.items)
        self.items.append(x)

    def discard(self, x):
        i = self.pos.pop(x, None)
        if i is None:
            return
        last = self.items.pop()
        if i < len(self.items):
            self.items[i] = last
            self.pos[last] = i

    def pick(self, rng):
        return self.items[rng.integers(len(self.items))]

    def __len__(self):
        return len(self.items)


@dataclass
class Triangulation:
    """A 4D CDT configuration in the gluing (object-identity) representation."""
    vtime: dict = field(default_factory=dict)        # vid -> time
    pent: dict = field(default_factory=dict)         # pid -> (5 sorted vids)
    nbr: dict = field(default_factory=dict)          # pid -> [5 neighbor pids]
    vinc: dict = field(default_factory=dict)         # vid -> set of pids containing it
    K: Optional[int] = None                          # periodic time period
    pent_bag: "Bag" = field(default_factory=Bag)     # O(1) random pentachoron
    vert_bag: "Bag" = field(default_factory=Bag)     # O(1) random vertex (in complex)
    n41: int = 0                                     # running count of (4,1)/(1,4)
    n32: int = 0                                     # running count of (3,2)/(2,3)
    nother: int = 0                                  # running count of non-CDT (should stay 0)
    _next_pid: int = 0
    _next_vid: int = 0
    # --- optional change-log (off by default; for local incremental actions) --
    _rec: bool = False
    _added: list = field(default_factory=list)
    _removed: list = field(default_factory=list)   # (pid, vs, nbr_at_removal)

    def begin_record(self):
        """Start/clear the add/remove change-log for the next surgery bracket."""
        self._rec = True
        self._added = []
        self._removed = []

    def take_record(self):
        """Return (added_pids, removed_records) since begin_record and clear."""
        a, r = self._added, self._removed
        self._added = []
        self._removed = []
        return a, r

    def _classify(self, vs):
        ts = sorted({self.vtime[v] for v in vs})
        if len(ts) != 2:
            return "nother"
        a = sum(1 for v in vs if self.vtime[v] == ts[0])
        if {a, 5 - a} == {4, 1}:
            return "n41"
        if {a, 5 - a} == {3, 2}:
            return "n32"
        return "nother"

    # --- construction ------------------------------------------------------
    def new_vertex(self, time: int) -> int:
        v = self._next_vid
        self._next_vid += 1
        self.vtime[v] = time
        # vinc[v] + vert_bag are populated when the first pentachoron is added
        return v

    def add_pentachoron(self, vs) -> int:
        """Add a pentachoron from 5 vertex ids; neighbors set later/by surgery."""
        pid = self._next_pid
        self._next_pid += 1
        vs = tuple(sorted(vs))
        self.pent[pid] = vs
        self.nbr[pid] = [-1, -1, -1, -1, -1]
        for v in vs:
            s = self.vinc.get(v)
            if s is None:
                self.vinc[v] = {pid}
                self.vert_bag.add(v)
            else:
                s.add(pid)
        self.pent_bag.add(pid)
        setattr(self, self._classify(vs), getattr(self, self._classify(vs)) + 1)
        if self._rec:
            self._added.append(pid)
        return pid

    def remove_pentachoron(self, pid):
        if self._rec:
            self._removed.append((pid, self.pent[pid], list(self.nbr[pid])))
        cls = self._classify(self.pent[pid])
        setattr(self, cls, getattr(self, cls) - 1)
        self.pent_bag.discard(pid)
        for v in self.pent[pid]:
            self.vinc[v].discard(pid)
            if not self.vinc[v]:
                del self.vinc[v]
                self.vert_bag.discard(v)
        del self.pent[pid]
        del self.nbr[pid]

    def has_edge(self, u, w) -> bool:
        """True if some pentachoron contains both vertices u and w."""
        su = self.vinc.get(u)
        sw = self.vinc.get(w)
        if not su or not sw:
            return False
        return not su.isdisjoint(sw)

    # --- typing / foliation ------------------------------------------------
    def ptype(self, pid):
        """(a,b) time-split of a pentachoron's 5 vertices across 2 slices.
        Valid CDT types are {4,1} and {3,2} (as sets). Periodic-time aware:
        a simplex spanning slices {0, K-1} is adjacent under S^1."""
        ts = sorted({self.vtime[v] for v in self.pent[pid]})
        if len(ts) == 1:
            return (5, 0)
        if len(ts) == 2:
            lo, hi = ts
            # treat the wraparound pair {0,K-1} as adjacent (handled by caller)
            a = sum(1 for v in self.pent[pid] if self.vtime[v] == lo)
            return (a, 5 - a)
        return (-1, -1)  # spans >2 slices -> non-causal

    def type_counts(self):
        """O(1): the running (N_41, N_32, N_other) maintained on add/remove."""
        return self.n41, self.n32, self.nother

    def _type_counts_recompute(self):
        """O(N) independent recompute (used by verify to cross-check)."""
        n41 = n32 = nother = 0
        for pid in self.pent:
            a, b = self.ptype(pid)
            if {a, b} == {4, 1}:
                n41 += 1
            elif {a, b} == {3, 2}:
                n32 += 1
            else:
                nother += 1
        return n41, n32, nother

    def n_pent(self):
        return len(self.pent)

    def f_vector(self):
        """(N0,N1,N2,N3,N4) computed from the pentachora (vertex labels used to
        enumerate sub-simplices -- valid for the init simplicial complex; for
        generalized configs N1,N2,N3 are counted as distinct gluing faces)."""
        N4 = len(self.pent)
        tets = set(); tris = set(); edges = set(); verts = set()
        for vs in self.pent.values():
            verts.update(vs)
            for f in combinations(vs, 4):
                tets.add(f)
            for f in combinations(vs, 3):
                tris.add(f)
            for f in combinations(vs, 2):
                edges.add(f)
        return (len(verts), len(edges), len(tris), len(tets), N4)


# ---------------------------------------------------------------------------
# Init: S^1 x S^3 staircase, bootstrapped from a valid vertex-tuple complex
# ---------------------------------------------------------------------------

def build_s1xs3(K: int = 6) -> Triangulation:
    """Build a valid S^1xS^3 CDT in the gluing representation. We reuse the
    (validated) v5 staircase vertex-tuple triangulation to get correct topology,
    then DERIVE the gluing pointers by matching each tetrahedral face to its two
    incident pentachora. After this, the structure is maintained purely by
    neighbor-pointer surgery (moves), never re-derived from vertex sets -- so it
    can become a generalized triangulation."""
    from v5_s1xs3_init import build_s1xs3_staircase
    T5 = build_s1xs3_staircase(K=K)

    T = Triangulation(K=K)
    # carry over vertex times
    for v, t in T5.vertex_time.items():
        T.vtime[v] = t
    T._next_vid = max(T5.vertex_time) + 1
    # add pentachora (preserve vertex-tuples)
    tuple_to_pid = {}
    for s in T5.simplices:
        pid = T.add_pentachoron(s)
        tuple_to_pid[tuple(sorted(s))] = pid
    # derive gluings: each tet face (4-vertex set) is shared by exactly 2 pents
    face_map = {}
    for pid, vs in T.pent.items():
        for i, f in enumerate(_faces(vs)):
            face_map.setdefault(f, []).append((pid, i))
    for f, inc in face_map.items():
        if len(inc) != 2:
            raise ValueError(f"init not a closed manifold: face {f} in {len(inc)} pentachora")
        (p1, i1), (p2, i2) = inc
        T.nbr[p1][i1] = p2
        T.nbr[p2][i2] = p1
    return T


# ---------------------------------------------------------------------------
# Independent verification
# ---------------------------------------------------------------------------

def gluing_audit(T: Triangulation):
    """Structural audit from the GLUING POINTERS ALONE -- valid even for
    generalized (non-simplicial) triangulations, because it never identifies a
    sub-simplex by its vertex set. Returns a dict of vertex-tuple-INDEPENDENT
    facts + a verdict of whether the current config happens to be simplicial
    (which is the only regime in which the vertex-tuple f_vector / link checks
    below are trustworthy)."""
    from collections import Counter
    N4 = len(T.pent)
    bad_nbr = tadpoles = self_energy = glued_slots = 0
    for pid, ns in T.nbr.items():
        seen = {}
        for q in ns:
            if q == -1 or q not in T.pent:
                bad_nbr += 1
                continue
            glued_slots += 1
            if q == pid:
                tadpoles += 1
            if pid not in T.nbr.get(q, ()):       # back-pointer must exist
                bad_nbr += 1
            seen[q] = seen.get(q, 0) + 1
        self_energy += sum(c - 1 for c in seen.values() if c > 1)
    # dual-graph connectivity
    connected = True
    if N4:
        start = next(iter(T.pent)); seenp = {start}; stack = [start]
        while stack:
            u = stack.pop()
            for w in T.nbr[u]:
                if w in T.pent and w not in seenp:
                    seenp.add(w); stack.append(w)
        connected = (len(seenp) == N4)
    # is the gluing consistent with vertex labels (=> simplicial)?
    fv_mismatch = 0
    for pid, vs in T.pent.items():
        for i, q in enumerate(T.nbr[pid]):
            if q == -1 or q not in T.pent:
                continue
            face_p = frozenset(vs[j] for j in range(5) if j != i)
            try:
                j = T.nbr[q].index(pid)
            except ValueError:
                continue
            qs = T.pent[q]
            if face_p != frozenset(qs[k] for k in range(5) if k != j):
                fv_mismatch += 1
    setc = Counter(frozenset(vs) for vs in T.pent.values())
    dup_pent = sum(c - 1 for c in setc.values() if c > 1)
    degen_pent = sum(1 for vs in T.pent.values() if len(set(vs)) != 5)
    is_simplicial = (fv_mismatch == 0 and dup_pent == 0 and degen_pent == 0)
    return {
        "N4": N4, "glued_slots": glued_slots, "n_tets": glued_slots // 2,
        "all_faces_glued": (glued_slots == 5 * N4),
        "bad_neighbor_pointers": bad_nbr, "tadpoles": tadpoles,
        "self_energy_pairs": self_energy, "connected": connected,
        "face_vertex_mismatch": fv_mismatch, "dup_pentachora": dup_pent,
        "degenerate_pentachora": degen_pent, "is_simplicial": is_simplicial,
    }


def verify(T: Triangulation, check_links: bool = True):
    """Independently verify T is a valid closed CDT 4-manifold. Returns
    (ok, report). The PRIMARY checks are gluing-based (always valid); the
    vertex-tuple f_vector / link checks are included only as a cross-check and
    are trusted in the verdict ONLY when the config is provably simplicial."""
    rep = {}
    ga = gluing_audit(T)
    rep["gluing"] = ga
    # gluing-based necessary manifold conditions (ALWAYS valid)
    gluing_ok = (ga["bad_neighbor_pointers"] == 0 and ga["all_faces_glued"]
                 and ga["tadpoles"] == 0 and ga["self_energy_pairs"] == 0
                 and ga["connected"])
    rep["ok_gluing_only"] = gluing_ok
    # CDT types (independent recompute vs running counts) -- gluing-independent
    n41, n32, nother = T._type_counts_recompute()
    rep["type_counts"] = (n41, n32, nother)
    rep["running_counts_match"] = ((n41, n32, nother) == T.type_counts())
    # vertex-tuple extras (RELIABLE ONLY WHEN SIMPLICIAL)
    fv = T.f_vector()
    chi = fv[0] - fv[1] + fv[2] - fv[3] + fv[4]
    rep["f_vector"] = fv
    rep["euler_char"] = chi
    rep["vertex_tuple_checks_reliable"] = ga["is_simplicial"]
    rep["link_failures"] = None
    if check_links:
        rep["link_failures"] = (_count_link_failures(T) if ga["is_simplicial"]
                                else "unreliable(generalized)")
    ok = (gluing_ok and nother == 0 and rep["running_counts_match"])
    if ga["is_simplicial"]:
        # vertex-tuple invariants are valid here, so demand them too
        ok = ok and (chi == 0) and (rep["link_failures"] in (0, None))
    return ok, rep


def _count_link_failures(T: Triangulation):
    """Build each vertex's link (3-complex of opposite tetrahedra) and check the
    necessary S^3 conditions (closed 3-manifold + chi 0), reusing v5_cdt_link
    via an adapter."""
    from v5_cdt_link import check_link_is_s3_necessary_conditions, CausalTriangulation
    # adapter: build a CausalTriangulation view (vertex_time + simplices set)
    adapter = CausalTriangulation()
    adapter.vertex_time = dict(T.vtime)
    for vs in T.pent.values():
        adapter.add_4simplex(vs)
    bad = 0
    for v in list(T.vtime.keys()):
        if v not in adapter.vertex_time:
            continue
        ok, _ = check_link_is_s3_necessary_conditions(adapter, v)
        if not ok:
            bad += 1
    return bad


# ---------------------------------------------------------------------------
# Pluggable action (Regge default; EPRL amplitude layers here later)
# ---------------------------------------------------------------------------

def regge_action(T: Triangulation, kappa_0: float, Delta: float, kappa_4: float):
    """Bare 4D CDT Regge action (AJL form):
    S = -(kappa_0 + 6 Delta) N_0 + kappa_4 (N_41 + N_32) + Delta N_41.

    O(1): N_0 = len(vinc), and N_41/N_32 are the running counts."""
    N0 = len(T.vinc)
    n41, n32, _ = T.type_counts()
    return -(kappa_0 + 6 * Delta) * N0 + kappa_4 * (n41 + n32) + Delta * n41


if __name__ == "__main__":
    print("=" * 72)
    print("  v6 gluing-based 4D CDT -- foundation self-test")
    print("=" * 72)
    for K in (4, 6, 10):
        T = build_s1xs3(K=K)
        ok, rep = verify(T)
        print(f"\nK={K}: N4={T.n_pent()}  verify={'PASS' if ok else 'FAIL'}")
        for k, v in rep.items():
            print(f"    {k}: {v}")
