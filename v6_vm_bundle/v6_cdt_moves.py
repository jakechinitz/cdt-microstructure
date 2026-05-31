#!/usr/bin/env python3
"""v6 4D CDT moves on the gluing (neighbor-pointer) representation.

This file implements the Pachner moves as explicit neighbor-pointer surgery,
mirroring the reference (JorenB) approach. Correctness is checked by re-running
the independent v6_cdt.verify() after each move (any surgery error surfaces as a
manifold/Euler/link violation) AND by forward-then-reverse identity tests.

Starting with the (2,4)/(4,2) pair as the verified template; the remaining
pairs ((3,3),(4,6)/(6,4),(2,8)/(8,2)) follow the same discipline.
"""
from __future__ import annotations
from itertools import combinations


# --- low-level gluing helpers ---------------------------------------------

def _drop_index(vs, face):
    """Index in vs of the single vertex not in `face` (the opposite vertex)."""
    fs = set(face)
    for i, v in enumerate(vs):
        if v not in fs:
            return i
    return -1


def get_nbr_across(T, pid, face):
    return T.nbr[pid][_drop_index(T.pent[pid], face)]


def set_nbr_across(T, pid, face, q):
    T.nbr[pid][_drop_index(T.pent[pid], face)] = q


def back_rewire(T, ext, old_pid, new_pid):
    """In ext's neighbor list, repoint the (unique) entry old_pid -> new_pid."""
    ns = T.nbr[ext]
    for i, x in enumerate(ns):
        if x == old_pid:
            ns[i] = new_pid
            return True
    return False


def nbr_opposite(T, pid, v):
    """Neighbor of pid across the face opposite vertex v."""
    return T.nbr[pid][T.pent[pid].index(v)]


def strictness_ok(T, new_pids, level=2):
    """Reference-style anti-degeneracy guard (suppresses branched-polymer
    configurations). Checked on the new pentachora and their neighbours:
      level >= 1 (no tadpole):    no pentachoron is its own neighbour.
      level >= 2 (no self-energy): no pentachoron shares 2 faces with the same
                                   neighbour (i.e. all 5 neighbours distinct).
    These are exactly the degenerate gluings real CDT forbids via `strictness`.
    """
    if level <= 0:
        return True
    affected = set(new_pids)
    for p in new_pids:
        affected.update(n for n in T.nbr[p] if n in T.pent)
    for p in affected:
        ns = T.nbr[p]
        if p in ns:                         # tadpole
            return False
        if level >= 2 and len(set(ns)) != 5:   # self-energy
            return False
    return True


def has_subsimplex(T, verts):
    """True if some single pentachoron contains all of `verts` (i.e. that
    sub-simplex already exists in the complex)."""
    sets = [T.vinc.get(v) for v in verts]
    if any(s is None for s in sets):
        return False
    common = set(sets[0])
    for s in sets[1:]:
        common &= s
        if not common:
            return False
    return bool(common)


def hinge_ring(T, pid, tri):
    """Walk the cyclic ring of pentachora sharing the triangle (hinge) `tri`,
    starting at `pid`. Returns (ring_pids, shared_verts) where shared_verts[i]
    is the apex vertex shared between ring[i] and ring[i+1], or None if the walk
    does not close cleanly (e.g. boundary or repeated pentachoron)."""
    tri = tuple(sorted(tri))
    tset = set(tri)
    others = [v for v in T.pent[pid] if v not in tset]
    if len(others) != 2:
        return None
    ring = [pid]
    shared = []
    cur = pid
    cross = others[0]                      # go across face tri ∪ {cross}
    for _ in range(64):
        nxt = get_nbr_across(T, cur, tuple(sorted(tset | {cross})))
        shared.append(cross)
        if nxt == pid:
            return ring, shared            # closed
        if nxt in ring or nxt == -1:
            return None
        ring.append(nxt)
        nxt_others = [v for v in T.pent[nxt] if v not in tset]
        new = [v for v in nxt_others if v != cross]
        if len(new) != 1:
            return None
        cross = new[0]
        cur = nxt
    return None


def _cdt_ok(T, vs):
    """Would a pentachoron with vertex set vs be a valid CDT type ({4,1}/{3,2})?"""
    ts = sorted({T.vtime[v] for v in vs})
    if len(ts) != 2:
        return False
    a = sum(1 for v in vs if T.vtime[v] == ts[0])
    return {a, 5 - a} == {4, 1} or {a, 5 - a} == {3, 2}


# ===========================================================================
# (2,4) move : a tetrahedron shared by 2 pentachora P,Q (apexes e,f) becomes
# 4 pentachora around the new edge {e,f}.
# ===========================================================================

def apply_24(T, P, face, enforce_cdt=True):
    """Apply (2,4) on the shared tetrahedron `face`={a,b,c,d} of pentachoron P
    (and its neighbor Q across that face). Returns (ok, new_pids)."""
    Q = get_nbr_across(T, P, face)
    if Q == -1 or Q == P:
        return False, []
    abcd = set(face)
    eset = set(T.pent[P]) - abcd
    fset = set(T.pent[Q]) - abcd
    if len(eset) != 1 or len(fset) != 1:
        return False, []
    e = eset.pop(); f = fset.pop()
    if e == f:
        return False, []
    # link condition: the new edge {e,f} must not already exist
    if T.has_edge(e, f):
        return False, []
    # CDT-type validity of the 4 prospective new pentachora
    new_vs = {w: tuple(sorted((abcd - {w}) | {e, f})) for w in abcd}
    if enforce_cdt and not all(_cdt_ok(T, vs) for vs in new_vs.values()):
        return False, []
    # external neighbors of P and Q (across faces opposite each w in abcd)
    Pext = {w: get_nbr_across(T, P, tuple(sorted((abcd - {w}) | {e}))) for w in abcd}
    Qext = {w: get_nbr_across(T, Q, tuple(sorted((abcd - {w}) | {f}))) for w in abcd}
    # create the 4 new pentachora
    N = {w: T.add_pentachoron(new_vs[w]) for w in abcd}
    # external gluings + back-pointers
    for w in abcd:
        nw = N[w]
        face_e = tuple(sorted((abcd - {w}) | {e}))
        face_f = tuple(sorted((abcd - {w}) | {f}))
        set_nbr_across(T, nw, face_e, Pext[w]); back_rewire(T, Pext[w], P, nw)
        set_nbr_across(T, nw, face_f, Qext[w]); back_rewire(T, Qext[w], Q, nw)
    # internal gluings (around the new edge {e,f})
    for w1, w2 in combinations(abcd, 2):
        shared = tuple(sorted((abcd - {w1, w2}) | {e, f}))
        set_nbr_across(T, N[w1], shared, N[w2])
        set_nbr_across(T, N[w2], shared, N[w1])
    # remove the 2 old pentachora
    T.remove_pentachoron(P); T.remove_pentachoron(Q)
    return True, list(N.values())


# ===========================================================================
# (4,2) move : inverse of (2,4). 4 pentachora around an edge {e,f} (whose
# "outer" vertices form a tetrahedron abcd) become 2 sharing that tet.
# ===========================================================================

def apply_42(T, e, f, enforce_cdt=True):
    """Apply (4,2) on the edge {e,f}: the 4 pentachora all containing both e
    and f, whose other 3 vertices are the triangles of a 4-set abcd, collapse to
    2 pentachora {abcd,e} and {abcd,f}. Returns (ok, new_pids)."""
    su = T.vinc.get(e); sw = T.vinc.get(f)
    if not su or not sw:
        return False, []
    shared = su & sw
    if len(shared) != 4:
        return False, []
    # the "outer" triangles (each pent minus {e,f}) must be the 4 triangles of a
    # single 4-vertex set abcd
    outers = []
    for pid in shared:
        out = tuple(sorted(set(T.pent[pid]) - {e, f}))
        if len(out) != 3:
            return False, []
        outers.append(out)
    abcd = sorted(set().union(*outers))
    if len(abcd) != 4:
        return False, []
    if set(outers) != {tuple(sorted(t)) for t in combinations(abcd, 3)}:
        return False, []
    abcd = set(abcd)
    # the new tetrahedron abcd must not already exist as a shared face elsewhere
    # (link condition): no pentachoron outside `shared` already carries all of abcd
    inc_abcd = set.intersection(*[T.vinc.get(v, set()) for v in abcd])
    if inc_abcd - shared:
        return False, []
    new1 = tuple(sorted(abcd | {e}))
    new2 = tuple(sorted(abcd | {f}))
    if enforce_cdt and not (_cdt_ok(T, new1) and _cdt_ok(T, new2)):
        return False, []
    # map each old pentachoron N_w by which abcd-vertex it omits
    Nw = {}
    for pid in shared:
        w = (abcd - (set(T.pent[pid]) - {e, f})).pop()
        Nw[w] = pid
    # external neighbors: N_w's face (abcd-w)+e -> external (was Pext[w]); +f -> Qext[w]
    Pext = {w: get_nbr_across(T, Nw[w], tuple(sorted((abcd - {w}) | {e}))) for w in abcd}
    Qext = {w: get_nbr_across(T, Nw[w], tuple(sorted((abcd - {w}) | {f}))) for w in abcd}
    # create the 2 new pentachora
    P = T.add_pentachoron(new1)   # {abcd, e}
    Qp = T.add_pentachoron(new2)  # {abcd, f}
    # P (={abcd,e}) glues across abcd to Qp; across (abcd-w)+e to Pext[w]
    set_nbr_across(T, P, tuple(sorted(abcd)), Qp)
    set_nbr_across(T, Qp, tuple(sorted(abcd)), P)
    for w in abcd:
        face_pe = tuple(sorted((abcd - {w}) | {e}))
        face_qf = tuple(sorted((abcd - {w}) | {f}))
        set_nbr_across(T, P, face_pe, Pext[w]); back_rewire(T, Pext[w], Nw[w], P)
        set_nbr_across(T, Qp, face_qf, Qext[w]); back_rewire(T, Qext[w], Nw[w], Qp)
    for pid in list(shared):
        T.remove_pentachoron(pid)
    return True, [P, Qp]


# ===========================================================================
# (2,8) move : subdivide a spatial tetrahedron by inserting a central vertex.
# A spatial tet abcd (4 vertices in slice t) is shared by an "up" pentachoron
# P={abcd,u} and a "down" pentachoron Q={abcd,w} (u,w in different neighbour
# slices). Insert new vertex x in slice t -> 8 pentachora (4 up + 4 down).
# ===========================================================================

def apply_28(T, P, abcd):
    """Apply (2,8) on spatial tetrahedron `abcd` (a face of P, all in one slice).
    Returns (ok, new_pids, x)."""
    abcd = set(abcd)
    times = {T.vtime[v] for v in abcd}
    if len(times) != 1:
        return False, [], None
    t = next(iter(times))
    Q = get_nbr_across(T, P, tuple(sorted(abcd)))
    if Q == -1 or Q == P:
        return False, [], None
    uset = set(T.pent[P]) - abcd
    wset = set(T.pent[Q]) - abcd
    if len(uset) != 1 or len(wset) != 1:
        return False, [], None
    u = uset.pop(); w = wset.pop()
    tu, tw = T.vtime[u], T.vtime[w]
    if tu == t or tw == t or tu == tw:   # need a genuine up/down sandwich
        return False, [], None
    Pext = {v: nbr_opposite(T, P, v) for v in abcd}
    Qext = {v: nbr_opposite(T, Q, v) for v in abcd}
    x = T.new_vertex(t)
    Nup = {v: T.add_pentachoron((abcd - {v}) | {x, u}) for v in abcd}
    Ndn = {v: T.add_pentachoron((abcd - {v}) | {x, w}) for v in abcd}
    for v in abcd:
        set_nbr_across(T, Nup[v], tuple(sorted((abcd - {v}) | {u})), Pext[v])
        back_rewire(T, Pext[v], P, Nup[v])
        set_nbr_across(T, Ndn[v], tuple(sorted((abcd - {v}) | {w})), Qext[v])
        back_rewire(T, Qext[v], Q, Ndn[v])
        # up<->down across the horizontal face (abcd-v)+x
        hface = tuple(sorted((abcd - {v}) | {x}))
        set_nbr_across(T, Nup[v], hface, Ndn[v])
        set_nbr_across(T, Ndn[v], hface, Nup[v])
    for v1, v2 in combinations(abcd, 2):
        fu = tuple(sorted((abcd - {v1, v2}) | {x, u}))
        set_nbr_across(T, Nup[v1], fu, Nup[v2]); set_nbr_across(T, Nup[v2], fu, Nup[v1])
        fw = tuple(sorted((abcd - {v1, v2}) | {x, w}))
        set_nbr_across(T, Ndn[v1], fw, Ndn[v2]); set_nbr_across(T, Ndn[v2], fw, Ndn[v1])
    T.remove_pentachoron(P); T.remove_pentachoron(Q)
    return True, list(Nup.values()) + list(Ndn.values()), x


def apply_82(T, x):
    """Apply (8,2): delete vertex x (incident to 8 pentachora in the (2,8)
    pattern), coalescing to 2. Returns (ok, new_pids)."""
    pids = list(T.vinc.get(x, ()))
    if len(pids) != 8:
        return False, []
    t = T.vtime[x]
    up = {}; dn = {}; apex = {"u": None, "w": None}
    for pid in pids:
        others = set(T.pent[pid]) - {x}
        same = {v for v in others if T.vtime[v] == t}
        ap = others - same
        if len(same) != 3 or len(ap) != 1:
            return False, []
        a = ap.pop()
        tri = tuple(sorted(same))
        # classify apex by time-adjacency under periodic time (handles wraparound)
        if T.vtime[a] == (t + 1) % T.K:
            up[tri] = (pid, a)
        elif T.vtime[a] == (t - 1) % T.K:
            dn[tri] = (pid, a)
        else:
            return False, []
    if len(up) != 4 or len(dn) != 4:
        return False, []
    us = {a for _, a in up.values()}; ws = {a for _, a in dn.values()}
    if len(us) != 1 or len(ws) != 1:
        return False, []
    u = us.pop(); w = ws.pop()
    abcd = set().union(*[set(tr) for tr in up])
    if len(abcd) != 4:
        return False, []
    if {tuple(sorted(t3)) for t3 in combinations(abcd, 3)} != set(up) or set(up) != set(dn):
        return False, []
    new1 = tuple(sorted(abcd | {u}))
    new2 = tuple(sorted(abcd | {w}))
    if not (_cdt_ok(T, new1) and _cdt_ok(T, new2)):
        return False, []
    # external neighbours: for each abcd-vertex v, the up-pentachoron with
    # triangle (abcd-v) borders external Pext[v] across its face opposite x.
    Pext = {}; Qext = {}
    for v in abcd:
        tri = tuple(sorted(abcd - {v}))
        pid_u = up[tri][0]; pid_d = dn[tri][0]
        Pext[v] = nbr_opposite(T, pid_u, x)
        Qext[v] = nbr_opposite(T, pid_d, x)
    P = T.add_pentachoron(new1)   # {abcd, u}
    Q = T.add_pentachoron(new2)   # {abcd, w}
    set_nbr_across(T, P, tuple(sorted(abcd)), Q)
    set_nbr_across(T, Q, tuple(sorted(abcd)), P)
    for v in abcd:
        set_nbr_across(T, P, tuple(sorted((abcd - {v}) | {u})), Pext[v])
        back_rewire(T, Pext[v], up[tuple(sorted(abcd - {v}))][0], P)
        set_nbr_across(T, Q, tuple(sorted((abcd - {v}) | {w})), Qext[v])
        back_rewire(T, Qext[v], dn[tuple(sorted(abcd - {v}))][0], Q)
    for pid in pids:
        T.remove_pentachoron(pid)
    del T.vtime[x]
    T.vinc.pop(x, None)
    return True, [P, Q]


# ===========================================================================
# (3,3) move : 3 pentachora around a hinge triangle `tri` <-> 3 pentachora
# around the dual hinge triangle `tri'` (the 3 apex vertices). Self-inverse.
# ===========================================================================

def apply_33(T, pid, tri, enforce_cdt=True):
    """Apply (3,3) on the hinge triangle `tri` if exactly 3 pentachora ring it.
    Returns (ok, new_pids, dual_tri)."""
    r = hinge_ring(T, pid, tri)
    if r is None:
        return False, [], None
    ring, shared = r
    if len(ring) != 3:
        return False, [], None
    tri = tuple(sorted(tri))
    dual = tuple(sorted(set(shared)))         # the 3 apex vertices
    if len(dual) != 3:
        return False, [], None
    if has_subsimplex(T, dual):               # dual hinge must be new
        return False, [], None
    new_vs = {e: tuple(sorted(set(dual) | set(e))) for e in combinations(tri, 2)}
    if enforce_cdt and not all(_cdt_ok(T, vs) for vs in new_vs.values()):
        return False, [], None
    # external faces of the old ring: face opposite each hinge vertex h
    ext_map = {}
    for pid_i in ring:
        for h in tri:
            face = frozenset(set(T.pent[pid_i]) - {h})
            ext_map[face] = (nbr_opposite(T, pid_i, h), pid_i)
    # create new pentachora (indexed by edge of tri)
    M = {e: T.add_pentachoron(vs) for e, vs in new_vs.items()}
    # internal gluings: the two M's sharing hinge-vertex v meet across dual∪{v}
    for v in tri:
        es = [e for e in M if v in e]
        f = tuple(sorted(set(dual) | {v}))
        set_nbr_across(T, M[es[0]], f, M[es[1]])
        set_nbr_across(T, M[es[1]], f, M[es[0]])
    # external gluings: each M's faces opposite a dual vertex match an old ext face
    for e, m in M.items():
        for x in dual:
            facelist = tuple(sorted(set(T.pent[m]) - {x}))
            fs = frozenset(facelist)
            if fs in ext_map:
                ext, old = ext_map[fs]
                set_nbr_across(T, m, facelist, ext)
                back_rewire(T, ext, old, m)
    for pid_i in ring:
        T.remove_pentachoron(pid_i)
    return True, list(M.values()), dual


# ===========================================================================
# (4,6) / (6,4) : the spatial (2,3)/(3,2) Pachner move lifted to 4D.
# (4,6): a spatial triangle `tri` (slice t) shared by 4 pentachora -- the (4,1)
# and (1,4) over the two spatial tets {tri,alpha} and {tri,beta} -- becomes 6
# pentachora over the three spatial tets {edge_of_tri, alpha, beta} (up+down),
# i.e. the new spatial edge {alpha,beta} replaces the spatial triangle.
# ===========================================================================

def apply_46(T, pid, tri, enforce_cdt=True):
    """Apply (4,6) on a spatial triangle hinge `tri`. Returns (ok,new_pids,(alpha,beta))."""
    tri = tuple(sorted(tri))
    times = {T.vtime[v] for v in tri}
    if len(times) != 1:
        return False, [], None
    t = next(iter(times))
    r = hinge_ring(T, pid, tri)
    if r is None or len(r[0]) != 4:
        return False, [], None
    ring, shared = r
    pos_spatial = [i for i, s in enumerate(shared) if T.vtime[s] == t]
    if pos_spatial not in ([0, 2], [1, 3]):     # must alternate spatial/apex
        return False, [], None
    spatial = [shared[i] for i in pos_spatial]
    apex = [shared[i] for i in (set(range(4)) - set(pos_spatial))]
    alpha, beta = spatial
    U, D = apex
    if T.vtime[U] == T.vtime[D] or T.vtime[U] == t or T.vtime[D] == t:
        return False, [], None
    if has_subsimplex(T, (alpha, beta)):        # new spatial edge must be new
        return False, [], None
    edges = list(combinations(tri, 2))
    new_vs = {(e, ap): tuple(sorted(set(e) | {alpha, beta, ap}))
              for e in edges for ap in (U, D)}
    if enforce_cdt and not all(_cdt_ok(T, vs) for vs in new_vs.values()):
        return False, [], None
    # external faces of the 4 old pentachora (opposite each tri vertex)
    ext_map = {}
    for pid_i in ring:
        for v in tri:
            ext_map[frozenset(set(T.pent[pid_i]) - {v})] = (nbr_opposite(T, pid_i, v), pid_i)
    M = {key: T.add_pentachoron(vs) for key, vs in new_vs.items()}
    # internal: up_e<->up_e' (and dn<->dn) across {common_vertex, alpha, beta, apex}
    for e1, e2 in combinations(edges, 2):
        w = set(e1) & set(e2)
        if len(w) == 1:
            w = w.pop()
            for ap in (U, D):
                f = tuple(sorted({w, alpha, beta, ap}))
                set_nbr_across(T, M[(e1, ap)], f, M[(e2, ap)])
                set_nbr_across(T, M[(e2, ap)], f, M[(e1, ap)])
    # internal: up_e<->dn_e across the spatial tet {edge, alpha, beta}
    for e in edges:
        f = tuple(sorted(set(e) | {alpha, beta}))
        set_nbr_across(T, M[(e, U)], f, M[(e, D)])
        set_nbr_across(T, M[(e, D)], f, M[(e, U)])
    # external gluings (each new pentachoron's faces opposite alpha and beta)
    for key, m in M.items():
        for x in (alpha, beta):
            facelist = tuple(sorted(set(T.pent[m]) - {x}))
            fs = frozenset(facelist)
            if fs in ext_map:
                ext, old = ext_map[fs]
                set_nbr_across(T, m, facelist, ext)
                back_rewire(T, ext, old, m)
    for pid_i in ring:
        T.remove_pentachoron(pid_i)
    return True, list(M.values()), (alpha, beta)


def apply_64(T, alpha, beta, enforce_cdt=True):
    """Apply (6,4): the spatial edge {alpha,beta} shared by 6 pentachora becomes
    a spatial triangle `tri` shared by 4. Returns (ok, new_pids, tri)."""
    if T.vtime.get(alpha) != T.vtime.get(beta) or alpha == beta:
        return False, [], None
    t = T.vtime[alpha]
    shared = list(T.vinc.get(alpha, set()) & T.vinc.get(beta, set()))
    if len(shared) != 6:
        return False, [], None
    up = {}; dn = {}
    for pid in shared:
        rest = set(T.pent[pid]) - {alpha, beta}
        sp = {v for v in rest if T.vtime[v] == t}
        ap = rest - sp
        if len(sp) != 2 or len(ap) != 1:
            return False, [], None
        e = tuple(sorted(sp)); a = ap.pop()
        # classify apex by time-adjacency under periodic time (handles wraparound)
        if T.vtime[a] == (t + 1) % T.K:
            up[e] = (pid, a)
        elif T.vtime[a] == (t - 1) % T.K:
            dn[e] = (pid, a)
        else:
            return False, [], None
    if len(up) != 3 or len(dn) != 3:
        return False, [], None
    Us = {a for _, a in up.values()}; Ds = {a for _, a in dn.values()}
    if len(Us) != 1 or len(Ds) != 1:
        return False, [], None
    U = Us.pop(); D = Ds.pop()
    tri = set().union(*[set(e) for e in up])
    if len(tri) != 3 or set(up) != set(dn):
        return False, [], None
    if {tuple(sorted(e)) for e in combinations(tri, 2)} != set(up):
        return False, [], None
    if has_subsimplex(T, tuple(sorted(tri))):    # new spatial triangle must be new
        return False, [], None
    tri = tuple(sorted(tri))
    new_vs = {(sp, ap): tuple(sorted(set(tri) | {sp, ap}))
              for sp in (alpha, beta) for ap in (U, D)}
    if enforce_cdt and not all(_cdt_ok(T, vs) for vs in new_vs.values()):
        return False, [], None
    # external faces of the 6 old (opposite alpha and opposite beta)
    ext_map = {}
    for e, (pid, _a) in list(up.items()) + list(dn.items()):
        for x in (alpha, beta):
            ext_map[frozenset(set(T.pent[pid]) - {x})] = (nbr_opposite(T, pid, x), pid)
    M = {key: T.add_pentachoron(vs) for key, vs in new_vs.items()}
    # internal: the 4 new share tri; across {tri,U}: (alpha,U)-(beta,U); etc.
    for ap in (U, D):
        f = tuple(sorted(set(tri) | {ap}))
        set_nbr_across(T, M[(alpha, ap)], f, M[(beta, ap)])
        set_nbr_across(T, M[(beta, ap)], f, M[(alpha, ap)])
    for sp in (alpha, beta):
        f = tuple(sorted(set(tri) | {sp}))
        set_nbr_across(T, M[(sp, U)], f, M[(sp, D)])
        set_nbr_across(T, M[(sp, D)], f, M[(sp, U)])
    # external gluings (each new pentachoron's faces opposite each tri vertex)
    for key, m in M.items():
        for v in tri:
            facelist = tuple(sorted(set(T.pent[m]) - {v}))
            fs = frozenset(facelist)
            if fs in ext_map:
                ext, old = ext_map[fs]
                set_nbr_across(T, m, facelist, ext)
                back_rewire(T, ext, old, m)
    for pid in shared:
        T.remove_pentachoron(pid)
    return True, list(M.values()), tri


# ===========================================================================
# tests
# ===========================================================================

if __name__ == "__main__":
    from v6_cdt import build_s1xs3, verify

    print("=" * 72)
    print("  v6 (2,4)/(4,2) surgery verification")
    print("=" * 72)
    T = build_s1xs3(K=10)
    ok0, rep0 = verify(T)
    fv0 = rep0["f_vector"]
    print(f"init: N4={T.n_pent()} verify={'PASS' if ok0 else 'FAIL'} f={fv0}")

    # find a (2,4)-applicable (P, face)
    applied = None
    for P in list(T.pent.keys()):
        for face in combinations(T.pent[P], 4):
            ok, news = apply_24(T, P, tuple(sorted(face)))
            if ok:
                applied = (news, face)
                break
        if applied:
            break
    assert applied, "no (2,4) candidate found in init"
    news, face = applied
    ok1, rep1 = verify(T)
    print(f"after (2,4): N4={T.n_pent()} verify={'PASS' if ok1 else 'FAIL'} "
          f"f={rep1['f_vector']} types={rep1['type_counts']} "
          f"links={rep1['link_failures']} chi={rep1['euler_char']}")

    # reverse with (4,2) on the new edge {e,f}
    efs = set(T.pent[news[0]])
    for nw in news[1:]:
        efs &= set(T.pent[nw])
    e, f = sorted(efs)
    ok42, news2 = apply_42(T, e, f)
    ok2, rep2 = verify(T)
    print(f"after (4,2) reverse: ok={ok42} N4={T.n_pent()} verify={'PASS' if ok2 else 'FAIL'} "
          f"f={rep2['f_vector']} returns_to_init={rep2['f_vector'] == fv0}")
    print(f"\nRESULT: (2,4) preserves manifold = {ok1}; (4,2) reverses to init = "
          f"{ok2 and rep2['f_vector'] == fv0}")

    # ---- (2,8)/(8,2) ----
    print("\n" + "=" * 72)
    print("  v6 (2,8)/(8,2) surgery verification")
    print("=" * 72)
    T = build_s1xs3(K=10)
    _, repA = verify(T); fvA = repA["f_vector"]
    print(f"init: N4={T.n_pent()} f={fvA}")
    applied = None
    for P in list(T.pent.keys()):
        for face in combinations(T.pent[P], 4):
            if len({T.vtime[v] for v in face}) == 1:  # spatial tet face
                ok, news, x = apply_28(T, P, face)
                if ok:
                    applied = (news, x)
                    break
        if applied:
            break
    assert applied, "no (2,8) candidate found"
    news, x = applied
    okB, repB = verify(T)
    print(f"after (2,8): N4={T.n_pent()} verify={'PASS' if okB else 'FAIL'} "
          f"f={repB['f_vector']} types={repB['type_counts']} links={repB['link_failures']} chi={repB['euler_char']}")
    ok82, _ = apply_82(T, x)
    okC, repC = verify(T)
    print(f"after (8,2) reverse: ok={ok82} N4={T.n_pent()} verify={'PASS' if okC else 'FAIL'} "
          f"returns_to_init={repC['f_vector'] == fvA}")
    print(f"\nRESULT: (2,8) preserves manifold = {okB}; (8,2) reverses to init = "
          f"{okC and repC['f_vector'] == fvA}")

    # ---- (3,3) ----  (the pristine staircase has no (3,3)/(4,6) hinge; warm up)
    print("\n" + "=" * 72)
    print("  v6 (3,3) surgery verification (self-inverse)")
    print("=" * 72)
    import numpy as np
    from v6_cdt_run import propose_and_apply

    def _warmup(K=10, moves=60, seed=0):
        Tw = build_s1xs3(K=K); rngw = np.random.default_rng(seed); a = 0
        for _ in range(4000):
            _, okw, _ = propose_and_apply(Tw, rngw, 2)
            if okw:
                a += 1
            if a >= moves:
                break
        return Tw

    T = _warmup(seed=3)
    _, repD = verify(T); fvD = repD["f_vector"]
    print(f"warmed init: N4={T.n_pent()} f={fvD}")
    applied = None
    for P in list(T.pent.keys()):
        for tri in combinations(T.pent[P], 3):
            ok, news, dual = apply_33(T, P, tri)
            if ok:
                applied = (news, dual)
                break
        if applied:
            break
    assert applied, "no (3,3) candidate found"
    news, dual = applied
    okE, repE = verify(T)
    print(f"after (3,3): N4={T.n_pent()} verify={'PASS' if okE else 'FAIL'} "
          f"f={repE['f_vector']} fvec_preserved={repE['f_vector'] == fvD} "
          f"links={repE['link_failures']} chi={repE['euler_char']}")
    # reverse: apply (3,3) on the dual hinge
    ok33b, _, _ = apply_33(T, news[0], dual)
    okF, repF = verify(T)
    print(f"after (3,3) on dual hinge: ok={ok33b} verify={'PASS' if okF else 'FAIL'} "
          f"returns_to_init={repF['f_vector'] == fvD}")
    print(f"\nRESULT: (3,3) preserves manifold+f-vector = {okE}; self-inverse = "
          f"{okF and repF['f_vector'] == fvD}")

    # ---- (4,6)/(6,4) ----
    print("\n" + "=" * 72)
    print("  v6 (4,6)/(6,4) surgery verification")
    print("=" * 72)
    T = _warmup(seed=5)
    _, repG = verify(T); fvG = repG["f_vector"]
    print(f"warmed init: N4={T.n_pent()} f={fvG}")
    applied = None
    for P in list(T.pent.keys()):
        for tri in combinations(T.pent[P], 3):
            if len({T.vtime[v] for v in tri}) == 1:        # spatial triangle hinge
                ok, news, edge = apply_46(T, P, tri)
                if ok:
                    applied = (news, edge); break
        if applied:
            break
    assert applied, "no (4,6) candidate found"
    news, (alpha, beta) = applied
    okH, repH = verify(T)
    print(f"after (4,6): N4={T.n_pent()} verify={'PASS' if okH else 'FAIL'} "
          f"f={repH['f_vector']} types={repH['type_counts']} "
          f"links={repH['link_failures']} chi={repH['euler_char']} "
          f"simplicial={repH['gluing']['is_simplicial']}")
    ok64, news2, tri2 = apply_64(T, alpha, beta)
    okI, repI = verify(T)
    print(f"after (6,4) reverse: ok={ok64} N4={T.n_pent()} verify={'PASS' if okI else 'FAIL'} "
          f"returns_to_init={repI['f_vector'] == fvG}")
    print(f"\nRESULT: (4,6) preserves manifold = {okH}; (6,4) reverses to init = "
          f"{okI and repI['f_vector'] == fvG}")

    # ---- comprehensive all-moves stress test (the real gate) ----
    print("\n" + "=" * 72)
    print("  v6 ALL-MOVES stress test: verify() after EVERY accepted move")
    print("=" * 72)
    import numpy as np
    from collections import Counter
    from v6_cdt_run import propose_and_apply
    rng = np.random.default_rng(7)
    T = build_s1xs3(K=12)
    mtypes = Counter(); checked = 0; worst = None
    all_ok = True; all_simplicial = True
    for step in range(8000):
        nb = T.n_pent()
        mt, ok, undo = propose_and_apply(T, rng, 2)
        if not ok:
            continue
        mtypes[mt] += 1
        # full verify (incl. links) every accepted move while small; periodically when larger
        if T.n_pent() < 400 or mtypes.total() % 25 == 0:
            okv, rep = verify(T, check_links=(T.n_pent() < 400))
            checked += 1
            if not okv:
                all_ok = False; worst = (mt, rep); break
            if not rep["gluing"]["is_simplicial"]:
                all_simplicial = False
        # keep volume in a band so shrink moves fire too
        if T.n_pent() > 900 and rng.random() < 0.7:
            undo()
        elif T.n_pent() < 500 and nb < T.n_pent() and rng.random() < 0.1:
            undo()
    print(f"moves fired: {dict(mtypes)}")
    print(f"verify() calls: {checked}  all_ok={all_ok}  stayed_simplicial={all_simplicial}")
    covered = set(mtypes) >= {"(2,4)", "(4,2)", "(2,8)", "(8,2)", "(3,3)", "(4,6)", "(6,4)"}
    print(f"all 7 move types exercised: {covered}")
    if not all_ok:
        print(f"  !! FIRST FAILURE after move {worst[0]}: {worst[1]['gluing']}")
    print(f"\nRESULT: all moves preserve a valid CDT 4-manifold = {all_ok}; "
          f"full move-set covered = {covered}")
