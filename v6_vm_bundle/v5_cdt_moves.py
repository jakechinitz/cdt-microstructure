# v5_cdt_moves.py -- the 4D CDT Pachner move set (Stage 1.B).
#
# Implements all 10 4D CDT moves (paired as forward/reverse), with the
# Hastings (asymmetric proposal) factors needed for proper detailed balance.
#
# Move pairs and what they do (AJL Nucl Phys B 610 (2001) Sec 7-8):
#   (3,3)             : exchange 3 simplices sharing a time-like triangle
#                       for 3 sharing the dual triangle. Self-inverse, f-vector
#                       preserving. (already in v5_cdt_lib.py)
#   (2,4) / (4,2)     : 2 simplices sharing a tet  <->  4 simplices sharing
#                       the dual edge. Df-vec = (0, +1, +4, +5, +2).
#   (4,6) / (6,4)     : 4 simplices sharing a spatial triangle  <->  6
#                       simplices sharing a spatial edge. Spatial 2-3 lift.
#   (2,8) / (8,2)     : 2 simplices sharing a spatial tetrahedron  <->  8
#                       simplices sharing a new spatial vertex at the
#                       tetrahedron's center. Spatial barycentric subdivision.
#
# Verification status:
#   - (2,4)/(4,2): geometry derived from first principles + tested via the
#                  invariant "forward then reverse is identity" + f-vector
#                  consistency check at every move. Hastings ratio = N_3 / N_1.
#   - (4,6)/(6,4): IMPLEMENTED but the move-pattern matching is the most
#                  fragile part; flagged for cross-check against CDT-plusplus
#                  / AJL §7.3.
#   - (2,8)/(8,2): IMPLEMENTED. The (2,8) is a 4D analog of barycentric
#                  subdivision of a spatial tetrahedron. Hastings includes
#                  the choice of which spatial tetrahedron to subdivide and
#                  which to coalesce.
#
# What this module does NOT do (deferred to Stage 1.C, 1.D, 2):
#   - Full simplicial-manifold link-condition check. The collision check
#     ("new simplex already exists" / "would-add edge already exists") is
#     a NECESSARY condition for link-S^3 preservation, not sufficient.
#     Full verification requires walking the local 1-skeleton and confirming
#     the link is combinatorially a 3-sphere (Stage 1.C).
#   - Volume fixing penalty (added in Stage 1.D).
#   - EPRL action (Stage 2).

from __future__ import annotations
from itertools import combinations
from collections import Counter
from typing import Optional
import numpy as np

from v5_cdt_lib import (
    CausalTriangulation,
    FourSimplex, Tetrahedron, Triangle,
    faces_of_4simplex, triangles_of_4simplex, edges_of_4simplex,
    find_three_three_partners, apply_three_three,
)


# ===========================================================================
# (2,4) MOVE
# ===========================================================================

def find_two_four_candidate(T: CausalTriangulation, tet: Tetrahedron):
    """For tetrahedron tet = {a,b,c,d} shared by exactly 2 4-simplices
    σ_1 = {a,b,c,d,e}, σ_2 = {a,b,c,d,f}, return (e, f, sigmas).
    Returns None if (a) tet not shared by exactly 2, or (b) the link
    condition is violated (edge {e,f} already exists)."""
    sset = T.tet_to_simplices.get(tet)
    if sset is None or len(sset) != 2:
        return None
    s1, s2 = tuple(sset)
    e_set = set(s1) - set(tet)
    f_set = set(s2) - set(tet)
    if len(e_set) != 1 or len(f_set) != 1:
        return None
    e = next(iter(e_set))
    f = next(iter(f_set))
    if e == f:
        return None
    # Link condition: {e,f} must not be an existing edge
    if T.has_edge(e, f):
        return None
    return (e, f, (s1, s2))


def cdt_types_valid(T: CausalTriangulation, simplices) -> bool:
    """Verify each simplex has a valid CDT type ((4,1), (3,2), (2,3), (1,4))."""
    for s in simplices:
        a, b = T.simplex_type(s)
        if {a, b} != {4, 1} and {a, b} != {3, 2}:
            return False
    return True


def apply_two_four(T: CausalTriangulation, tet: Tetrahedron, enforce_cdt: bool = True):
    """Apply (2,4): remove the 2 simplices sharing `tet`, add 4 sharing the
    dual edge. Returns (success, removed, added) tuple."""
    cand = find_two_four_candidate(T, tet)
    if cand is None:
        return False, [], []
    e, f, (s1, s2) = cand
    tet_v = list(tet)
    # New simplices: each triangle of tet plus the new edge {e, f}
    new_simplices = []
    for tri in combinations(tet_v, 3):
        new_s = tuple(sorted(tri + (e, f)))
        if new_s in T.simplices:
            return False, [], []  # collision
        new_simplices.append(new_s)
    # CDT type validity check on the new simplices
    if enforce_cdt:
        # build a hypothetical T with the new simplices for type checking
        if not all(
            _cdt_type_ok(T, s) for s in new_simplices
        ):
            return False, [], []
    # Apply
    T.remove_4simplex(s1)
    T.remove_4simplex(s2)
    for s in new_simplices:
        T.add_4simplex(s)
    return True, [s1, s2], new_simplices


def _cdt_type_ok(T: CausalTriangulation, s):
    """Check a (proposed) simplex would have a valid CDT type."""
    a, b = T.simplex_type(s)
    return {a, b} == {4, 1} or {a, b} == {3, 2}


# ===========================================================================
# (4,2) MOVE -- reverse of (2,4)
# ===========================================================================

def find_four_two_candidate(T: CausalTriangulation, edge):
    """For edge {e,f} shared by exactly 4 simplices forming the (4,2)
    pattern {e,f,xyz} for 4 choices of xyz from a 4-vertex set, return
    (e, f, abcd, sigmas). Returns None if not the right pattern.

    Recognition: edges shared by exactly 4 simplices whose "outside-edge"
    vertices form the C(4,3)=4 triangles of a single 4-vertex set, AND
    that set is not already a tetrahedron in T."""
    edge = (min(edge), max(edge))
    sset = T.edge_to_simplices.get(edge)
    if sset is None or len(sset) != 4:
        return None
    e, f = edge
    # For each simplex, its "outside" vertices (not in {e,f}) form a triangle
    outside_tris = []
    for s in sset:
        out = tuple(sorted(set(s) - {e, f}))
        if len(out) != 3:
            return None
        outside_tris.append(out)
    # The 4 outside-triangles should be exactly the 4 triangles of a 4-vertex
    # set {a,b,c,d}.
    abcd = sorted(set().union(*outside_tris))
    if len(abcd) != 4:
        return None
    expected = set(tuple(sorted(t)) for t in combinations(abcd, 3))
    if set(outside_tris) != expected:
        return None
    abcd = tuple(abcd)
    # Link condition: tetrahedron {a,b,c,d} must NOT exist
    if T.has_tetrahedron(abcd):
        return None
    return (e, f, abcd, tuple(sset))


def apply_four_two(T: CausalTriangulation, edge, enforce_cdt: bool = True):
    """Apply (4,2): remove 4 simplices sharing `edge`, add 2 sharing the
    dual tetrahedron. Returns (success, removed, added)."""
    cand = find_four_two_candidate(T, edge)
    if cand is None:
        return False, [], []
    e, f, abcd, removed = cand
    new1 = tuple(sorted(abcd + (e,)))
    new2 = tuple(sorted(abcd + (f,)))
    if new1 in T.simplices or new2 in T.simplices:
        return False, [], []
    if enforce_cdt:
        if not (_cdt_type_ok(T, new1) and _cdt_type_ok(T, new2)):
            return False, [], []
    for s in removed:
        T.remove_4simplex(s)
    T.add_4simplex(new1)
    T.add_4simplex(new2)
    return True, list(removed), [new1, new2]


# ===========================================================================
# (4,6) MOVE -- spatial triangle <-> spatial edge swap
# ===========================================================================
# A spatial triangle T = {p,q,r} (all 3 vertices in same time slice) shared
# by 4 four-simplices: 2 going "up" to t+1, 2 going "down" to t-1.
# After (4,6): T is replaced by 6 simplices sharing a spatial edge in some
# slice. The exact pattern is fragile and merits cross-check against AJL.
#
# Implementation strategy: identify the 4 sharing simplices, classify them
# into "up" pair and "down" pair, then construct the 6 dual simplices.

def find_four_six_candidate(T: CausalTriangulation, tri: Triangle):
    """For SPATIAL triangle tri (all 3 vertices same time slice) shared by
    exactly 4 simplices in the (4,6) "clean stacking" pattern: 2 sigmas
    with apexes in slice A and 2 sigmas with apexes in slice B, where A != B
    (and both != t). Within each group of 2, the apexes must be the SAME
    vertex, and the 2 spatial neighbors (alpha, beta) must match across the
    two groups. Link condition: spatial edge {alpha, beta} must NOT exist.

    Periodic-time-safe: A and B can be any two slices, not necessarily t+1
    and t-1 in integer arithmetic.

    Returns (tri, alpha, beta, U, D, up_sims, dn_sims) or None.
    """
    sset = T.tri_to_simplices.get(tri)
    if sset is None or len(sset) != 4:
        return None
    times = {T.vertex_time[v] for v in tri}
    if len(times) != 1:
        return None
    t = next(iter(times))
    # Decompose each sigma into (spatial_neighbor, apex), group by apex-slice
    by_slice = {}
    for s in sset:
        extras = [v for v in s if v not in tri]
        if len(extras) != 2:
            return None
        spatial = [v for v in extras if T.vertex_time[v] == t]
        timed = [v for v in extras if T.vertex_time[v] != t]
        if len(spatial) != 1 or len(timed) != 1:
            return None
        sp, ap = spatial[0], timed[0]
        ap_t = T.vertex_time[ap]
        by_slice.setdefault(ap_t, []).append((s, sp, ap))
    # Need exactly 2 distinct apex-slices, each with 2 sigmas
    if len(by_slice) != 2 or not all(len(v) == 2 for v in by_slice.values()):
        return None
    apex_slices = sorted(by_slice.keys())
    # Arbitrary labels: lower slice = "down", upper = "up"
    above_grp = by_slice[apex_slices[1]]
    below_grp = by_slice[apex_slices[0]]
    # Within each group, the 2 apexes must be the same vertex
    if above_grp[0][2] != above_grp[1][2]:
        return None
    U = above_grp[0][2]
    if below_grp[0][2] != below_grp[1][2]:
        return None
    D = below_grp[0][2]
    # The 2 spatial neighbors alpha, beta must be the same in both groups
    alpha_a, beta_a = above_grp[0][1], above_grp[1][1]
    if alpha_a == beta_a:
        return None
    if {alpha_a, beta_a} != {below_grp[0][1], below_grp[1][1]}:
        return None
    alpha, beta = (min(alpha_a, beta_a), max(alpha_a, beta_a))
    # Link condition: spatial edge {alpha, beta} must NOT already exist
    if T.has_edge(alpha, beta):
        return None
    return (tri, alpha, beta, U, D,
            tuple(s for s, _, _ in above_grp),
            tuple(s for s, _, _ in below_grp))


def apply_four_six(T: CausalTriangulation, tri: Triangle, enforce_cdt: bool = True):
    """Apply (4,6): replace 4 simplices sharing spatial triangle tri with
    6 simplices sharing the new spatial edge {alpha, beta}. The 6 are:
    for each of the 3 edges of tri, a "wedge" simplex = tri_edge + {alpha,
    beta, U} (up) and tri_edge + {alpha, beta, D} (down).
    """
    cand = find_four_six_candidate(T, tri)
    if cand is None:
        return False, [], []
    tri, alpha, beta, U, D, ups, dns = cand
    tri_v = list(tri)
    new_simplices = []
    for tri_edge in combinations(tri_v, 2):
        new_up = tuple(sorted(tri_edge + (alpha, beta, U)))
        new_dn = tuple(sorted(tri_edge + (alpha, beta, D)))
        new_simplices.extend([new_up, new_dn])
    for s in new_simplices:
        if s in T.simplices:
            return False, [], []
    if enforce_cdt:
        for s in new_simplices:
            if not _cdt_type_ok(T, s):
                return False, [], []
    removed = list(ups) + list(dns)
    for s in removed:
        T.remove_4simplex(s)
    for s in new_simplices:
        T.add_4simplex(s)
    return True, removed, new_simplices


def find_six_four_candidate(T: CausalTriangulation, edge):
    """For spatial edge {alpha, beta} (both endpoints same slice t) shared
    by exactly 6 simplices in the (6,4)-reverse pattern (3 up + 3 down,
    each triple sharing the same apex above/below). Returns
    (tri, alpha, beta, U, D, up_sims, dn_sims) where tri is the spatial
    triangle that the move will create."""
    edge = (min(edge), max(edge))
    sset = T.edge_to_simplices.get(edge)
    if sset is None or len(sset) != 6:
        return None
    a, b = edge
    t_a = T.vertex_time.get(a); t_b = T.vertex_time.get(b)
    if t_a != t_b:
        return None  # edge not spatial
    t = t_a
    by_slice = {}
    for s in sset:
        others = [v for v in s if v != a and v != b]
        if len(others) != 3:
            return None
        spatial_others = [v for v in others if T.vertex_time[v] == t]
        timed_others = [v for v in others if T.vertex_time[v] != t]
        if len(spatial_others) != 2 or len(timed_others) != 1:
            return None
        apex = timed_others[0]
        spatial_pair = tuple(sorted(spatial_others))
        ap_t = T.vertex_time[apex]
        by_slice.setdefault(ap_t, []).append((s, spatial_pair, apex))
    # Need exactly 2 distinct apex-slices, each with 3 sigmas
    if len(by_slice) != 2 or not all(len(v) == 3 for v in by_slice.values()):
        return None
    apex_slices = sorted(by_slice.keys())
    above_grp = by_slice[apex_slices[1]]
    below_grp = by_slice[apex_slices[0]]
    # All 3 in each group must share the same apex
    U_set = {x[2] for x in above_grp}
    if len(U_set) != 1: return None
    U = U_set.pop()
    D_set = {x[2] for x in below_grp}
    if len(D_set) != 1: return None
    D = D_set.pop()
    # The 3 spatial-pairs in above_grp should be the 3 edges of a triangle
    up_pairs = {x[1] for x in above_grp}
    tri_vertices = set().union(*up_pairs)
    if len(tri_vertices) != 3:
        return None
    expected_pairs = {tuple(sorted(p)) for p in combinations(tri_vertices, 2)}
    if up_pairs != expected_pairs:
        return None
    dn_pairs = {x[1] for x in below_grp}
    if dn_pairs != expected_pairs:
        return None
    tri = tuple(sorted(tri_vertices))
    # Link condition: the new spatial triangle tri must NOT already exist
    if T.has_triangle(tri):
        return None
    return (tri, a, b, U, D,
            tuple(s for s, _, _ in above_grp),
            tuple(s for s, _, _ in below_grp))


def apply_six_four(T: CausalTriangulation, edge, enforce_cdt: bool = True):
    """Apply (6,4): replace 6 simplices sharing spatial edge with 4 simplices
    sharing the spatially dual triangle tri."""
    cand = find_six_four_candidate(T, edge)
    if cand is None:
        return False, [], []
    tri, a, b, U, D, ups, dns = cand
    new_simplices = []
    # New (up) simplices: tri + {a, U}, tri + {b, U}
    new_simplices.append(tuple(sorted(tri + (a, U))))
    new_simplices.append(tuple(sorted(tri + (b, U))))
    new_simplices.append(tuple(sorted(tri + (a, D))))
    new_simplices.append(tuple(sorted(tri + (b, D))))
    for s in new_simplices:
        if s in T.simplices:
            return False, [], []
    if enforce_cdt:
        for s in new_simplices:
            if not _cdt_type_ok(T, s):
                return False, [], []
    removed = list(ups) + list(dns)
    for s in removed:
        T.remove_4simplex(s)
    for s in new_simplices:
        T.add_4simplex(s)
    return True, removed, new_simplices


# ===========================================================================
# (2,8) MOVE -- insert spatial vertex
# ===========================================================================
# Insert a new vertex v_new at the "center" of a spatial tetrahedron shared
# by a (1,4)+(4,1) pair (sandwich: one simplex going up, one going down,
# sharing a spatial tetrahedron).

def find_two_eight_candidate(T: CausalTriangulation, spatial_tet: Tetrahedron):
    """For a spatial tetrahedron (all 4 vertices same slice) shared by
    exactly 2 simplices σ_up (extending one apex in one neighbor slice) and
    σ_dn (apex in a DIFFERENT neighbor slice), return (spatial_tet,
    apex_up, apex_dn, σ_up, σ_dn).

    Loose check (periodic-time-safe): apex times must be different from each
    other and different from t. The naming of which apex is "up" vs "down"
    is by integer-time order (arbitrary, but consistent)."""
    sset = T.tet_to_simplices.get(spatial_tet)
    if sset is None or len(sset) != 2:
        return None
    times = {T.vertex_time[v] for v in spatial_tet}
    if len(times) != 1:
        return None
    t = next(iter(times))
    s1, s2 = tuple(sset)
    apex1 = next(iter(set(s1) - set(spatial_tet)))
    apex2 = next(iter(set(s2) - set(spatial_tet)))
    t1 = T.vertex_time[apex1]
    t2 = T.vertex_time[apex2]
    # Loose check: apexes in DIFFERENT slices, both different from t
    if t1 == t2 or t1 == t or t2 == t:
        return None
    # Arbitrary up/down assignment: smaller integer time is "down"
    if t1 < t2:
        apex_dn, apex_up = apex1, apex2
        s_dn, s_up = s1, s2
    else:
        apex_dn, apex_up = apex2, apex1
        s_dn, s_up = s2, s1
    return (spatial_tet, apex_up, apex_dn, s_up, s_dn)


def apply_two_eight(T: CausalTriangulation, spatial_tet: Tetrahedron, enforce_cdt: bool = True):
    """Insert vertex v_new at center of spatial_tet, splitting the (1,4)+(4,1)
    pair into 8 simplices (4 going up + 4 going down).
    Df-vec = (+1, +4, +6, +4, +6) approximately; exact values verified by
    invariant check at apply time."""
    cand = find_two_eight_candidate(T, spatial_tet)
    if cand is None:
        return False, [], []
    tet_v, apex_up, apex_dn, s_up, s_dn = cand
    # New vertex
    v_new = T.next_vertex_id()
    t = T.vertex_time[tet_v[0]]
    # The 8 new simplices: each takes 3 of the 4 spatial-tet vertices plus
    # v_new plus the apex (up or down)
    new_simplices = []
    for v_drop in tet_v:
        face = tuple(sorted(set(tet_v) - {v_drop}))  # 3 vertices of a face of spatial_tet
        new_up = tuple(sorted(face + (v_new, apex_up)))
        new_dn = tuple(sorted(face + (v_new, apex_dn)))
        new_simplices.extend([new_up, new_dn])
    # Collision checks
    for s in new_simplices:
        if s in T.simplices:
            return False, [], []
    # Add new vertex first (so type checks work)
    T.vertex_time[v_new] = t
    if enforce_cdt:
        for s in new_simplices:
            if not _cdt_type_ok(T, s):
                del T.vertex_time[v_new]
                return False, [], []
    # Apply
    T.remove_4simplex(s_up)
    T.remove_4simplex(s_dn)
    for s in new_simplices:
        T.add_4simplex(s)
    return True, [s_up, s_dn], new_simplices


# ===========================================================================
# (8,2) MOVE -- delete a spatial vertex of degree 4 in its slice
# ===========================================================================

def find_eight_two_candidate(T: CausalTriangulation, v: int):
    """For a vertex v (in some slice t) of "local degree 8" in the (8,2)
    sense: v is the apex of 8 simplices forming the (2,8) pattern dual."""
    if v not in T.vertex_time:
        return None
    t = T.vertex_time[v]
    # Find all simplices containing v (O(1) via maintained incidence map)
    sset = list(T.vertex_to_simplices.get(v, ()))
    if len(sset) != 8:
        return None
    # Each simplex has 4 other vertices. The 8 simplices should pair into
    # 4 up + 4 down (each pair sharing 3 of the 4 spatial-tet vertices around v).
    by_slice = {}
    for s in sset:
        others = set(s) - {v}
        same_slice = [u for u in others if T.vertex_time[u] == t]
        other_slice = [u for u in others if T.vertex_time[u] != t]
        if len(same_slice) != 3 or len(other_slice) != 1:
            return None
        u_apex = other_slice[0]
        ap_t = T.vertex_time[u_apex]
        by_slice.setdefault(ap_t, []).append((s, tuple(sorted(same_slice)), u_apex))
    # Need exactly 2 distinct apex-slices, each with 4 sigmas
    if len(by_slice) != 2 or not all(len(x) == 4 for x in by_slice.values()):
        return None
    apex_slices = sorted(by_slice.keys())
    up_sim = by_slice[apex_slices[1]]
    dn_sim = by_slice[apex_slices[0]]
    up_apexes = {x[2] for x in up_sim}
    dn_apexes = {x[2] for x in dn_sim}
    if len(up_apexes) != 1 or len(dn_apexes) != 1:
        return None
    apex_up = up_apexes.pop()
    apex_dn = dn_apexes.pop()
    # The 4 "spatial faces" around v in up_sim should be the 4 triangles
    # of some spatial tetrahedron {a,b,c,d}
    up_faces = set(x[1] for x in up_sim)
    abcd = sorted(set().union(*up_faces))
    if len(abcd) != 4:
        return None
    expected = {tuple(sorted(t)) for t in combinations(abcd, 3)}
    if up_faces != expected:
        return None
    spatial_tet = tuple(abcd)
    # Link condition: after deletion, spatial_tet must not already exist
    if T.has_tetrahedron(spatial_tet):
        return None
    return (v, spatial_tet, apex_up, apex_dn, [x[0] for x in up_sim], [x[0] for x in dn_sim])


def apply_eight_two(T: CausalTriangulation, v: int, enforce_cdt: bool = True):
    """Apply (8,2): delete vertex v, coalesce 8 simplices into 2."""
    cand = find_eight_two_candidate(T, v)
    if cand is None:
        return False, [], []
    v, spatial_tet, apex_up, apex_dn, up_sim, dn_sim = cand
    new_up = tuple(sorted(spatial_tet + (apex_up,)))
    new_dn = tuple(sorted(spatial_tet + (apex_dn,)))
    if new_up in T.simplices or new_dn in T.simplices:
        return False, [], []
    if enforce_cdt:
        if not (_cdt_type_ok(T, new_up) and _cdt_type_ok(T, new_dn)):
            return False, [], []
    for s in up_sim + dn_sim:
        T.remove_4simplex(s)
    del T.vertex_time[v]
    T.add_4simplex(new_up)
    T.add_4simplex(new_dn)
    return True, list(up_sim) + list(dn_sim), [new_up, new_dn]


# ===========================================================================
# Move dispatcher with proper Hastings ratios
# ===========================================================================

MOVE_NAMES = ["(3,3)", "(2,4)", "(4,2)", "(4,6)", "(6,4)", "(2,8)", "(8,2)"]


def propose_move(T: CausalTriangulation, rng):
    """Pick a random move type uniformly, propose a move of that type by
    picking a target sub-simplex uniformly. Returns (move_type, success,
    removed, added, hastings_log_ratio).

    The Hastings ratio q(x'->x)/q(x->x') is computed as
        (1/n_targets_reverse) / (1/n_targets_forward)
        = n_targets_forward / n_targets_reverse
    where n_targets_forward = number of valid target sub-simplices for the
    forward move in state x, and n_targets_reverse = same for the reverse
    move in state x' (post-move state).

    Conservative: includes only the sub-simplex selection step in the
    Hastings ratio. The (1/n_move_types) selection cancels because forward
    and reverse use the same n_move_types.
    """
    move_type = MOVE_NAMES[rng.integers(0, len(MOVE_NAMES))]

    # Uniform sampling + candidate counts now come from O(1) bags maintained in
    # CausalTriangulation.add/remove_4simplex (was O(N) list builds + filters).
    # Proposal distribution and Hastings counts are identical to the previous
    # list-based version; only the per-attempt cost changed (O(N) -> O(1)).
    # n_fwd is read BEFORE the move, n_rev (reverse-candidate count) AFTER it.

    if move_type == "(3,3)":
        if len(T.tl3_bag) == 0:
            return move_type, False, [], [], 0.0
        n_fwd = len(T.tl3_bag)
        tri = T.tl3_bag.pick(rng)
        ok, removed, added = apply_three_three(T, tri)
        if not ok:
            return move_type, False, [], [], 0.0
        n_rev = len(T.tl3_bag)
        hastings = np.log(max(n_fwd, 1)) - np.log(max(n_rev, 1))
        return move_type, True, removed, added, hastings

    if move_type == "(2,4)":
        if len(T.tet_bag) == 0:
            return move_type, False, [], [], 0.0
        n_fwd_total = len(T.tet_bag)
        tet = T.tet_bag.pick(rng)
        ok, removed, added = apply_two_four(T, tet)
        if not ok:
            return move_type, False, [], [], 0.0
        n_rev_total = len(T.edge_bag)
        hastings = np.log(n_fwd_total) - np.log(max(n_rev_total, 1))
        return move_type, True, removed, added, hastings

    if move_type == "(4,2)":
        if len(T.edge_bag) == 0:
            return move_type, False, [], [], 0.0
        n_fwd_total = len(T.edge_bag)
        edge = T.edge_bag.pick(rng)
        ok, removed, added = apply_four_two(T, edge)
        if not ok:
            return move_type, False, [], [], 0.0
        n_rev_total = len(T.tet_bag)
        hastings = np.log(n_fwd_total) - np.log(max(n_rev_total, 1))
        return move_type, True, removed, added, hastings

    if move_type == "(4,6)":
        if len(T.spatial_tri_bag) == 0:
            return move_type, False, [], [], 0.0
        n_fwd_total = len(T.spatial_tri_bag)
        tri = T.spatial_tri_bag.pick(rng)
        ok, removed, added = apply_four_six(T, tri)
        if not ok:
            return move_type, False, [], [], 0.0
        n_rev_total = len(T.spatial_edge_bag)
        hastings = np.log(n_fwd_total) - np.log(max(n_rev_total, 1))
        return move_type, True, removed, added, hastings

    if move_type == "(6,4)":
        if len(T.spatial_edge_bag) == 0:
            return move_type, False, [], [], 0.0
        n_fwd_total = len(T.spatial_edge_bag)
        edge = T.spatial_edge_bag.pick(rng)
        ok, removed, added = apply_six_four(T, edge)
        if not ok:
            return move_type, False, [], [], 0.0
        n_rev_total = len(T.spatial_tri_bag)
        hastings = np.log(n_fwd_total) - np.log(max(n_rev_total, 1))
        return move_type, True, removed, added, hastings

    if move_type == "(2,8)":
        if len(T.spatial_tet_bag) == 0:
            return move_type, False, [], [], 0.0
        n_fwd_total = len(T.spatial_tet_bag)
        tet = T.spatial_tet_bag.pick(rng)
        ok, removed, added = apply_two_eight(T, tet)
        if not ok:
            return move_type, False, [], [], 0.0
        # Reverse (8,2) picks any vertex (most rejected by the structural check);
        # the proposal denominator is the total vertex count, post-move.
        n_rev_total = len(T.vertex_bag)
        hastings = np.log(n_fwd_total) - np.log(max(n_rev_total, 1))
        return move_type, True, removed, added, hastings

    if move_type == "(8,2)":
        if len(T.vertex_bag) == 0:
            return move_type, False, [], [], 0.0
        n_fwd_total = len(T.vertex_bag)
        v = T.vertex_bag.pick(rng)
        ok, removed, added = apply_eight_two(T, v)
        if not ok:
            return move_type, False, [], [], 0.0
        n_rev_total = len(T.spatial_tet_bag)
        hastings = np.log(n_fwd_total) - np.log(max(n_rev_total, 1))
        return move_type, True, removed, added, hastings

    return move_type, False, [], [], 0.0


# ===========================================================================
# Move reversal (for rejected proposals)
# ===========================================================================

def undo_move(T: CausalTriangulation, move_type: str, removed, added,
              vertex_change: Optional[tuple] = None):
    """Reverse a successfully-applied move so we can revert on rejection.
    Handles vertex addition/deletion (for (2,8)/(8,2)) by tracking changes
    in vertex_time.

    vertex_change is one of:
        None                      -- no vertex change
        ("added", v_id, time)     -- a (2,8) added vertex v_id; undo by deleting it
        ("removed", v_id, time)   -- an (8,2) removed vertex v_id; undo by restoring
    """
    # remove the added simplices first (must happen before vertex deletion
    # because remove_4simplex needs the vertex to exist in the lookup tables)
    for s in added:
        T.remove_4simplex(s)
    # Handle vertex change
    if vertex_change is not None:
        kind, v_id, t = vertex_change
        if kind == "added":
            # The move added v_id; undo by deleting (the added simplices are
            # already removed above, so the vertex is now orphaned)
            if v_id in T.vertex_time:
                del T.vertex_time[v_id]
        elif kind == "removed":
            # The move removed v_id; restore it before re-adding old simplices
            T.vertex_time[v_id] = t
    # Restore removed simplices
    for s in removed:
        T.add_4simplex(s)


def propose_move_tracked(T: CausalTriangulation, rng):
    """Wrapper around propose_move that also returns vertex_change info so the
    Metropolis loop can call undo_move correctly on rejection.

    Only (2,8) adds a vertex and (8,2) removes one, so the changed vertex is
    found in O(1) as the set-difference of the added vs removed simplex
    vertices -- no O(N) snapshot of vertex_time. For (8,2) the deleted vertex's
    time is reconstructed as the modal slice of the other vertices of a removed
    simplex (3 spatial neighbors at that slice + 1 apex), which equals it."""
    move_type, ok, removed, added, hastings = propose_move(T, rng)
    if not ok:
        return move_type, False, [], [], 0.0, None
    vertex_change = None
    if move_type == "(2,8)":
        v = next(iter(set().union(*added) - set().union(*removed)))
        vertex_change = ("added", v, T.vertex_time[v])
    elif move_type == "(8,2)":
        v = next(iter(set().union(*removed) - set().union(*added)))
        times = [T.vertex_time[x] for x in removed[0] if x != v]
        t = Counter(times).most_common(1)[0][0]
        vertex_change = ("removed", v, t)
    return move_type, True, removed, added, hastings, vertex_change
