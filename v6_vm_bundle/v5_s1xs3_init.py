# v5_s1xs3_init.py -- proper S^1 x S^3 stacked initialization for CDT v5.
#
# Per spec:
#   - K >= 3 time slices, periodic identification (t+K ≡ t)
#   - Each spatial slice = boundary-of-4-simplex (∂Δ^4): 5 vertices, 5 spatial
#     tetrahedra forming a triangulation of S^3
#   - For each spatial tet (a,b,c,d) in slice t (ordered a<b<c<d), the prism
#     to slice t+1 is decomposed via STAIRCASE into 4 4-simplices:
#         (4,1): [a_t, b_t, c_t, d_t,     d_{t+1}]
#         (3,2): [a_t, b_t, c_t, c_{t+1}, d_{t+1}]
#         (2,3): [a_t, b_t, b_{t+1}, c_{t+1}, d_{t+1}]
#         (1,4): [a_t, a_{t+1}, b_{t+1}, c_{t+1}, d_{t+1}]
#     Total per init: 5 tets x 4 stair pieces x K prisms = 20K 4-simplices
#
# Why the staircase works:
#   For two spatial tets T = {a,b,c,d} and T' = {a,b,c,e} that share the
#   triangle {a,b,c} (so e is the "fifth" vertex of ∂Δ^4):
#     - In T's prism, the (3,2) is {a_t, b_t, c_t, c_{t+1}, d_{t+1}}.
#       Its face dropping d_{t+1} is the tetrahedron {a_t, b_t, c_t, c_{t+1}}.
#     - In T's prism, the (3,2) is {a_t, b_t, c_t, c_{t+1}, e_{t+1}}.
#       Its face dropping e_{t+1} is the SAME tetrahedron {a_t, b_t, c_t, c_{t+1}}.
#   So both prisms' (3,2) simplices share this tet -> manifold ok across
#   the shared triangle. The same holds for (4,1)/(1,4)/(2,3) walls.
#
# Result: each spatial tetrahedron T in slice t is shared by EXACTLY 2
# 4-simplices: the (4,1) of prism(T, t->t+1) and the (1,4) of prism(T,
# t-1->t). That is a real (4,1)+(1,4) SANDWICH -> the (2,8) Pachner move
# finally has candidates.
#
# Note on periodic time and simplex_type:
#   At the wraparound prism (slice K-1 -> slice 0), simplex_type sees times
#   {0, K-1}. By integer order, "lower" is 0 and "upper" is K-1, so a (4,1)
#   staircase piece (4 in slice K-1, 1 in slice 0) is CLASSIFIED as (1,4) by
#   the simplex_type function -- the type LABEL flips at the wraparound but
#   the SET-based CDT-validity check ({a,b} in {{4,1},{3,2}}) still passes.
#   The Pachner move finders have been loosened to accept any 2 apex slices
#   that are both DIFFERENT from t and from each other (the natural
#   CDT-validity criterion under periodicity).

from __future__ import annotations
from itertools import combinations

from v5_cdt_lib import CausalTriangulation


def build_s1xs3_staircase(K: int = 3, verbose: bool = False) -> CausalTriangulation:
    """Build a proper S^1 x S^3 stacked initial triangulation.

    Args:
        K: number of time slices (K >= 3 for periodic adjacency to give
           distinct (t+1) mod K and (t-1) mod K)

    Returns:
        CausalTriangulation with 5K vertices, 20K 4-simplices, K_periodic=K.
    """
    if K < 3:
        raise ValueError(f"Need K >= 3 for proper periodic time; got K={K}")
    T = CausalTriangulation()
    T.K_periodic = K

    # Vertex labeling: global ID = 5*t + i for slice t, intra-slice index i
    def v(t, i):
        return 5 * (t % K) + i

    # Assign time labels
    for t in range(K):
        for i in range(5):
            T.vertex_time[v(t, i)] = t

    # 5 spatial tetrahedra per slice = 4-subsets of {0..4}
    spatial_tets = list(combinations(range(5), 4))

    # For each pair of adjacent slices, build the 4-simplex prism
    for t in range(K):
        for tet in spatial_tets:
            a, b, c, d = tet  # ordered a<b<c<d
            # The 4 staircase pieces (vertex tuples auto-sorted by add_4simplex)
            T.add_4simplex((v(t, a), v(t, b), v(t, c), v(t, d), v(t+1, d)))           # (4,1)
            T.add_4simplex((v(t, a), v(t, b), v(t, c), v(t+1, c), v(t+1, d)))         # (3,2)
            T.add_4simplex((v(t, a), v(t, b), v(t+1, b), v(t+1, c), v(t+1, d)))       # (2,3)
            T.add_4simplex((v(t, a), v(t+1, a), v(t+1, b), v(t+1, c), v(t+1, d)))     # (1,4)

    if verbose:
        print(f"  K = {K} slices, vertices = {5*K}, simplices = {20*K}")
        print(f"  f-vector: {T.f_vector()}")
        N_41, N_32, N_other = T.type_counts()
        print(f"  type counts (by set): N_(4,1)={N_41}, N_(3,2)={N_32}, "
              f"N_other(non-CDT)={N_other}")

    return T


# ===========================================================================
# Validation
# ===========================================================================

def check_manifold_explicit(T: CausalTriangulation):
    """Every tetrahedral face must be incident to exactly 2 4-simplices
    (closed-manifold condition for the 3-face/4-cell pairing)."""
    bad = []
    for tet, sset in T.tet_to_simplices.items():
        if len(sset) != 2:
            bad.append((tet, len(sset)))
    return bad


def check_cdt_types(T: CausalTriangulation):
    """Every 4-simplex must have CDT-valid type ({4,1} or {3,2} as a set)."""
    bad = []
    for s in T.simplices:
        a, b = T.simplex_type(s)
        if {a, b} not in [{4, 1}, {3, 2}]:
            bad.append((s, (a, b)))
    return bad


def count_move_candidates(T: CausalTriangulation):
    """Count link-legal candidates for each Pachner move type."""
    from v5_cdt_moves import (
        find_three_three_partners,
        find_two_four_candidate,
        find_four_two_candidate,
        find_four_six_candidate,
        find_six_four_candidate,
        find_two_eight_candidate,
        find_eight_two_candidate,
    )

    counts = {"(3,3)": 0, "(2,4)": 0, "(4,2)": 0,
              "(4,6)": 0, "(6,4)": 0, "(2,8)": 0, "(8,2)": 0}

    # (3,3) over time-like triangles with exactly 3 incident simplices
    for tri in T.timelike_triangles_with_three_simplices():
        if find_three_three_partners(T, tri) is not None:
            counts["(3,3)"] += 1

    # (2,4) over all tetrahedra
    for tet in T.tet_to_simplices:
        if find_two_four_candidate(T, tet) is not None:
            counts["(2,4)"] += 1

    # (4,2) over all edges
    for edge in T.edge_to_simplices:
        if find_four_two_candidate(T, edge) is not None:
            counts["(4,2)"] += 1

    # (4,6) over spatial triangles
    for tri in T.tri_to_simplices:
        times = {T.vertex_time[u] for u in tri}
        if len(times) == 1:
            if find_four_six_candidate(T, tri) is not None:
                counts["(4,6)"] += 1

    # (6,4) over spatial edges
    for edge in T.edge_to_simplices:
        if T.vertex_time[edge[0]] == T.vertex_time[edge[1]]:
            if find_six_four_candidate(T, edge) is not None:
                counts["(6,4)"] += 1

    # (2,8) over spatial tetrahedra
    for tet in T.tet_to_simplices:
        times = {T.vertex_time[u] for u in tet}
        if len(times) == 1:
            if find_two_eight_candidate(T, tet) is not None:
                counts["(2,8)"] += 1

    # (8,2) over vertices
    for v_id in T.vertex_time:
        if find_eight_two_candidate(T, v_id) is not None:
            counts["(8,2)"] += 1

    return counts


if __name__ == "__main__":
    print("=" * 76)
    print("  v5 S^1 x S^3 staircase init -- validation")
    print("=" * 76)
    print()

    from v5_cdt_link import verify_manifold

    for K in [3, 4, 5]:
        print(f"--- K = {K} ---")
        T = build_s1xs3_staircase(K=K, verbose=True)

        bad_mfd = check_manifold_explicit(T)
        bad_types = check_cdt_types(T)
        ok_link, link_failures = verify_manifold(T)
        counts = count_move_candidates(T)

        print(f"  every tet shared by 2 simplices: "
              f"{'PASS' if not bad_mfd else 'FAIL ('+str(len(bad_mfd))+' bad tets)'}")
        if bad_mfd[:3]:
            for tet, cnt in bad_mfd[:3]:
                print(f"    tet {tet} -> {cnt} simplices")

        print(f"  all simplex types CDT-valid (set in [{{4,1}},{{3,2}}]): "
              f"{'PASS' if not bad_types else 'FAIL ('+str(len(bad_types))+' bad)'}")
        if bad_types[:3]:
            for s, ty in bad_types[:3]:
                print(f"    simplex {s} has type {ty}")

        print(f"  Stage 1.C link-as-S^3 check: "
              f"{'PASS' if ok_link else f'FAIL ({len(link_failures)} vertices)'}")
        if link_failures and not ok_link:
            for v_id, fs in list(link_failures.items())[:3]:
                print(f"    vertex {v_id}: {fs}")

        print(f"  Move candidate counts (zero -> move cannot fire):")
        for move, cnt in counts.items():
            mark = "OK" if cnt > 0 else "ZERO"
            print(f"    {move:>6}: {cnt:>4}  {mark}")
        print()
