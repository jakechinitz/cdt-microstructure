# v5_cdt_link.py -- Stage 1.C: link condition checker.
#
# For a 4D simplicial complex T to be a PL 4-manifold, the link of every
# vertex must be a triangulation of S^3. Verifying "is this 3-complex S^3?"
# rigorously is undecidable in general (Markov), but practical CDT codes
# use a chain of necessary conditions that catch the vast majority of
# violations:
#
#   (1) The link is a pure 3-complex (every 2-simplex is in some 3-simplex)
#   (2) Every 2-simplex in the link is shared by exactly 2 3-simplices
#       (3-manifold condition)
#   (3) The link is connected
#   (4) Euler characteristic chi(link) = 0  (chi(S^3) = 0)
#   (5) The link of every edge in the link is a circle (1-manifold,
#       chi=0, every vertex shared by exactly 2 edges)
#
# Conditions (1)-(5) are NECESSARY for S^3 but not sufficient: a non-trivial
# homotopy 3-sphere (Poincare-like manifold) would pass all five. Catching
# such things requires computing the fundamental group or normal-surface
# analysis (real combinatorial topology), which is genuinely beyond scope
# for a single session. For CDT-style simulations these conditions are
# known to be effective in practice; CDT-plusplus uses an equivalent
# chain. See Ambjorn-Jurkiewicz-Loll hep-th/0105267 Sec 7.

from __future__ import annotations
from itertools import combinations

from v5_cdt_lib import CausalTriangulation


def vertex_link(T: CausalTriangulation, v):
    """Compute the link of vertex v as a 3-dim simplicial complex.

    Returns a dict with:
        verts: set of vertices in the link
        edges: set of edge tuples in the link
        tris: set of triangle tuples in the link
        tets: set of tetrahedron tuples in the link
    """
    # Find all 4-simplices containing v; their tetrahedral faces opposite to
    # v form the 3-simplices of link(v).
    link_tets = set()
    for s in T.simplices:
        if v in s:
            opp = tuple(sorted(x for x in s if x != v))
            link_tets.add(opp)
    # Derive lower-dim faces of the link
    link_tris = set()
    link_edges = set()
    link_verts = set()
    for tet in link_tets:
        link_verts.update(tet)
        for tri in combinations(tet, 3):
            link_tris.add(tuple(sorted(tri)))
        for edge in combinations(tet, 2):
            link_edges.add(tuple(sorted(edge)))
    return {
        "verts": link_verts,
        "edges": link_edges,
        "tris": link_tris,
        "tets": link_tets,
    }


def check_3manifold(link):
    """Check that the link is a closed 3-manifold:
        - every 2-simplex shared by exactly 2 3-simplices
        - every 1-simplex's "link in the link" is a circle
    Returns (ok, reason_or_None).
    """
    # Triangle-tet sharing: each triangle should be in exactly 2 tets
    tri_count = {}
    for tet in link["tets"]:
        for tri in combinations(tet, 3):
            t = tuple(sorted(tri))
            tri_count[t] = tri_count.get(t, 0) + 1
    for tri, cnt in tri_count.items():
        if cnt != 2:
            return False, f"triangle {tri} shared by {cnt} tetrahedra (need 2)"
    # Edge-link-as-circle check: for each edge, the set of triangles
    # containing it forms a 1-complex; the vertices of that 1-complex (the
    # third-vertex-of-each-triangle) should pair up to form a single cycle.
    for edge in link["edges"]:
        u, w = edge
        third_vertices = []
        for tri in link["tris"]:
            if u in tri and w in tri:
                third_vertices.append(next(x for x in tri if x != u and x != w))
        if len(third_vertices) < 3:
            return False, f"edge {edge} in only {len(third_vertices)} triangles"
        # The third vertices, paired by which tetrahedron pairs them, must
        # form a single cycle. Build adjacency among third vertices: two
        # third vertices a, b are adjacent iff there's a tet containing
        # {edge, a, b}.
        adj = {x: set() for x in third_vertices}
        for tet in link["tets"]:
            if u in tet and w in tet:
                rest = [x for x in tet if x != u and x != w]
                if len(rest) == 2:
                    a, b = rest
                    adj[a].add(b); adj[b].add(a)
        # Every vertex of the adjacency must have degree exactly 2 (for a
        # circle).
        for x, nbrs in adj.items():
            if len(nbrs) != 2:
                return False, f"edge {edge} link not a circle: vertex {x} has degree {len(nbrs)}"
        # And it should be connected (single cycle, not multiple)
        if not _connected(adj):
            return False, f"edge {edge} link not connected"
    return True, None


def _connected(adj):
    """Check if the graph given by adjacency dict is connected."""
    if not adj:
        return True
    start = next(iter(adj))
    visited = {start}
    stack = [start]
    while stack:
        x = stack.pop()
        for n in adj[x]:
            if n not in visited:
                visited.add(n); stack.append(n)
    return len(visited) == len(adj)


def euler_char_3d(link):
    """Euler char of the 3-complex: V - E + F - C."""
    return len(link["verts"]) - len(link["edges"]) + len(link["tris"]) - len(link["tets"])


def check_link_is_s3_necessary_conditions(T: CausalTriangulation, v):
    """Check necessary conditions for link(v) to be a triangulation of S^3.
    Returns (ok, list_of_failures)."""
    link = vertex_link(T, v)
    failures = []
    # (1) Non-empty
    if not link["tets"]:
        failures.append("link is empty")
        return False, failures
    # (4) Euler characteristic
    chi = euler_char_3d(link)
    if chi != 0:
        failures.append(f"euler char = {chi} != 0")
    # (3) Connected (on the tets, via shared triangles)
    tet_adj = {tet: set() for tet in link["tets"]}
    for tri, cnt in _tri_tet_map(link).items():
        if cnt[0] > 1:
            for i in range(len(cnt[1])):
                for j in range(i+1, len(cnt[1])):
                    tet_adj[cnt[1][i]].add(cnt[1][j])
                    tet_adj[cnt[1][j]].add(cnt[1][i])
    if not _connected(tet_adj):
        failures.append("link tets not connected")
    # (2,5) 3-manifold conditions
    is_3man, reason = check_3manifold(link)
    if not is_3man:
        failures.append(f"not a closed 3-manifold: {reason}")
    return len(failures) == 0, failures


def _tri_tet_map(link):
    """Map each triangle to (count, list of tets containing it)."""
    m = {}
    for tet in link["tets"]:
        for tri in combinations(tet, 3):
            t = tuple(sorted(tri))
            if t not in m:
                m[t] = [0, []]
            m[t][0] += 1
            m[t][1].append(tet)
    return m


def verify_manifold(T: CausalTriangulation, verbose: bool = False):
    """Run link-condition check on every vertex. Returns (ok, dict_of_failures)."""
    failures = {}
    for v in list(T.vertex_time.keys()):
        ok, fails = check_link_is_s3_necessary_conditions(T, v)
        if not ok:
            failures[v] = fails
            if verbose:
                print(f"  vertex {v}: link condition FAILS: {fails}")
    return len(failures) == 0, failures


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from v5_cdt_lib import initial_boundary_of_5simplex, initial_grown

    print("=" * 76)
    print("  Stage 1.C: link condition checker -- self-test")
    print("=" * 76)
    print()

    # Test 1: boundary of 5-simplex is a known PL S^4. Every vertex link
    # should be a PL S^3.
    print("Test 1: boundary of 5-simplex (known S^4)")
    T = initial_boundary_of_5simplex(time_split=(3, 3))
    ok, fails = verify_manifold(T)
    print(f"  manifold check: {'PASS' if ok else 'FAIL'}")
    if not ok:
        for v, fs in fails.items():
            print(f"    vertex {v}: {fs}")

    # Test 2: an (1,5)-grown triangulation -- still a PL S^4, should pass.
    print()
    print("Test 2: initial_grown(target_N4=30) -- (1,5)-grown S^4")
    T = initial_grown(target_N4=30, seed=0)
    ok, fails = verify_manifold(T)
    print(f"  manifold check: {'PASS' if ok else 'FAIL'}")
    if not ok:
        print(f"  {len(fails)} vertices fail link condition")
        for v, fs in list(fails.items())[:3]:
            print(f"    vertex {v}: {fs}")

    # Test 3: verify that after running the full Stage 1.B Metropolis loop,
    # the manifold property is preserved (= moves don't break manifoldness).
    print()
    print("Test 3: post-Stage-1.B-smoke manifold check")
    import numpy as np
    from v5_cdt_moves import propose_move_tracked, undo_move
    from v5_cdt_lib import total_action
    rng = np.random.default_rng(0)
    T = initial_grown(target_N4=30, seed=0)
    S = total_action(T, 3.0, 0.4, 1.0)
    for _ in range(500):
        move_type, ok, removed, added, hastings, vch = propose_move_tracked(T, rng)
        if not ok:
            continue
        S_new = total_action(T, 3.0, 0.4, 1.0)
        if rng.random() < np.exp(-(S_new - S) + hastings):
            S = S_new
        else:
            undo_move(T, move_type, removed, added, vch)
    ok, fails = verify_manifold(T)
    print(f"  manifold check after 500 attempted moves: {'PASS' if ok else 'FAIL'}")
    if not ok:
        print(f"  {len(fails)} vertices broke manifoldness (moves are NOT manifold-preserving)")
        for v, fs in list(fails.items())[:3]:
            print(f"    vertex {v}: {fs}")
        print(f"  This is the diagnostic Stage 1.C exists to surface --")
        print(f"  flags moves that pass collision checks but violate link-S^3.")
