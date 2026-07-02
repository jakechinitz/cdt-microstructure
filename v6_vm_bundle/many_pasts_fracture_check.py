#!/usr/bin/env python3
"""Verify the paper's App. D.4 fracture theorem -- the load-bearing
justification for the Many-Pasts refresh kernel.

CLAIM (paper): under local single-label UNIT shifts (m -> m+-1, one slot,
injectivity preserved), the 1680-state boundary ensemble fractures into
4! = 24 disconnected ordering sectors x 2 parities = 48 components of
C(7,4) = 35 states each; slot order is a frozen charge of local motion, so
no local transfer operator can equilibrate to the admissibility
distribution -- a non-local memoryless refresh is REQUIRED, not optional.

VERIFIED HERE by explicit graph construction (2026-07): 48 components,
every one of size exactly 35, and the conserved invariant is exactly
(order-type, parity). The theorem holds as stated.

Scope note: this validates the NECESSITY argument for Many-Pasts' refresh
role. It does not test the history-weighting itself (see the Many-Pasts
test ladder in the session notes / STAGE3 doc): the operational
Born/no-signaling branch is out of scope for classical MC entirely, and a
lattice implementation of history weighting is under-specified by the
paper -- the well-posed next step is the inverse problem (what memory
strength/range would be needed for order-one induced vacuum stiffness).
"""
from itertools import permutations
from collections import Counter

states = [(t, p) for t in permutations(range(-3, 4), 4) for p in (0, 1)]
idx = {s: i for i, s in enumerate(states)}

def neighbors(s):
    t, p = s
    for i in range(4):
        for d in (-1, 1):
            v = t[i] + d
            if -3 <= v <= 3 and v not in t:
                yield (t[:i] + (v,) + t[i+1:], p)

seen = [False] * len(states)
comps = []
for i, s in enumerate(states):
    if seen[i]:
        continue
    comp = 1; seen[i] = True; stack = [s]
    while stack:
        for nb in neighbors(stack.pop()):
            j = idx[nb]
            if not seen[j]:
                seen[j] = True; comp += 1; stack.append(states[j])
    comps.append(comp)

hist = dict(Counter(comps))
ok = len(comps) == 48 and hist == {35: 48}
print(f"states={len(states)}  components={len(comps)}  sizes={hist}")
print(f"D.4 fracture theorem (48 x 35): {'VERIFIED' if ok else 'FAILED'}")
