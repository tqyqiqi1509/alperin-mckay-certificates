#!/usr/bin/env python3
"""Exact checks for the integral G2 torus action in the accompanying paper.

This program checks the finite matrix and root identities only. It does not
replace the cited block classification or the proof of the inductive condition.
It uses only the Python standard library.
"""
import json
from collections import deque

I = (1, 0, 0, 1)
r = (-1, -1, 1, 0)
u = (0, 1, 1, 0)
z = (-1, 0, 0, -1)


def mul(a, b):
    return (a[0]*b[0] + a[1]*b[2],
            a[0]*b[1] + a[1]*b[3],
            a[2]*b[0] + a[3]*b[2],
            a[2]*b[1] + a[3]*b[3])


assert mul(mul(r, r), r) == I
assert mul(u, u) == I
assert mul(mul(u, r), u) == mul(r, r)
assert mul(z, z) == I
assert mul(z, r) == mul(r, z)
assert mul(z, u) == mul(u, z)
group = {I}
queue = deque([I])
while queue:
    x = queue.popleft()
    for g in [r, u, z]:
        y = mul(x, g)
        if y not in group:
            group.add(y)
            queue.append(y)
assert len(group) == 12
assert all(entry in {-1, 0, 1} for g in group for entry in g)
for prime in [3, 5]:
    assert len({tuple(v % prime for v in g) for g in group}) == 12

# The six positive/negative pairs of root characters in z1,z2 coordinates.
positive = [(1, 0), (0, 1), (-1, -1), (1, -1), (2, 1), (1, 2)]
roots = sorted(set(positive + [(-a, -b) for a, b in positive]))
assert len(roots) == 12
# a1 = (zeta,zeta^-1), a2 = (zeta,1).
values = [{"root": [a, b], "on_a1": a-b, "on_a2": a}
          for a, b in roots]
assert sum(item["on_a1"] == 0 for item in values) == 2
assert all(item["on_a2"] != 0
           for item in values if item["on_a1"] == 0)
# No root vanishes on both generators modulo an odd prime: at least one
# evaluation is +1 or -1 for each root, except (1,-1) and its negative,
# whose second evaluation is still +1 or -1.
assert all(abs(item["on_a1"]) == 1 or abs(item["on_a2"]) == 1
           for item in values)

print(json.dumps({
    "weyl_group_order": len(group),
    "matrices": [list(g) for g in sorted(group)],
    "matrix_entries_in_zero_plus_minus_one": True,
    "faithful_reduction_mod_3_and_5": True,
    "all_odd_primes_follow_from_entry_bound": True,
    "root_count": len(roots),
    "root_evaluations": values,
    "root_generator_chain_verified": True
}, indent=2))
