#!/usr/bin/env python3
"""Independent D4 root check for the third triality review.

This verifies the lattice identities and their reduction to the unequal
3-power module. It does not test block membership or prove iAM.
"""
import json
from collections import deque
from math import gcd
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

C = sp.Matrix([[2, -1, 0, 0], [-1, 2, -1, -1],
               [0, -1, 2, 0], [0, -1, 0, 2]])
eye = sp.eye(4)
simple = [eye[:, i] for i in range(4)]
reflections = [eye - r * (r.T * C) for r in simple]
roots = set(tuple(r) for r in simple)
queue = deque(roots)
while queue:
    r = sp.Matrix(queue.popleft())
    for u in reflections:
        t = tuple(u * r)
        if t not in roots:
            roots.add(t)
            queue.append(t)
assert len(roots) == 24

alpha = sp.Matrix([0, 1, 0, 0])
beta = sp.ones(4, 1)
r_alpha = eye - alpha * (alpha.T * C)
r_beta = eye - beta * (beta.T * C)
w = r_alpha * r_beta
assert w**3 == eye
assert w - eye == sp.Matrix.hstack(-beta, beta-alpha, -beta, -beta)
assert abs(int(sp.Matrix.hstack(alpha, beta, simple[0], simple[2]).det())) == 1

# After quotienting by alpha,beta, take x=r1,y=r3 and r4=-x-y.
projection = sp.Matrix([[1, 0, 0, -1], [0, 0, 1, -1]])
rho = sp.zeros(4)
for source, target in [(0, 2), (1, 1), (2, 3), (3, 0)]:
    rho[target, source] = 1
R = sp.Matrix([[0, -1], [1, -1]])
assert projection * rho == R * projection
restricted = {}
for r in roots:
    t = tuple(projection * sp.Matrix(r))
    restricted[t] = restricted.get(t, 0) + 1
assert restricted[(0, 0)] == 6
assert sorted(v for k, v in restricted.items() if k != (0, 0)) == [3]*6

# All roots are nontrivial on the full unequal module A.
mod3 = sp.Matrix([[1, 0, 1, 1], [0, 1, 0, 0]])
assert all(any(int(v) % 3 for v in mod3*sp.Matrix(r)) for r in roots)

fields = []
for q in [2, 4, 5, 7, 8, 11, 13, 16, 19, 25, 31, 37, 49, 64, 73]:
    if q % 3 == 0:
        continue
    eps = 1 if q % 3 == 1 else -1
    value = q-eps
    a = 0
    while value % 3 == 0:
        value //= 3
        a += 1
    small, large = 3**a, 3**(a+1)
    m = q*q + eps*q + 1
    assert m % 3 == 0 and gcd(m//3, 6) == 1
    assert (m//3 - 1) % small == 0
    lattice = sp.Matrix([[1, 0], [0, 3], [eps*q, 0], [q*q, 0]])
    sigma = sp.Matrix([[-2, 3], [-1, 1]])
    tau = sp.Matrix([[1, 0], [1, -1]])
    assert all(int(v) % large == 0 for v in w*lattice-lattice*sigma)
    assert all(int(v) % large == 0 for v in r_alpha*lattice-lattice*tau)
    smith = smith_normal_form(q*R-eps*sp.eye(2), domain=sp.ZZ)
    assert abs(int(smith[0, 0])) == 1
    assert abs(int(smith[1, 1])) == m
    fields.append({"q": q, "epsilon": eps, "a": a,
                   "fixed_torus_order": m,
                   "fixed_torus_smith_invariants": [1, m]})

print(json.dumps({
    "root_count": len(roots),
    "cartan_determinant": int(C.det()),
    "w": [[int(v) for v in row] for row in w.tolist()],
    "primitive_fixed_torus_quotient": True,
    "roots_zero_on_fixed_torus": 6,
    "roots_zero_on_A": 0,
    "integer_specialization_checks": fields
}, indent=2))
