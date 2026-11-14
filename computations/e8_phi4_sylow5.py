"""The E8 Phi4 torus and its elliptic Sylow-5 action, for all field sizes.

Roots use doubled orthonormal coordinates; the final certificates use the
integral Bourbaki E8 lattice. GAP supplies an actual order-five element in the
240-point root permutation representation of the relative group G31.
"""
from itertools import combinations, product
from pathlib import Path
import ast
import json
import re
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

OUT = Path(__file__).resolve().parent
roots = []
for i, j in combinations(range(8), 2):
    for si, sj in product([-2, 2], repeat=2):
        a = [0] * 8
        a[i], a[j] = si, sj
        roots.append(tuple(a))
roots += [a for a in product([-1, 1], repeat=8) if sum(x < 0 for x in a) % 2 == 0]
roots = sorted(roots)
assert len(roots) == len(set(roots)) == 240
idx = {a: i + 1 for i, a in enumerate(roots)}

def w_action(a):
    ans = []
    for i in range(0, 8, 2):
        ans.extend([-a[i + 1], a[i]])
    return tuple(ans)

def reflect(a, x):
    dot = sum(ai * xi for ai, xi in zip(a, x))
    assert dot % 4 == 0
    return tuple(xi - dot // 4 * ai for ai, xi in zip(a, x))

generators = set()
for a in roots:
    wa = w_action(a)
    assert sum(ai * bi for ai, bi in zip(a, wa)) == 0
    permutation = tuple(idx[reflect(a, reflect(wa, x))] for x in roots)
    assert all(
        roots[permutation[idx[w_action(x)] - 1] - 1] ==
        w_action(roots[permutation[idx[x] - 1] - 1])
        for x in roots
    )
    generators.add(permutation)
assert len(generators) == 60

simple = [
    (1, -1, -1, -1, -1, -1, -1, 1),
    (2, 2, 0, 0, 0, 0, 0, 0),
    (-2, 2, 0, 0, 0, 0, 0, 0),
    (0, -2, 2, 0, 0, 0, 0, 0),
    (0, 0, -2, 2, 0, 0, 0, 0),
    (0, 0, 0, -2, 2, 0, 0, 0),
    (0, 0, 0, 0, -2, 2, 0, 0),
    (0, 0, 0, 0, 0, -2, 2, 0),
]
B = sp.Matrix.hstack(*(sp.Matrix(a) for a in simple))
C = B.T * B / 4
assert C.det() == 1
w = B.inv() * sp.Matrix.hstack(*(sp.Matrix(w_action(a)) for a in simple))
assert all(x.q == 1 for x in w)
assert w**2 == -sp.eye(8)
assert w.charpoly().as_expr().expand() == ((sp.Symbol("lambda")**2 + 1)**4).expand()

def gaplist(xs):
    return "[" + ",".join(map(str, xs)) + "]"

gap = "W4:=Group(" + ",".join(
    "PermList(" + gaplist(g) + ")" for g in sorted(generators)
) + ");;\n"
gap += "S:=SylowSubgroup(W4,5);; u:=GeneratorsOfGroup(S)[1];;\n"
gap += 'Print("ORDER|",Size(W4),"\\nSYLOW|",Size(S),"\\n");\n'
gap += 'Print("ACTION|",List(' + gaplist([idx[a] for a in simple]) + ',i->i^u),"\\n");\n'
gap += "QUIT;\n"
(OUT / "e8_phi4_sylow5.g").write_text(gap)

def nullspace_mod(matrix, p):
    a = [[int(x) % p for x in row] for row in matrix.tolist()]
    rows, cols = len(a), len(a[0])
    pivots = []
    k = 0
    for j in range(cols):
        pivot = next((i for i in range(k, rows) if a[i][j]), None)
        if pivot is None:
            continue
        a[k], a[pivot] = a[pivot], a[k]
        fac = pow(a[k][j], -1, p)
        a[k] = [fac * x % p for x in a[k]]
        for i in range(rows):
            if i != k:
                fac = a[i][j]
                a[i] = [(x - fac * y) % p for x, y in zip(a[i], a[k])]
        pivots.append(j)
        k += 1
    basis = []
    for j in range(cols):
        if j in pivots:
            continue
        v = [0] * cols
        v[j] = 1
        for i, pivot in enumerate(pivots):
            v[pivot] = -a[i][j] % p
        basis.append(v)
    return sp.Matrix.hstack(*(sp.Matrix(v) for v in basis))

root_coordinates = [B.inv() * sp.Matrix(a) for a in roots]
assert all(x.q == 1 for r in root_coordinates for x in r)
torus_data = []
for qmod in [2, 3]:
    base = nullspace_mod(qmod * w - sp.eye(8), 5)
    assert base.cols == 4
    vanished = [
        roots[i] for i, r in enumerate(root_coordinates)
        if all(int(x) % 5 == 0 for x in r.T * C * base)
    ]
    assert vanished == []
    torus_data.append({
        "q_mod5": qmod, "dimension": base.cols,
        "cocharacter_basis_mod5": base.tolist(),
        "roots_centralizing_whole_A5": len(vanished),
    })

record = {
    "scope": "all q with ord_5(q)=4 and every a=v5(q^2+1)>=1",
    "roots": 240, "commuting_reflections": 60,
    "E8_gram": C.tolist(), "phi4_operator": w.tolist(),
    "torus_mod5": torus_data,
}
observed = OUT / "e8_phi4_sylow5_gap_observed.txt"
if observed.exists():
    s = observed.read_text()
    assert "Error," not in s
    order = int(re.search(r"ORDER\|(\d+)", s)[1])
    sylow = int(re.search(r"SYLOW\|(\d+)", s)[1])
    images = ast.literal_eval(re.search(r"ACTION\|(\[[\s\d,]*\])", s)[1])
    assert order == 46080 and sylow == 5
    u = B.inv() * sp.Matrix.hstack(*(sp.Matrix(roots[i - 1]) for i in images))
    assert all(x.q == 1 for x in u)
    assert u**5 == sp.eye(8) and u != sp.eye(8)
    assert u * w == w * u and u.T * C * u == C
    assert (u - sp.eye(8)).rank() == 8
    d = smith_normal_form(u - sp.eye(8), domain=ZZ)
    invariants = [abs(int(d[i, i])) for i in range(8)]
    assert all(x in [1, 5] for x in invariants)
    assert abs((u - sp.eye(8)).det()) == 25
    record.update({
        "relative_group_order": order, "relative_sylow5_order": sylow,
        "sylow_generator_integral": u.tolist(),
        "fixed_space_dimension_Q": 0, "coinvariant_smith": invariants,
        "coinvariant_order": 25,
    })

record = json.loads(json.dumps(record, default=int))
(OUT / "e8_phi4_sylow5.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps({k: record[k] for k in record if k not in
    ["E8_gram", "phi4_operator", "torus_mod5", "sylow_generator_integral"]}))
