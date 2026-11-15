"""Independent integer certificate for E7 source lattices and signed carrier.

No earlier certificate is imported.  The minuscule Chevalley generators are
built afresh from the 56 weights and the seven displayed simple roots.
"""
from collections import deque
from functools import reduce
from itertools import combinations
from math import gcd
from pathlib import Path
import json
import numpy as np
import sympy as sy
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

C = 2 * sy.eye(7)
for i, j in [(1, 3), (3, 4), (4, 5), (5, 6), (6, 7), (2, 4)]:
    C[i - 1, j - 1] = C[j - 1, i - 1] = -1
K = sy.Matrix.hstack(*[sy.Matrix(v) for v in [
    (0, 0, 1, 0, 0, 0, 0), (0, 1, 1, 2, 2, 2, 1),
    (0, 1, 1, 2, 1, 0, 0), (1, 0, 0, 0, 0, 0, 0)]])
k = sy.Matrix([0, 1, 0, 0, 1, 0, 1])
H = sy.eye(7)[:, [1, 4, 6]]
regroots = C.row_join(-C * k / 2)
KC = K.col_join(sy.zeros(1, 4)).row_join(k.col_join(sy.Matrix([2])))
def minor_gcd(M):
    return reduce(gcd, [abs(int(M[list(rows), :].det()))
                       for rows in combinations(range(M.rows), M.cols)], 0)
def smith_entries(M):
    D = smith_normal_form(M, domain=ZZ)
    return [abs(int(D[i, i])) for i in range(min(M.shape))]

out = {"K_Gram": (K.T * C * K).tolist(), "K_primitive": minor_gcd(K),
       "HK_index": int(abs(sy.Matrix.hstack(H, K).det())),
       "natural_H_root_Smith": smith_entries(C[[1, 4, 6], :]),
       "regular_H_root_Smith": smith_entries(regroots[[1, 4, 6], :]),
       "K_plus_regular_center_primitive": minor_gcd(KC)}
K0 = K.row_join(sy.eye(7)[:, 6])
K0C = K0.col_join(sy.zeros(1, 5)).row_join(k.col_join(sy.Matrix([2])))
out.update({"A1_square_perpendicular_primitive": minor_gcd(K0),
            "A1_square_plus_regular_center_index": minor_gcd(K0C),
            "natural_A1_square_root_Smith": smith_entries(C[[1, 4], :]),
            "regular_A1_square_root_Smith": smith_entries(regroots[[1, 4], :])})

# The eight-coordinate realization; all vectors are scaled by four.
ee = np.eye(8, dtype=np.int64)
simples = np.array([4 * (ee[5] - ee[6]),
                   2 * np.array([-1, -1, -1, -1, 1, 1, 1, 1]),
                   4 * (ee[4] - ee[5]), 4 * (ee[3] - ee[4]),
                   4 * (ee[2] - ee[3]), 4 * (ee[1] - ee[2]),
                   4 * (ee[0] - ee[1])])
weights = sorted({tuple(sign * (4 * ee[i] + 4 * ee[j] - np.ones(8, dtype=np.int64)))
                  for i, j in combinations(range(8), 2) for sign in [1, -1]})
idx = {w: i for i, w in enumerate(weights)}
weights = np.array(weights, dtype=np.int64)
dim = len(weights)
ident = np.eye(dim, dtype=np.int64)
E, F, N = [], [], []
for a in simples:
    ep = np.zeros((dim, dim), dtype=np.int64)
    fm = ep.copy()
    for j, mu in enumerate(weights):
        val = int(mu @ a)
        assert val % 16 == 0
        m = val // 16
        assert m in [-1, 0, 1]
        if m == -1:
            ep[idx[tuple(mu + a)], j] = 1
        if m == 1:
            fm[idx[tuple(mu - a)], j] = 1
    E.append(ep); F.append(fm)
    N.append((ident + ep) @ (ident - fm) @ (ident + ep))

checks = 0
for i in range(7):
    for j in range(7):
        expected = np.diag((weights @ simples[i]) // 16) if i == j else np.zeros_like(ident)
        assert np.array_equal(E[i] @ F[j] - F[j] @ E[i], expected)
        checks += 1
        if i != j:
            if C[i, j] == 0:
                assert np.array_equal(E[i] @ E[j], E[j] @ E[i])
                assert np.array_equal(F[i] @ F[j], F[j] @ F[i])
            else:
                assert not np.any(E[i] @ E[j] @ E[i])
                assert not np.any(F[i] @ F[j] @ F[i])
            checks += 2
for i in range(7):
    assert np.array_equal(N[i].T @ N[i], ident)
    for j in range(i):
        if C[i, j] == 0:
            assert np.array_equal(N[i] @ N[j], N[j] @ N[i])
        else:
            assert np.array_equal(N[i] @ N[j] @ N[i], N[j] @ N[i] @ N[j])
        checks += 1

def product(mats):
    result = ident.copy()
    for a in mats:
        result = result @ a
    return result
cox = product(N)
n0 = np.linalg.matrix_power(cox, 9)
nH = N[0] @ N[2] @ N[0]
z0 = n0 @ nH.T
assert np.array_equal(n0 @ n0, -ident)
assert np.array_equal(nH @ nH, ident)
assert np.array_equal(z0 @ z0, -ident)
assert all(np.array_equal(n0 @ n, n @ n0) for n in N)
assert np.array_equal(z0 @ E[0] @ z0.T, E[2])
assert np.array_equal(z0 @ E[2] @ z0.T, E[0])

# A5 simple roots perpendicular to A2=<alpha1,alpha3>.
perp = [4 * (ee[0] - ee[1]), 4 * (ee[1] - ee[2]),
        4 * (ee[2] - ee[3]), 4 * (ee[3] - ee[7]), simples[1]]
K5 = sy.Matrix.hstack(*[sy.Matrix(simples).T.gauss_jordan_solve(sy.Matrix(b))[0]
                       for b in perp])
H3 = sy.eye(7)[:, [0, 2]]
out.update({"A2_perpendicular_coroot_basis": K5.tolist(),
            "A2_perpendicular_primitive": minor_gcd(K5),
            "A2_A5_index": int(abs(sy.Matrix.hstack(H3, K5).det())),
            "natural_A2_root_Smith": smith_entries(C[[0, 2], :])})
def reflect(root, a):
    dot = int(root @ a)
    assert dot % 16 == 0
    return root - (dot // 16) * a
# Transport alpha1 through the root orbit, retaining its actual Tits lift.
start = tuple(simples[0])
words = {start: ()}; queue = deque([start])
while queue:
    root = queue.popleft()
    for j, a in enumerate(simples):
        new = tuple(reflect(np.array(root), a))
        if new not in words:
            words[new] = words[root] + (j,)
            queue.append(new)
assert len(words) == 126
P = []
for beta in perp:
    word = words[tuple(beta)]
    transport = product([N[j] for j in reversed(word)])
    nb = transport @ N[0] @ transport.T
    P.append(nb)
    assert np.array_equal(z0 @ nb, nb @ z0)
    assert all(np.array_equal(nb @ E[j], E[j] @ nb) for j in [0, 2])
    assert all(np.array_equal(nb @ F[j], F[j] @ nb) for j in [0, 2])
# The choices may differ by signs on A5 simple root generators. These are a
# torus-conjugate A5 pinning and therefore still satisfy the A5 relations.
for i in range(5):
    for j in range(i):
        if abs(i - j) == 1:
            assert np.array_equal(P[i] @ P[j] @ P[i], P[j] @ P[i] @ P[j])
        else:
            assert np.array_equal(P[i] @ P[j], P[j] @ P[i])
a1 = P[0] @ P[1]; a2 = P[3] @ P[4]
target = [3, 4, 5, 0, 1, 2]; perm = list(range(6)); exchange_word = []
for i, value in enumerate(target):
    j = perm.index(value)
    while j > i:
        perm[j], perm[j-1] = perm[j-1], perm[j]
        exchange_word.append(j-1); j -= 1
n = product([P[j] for j in exchange_word])
assert np.array_equal(np.linalg.matrix_power(a1, 3), ident)
assert np.array_equal(np.linalg.matrix_power(a2, 3), ident)
assert np.array_equal(a1 @ a2, a2 @ a1)
assert np.array_equal(n @ n, -ident)
assert np.array_equal(n @ a1 @ n.T, a2)
assert np.array_equal(z0 @ n, n @ z0)
u3 = a1 @ a2; u6 = u3 @ z0
assert np.array_equal(np.linalg.matrix_power(u3, 3), ident)
assert np.array_equal(np.linalg.matrix_power(u6, 6), -ident)
assert np.array_equal(np.linalg.matrix_power(u6, 12), ident)
out.update({"weights": dim, "roots": len(words), "simple_Chevalley_and_braid_checks": checks,
            "longest_square": "-I56", "A2_longest_square": "I56",
            "z0_square": "-I56", "z0_positive_A2_diagram": True,
            "A5_exchange_word_1_based": [j+1 for j in exchange_word],
            "A5_perpendicular_Tits_words_1_based": [[j+1 for j in words[tuple(b)]] for b in perp],
            "signed_carrier_relations": True, "u3_order": 3, "u6_order": 12,
            "u6_sixth_power": "-I56"})

# Convert exact SymPy integers before serializing.
out = json.loads(json.dumps(out, default=int))
path = Path(__file__).with_suffix('.json')
path.write_text(json.dumps(out, indent=2))
print(json.dumps(out))
