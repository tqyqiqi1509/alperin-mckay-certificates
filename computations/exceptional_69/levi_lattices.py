"""Independent integral Levi data for construction 69.

Input: the Bourbaki E7/E8 Cartan graphs only.  No predecessor's
root, Weyl, Frobenius or centre arrays are read.  This is a uniform
root-datum certificate, not a finite-q character calculation.
"""
from collections import deque
from fractions import Fraction
from pathlib import Path
import json
import math
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

OUT = Path(__file__).resolve().parent


def cartan(n):
    C = sp.eye(n) * 2
    edges = [(0, 2), (2, 3), (3, 4), (4, 5), (1, 3)]
    edges += [(i, i + 1) for i in range(5, n - 1)]
    for a, b in edges:
        C[a, b] = C[b, a] = -1
    return C


def roots_from_reflections(C):
    n = C.rows
    seeds = [tuple(int(i == j) for i in range(n)) for j in range(n)]
    roots = set(seeds)
    queue = deque(seeds)
    while queue:
        v = queue.popleft()
        for i in range(n):
            w = list(v)
            w[i] -= sum(int(C[i, j]) * v[j] for j in range(n))
            w = tuple(w)
            if w not in roots:
                roots.add(w)
                queue.append(w)
    assert all((sp.Matrix(v).T * C * sp.Matrix(v))[0] == 2 for v in roots)
    return roots


def simple_reflections(C):
    n = C.rows
    ans = []
    for i in range(n):
        s = sp.eye(n)
        for j in range(n):
            s[i, j] -= C[i, j]
        assert s * s == sp.eye(n)
        assert s.T * C * s == C
        ans.append(s)
    return ans


def longest(C, I):
    S = simple_reflections(C)
    w = sp.eye(C.rows)
    word = []
    while True:
        choices = [i for i in I if all(a >= 0 for a in w[:, i])]
        if not choices:
            break
        i = choices[0]
        word.append(i + 1)
        w = w * S[i]
        assert len(word) <= 120
    assert w * w == sp.eye(C.rows)
    return w, word


def integer_nullbasis(A):
    # For the four displayed subsystems the rational nullspace basis
    # normalizes to a saturated integral basis.  Saturation is checked
    # below by gcd of all maximal minors, not assumed.
    basis = []
    for v in A.nullspace():
        den = sp.ilcm(*[a.q for a in v])
        z = v * den
        content = math.gcd(*[int(a) for a in z])
        basis.append(z / content)
    V = sp.Matrix.hstack(*basis)
    import itertools
    minors = [abs(int(V[list(J), :].det()))
              for J in itertools.combinations(range(V.rows), V.cols)]
    assert math.gcd(*minors) == 1
    return V


def data(n, I, name):
    C = cartan(n)
    roots = roots_from_reflections(C)
    Hroots = [v for v in roots if all(v[j] == 0 for j in range(n) if j not in I)]
    A = C[I, :]
    V = integer_nullbasis(A)
    outside = [i for i in range(n) if i not in I]
    HplusZ = sp.Matrix.hstack(*[sp.eye(n)[:, i] for i in I], V)
    S = smith_normal_form(A, domain=ZZ)
    quotient = smith_normal_form(V[outside, :], domain=ZZ)
    wH, wordH = longest(C, I)
    wG, wordG = longest(C, list(range(n)))
    assert wG == -sp.eye(n)
    r = wG * wH
    assert r * r == sp.eye(n)
    assert r * V == -V
    permutation = {}
    for i in I:
        image = r[:, i]
        j = next(j for j in I if image == sp.eye(n)[:, j])
        permutation[i + 1] = j + 1
    answer = {
        "name": name,
        "ambient_rank": n,
        "levi_nodes": [i + 1 for i in I],
        "ambient_roots": len(roots),
        "derived_roots": len(Hroots),
        "cartan_determinant": int(C.det()),
        "centre_character_smith": [abs(int(S[i, i])) for i in range(len(I))],
        "connected_centre": all(abs(int(S[i, i])) == 1 for i in range(len(I))),
        "primitive_centre_cocharacters": [[int(a) for a in V[:, j]] for j in range(V.cols)],
        "quotient_torus_inclusion": [[int(a) for a in row] for row in V[outside, :].tolist()],
        "central_isogeny_smith": [abs(int(quotient[i, i])) for i in range(V.cols)],
        "central_isogeny_index": abs(int(HplusZ.det())),
        "opposition_return": [[int(a) for a in row] for row in r.tolist()],
        "derived_pinning_permutation": permutation,
        "longest_derived_word": wordH,
        "longest_ambient_length": len(wordG),
        "return_on_centre": [[int(a) for a in row] for row in (-sp.eye(V.cols)).tolist()],
    }
    if V.cols == 1:
        m = int(V[outside[0], 0])
        answer["central_kernel_generator_in_derived_coroots"] = [
            str(Fraction(int(V[i, 0]), m) % 1) for i in I]
        # pi(z_v(t))=t^m is literal: pi has character omega_j.
        answer["power_map_exponent"] = m
    return answer


if __name__ == "__main__":
    values = [
        data(7, list(range(6)), "E7_E6T1"),
        data(8, list(range(6)) + [7], "E8_E6A1T1"),
        data(8, list(range(7)), "E8_E7T1"),
        data(8, list(range(6)), "E8_E6T2"),
    ]
    certificate = {
        "main_theorem_closed": False,
        "input": "Bourbaki Cartan graph only",
        "meaning": "Uniform integral root-datum identities; no character theorem is inferred by computation.",
        "levis": values,
    }
    (OUT / "levi_lattices.json").write_text(json.dumps(certificate, indent=2) + "\n")
    for v in values:
        print(v["name"], "roots", v["ambient_roots"], v["derived_roots"],
              "connected centre", v["connected_centre"], "kernel index", v["central_isogeny_index"],
              "cochars", v["primitive_centre_cocharacters"],
              "pinning", v["derived_pinning_permutation"])
