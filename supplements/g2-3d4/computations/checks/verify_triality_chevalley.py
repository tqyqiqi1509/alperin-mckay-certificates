#!/usr/bin/env python3
"""Independent triality arithmetic from integral SO(8) root matrices.

The integral root-lift tests hold in every characteristic after reduction.
Their lift from SO(8) to Spin(8) uses triality-fixedness and the trivial
triality-fixed subgroup of the central kernel, as explained in the report.
Finite torus tests independently derive the formulas from the root action.
They do not certify the global ordinary-character or block-classification
inputs. Python standard library only.
"""
from collections import Counter, deque
from fractions import Fraction
from itertools import product
import json


def eye(n):
    return tuple(tuple(int(i == j) for j in range(n)) for i in range(n))


def mul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def power(a, n):
    out = eye(len(a))
    while n:
        if n & 1:
            out = mul(out, a)
        a = mul(a, a)
        n //= 2
    return out


def inv(a):
    n = len(a)
    rows = [[Fraction(x) for x in a[i]] +
            [Fraction(i == j) for j in range(n)] for i in range(n)]
    for j in range(n):
        i = next(i for i in range(j, n) if rows[i][j])
        rows[i], rows[j] = rows[j], rows[i]
        value = rows[j][j]
        rows[j] = [x / value for x in rows[j]]
        for i in range(n):
            if i != j:
                value = rows[i][j]
                rows[i] = [x - value * y for x, y in zip(rows[i], rows[j])]
    out = tuple(tuple(row[n:]) for row in rows)
    assert mul(a, out) == eye(n) == mul(out, a)
    return out


def integral(a):
    assert all(Fraction(x).denominator == 1 for row in a for x in row)
    return tuple(tuple(int(x) for x in row) for row in a)


def apply(a, v):
    return tuple(sum(x * y for x, y in zip(row, v)) for row in a)


def root_lift(i, j, plus):
    """Chevalley n_r(1) for r=e_i-e_j or e_i+e_j, indices 0..3."""
    pos = [list(r) for r in eye(8)]
    neg = [list(r) for r in eye(8)]
    if plus:
        pos[i][j + 4] += 1
        pos[j][i + 4] -= 1
        # x_{-r}(-1)
        neg[j + 4][i] -= 1
        neg[i + 4][j] += 1
    else:
        pos[i][j] += 1
        pos[j + 4][i + 4] -= 1
        neg[j][i] -= 1
        neg[i + 4][j + 4] += 1
    return mul(mul(tuple(map(tuple, pos)), tuple(map(tuple, neg))),
               tuple(map(tuple, pos)))


def hminus(i, j):
    return tuple(tuple((-1 if a in (i, j, i + 4, j + 4) else 1)
                       if a == b else 0 for b in range(8)) for a in range(8))


N = [root_lift(0, 1, False), root_lift(1, 2, False),
     root_lift(2, 3, False), root_lift(2, 3, True)]
H = [hminus(0, 1), hminus(1, 2), hminus(2, 3), hminus(2, 3)]
assert all(h == power(n, 2) for h, n in zip(H, N))
u = mul(mul(mul(H[0], H[2]), H[3]), N[1])
v = mul(mul(mul(H[1], N[0]), N[2]), N[3])
n0 = mul(mul(mul(N[0], N[2]), N[3]), root_lift(0, 1, True))
sigma = power(mul(u, v), 2)
iota = power(mul(u, v), 3)
assert power(u, 2) == power(v, 2) == power(n0, 2) == eye(8)
assert power(mul(u, v), 6) == eye(8)
assert mul(n0, u) == mul(u, n0) and mul(n0, v) == mul(v, n0)
assert iota == n0
assert power(sigma, 3) == eye(8)
assert mul(mul(u, sigma), u) == power(sigma, 2)


def signed_basis_images(matrix):
    """Encode the image of the jth basis vector by its signed 1-based index."""
    images = []
    for j in range(8):
        nonzero = [(i, matrix[i][j]) for i in range(8) if matrix[i][j]]
        assert len(nonzero) == 1
        i, sign = nonzero[0]
        assert sign in (-1, 1)
        images.append(sign * (i + 1))
    return tuple(images)


# These are the three printed columns of Appendix A, in the basis
# e1,e2,e3,e4,f1,f2,f3,f4. Compare them with the independently multiplied
# Chevalley root matrices, not with matrices reconstructed from this table.
SIGNED_TABLE = {
    "u": (-1, -3, -2, 4, -5, -7, -6, 8),
    "v": (2, 1, 7, -8, 6, 5, 3, -4),
    "n0": (-5, -6, -7, -8, -1, -2, -3, -4),
}
for name, matrix in (("u", u), ("v", v), ("n0", n0)):
    assert signed_basis_images(matrix) == SIGNED_TABLE[name]

K = {eye(8)}
queue = deque(K)
while queue:
    a = queue.popleft()
    for b in (u, v):
        c = mul(a, b)
        if c not in K:
            K.add(c)
            queue.append(c)
assert len(K) == 12
assert len({tuple(tuple(z % 2 for z in row) for row in k) for k in K}) == 12


def class_count(group):
    remaining = set(group)
    count = 0
    while remaining:
        g = next(iter(remaining))
        conjugates = {integral(mul(mul(h, g), inv(h))) for h in group}
        assert conjugates <= group
        remaining -= conjugates
        count += 1
    return count


K0 = {mul(power(sigma, j), power(u, k)) for j in range(3) for k in range(2)}
assert len(K0) == 6
K_CLASSES = class_count(K)
S3_CLASSES = class_count(K0)
assert (K_CLASSES, S3_CLASSES) == (6, 3)

# Columns are the four simple coroots in the orthonormal D4 model.
B = ((1, 0, 0, 0), (-1, 1, 0, 0), (0, -1, 1, 1), (0, 0, -1, 1))
Binv = inv(B)


def cocharacter_action(n):
    action = []
    for row in n[:4]:
        positions = [j for j, z in enumerate(row) if z]
        assert len(positions) == 1
        j = positions[0]
        action.append(tuple((1 if j < 4 else -1) * int(k == j % 4)
                            for k in range(4)))
    return integral(mul(mul(Binv, tuple(action)), B))


U = cocharacter_action(u)
V = cocharacter_action(v)
W = cocharacter_action(sigma)
J = cocharacter_action(iota)
assert J == tuple(tuple(-z for z in row) for row in eye(4))
assert W == power(mul(U, V), 2)
CARTAN = ((2, -1, 0, 0), (-1, 2, -1, -1),
          (0, -1, 2, 0), (0, -1, 0, 2))
simple = tuple(tuple(int(i == j) for i in range(4)) for j in range(4))
reflections = []
for i in range(4):
    reflections.append(tuple(tuple(int(j == k) - int(j == i) * CARTAN[i][k]
                                   for k in range(4)) for j in range(4)))
assert U == reflections[1]
assert V == mul(mul(reflections[0], reflections[2]), reflections[3])
roots = set(simple)
queue = deque(simple)
while queue:
    root = queue.popleft()
    for s in reflections:
        new = apply(s, root)
        if new not in roots:
            roots.add(new)
            queue.append(new)
assert len(roots) == 24
positive = sorted(r for r in roots if all(x >= 0 for x in r))
assert len(positive) == 12
restricted_positive = Counter(((r[0] + r[2] + r[3]) % 3, r[1] % 3)
                              for r in positive)
assert (0, 0) not in restricted_positive


def valuation3(n):
    answer = 0
    while n % 3 == 0:
        n //= 3
        answer += 1
    return answer


def actual_torus_check(q):
    eps = 1 if q % 3 == 1 else -1
    Q = eps * q
    a = valuation3(q - eps)
    m = q * q + eps * q + 1
    c = m // 3
    n1, n2 = q ** 3 - eps, q - eps
    # Common denominator n1: t1 has coordinate x, t2 has m*y.
    L = ((1, 0), (0, m), (Q, 0), (Q * Q, 0))
    assert n1 == m * n2

    def action(w, point):
        full = apply(L, point)
        image = tuple(x % n1 for x in apply(w, full))
        assert image[1] % m == 0
        answer = (image[0], (image[1] // m) % n2)
        assert image == tuple(x % n1 for x in apply(L, answer))
        return answer

    # Derive and compare images on torus generators, not on the alleged formula.
    for point in ((1, 0), (0, 1)):
        x, y = point
        assert action(U, point) == (x % n1, (x-y) % n2)
        assert action(V, point) == ((m*y-x) % n1, y % n2)
        assert action(W, point) == (((1-m)*x+m*y) % n1, (-x+y) % n2)
        assert action(J, point) == ((-x) % n1, (-y) % n2)

    # On A use denominator 3^(a+1), with t2 exponent 3*y.
    big, small = 3 ** (a + 1), 3 ** a
    LA = ((1, 0), (0, 3), (Q, 0), (Q * Q, 0))
    target_sigma = ((-2, 3), (-1, 1))
    target_tau = ((1, 0), (1, -1))
    for w, target in ((W, target_sigma), (U, target_tau)):
        one = mul(w, LA)
        other = mul(LA, target)
        assert all((one[i][j]-other[i][j]) % big == 0
                   for i in range(4) for j in range(2))
    assert (m // 3 - 1) % small == 0
    assert m % 2 == 1 and valuation3(m) == 1

    # Root restriction uses different denominators in the two target factors.
    root_restrictions = []
    for root in roots:
        first = 2*root[0]-m*root[1]+2*Q*root[2]+2*Q*Q*root[3]
        second = -root[0]+2*root[1]-root[2]-root[3]
        assert (first % big, second % small) != (0, 0)
        # Target change: (a,b) -> (2*a, a+2*b), computed modulo 3.
        assert ((2*first) % 3, (first+2*second) % 3) == (
            (root[0]+root[2]+root[3]) % 3, root[1] % 3)
        root_restrictions.append((first % big, second % small))

    # Exact enumeration for several small finite tori, including q=2.
    fixed = None
    if n1*n2 <= 200000:
        fixed_points = [(x, y) for x in range(n1) for y in range(n2)
                        if action(W, (x, y)) == (x, y)]
        expected = [(j*n2, 0) for j in range(m)]
        assert fixed_points == expected
        assert all(action(U, point) == point for point in fixed_points)
        fixed = len(fixed_points)
    return {"q": q, "epsilon": eps, "a": a,
            "root_derived_module_action": True,
            "finite_torus_fixed_order": fixed,
            "uniform_fixed_order_formula": m,
            "root_restrictions_nontrivial": len(root_restrictions)}


def dual_module_check(a):
    big, small = 3**(a+1), 3**a
    normalize = lambda x, y: (x % big, y % small)
    sigma_a = lambda z: normalize(-2*z[0]+3*z[1], -z[0]+z[1])
    tau_a = lambda z: normalize(z[0], z[0]-z[1])
    # A character with label (b,c) has exponent b*x/big + c*y/small.
    sigma_dual = lambda z: normalize(-2*z[0]-3*z[1], z[0]+z[1])
    tau_dual = lambda z: normalize(z[0]+3*z[1], -z[1])
    points = list(product(range(big), range(small)))
    fixed_a = [z for z in points if sigma_a(z) == z]
    fixed_dual = [z for z in points if sigma_dual(z) == z]
    expected = [(small*j, 0) for j in range(3)]
    assert fixed_a == expected == fixed_dual
    assert all(tau_dual(z) == z for z in fixed_dual)
    image = {normalize(sigma_a(z)[0]-z[0], sigma_a(z)[1]-z[1]) for z in points}
    assert image == set(product(range(0, big, 3), range(small)))
    actions = set()
    for h, k, j in product(range(2), range(2), range(3)):
        images = []
        for point in ((1, 0), (0, 1)):
            for _ in range(j):
                point = sigma_a(point)
            if k:
                point = tau_a(point)
            if h:
                point = normalize(-point[0], -point[1])
            images.append(point)
        actions.add(tuple(images))
    assert len(actions) == 12
    return {"a": a, "module_order": big*small,
            "sigma_fixed": fixed_a, "sigma_fixed_dual": fixed_dual,
            "commutator_order": len(image), "faithful_complement_order": 12,
            "principal_height_zero": K_CLASSES+S3_CLASSES,
            "nonprincipal_height_zero": 3*S3_CLASSES}


def return_checks():
    records = []
    for p in (2, 5, 7, 11, 13, 17, 19):
        for f in range(1, 9):
            q = p**f
            eps = 1 if q % 3 == 1 else -1
            m = q*q+eps*q+1
            assert pow(p, 3*f, m) == eps % m
            assert pow(p, 3*f, 3) == eps % 3
            # The manuscript uses chi^a = chi composed with a^{-1}.
            assert pow(p, -3*f, m) == eps % m
            assert pow(p, -3*f, 3) == eps % 3
            records.append({"p": p, "f": f, "epsilon": eps})
    return records


def q4_nonprincipal_action():
    # c=7; F_2^3 fixes lambda of order 7, but sends j to -j modulo 3.
    c, p, power_index = 7, 2, 3
    assert pow(p, -power_index, c) == 1
    labels = list(product(range(3), range(S3_CLASSES)))
    action = lambda pair: ((pow(p, -power_index, 3)*pair[0]) % 3, pair[1])
    fixed = [pair for pair in labels if action(pair) == pair]
    assert len(fixed) == 3
    assert all(action(action(pair)) == pair for pair in labels)
    return {"q": 4, "c": c, "automorphism": "F_2^3",
            "block_lambda_order": 7, "height_zero_count": len(labels),
            "fixed_height_zero_count": len(fixed),
            "two_cycles": (len(labels)-len(fixed))//2}


if __name__ == "__main__":
    result = {
        "integral_chevalley_lifts": {
            "u_squared": 1, "v_squared": 1, "uv_order": 6,
            "complement_order": len(K), "iota_equals_n0": iota == n0,
            "complement_order_after_mod_2_reduction": 12,
            "n0_centralizes_complement": True,
            "signed_basis_table": SIGNED_TABLE,
            "signed_basis_table_verified": True,
            "cocharacter_u": U, "cocharacter_v": V,
            "cocharacter_sigma": W,
            "K_conjugacy_classes": K_CLASSES,
            "S3_conjugacy_classes": S3_CLASSES},
        "d4_roots": len(roots),
        "positive_root_restriction_counts": [
            {"image": k, "multiplicity": v}
            for k, v in sorted(restricted_positive.items())],
        "actual_torus_tests": [actual_torus_check(q)
            for q in (2,4,5,7,8,11,13,16,19,25,31,37,49,64,73)],
        "unequal_module_tests": [dual_module_check(a) for a in range(1,5)],
        "field_return_tests": return_checks(),
        "nontrivial_individual_block_stabilizer": q4_nonprincipal_action()}
    print(json.dumps(result, indent=2))
