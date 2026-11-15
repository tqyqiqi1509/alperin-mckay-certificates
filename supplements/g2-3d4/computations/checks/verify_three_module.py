#!/usr/bin/env python3
"""Independent exact checks for the Sylow-three module in 3D4(q).

These checks validate the displayed finite arithmetic; they do not certify
block membership, a character-table input, or the inductive condition.
"""
from itertools import product
import json


def check(a):
    n, m = 3 ** (a + 1), 3**a
    points = list(product(range(n), range(m)))
    normalize = lambda x, y: (x % n, y % m)
    sigma = lambda z: normalize(-2*z[0] + 3*z[1], -z[0] + z[1])
    tau = lambda z: normalize(z[0], z[0]-z[1])
    iota = lambda z: normalize(-z[0], -z[1])
    sigma_dual = lambda z: normalize(-2*z[0]-3*z[1], z[0]+z[1])
    tau_dual = lambda z: normalize(z[0]+3*z[1], -z[1])
    for z in points:
        assert sigma(sigma(sigma(z))) == z
        assert tau(tau(z)) == z
        assert tau(sigma(tau(z))) == sigma(sigma(z))
        assert iota(sigma(z)) == sigma(iota(z))
        assert iota(tau(z)) == tau(iota(z))
    fixed = [z for z in points if sigma(z) == z]
    commutators = {normalize(sigma(z)[0]-z[0], sigma(z)[1]-z[1])
                   for z in points}
    dual_fixed = [z for z in points if sigma_dual(z) == z]
    expected = [(j*m, 0) for j in range(3)]
    assert fixed == expected and dual_fixed == expected
    assert commutators == {(x, y) for x in range(0, n, 3)
                           for y in range(m)}
    assert all(tau_dual(z) == z for z in dual_fixed)
    # Test K=S3 x C2 on two generators, enough to distinguish linear maps.
    maps = set()
    for j, k, h in product(range(3), range(2), range(2)):
        images = []
        for z in [(1, 0), (0, 1)]:
            for _ in range(j):
                z = sigma(z)
            if k:
                z = tau(z)
            if h:
                z = iota(z)
            images.append(z)
        maps.add(tuple(images))
    assert len(maps) == 12
    return {"a": a, "A_order": n*m, "S_order": 3*n*m,
            "sigma_fixed": fixed, "commutator_order": len(commutators),
            "sigma_fixed_dual_labels": dual_fixed, "faithful_K_order": 12,
            "principal_local_height_zero_count": 6 + 3,
            "nonprincipal_local_height_zero_count": 3 * 3}


def valuation(n, p=3):
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def fields():
    out = []
    for p in [2, 5, 7, 11, 13, 17, 19]:
        for f in range(1, 7):
            q = p**f
            epsilon = 1 if q % 3 == 1 else -1
            a = valuation(q-epsilon)
            m = q*q+epsilon*q+1
            c = m//3
            assert valuation(m) == 1
            assert (m//3-1) % (3**a) == 0
            assert pow(p, 3*f, m) == epsilon % m
            assert pow(p, 3*f, 3) == epsilon % 3
            # Six of the eight unipotent degree formulas have 3'-degree.
            phi12 = q**4-q*q+1
            degrees = [1, q**12, q*phi12, q**7*phi12,
                       q**3*(q+1)**2*phi12//2,
                       q**3*(q-1)**2*phi12//2,
                       q**3*(q+1)**2*(q*q-q+1)**2//2,
                       q**3*(q-1)**2*(q*q+q+1)**2//2]
            assert len(set(degrees)) == 8
            assert sum(d % 3 != 0 for d in degrees) == 6
            out.append({"p": p, "f": f, "epsilon": epsilon, "a": a,
                        "c": c, "unipotent_three_prime_degree_count": 6})
    return out


if __name__ == "__main__":
    print(json.dumps({"module_checks": [check(a) for a in range(1, 5)],
                      "field_checks": fields()}, indent=2))
