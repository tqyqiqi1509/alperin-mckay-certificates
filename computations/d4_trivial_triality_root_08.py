"""Exact D4 root-lattice certificate for the trivial-centre triality route.

The variables q below are residue classes modulo the order of the torus
point.  No finite-field examples or ordinary character tables are sampled.
"""
from itertools import permutations, product
from pathlib import Path
import json
import sympy as sp

OUT = Path(__file__).with_suffix('.json')
E = sp.eye(4)
B = sp.Matrix.hstack(E[:, 0]-E[:, 1], E[:, 1]-E[:, 2],
                     E[:, 2]-E[:, 3], E[:, 2]+E[:, 3])
roots = []
for i in range(4):
    for j in range(i+1, 4):
        for a, b in product((-1, 1), repeat=2):
            roots.append(a*E[:, i]+b*E[:, j])

W = []
for p in permutations(range(4)):
    for signs in product((-1, 1), repeat=4):
        if sp.prod(signs) == 1:
            W.append(sp.Matrix.hstack(*(signs[i]*E[:, p[i]] for i in range(4))))
assert len(W) == 192

def in_Q(v):
    return all(x.q == 1 for x in v) and sum(v) % 2 == 0

def in_P(v):
    return all(x.q == 1 for x in v) or all(x.q == 2 for x in v)

def colour(v):
    assert in_P(v)
    if all(x.q == 1 for x in v):
        return '0' if sum(v) % 2 == 0 else 'vector'
    return 'spin_even' if sum(v) % 2 == 0 else 'spin_odd'

Gamma = []
for p in permutations((0, 2, 3)):
    b = [None]*4
    b[1] = B[:, 1]
    for i, j in zip((0, 2, 3), p):
        b[i] = B[:, j]
    g = sp.Matrix.hstack(*b)*B.inv()
    assert g.T*g == E
    Gamma.append(g)
assert len({tuple(g) for g in Gamma}) == 6
Aut = [w*g for w in W for g in Gamma]
assert len({tuple(g) for g in Aut}) == 1152
assert {x for g in Aut for x in g} <= {0, 1, -1, sp.Rational(1, 2), sp.Rational(-1, 2)}

points = {
    'isolated_D2_D2': sp.Matrix([sp.Rational(1, 2), sp.Rational(1, 2), 0, 0]),
    'strict_A1_T3': sp.Matrix([sp.Rational(1, 2), sp.Rational(1, 4), sp.Rational(1, 4), 0]),
}
result = {'W_order': len(W), 'root_automorphism_order': len(Aut), 'points': {}}
for name, x in points.items():
    rc = [r for r in roots if (r.dot(x)).q == 1]
    wc = [w for w in W if in_Q(w*x-x)]
    wa = [w for w in W if in_P(w*x-x)]
    aug = [w for w in Aut if in_P(w*x-x)]
    translations = {c:sum(colour(w*x-x) == c for w in wa)
                    for c in ('0', 'vector', 'spin_even', 'spin_odd')}
    assert len(wa) == 4*len(wc)
    assert all(n == len(wc) for n in translations.values())
    assert len(aug) == 6*len(wa)
    root_basis = sp.Matrix.hstack(*rc).columnspace()
    centre = sp.Matrix.hstack(*root_basis).T.nullspace()
    rational = []
    for q in (1, 3):
        candidates = [w for w in W if in_P(w*(q*x)-x)]
        zero = [w for w in candidates if in_Q(w*(q*x)-x)]
        assert len(zero) == len(wc)
        f = E if name == 'isolated_D2_D2' or q == 1 else -E
        assert in_Q(f*(q*x)-x)
        if centre:
            cb = sp.Matrix.hstack(*centre)
            fm = cb.gauss_jordan_solve(f*cb)[0]
            centre_polynomial = str(sp.factor(fm.charpoly().as_expr()))
        else:
            centre_polynomial = '1'
        # All zero-return choices differ by the actual connected Weyl group.
        assert all(any(w == v*f for v in wc) for w in zero)
        rational.append({'q_mod_4': q, 'zero_return_transporter_count': len(zero),
                         'representative_F': f.tolist(),
                         'central_torus_F_polynomial': centre_polynomial})
    result['points'][name] = {
        'point': [str(a) for a in x], 'centralizer_roots': len(rc),
        'connected_W_order': len(wc), 'adjoint_stabilizer_order': len(wa),
        'full_triality_stabilizer_order': len(aug),
        'central_translation_fibres': translations,
        'trivial_colour_rational_forms': rational,
    }

# Full outer graph exclusion for the toral defect: a root automorphism
# fixing the full torus lattice pointwise is identity.  For the optional
# A1-source interface, its perpendicular hyperplane has only the inner
# A1 reflection as a nontrivial pointwise stabilizer.
alpha = B[:, 1]
perp = sp.Matrix([list(alpha)]).nullspace()
fix = [a for a in Aut if all(a*v == v for v in perp)]
assert len(fix) == 2
assert all(any(a == w for w in W) for a in fix)
result['A1_perpendicular_pointwise_stabilizer'] = {'order': 2, 'all_inner': True}

# The three pairs of mutually perpendicular end roots have distinct
# W(D4)-orbits. Triality permutes all three transitively.
pairs = [(B[:, i], B[:, j]) for i, j in ((0, 2), (0, 3), (2, 3))]
def root_set(a, b):
    return {tuple(a), tuple(-a), tuple(b), tuple(-b)}
for i in range(3):
    for j in range(i):
        assert not any(root_set(w*pairs[i][0], w*pairs[i][1]) ==
                       root_set(*pairs[j]) for w in W)
assert all(any(root_set(g*pairs[0][0], g*pairs[0][1]) == root_set(*p)
               for g in Gamma) for p in pairs)
result['D2_Levi_classes'] = {'inner_orbits': 3, 'triality_transitive': True}

OUT.write_text(json.dumps(result, indent=2, default=int)+'\n')
print(json.dumps(result, indent=2, default=int))
