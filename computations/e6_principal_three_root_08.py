"""An integral Weyl witness for the matching E6 Sylow-three centralizer.

This verifies one root-lattice identity uniformly for all field sizes.
The use of the identity in finite groups is proved in the accompanying note.
"""
from pathlib import Path
import json
import sympy as S

# Upstream provenance only; the proof below uses the explicit word.
source = 'https://github.com/jmichel7/gap3-jm/blob/64365f30757f5374c86511b377d7d98eb15c1e02/pkg/chevie/tbl/weyle6.g'
recorded_source_sha256 = 'b622f3f9bdd658772b9848b579f40fb6f918acaf97a13ffdcdbee1cbed8d4790'
# The fourteenth ClassInfo word, named E_6(a_1), in the pinned source.
word = [1, 3, 4, 3, 2, 4, 5, 6]
C = 2*S.eye(6)
for i, j in [(1,3), (3,4), (4,5), (5,6), (2,4)]:
    C[i-1,j-1] = C[j-1,i-1] = -1
reflections = []
for i in range(6):
    r = S.eye(6)
    for j in range(6):
        r[i,j] -= C[i,j]
    assert r*r == S.eye(6) and r.T*C*r == C
    reflections.append(r)
w = S.eye(6)
for i in word:
    w *= reflections[i-1]
x = S.symbols('x')
assert C.det() == 3
assert w**9 == S.eye(6) and w**3 != S.eye(6)
assert w.charpoly(x).as_expr() == S.cyclotomic_poly(9,x)
assert (S.eye(6)-w).det() == 3
# The weight-lattice matrix is conjugate by the Cartan matrix.
weight_w = C*w*C.inv()
assert all(z.q == 1 for z in weight_w)
assert (S.eye(6)-weight_w).det() == 3
result = {
    'source': source,
    'recorded_source_sha256': recorded_source_sha256,
    'class': 'E_6(a_1)', 'word': word, 'order': 9,
    'root_matrix': w.tolist(), 'weight_matrix': weight_w.tolist(),
    'characteristic_polynomial': str(w.charpoly(x).as_expr()),
    'determinant_one_minus_w': 3,
    'scope': 'All matching E6^epsilon(q), 3 divides q-epsilon, not characteristic3',
    'mathematical_use': 'T^w is a diagonalizable group of order3; hence C_G(P)=Z(G) once A=T_3 and an actual Sylow P over the Weyl Sylow are chosen.'
}
Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2, default=int)+'\n')
print(json.dumps({'checked': True, 'order': 9, 'charpoly': result['characteristic_polynomial'], 'det': 3}))
