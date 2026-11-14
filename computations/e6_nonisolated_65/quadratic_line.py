"""Exact same-F central quadratic point and its actual signed S3 orbit.

Input: the E6 Cartan matrix, x6's fundamental-coweight coefficients,
and the independently reconstructed perpendicular A2^3 carrier roots.
There are no finite-field samples and no GAP table identifiers.
"""
from itertools import permutations
from pathlib import Path
import json
import sympy as s

C=s.Matrix([[2,0,-1,0,0,0],[0,2,0,-1,0,0],
            [-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],
            [0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
x6=C.inv()*s.Matrix([0,1,1,0,1,0])/6
e4=s.eye(6)[:,3]
diff=3*x6-e4/2
assert all(v.q==1 for v in diff)
H=s.Matrix.hstack(*map(s.Matrix,[(0,0,1,1,1,1),
                               (0,1,0,1,1,1),
                               (1,0,1,1,1,0),
                               (1,1,1,1,0,0)]))
assert H.T*C*e4==s.zeros(4,1)
assert H.T*C*(3*x6)==s.ones(4,1)

section={}
for p in permutations(range(3)):
    P=s.zeros(3)
    for j,i in enumerate(p): P[i,j]=1
    R=P.det()*P
    assert R.det()==1 and R.T*R==s.eye(3)
    section[p]=R
for p,R in section.items():
    for q,S in section.items():
        pq=tuple(p[q[j]] for j in range(3))
        assert R*S==section[pq]

z=s.diag(-1,-1,1)
orbit=set()
stabilizer=[]
for p,R in section.items():
    for tau_power in (0,1):
        z_tau=z if tau_power==0 else z.inv().T
        conjugate=R*z_tau*R.inv()
        orbit.add(tuple(conjugate))
        if conjugate==z: stabilizer.append([list(p),tau_power])
assert len(orbit)==3 and len(stabilizer)==4
assert all(s.Matrix(3,3,zv).inv().T==s.Matrix(3,3,zv) for zv in orbit)

out={
 "method":"exact characteristic-zero root and signed matrix identities",
 "three_x6":[str(v) for v in 3*x6],
 "three_x6_minus_half_alpha4_coroot":[str(v) for v in diff],
 "H_pairing_alpha4":[str(v) for v in H.T*C*e4],
 "H_pairing_three_x6":[str(v) for v in H.T*C*(3*x6)],
 "signed_section_order":len(section),
 "all_signed_matrices_determinant_one":True,
 "all_signed_matrices_orthogonal":True,
 "all_signed_section_multiplication_relations":True,
 "two_torsion_orbit":[list(map(int,zv)) for zv in sorted(orbit)],
 "orbit_size":len(orbit),
 "Q12_stabilizer_order":len(stabilizer),
 "Q12_stabilizer_permutation_and_tau":stabilizer,
 "tau_fixes_each_two_torsion_point":True,
 "source":"z2'=s6^3=h_alpha4(-1), no rational transporter required"
}
target=Path(__file__).with_suffix('.json')
target.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
