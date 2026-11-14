"""Integral D4 row-33 carrier and its whole-q degree, rebuilt at node 58."""
from pathlib import Path
from itertools import combinations, product
import json
import sympy as S

out=Path(__file__).parent
I=S.eye(8); J=S.zeros(8)
for i in range(4): J[i,4+i]=J[4+i,i]=1
roots=[]; root_matrices=[]
for i,j in combinations(range(4),2):
    for a,b in product((-1,1), repeat=2):
        r=S.zeros(4,1);r[i]=a;r[j]=b
        X=S.zeros(8)
        if a==1 and b==-1:X[i,j]=1;X[4+j,4+i]=-1
        elif a==-1 and b==1:X[j,i]=1;X[4+i,4+j]=-1
        elif a==1 and b==1:X[i,4+j]=1;X[j,4+i]=-1
        else:X[4+j,i]=1;X[4+i,j]=-1
        assert X.T*J+J*X==S.zeros(8)
        roots.append(r);root_matrices.append(X)
def X(r):return root_matrices[next(i for i,t in enumerate(roots) if t==S.Matrix(r))]
e=X([1,-1,0,0])+X([0,1,-1,0])+X([0,1,1,0])
chain=[I[:,4]]
for k in range(4):chain.append(e*chain[-1])
assert e*chain[-1]==S.zeros(8,1)
others=[I[:,2]-I[:,6], I[:,3], I[:,7]]
chain_matrix=S.Matrix.hstack(*chain,*others)
assert abs(chain_matrix.det())==8
assert all(e*v==S.zeros(8,1) for v in others)
eta=S.Matrix([4,2,0,0])
grades=[int((eta.T*r)[0]) for r in roots]
assert all(g%2==0 for g in grades)
hs=[]
for i in range(4):
    h=S.zeros(8);h[i,i]=1;h[4+i,4+i]=-1;hs.append(h)
lie=root_matrices+hs
trace_form=S.Matrix([[(a*b).trace() for b in lie] for a in lie])
assert abs(trace_form.det())==2**28
g0=[a for a,g in zip(root_matrices,grades) if g==0]+hs
g2=[a for a,g in zip(root_matrices,grades) if g==2]
flat=lambda a:S.Matrix(list(a))
B2=S.Matrix.hstack(*map(flat,g2))
ad=S.Matrix.hstack(*[B2.gauss_jordan_solve(flat(a*e-e*a))[0] for a in g0])
mins=[ad[:,cc].det() for cc in combinations(range(ad.cols),ad.rows)]
assert 2 in mins or -2 in mins
graph=I.copy();graph[:,3],graph[:,7]=I[:,7],I[:,3]
assert graph*graph==I and graph.T*J*graph==J
assert graph*e*graph==e
u=(I+e)*(I-e).inv()
assert u.T*J*u==J and (u-I)*(u+I).inv()==e
assert graph*u*graph==u
assert all(S.denom(t)==1 for t in u)
d=json.loads((out/'independent_geometry.json').read_text())['cases']['D5A3_quartic33']
K=S.Matrix.hstack(*[S.Matrix(r) for r in d['K_base']]);F=S.Matrix(d['F'])
Z=S.Matrix.hstack(*K.T.nullspace());FZ=(Z.T*Z).inv()*Z.T*F*Z
x,q=S.symbols('x q')
center=S.factor(FZ.charpoly(x).as_expr());assert center==(x+1)*(x*x+1)**2
phi=lambda k:S.cyclotomic_poly(k,q)
lord=phi(4)**2*(q*q-1)*(q**4-1)**2*(q**6-1)
cord=phi(2)*phi(4)**2*(q*q-1)*(q**3+1)*(q**4-1)
degree=S.cancel(lord/cord*q*q*phi(4))
assert S.expand(degree-q*q*phi(1)**2*phi(3)*phi(4)**2)==0
ennola=S.expand(degree.subs(q,-q))
assert S.expand(ennola-q*q*phi(2)**2*phi(6)*phi(4)**2)==0
assert S.Poly(degree,q).coeff_monomial(q*q)==S.Poly(ennola,q).coeff_monomial(q*q)==1
record={
 'chain_det':str(chain_matrix.det()),'partition':[5,1,1,1],
 'even_grades':sorted(set(grades)), 'grade0_to_grade2_rank':ad.rank(),
 'prime_to_all_odd_minor':2,'trace_form_det':str(trace_form.det()),
 'integral_Cayley':True,'actual_graph_fixes_e_and_u':True,
 'source_centralizer_torus':str(center),'source_degree':str(S.factor(degree)),
 'source_Ennola_degree':str(S.factor(ennola)),'lowest_coefficient':1,
 'chain':chain_matrix.tolist(), 'ad_grade0_grade2':ad.tolist(),
 'e':e.tolist(),'u':u.tolist()}
(out/'even_gggr.json').write_text(json.dumps(record,indent=2,default=str)+'\n')
print(json.dumps({k:v for k,v in record.items() if k not in ['chain','ad_grade0_grade2','e','u']},indent=2))
