"""Exact, symbolic D4 identities used in the all-q triality construction.

No finite field q or valuation a is sampled.  The assertions concern integral
root lattices, symbolic order polynomials, and literal finite carrier matrices.
"""
from itertools import combinations, product
from pathlib import Path
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
C = s.Matrix([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]])
I = s.eye(4)
alpha = s.Matrix([0,1,0,0])
beta = s.Matrix([1,1,1,1])
reflection = lambda v: I-v*(v.T*C)
tau = reflection(alpha)
sigma = reflection(alpha)*reflection(beta)
rho = s.Matrix([[0,0,0,1],[0,1,0,0],[1,0,0,0],[0,0,1,0]])
Q = s.Matrix([[1,0,0,-1],[0,0,1,-1]])
R = s.Matrix([[0,-1],[1,-1]])
assert sigma**3 == I and rho**3 == I
assert sigma*rho == rho*sigma
assert tau*sigma*tau == sigma**-1
assert Q*(sigma-I) == s.zeros(2,4)
assert Q*rho == R*Q
assert Q*tau == Q
assert smith_normal_form(sigma-I, domain=s.ZZ) == s.diag(1,1,0,0)
assert C.det() == 4  # the simply connected/adjoint comparison is 3-prime

q = s.symbols('q')
M = (q*rho-I).row_join(sigma-I)
minor_gcd = s.Integer(0)
for cols in combinations(range(8), 4):
    minor_gcd = s.gcd(minor_gcd, s.expand(M[:,cols].det()))
assert s.factor(minor_gcd) == q*q+q+1
assert s.factor((q*R-s.eye(2)).det()) == q*q+q+1

# All D4 roots, expressed in the simple-root basis.
simple = s.Matrix([[1,0,0,0],[-1,1,0,0],[0,-1,1,1],[0,0,-1,1]])
roots = []
for i,j in combinations(range(4),2):
    for x,y in product([-1,1], repeat=2):
        v=s.zeros(4,1);v[i]=x;v[j]=y
        z=simple.inv()*v
        assert all(t.q==1 for t in z)
        roots.append(z)
restricted={}
for z in roots:
    k=tuple(int(x) for x in Q*z)
    restricted[k]=restricted.get(k,0)+1
assert restricted == {(0,0):6,(1,0):3,(-1,0):3,(0,1):3,(0,-1):3,(1,1):3,(-1,-1):3}

# Exact Weyl normalizer of the A2 subsystem: the quotient is C2.
rootkey=lambda z:tuple(int(x) for x in z)
a2roots={rootkey(z) for z in roots if Q*z==s.zeros(2,1)}
generators=[reflection(s.eye(4)[:,i]) for i in range(4)]
matrixkey=lambda z:tuple(int(x) for x in z)
weyl={matrixkey(I):I}; todo_weyl=[I]
while todo_weyl:
    g=todo_weyl.pop()
    for t in generators:
        gt=g*t; key=matrixkey(gt)
        if key not in weyl:
            weyl[key]=gt;todo_weyl.append(gt)
assert len(weyl)==192
stabilizer=[g for g in weyl.values()
            if {rootkey(g*s.Matrix(z)) for z in a2roots}==a2roots]
assert len(stabilizer)==12
subsystems={tuple(sorted(rootkey(g*s.Matrix(z)) for z in a2roots))
            for g in weyl.values()}
assert len(subsystems)==16

# At every q == epsilon mod 3, this is the mod-3 quotient by F-1.
Q3=s.Matrix([[1,0,1,1],[0,1,0,0]])
assert all(any(int(x)%3 for x in Q3*z) for z in roots)
tau3=s.Matrix([[1,0],[1,2]])
sigma3=s.Matrix([[1,0],[2,1]])
red=lambda A:A.applyfunc(lambda x:int(x)%3)
assert red(Q3*tau)==red(tau3*Q3)
assert red(Q3*sigma)==red(sigma3*Q3)
flat=lambda A:tuple(int(x) for x in A)
carrier={flat(s.eye(2)):s.eye(2)}
todo=[s.eye(2)]
while todo:
    A=todo.pop()
    for B in [sigma3,tau3,-s.eye(2)]:
        AB=red(A*B);k=flat(AB)
        if k not in carrier:carrier[k]=AB;todo.append(AB)
assert len(carrier)==12

# Literal all-a action on Z/(3^(a+1)) x Z/(3^a).
Sm=s.Matrix([[-2,3],[-1,1]])
Tm=s.Matrix([[1,0],[1,-1]])
assert Sm**3==s.eye(2)
assert Tm**2==s.eye(2)
assert Tm*Sm*Tm==Sm**-1
assert Sm-s.eye(2)==s.Matrix([[-3,3],[-1,0]])
# Exact identity: c-1=(q-epsilon)(q+2 epsilon)/3.
for eps in [-1,1]:
    c=(q*q+eps*q+1)/3
    assert s.expand(c-1-(q-eps)*(q+2*eps)/3)==0

Phi=lambda n:s.cyclotomic_poly(n,q)
degrees={
 '1':s.Integer(1),'St':q**12,
 "phi13prime":q*Phi(12),"phi13doubleprime":q**7*Phi(12),
 'phi22':q**3*Phi(2)**2*Phi(12)/2,
 '3D4plus':q**3*Phi(1)**2*Phi(12)/2,
 'phi21':q**3*Phi(2)**2*Phi(6)**2/2,
 '3D4minus':q**3*Phi(1)**2*Phi(3)**2/2,
}
# Independent primary-data route: corrected CHEVIE weyl3d4.g has
# C2 family in positions [7,4,3,8] and these sparse fake degrees.
familyfake=s.Matrix([q**3+q**9,-q**4+2*q**6-q**8,q**4+q**8,0])
fourier=s.Matrix([[1,1,1,1],[1,1,-1,-1],
                  [1,-1,1,-1],[1,-1,-1,1]])/2
for got,key in zip(fourier*familyfake,
                   ['phi21','3D4plus','phi22','3D4minus']):
    assert s.expand(got-degrees[key])==0
assert s.factor(degrees['phi21']-degrees['3D4minus'])==2*q**6
assert s.factor(degrees['phi22']-degrees['phi21'])==q**4*(q-1)**2*(q+1)**2
# Unequal-parameter D12 one-dimensional Hecke labels, symbolically.
counts=[(0,0)]
for n in range(1,6):counts.extend([((n+1)//2,n//2),(n//2,(n+1)//2)])
counts.append((3,3))
P=sum(q**(x+3*y) for x,y in counts)
assert s.factor(P/sum(q**(-x+3*y) for x,y in counts))==degrees['phi13prime']
assert s.factor(P/sum(q**(x-3*y) for x,y in counts))==degrees['phi13doubleprime']

record={
 'scope':'all q, all a; no q samples',
 'cartan':C.tolist(),'triality':rho.tolist(),'sigma':sigma.tolist(),
 'quotient_map':Q.tolist(),'quotient_triality':R.tolist(),
 'smith_sigma_minus_identity':[1,1,0,0],
 'symbolic_full_minor_gcd':str(s.factor(minor_gcd)),
 'restricted_roots':{str(k):v for k,v in restricted.items()},
 'weyl_order':len(weyl),'a2_subsystems':len(subsystems),
 'a2_weyl_normalizer_order':len(stabilizer),'a2_relative_normalizer_order':2,
 'roots_checked':len(roots),'all_roots_nonzero_on_mod3_torus':True,
 'faithful_mod3_carrier_order':len(carrier),
 'literal_sylow_module_sigma':Sm.tolist(),'literal_sylow_module_tau':Tm.tolist(),
 'generic_unipotent_degrees':{k:str(s.factor(v)) for k,v in degrees.items()},
 'degrees_match_corrected_chevie_fourier_data':True,
 'all_assertions_passed':True,
}
Path(__file__).with_suffix('.json').write_text(json.dumps(record,indent=2,default=str)+'\n')
print(json.dumps({'all_assertions_passed':True,'roots':24,'carrier_order':12,'symbolic':True}))
