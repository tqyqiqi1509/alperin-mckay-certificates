"""Independent uniform witnesses for the matching E6 D4 source.

No finite defining field is sampled.  The E6 root computation is integral;
the 3-residue computation is the fixed coefficient module, and the complex
matrices give a representation of the actual GL2(3) inertia of Q8.
"""
from collections import Counter
from itertools import product, combinations
from pathlib import Path
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

C = s.Matrix([[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],
              [0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
simple = [s.eye(6)[:,i] for i in range(6)]
refs = []
for i in range(6):
    refs.append(s.eye(6)-simple[i]*C[i,:])
roots = {tuple(v) for v in simple}
todo = list(roots)
while todo:
    v = s.Matrix(todo.pop())
    for r in refs:
        u = tuple(r*v)
        if u not in roots:
            roots.add(u); todo.append(u)
assert len(roots) == 72

B = s.Matrix([[2,0],[1,1],[2,1],[2,2],[1,2],[0,2]])
source = [1,2,3,4]
assert C[source,:]*B == s.zeros(4,2)
assert s.gcd_list([int(B[list(I),:].det()) for I in combinations(range(6),2)]) == 1
counts = Counter(tuple((s.Matrix(v).T*C*B)) for v in roots)
assert counts == Counter({(0,0):24,(2,-1):8,(-2,1):8,(-1,2):8,
                          (1,-2):8,(1,1):8,(-1,-1):8})
# Root columns in the simply connected character lattice of E6.
source_columns = C[:,source]
minor_gcd = abs(s.gcd_list([int(source_columns[list(I),:].det())
                          for I in combinations(range(6),4)]))
assert minor_gcd == 1

# An independent sparse minuscule model checks actual Tits lifts, not just W.
weights={(1,0,0,0,0,0)}; frontier=list(weights)
while frontier:
    wt=frontier.pop()
    for i in range(6):
        nxt=tuple(wt[j]-wt[i]*int(C[j,i]) for j in range(6))
        if nxt not in weights: weights.add(nxt);frontier.append(nxt)
assert len(weights)==27
ident={wt:(wt,1) for wt in weights}
N=[]
for i in range(6):
    N.append({wt:(tuple(wt[j]-wt[i]*int(C[j,i]) for j in range(6)),
                    -1 if wt[i]==1 else 1) for wt in weights})
def compose(a,b):
    return {x:(a[y][0],sgn*a[y][1]) for x,(y,sgn) in b.items()}
def tits_word(word):
    out=ident
    for i in word: out=compose(out,N[i-1])
    return out
word1=[1,3,4,2,5,4,3,1];word2=[6,5,4,2,3,4,5,6]
q1=tits_word(word1);q2=tits_word(word2);sig=compose(q1,q2)
assert compose(q1,q1)==compose(q2,q2)==ident
assert compose(compose(sig,sig),sig)==ident
for q,perm in [(q1,{2:5,3:3,4:4,5:2}),(q2,{2:3,3:2,4:4,5:5})]:
    for j,k in perm.items():
        for sign in [-1,1]:
            for wt in weights:
                target,sgn=q[wt]
                assert (wt[j-1]==sign)==(target[k-1]==sign)
                if wt[j-1]==sign:
                    nxt=tuple(wt[l]-sign*int(C[l,j-1]) for l in range(6))
                    expected=tuple(target[l]-sign*int(C[l,k-1]) for l in range(6))
                    assert q[nxt]==(expected,sgn)

def weyl_word(word):
    out=s.eye(6)
    for i in word: out=out*refs[i-1]
    return out
qw1=weyl_word(word1);qw2=weyl_word(word2)
gamma=s.zeros(6)
for i,j in enumerate([5,1,4,3,2,0]):gamma[j,i]=1
u=qw1*qw2*qw1
assert u**2 == s.eye(6)
assert u**-1*gamma*B == -B

w = s.Matrix([[0,-1],[1,-1]])
r = s.Matrix([[-1,1],[0,1]])
assert qw1*B == B*r
assert qw1*qw2*B == B*w
M = w-s.eye(2)
rad = s.Matrix([[0,1],[1,0]])
assert w**3 == r**2 == s.eye(2)
assert r*w*r == w**-1
assert M**2 == -3*w and M*r == rad*M
assert abs(M.det()) == 3
SNF = smith_normal_form(M,domain=ZZ)
assert sorted(abs(int(SNF[i,i])) for i in range(2)) == [1,3]

# Formal, valuation-independent self-twist equations: solutions are
# t0+j*z and their negatives, with z=n*(1,2)/3 and t0=n*(0,-1)/3.
n = s.symbols('n',integer=True,positive=True)
z = s.Matrix([n/3,2*n/3]); t0 = s.Matrix([0,-n/3])
assert M*t0 == z
assert w*z-z == s.Matrix([-n,-n])
assert r*z-z == s.zeros(2,1)
assert (r+s.eye(2))*z == 2*z
assert s.expand(3*(1+1+4)+2*(n-3)*9+(n*n-3*n+6)/6*36-6*n*n) == 0

# The GL2(3) representation is checked over Q(i,sqrt(2)).
I = s.eye(2)
U = s.diag(s.I,-s.I)
V = s.Matrix([[0,1],[-1,0]])
W = U*V
AA = (-I+U+V-W)/2
BB = (V+W)/(s.I*s.sqrt(2))
normal = lambda A: A.applyfunc(lambda x:s.expand(s.simplify(x)))
assert U**2 == V**2 == W**2 == -I
assert normal(AA**3) == I and normal(BB**2) == I
assert normal(BB*AA*BB) == normal(AA**-1)
assert normal(AA*U*AA**-1) == V
assert normal(AA*V*AA**-1) == -W
assert normal(BB*U*BB**-1) == -U
assert normal(BB*V*BB**-1) == W

def key(A):
    return tuple(normal(A))
q8 = [I,-I,U,-U,V,-V,W,-W]
complement = [I,AA,normal(AA**2),BB,normal(AA*BB),normal(AA**2*BB)]
all_g = {key(x*y): normal(x*y) for x in q8 for y in complement}
assert len(all_g) == 48
for g in all_g.values():
    for h in [U,V,AA,BB]:
        assert key(g*h) in all_g

# Explicit conjugation on the three nontrivial quaternion lines.
line_basis = [U,V,W]
def conjugation(A):
    cols=[]
    for x in line_basis:
        v=normal(A*x*A**-1)
        coeff=[normal(v*y).trace()/(-2) for y in line_basis]
        assert normal(sum((a*y for a,y in zip(coeff,line_basis)),s.zeros(2))) == v
        cols.append(s.Matrix(coeff))
    return s.Matrix.hstack(*cols)
AU,AV,Aa,Ab = map(conjugation,[U,V,AA,BB])
assert AU == s.diag(1,-1,-1) and AV == s.diag(-1,1,-1)
assert Aa == s.Matrix([[0,0,-1],[1,0,0],[0,-1,0]])
assert Ab == s.Matrix([[-1,0,0],[0,0,1],[0,1,0]])

# Fixed characteristic-3 coefficient module: Q8 is regular on V^* - 0.
u = ((0,1),(2,0)); v = ((1,1),(1,2))
def mul(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2))%3
                       for j in range(2)) for i in range(2))
one=((1,0),(0,1)); Q={one}; frontier=[one]
while frontier:
    a=frontier.pop()
    for b in (u,v):
        c=mul(a,b)
        if c not in Q: Q.add(c);frontier.append(c)
assert len(Q)==8
vectors=[x for x in product(range(3),repeat=2) if x!=(0,0)]
for x in vectors:
    orbit={tuple(sum(x[k]*a[k][j] for k in range(2))%3 for j in range(2)) for a in Q}
    assert len(orbit)==8
GL=[((a,b),(c,d)) for a,b,c,d in product(range(3),repeat=4) if (a*d-b*c)%3]
assert len(GL)==48
for g in GL:
    gi=next(h for h in GL if mul(g,h)==one)
    assert {mul(mul(g,x),gi) for x in Q} == Q

out={
 'passed':True,'scope':'Integral E6 and fixed universal coefficient groups; no defining-field sample',
 'root_count':72,'root_restriction_multiplicities':{str(k):v for k,v in counts.items()},
 'sparse_minuscule_dimension':27,'Tits_actual_S3_relations':True,
 'positive_negative_D4_pinning_signs':'+1','corrected_graph_on_center':'-I',
 'D4_center_primitive_minor_gcd':int(minor_gcd),
 'center_cocharacters':B.tolist(),'w':w.tolist(),'reflection':r.tolist(),
 'diagonal_isogeny_M':M.tolist(),'M_Smith_invariant_factors':[1,3],
 'formal_self_twist_identity':True,'ordinary_sum_of_squares_identity':True,
 'quaternion_GL2_complex_representation_order':len(all_g),
 'traceless_A':Aa.tolist(),'traceless_B':Ab.tolist(),
 'Q8_on_nonzero_covectors':'regular','normalizer_Q8_in_GL2_3_order':48,
 'base_coefficient_degrees':[1,1,1,1,2,8],
 'upper_SL2_coefficient_degrees':[1,1,1,2,2,2,8,8,8,3]
}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2,default=int)+'\n')
print(json.dumps({'passed':True,'roots':72,'quaternion_extension_order':48,'normalizer_order':48}))
