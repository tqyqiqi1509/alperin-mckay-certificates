"""Independent Euclidean E8 witness, with the D4+D4 integral gluing retained.

This is a uniform root-lattice computation, with no finite field sampled.
Roots are stored after multiplication by two to avoid floating point.
"""
from itertools import combinations, product, permutations
from collections import deque, Counter
from pathlib import Path
import json
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

def key(m): return tuple(map(int, m))
def mat(k, n=4): return S.Matrix(n, n, k)
def mult(a,b,n=4):
    return tuple(sum(a[n*i+k]*b[n*k+j] for k in range(n)) for i in range(n) for j in range(n))
def act(a,b,n=4): return tuple(sum(a[n*i+j]*b[j] for j in range(n)) for i in range(n))
def closure(gens):
    unit=key(S.eye(4)); seen={unit}; queue=deque([unit])
    while queue:
        a=queue.popleft()
        for g in gens:
            b=mult(a,g)
            if b not in seen:seen.add(b);queue.append(b)
    return sorted(seen)
def lattice_e8(x):
    twice=[2*S.Rational(t) for t in x]
    if not all(t.q==1 for t in twice):return False
    parity={int(t)%2 for t in twice}
    return len(parity)==1 and sum(x)%2==0

roots=set()
for i,j in combinations(range(8),2):
    for a,b in product((-2,2),repeat=2):
        r=[0]*8;r[i]=a;r[j]=b;roots.add(tuple(r))
for signs in product((-1,1),repeat=8):
    if sum(t<0 for t in signs)%2==0:roots.add(signs)
assert len(roots)==240

# D4 is on coordinates 1..4; its central torus has D4 root lattice on 5..8.
R=S.Matrix([[1,0,0,0],[-1,1,0,0],[0,-1,1,1],[0,0,-1,1]])
Q=R.T*R
assert R.det()==2 and Q.det()==4
base=S.zeros(8,4);base[4:8,0:4]=R
restr=[]
for r in roots:
    row=(S.Matrix(1,8,r)*base)/2
    assert all(t.q==1 for t in row)
    restr.append(tuple(map(int,row)))
assert sum(not any(r) for r in restr)==24
assert sum(all(t%3==0 for t in r) for r in restr)==24
assert sum(all(t%5==0 for t in r) for r in restr)==24

gens=[]; lifts=[]
for i in range(4):
    w=S.eye(4)
    for j in range(4):w[i,j]-=Q[i,j]
    gens.append(key(w))
    lift=S.eye(8);lift[4:8,4:8]=R*w*R.inv();lifts.append(lift)
for perm in ((2,1,0,3),(2,1,3,0)):
    p=S.zeros(4)
    for j in range(4):p[perm[j],j]=1
    gens.append(key(p))
    trial=R*p*R.inv()
    lifts.append(S.diag(trial,trial))
for lift in lifts:
    assert lift.T*lift==S.eye(8)
    image={tuple(lift*S.Matrix(r)) for r in roots}
    assert image==roots  # paired triality preserves actual E8 gluing
    assert all(lattice_e8(lift*base[:,i]) for i in range(4))
W=closure(gens);assert len(W)==1152
unit=key(S.eye(4))
rows=[]
for z in product(range(2),repeat=4):
    if not any(z):continue
    nr=sum(sum(a*b for a,b in zip(r,z))%2==0 for r in restr)
    stabilizer=[w for w in W if tuple(t%2 for t in act(w,z))==z]
    o3=[w for w in stabilizer if w!=unit and mult(mult(w,w),w)==unit]
    snfs=Counter()
    for w in o3:
        d=smith_normal_form(mat(w)-S.eye(4),domain=ZZ)
        snfs[tuple(abs(int(d[i,i])) for i in range(4))]+=1
    assert all(f==(1,1,0,0) for f in snfs)
    rows.append({'z':z,'centralizing_roots':nr,'stabilizer_order':len(stabilizer),
                 'order3_count':len(o3),'smith':dict((str(k),v) for k,v in snfs.items())})
assert Counter((r['centralizing_roots'],r['stabilizer_order']) for r in rows)==Counter({(112,384):3,(128,96):12})

order4=[]
for z in product(range(4),repeat=4):
    if all(t%2==0 for t in z):continue
    nr=sum(sum(a*b for a,b in zip(r,z))%4==0 for r in restr)
    if nr!=52:continue
    stab=[w for w in W if tuple(t%4 for t in act(w,z))==z]
    assert len(stab)==48
    order4.append(z)
assert len(order4)==24

# D8 central involution: y=2e5. The four negative 2-cycles lie in W(D8).
y=S.zeros(8,1);y[4]=2
assert lattice_e8(y) and not lattice_e8(y/2)
assert sum((S.Matrix(1,8,r)*y)[0]%4==0 for r in roots)==112
J=S.Matrix([[0,-1],[1,0]])
w4=S.diag(J,J,J,J)
assert S.expand(w4.charpoly().as_expr()-(S.Symbol('lambda')**2+1)**4)==0
assert {tuple(w4*S.Matrix(r)) for r in roots}==roots
assert lattice_e8((w4*y-y)/2)  # fixes s, retaining the E8 torus lattice
assert S.expand((-S.eye(8)).charpoly().as_expr()-(S.Symbol('lambda')+1)**8)==0

out={'model':'Euclidean E8 with doubled roots; actual D4+D4 lattice gluing',
     'root_count':len(roots),'central_lattice_basis':base.tolist(),
     'central_gram':Q.tolist(),'ambient_lifts_preserve_all_roots':True,
     'full_relative_order':len(W),'involution_rows':rows,
     'D5A3_order4_parameters':order4,
     'D8_Phi4_word_matrix':w4.tolist(),
     'D8_involution_cocharacter':list(y),
     'pure_torus_witnesses':{'e1':'identity','e2':'minus identity','e4':'four negative two-cycles'},
     'scope':'All odd q; all ell=3,5 valuations. No numerical field sample.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2,default=int)+'\n')
print(json.dumps({'passed':True,'roots':240,'relative_order':1152,'D8_Phi4_preserves_involution':True,
                  'order4_D5A3_parameters':24}))
