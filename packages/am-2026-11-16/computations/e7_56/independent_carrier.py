"""Independent E7 natural carrier in the A7 Euclidean realization.

No earlier root, weight, matrix, or group output is imported.
Weights are the two exterior-square A7 orbits, not a Weyl orbit
generated in fundamental-weight coordinates. Integer actions only.
"""
from itertools import combinations, product
from collections import deque
from pathlib import Path
import json
import sympy as s

OUT=Path(__file__).resolve().parent
def diff(i,j):
    x=[0]*8;x[i]=2;x[j]=-2;return tuple(x)
roots=set(diff(i,j) for i in range(8) for j in range(8) if i!=j)
for pos in combinations(range(8),4):
    roots.add(tuple(1 if i in pos else -1 for i in range(8)))
assert len(roots)==126
alpha=[diff(5,6),(-1,-1,-1,-1,1,1,1,1),diff(4,5),diff(3,4),diff(2,3),diff(1,2),diff(0,1)]
B=s.Matrix.hstack(*(s.Matrix(x) for x in alpha)); C=B.T*B/4
assert C.det()==2
simple_coords={r:tuple(B.gauss_jordan_solve(s.Matrix(r))[0]) for r in roots}
assert all(all(z.q==1 for z in c) and (all(z>=0 for z in c) or all(z<=0 for z in c)) for c in simple_coords.values())
weights=sorted({tuple(sign*(3 if i in pos else -1) for i in range(8)) for pos in combinations(range(8),2) for sign in [-1,1]})
assert len(weights)==56
wi={w:i for i,w in enumerate(weights)}
def pairing(w,r):
    q=sum(x*y for x,y in zip(w,r));assert q%8==0
    return q//8
assert {pairing(w,r) for w in weights for r in roots}=={-1,0,1}
def move(w,r,sign=1):return tuple(x+sign*2*y for x,y in zip(w,r))
def raising(w,i):return move(w,alpha[i]) if pairing(w,alpha[i])==-1 else None
def lowering(w,i):return move(w,alpha[i],-1) if pairing(w,alpha[i])==1 else None
rel=0
for w in weights:
    for i,j in product(range(7),repeat=2):
        f=lowering(w,j);ef=raising(f,i) if f else None
        e=raising(w,i);fe=lowering(e,j) if e else None
        if i!=j:assert ef==fe
        else:
            m=pairing(w,alpha[i]);assert (int(ef==w)-int(fe==w))==m
        if i!=j and C[i,j]==0:
            for op in [raising,lowering]:
                u=op(w,j);v=op(w,i)
                assert (op(u,i) if u else None)==(op(v,j) if v else None)
        elif i!=j and C[i,j]==-1:
            for op in [raising,lowering]:
                u=op(w,i);u=op(u,j) if u else None;u=op(u,i) if u else None
                assert u is None
        rel+=1

# Signed actions use an independent +/- (index+1) coding.
ID=tuple(range(1,57))
def compose(a,b):return tuple((1 if j>0 else -1)*a[abs(j)-1] for j in b)
def inverse(a):
    ans=[0]*len(a)
    for i,j in enumerate(a):ans[abs(j)-1]=(i+1)*(1 if j>0 else -1)
    return tuple(ans)
def power(a,n):
    if n<0:return power(inverse(a),-n)
    ans=tuple(range(1,len(a)+1))
    while n:
        if n&1:ans=compose(ans,a)
        a=compose(a,a);n//=2
    return ans
def simple_operator(i):
    r=alpha[i];ans=[]
    for w in weights:
        m=pairing(w,r);z=move(w,r,-m)
        ans.append((wi[z]+1)*(-1 if m==1 else 1))
    return tuple(ans)
ns=[simple_operator(i) for i in range(7)]
for i,j in product(range(7),repeat=2):
    if i==j:
        assert power(ns[i],2)==tuple((k+1)*(-1 if pairing(w,alpha[i]) else 1) for k,w in enumerate(weights))
    elif C[i,j]==0:assert compose(ns[i],ns[j])==compose(ns[j],ns[i])
    else:assert compose(compose(ns[i],ns[j]),ns[i])==compose(compose(ns[j],ns[i]),ns[j])

def pinning_sign(a,i,j):
    vals=[]
    for w in weights:
        u=raising(w,i)
        if u is None:continue
        x=a[wi[u]];y=a[wi[w]]
        assert raising(weights[abs(y)-1],j)==weights[abs(x)-1]
        vals.append((1 if x>0 else -1)*(1 if y>0 else -1))
    assert len(set(vals))==1
    return vals[0]

# Perpendicular A5 recovered from Euclidean orthogonality.
H=[0,2]
kr={r for r in roots if all(sum(x*y for x,y in zip(r,alpha[i]))==0 for i in H)}
positive={r for r in kr if all(c>=0 for c in simple_coords[r])}
kb=sorted(r for r in positive if not any(tuple(x-y for x,y in zip(r,t)) in positive for t in positive))
assert len(kr)==30 and len(kb)==5
adj={i:[j for j in range(5) if sum(x*y for x,y in zip(kb[i],kb[j]))==-4] for i in range(5)}
chain=[next(i for i in adj if len(adj[i])==1)]
while len(chain)<5:chain.append(next(j for j in adj[chain[-1]] if j not in chain))
kb=[kb[i] for i in chain]

# BFS transports roots, with a fresh actual operator at every vertex.
def refl(r,i):
    m=sum(x*y for x,y in zip(r,alpha[i]))//4
    return tuple(x-m*y for x,y in zip(r,alpha[i]))
root_lifts={alpha[0]:(ID,0,[])};q=deque([alpha[0]])
while q:
    r=q.popleft();a,j,word=root_lifts[r]
    for i in range(7):
        t=refl(r,i)
        if t not in root_lifts:
            root_lifts[t]=(compose(ns[i],a),j,[i]+word);q.append(t)
assert len(root_lifts)==126
kn=[]
for r in kb:
    a,j,word=root_lifts[r]
    x=compose(compose(a,ns[j]),inverse(a));kn.append(x)
    assert all(pinning_sign(x,i,i)==1 for i in H)
for i,j in product(range(5),repeat=2):
    if i==j:continue
    if abs(i-j)==1:assert compose(compose(kn[i],kn[j]),kn[i])==compose(compose(kn[j],kn[i]),kn[j])
    else:assert compose(kn[i],kn[j])==compose(kn[j],kn[i])

# Choose simple pinning signs by comparing the conjugated operators.
# They already satisfy the Chevalley A5 braid and square relations.
I6=tuple(range(1,7));n6=[];s6=[]
for i in range(5):
    x=list(I6);x[i]=-(i+2);x[i+1]=i+1;n6.append(tuple(x));s6.append(tuple(abs(z) for z in x))
permutation_words={I6:[]};q=deque([I6])
while q:
    a=q.popleft()
    for i,b in enumerate(s6):
        c=compose(a,b)
        if c not in permutation_words:permutation_words[c]=permutation_words[a]+[i];q.append(c)
assert len(permutation_words)==720
def monomial_lift(target):
    word=permutation_words[tuple(abs(j) for j in target)]
    p6=I6;p56=ID
    for i in word:p6=compose(p6,n6[i]);p56=compose(p56,kn[i])
    diagonal=compose(target,inverse(p6))
    assert all(abs(j)==i+1 for i,j in enumerate(diagonal))
    bits=[int(j<0) for j in diagonal];assert sum(bits)%2==0
    tor56=ID;tor6=I6;prefix=0
    for i in range(5):
        prefix^=bits[i]
        if prefix:tor56=compose(tor56,power(kn[i],2));tor6=compose(tor6,power(n6[i],2))
    assert compose(tor6,p6)==target
    return compose(tor56,p56)
a1=monomial_lift((2,3,1,4,5,6));a2=monomial_lift((1,2,3,5,6,4))
n=monomial_lift((-4,-5,-6,1,2,3))
c=ID
for ni in ns:c=compose(c,ni)
n0=power(c,9)
assert all(weights[abs(n0[i])-1]==tuple(-v for v in w) for i,w in enumerate(weights))
nh=compose(compose(ns[0],ns[2]),ns[0]);z=compose(n0,inverse(nh))
zG=tuple(-i for i in ID)
assert power(n0,2)==zG and power(nh,2)==ID
assert pinning_sign(z,0,2)==pinning_sign(z,2,0)==1
assert all(compose(z,ni)==compose(ni,z) for ni in kn)
assert power(a1,3)==power(a2,3)==ID
assert compose(a1,a2)==compose(a2,a1)
assert compose(compose(n,a1),inverse(n))==a2
assert power(n,2)==power(z,2)==zG
v3=compose(a1,a2);v6=compose(v3,z)
assert power(v3,3)==ID and power(v6,6)==zG
gens=[a1,a2,n,z]
assert all(compose(v,x)==compose(x,v) for v in [v3,v6] for x in gens)

# Enumerate actual 72-element group while checking the C4 line on ALL edges.
elts={ID:0};q=deque([ID]);values=[0,0,1,1]
while q:
    a=q.popleft()
    for g,val in zip(gens,values):
        b=compose(a,g);ex=(elts[a]+val)%4
        if b in elts:assert elts[b]==ex
        else:elts[b]=ex;q.append(b)
assert len(elts)==72 and elts[zG]==2
unsigned={tuple(abs(j) for j in a) for a in elts}
kernel=[a for a in elts if tuple(abs(j) for j in a)==ID]
assert len(unsigned)==36 and set(kernel)=={ID,zG}
central=[a for a in elts if all(compose(a,g)==compose(g,a) for g in gens)]
assert len(central)==12

record={'realization':'A7 coordinate E7; weights +/- exterior-square weights',
 'roots':126,'weights':56,'simple_Cartan':[list(map(int,C.row(i))) for i in range(7)],
 'chevalley_weight_relations':rel,'simple_Tits_braids':True,
 'H_nodes':[1,3],'K_simple_roots_scaled2':kb,
 'K_lift_words':[root_lifts[r][2] for r in kb],
 'longest_word':'(n1 n2 n3 n4 n5 n6 n7)^9',
 'U_order':len(elts),'W_order':len(unsigned),'U_center_order':len(central),
 'toral_kernel_exact_central2':True,'z0_squared':'-I56','n_squared':'-I56',
 'u6_sixth_power':'-I56','faithful_line_c4_all_Cayley_edges':True,
 'field_return_exponent_preserved':12}
(OUT/'independent_carrier.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
