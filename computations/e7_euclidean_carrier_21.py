"""Independent all-q E7/D4 carrier calculation in the Euclidean E8 model.

Uses doubled eight-dimensional root/weight coordinates, not node 16's
Cartan-lattice orbit construction.  Every displayed relation is integral.
"""
from itertools import combinations, product
from collections import Counter, deque
from pathlib import Path
import json

def add(x, y): return tuple(a+b for a,b in zip(x,y))
def scale(a, x): return tuple(a*b for b in x)
def dot(x, y):
    a=sum(b*c for b,c in zip(x,y))
    assert a % 4 == 0
    return a//4
def unit(i, a=2): return tuple(a if j==i else 0 for j in range(8))

simple=[(1,-1,-1,-1,-1,-1,-1,1), add(unit(0),unit(1))]
simple += [add(unit(i),unit(i-1,-2)) for i in range(1,6)]
roots=set()
for i,j in combinations(range(6),2):
    for a,b in product((-2,2),repeat=2): roots.add(add(unit(i,a),unit(j,b)))
roots.update((add(unit(7),unit(6,-2)),add(unit(6),unit(7,-2))))
for signs in product((-1,1),repeat=6):
    if signs.count(-1)%2:
        for a in (-1,1): roots.add(signs+(a,-a))
assert len(roots)==126

weights=set()
for i in range(6):
    for a,b in product((-1,1),repeat=2):
        weights.add(add(unit(i,2*a),add(unit(7,b),unit(6,-b))))
for signs in product((-1,1),repeat=6):
    if signs.count(-1)%2==0: weights.add(signs+(0,0))
weights=sorted(weights); wi={w:i for i,w in enumerate(weights)}
assert len(weights)==56
def reflect(x,r): return add(x,scale(-dot(x,r),r))
for r in simple:
    assert {reflect(x,r) for x in roots}==roots
    assert {reflect(x,r) for x in weights}==set(weights)

ID=tuple((i,1) for i in range(56))
def mul(a,b): return tuple((a[j][0],s*a[j][1]) for j,s in b)
def inv(a):
    out=[None]*56
    for i,(j,s) in enumerate(a): out[j]=(i,s)
    return tuple(out)
def power(a,k):
    v=ID
    for _ in range(k): v=mul(v,a)
    return v
n=[]
for r in simple:
    n.append(tuple((wi[reflect(w,r)], -1 if dot(w,r)==1 else 1) for w in weights))
def word(indices):
    v=ID
    for i in indices: v=mul(v,n[i-1])
    return v
def h(r): return tuple((i,(-1)**dot(w,r)) for i,w in enumerate(weights))
def raise_at(w,i):
    return add(w,simple[i]) if dot(w,simple[i])==-1 else None
def lower_at(w,i):
    return add(w,scale(-1,simple[i])) if dot(w,simple[i])==1 else None
for w in weights:
    for i in range(7):
        for j in range(7):
            a=lower_at(w,j); ef=raise_at(a,i) if a is not None else None
            b=raise_at(w,i); fe=lower_at(b,j) if b is not None else None
            if i!=j: assert ef==fe
            else:
                c=dot(w,simple[i])
                assert (ef,fe)==((w,None) if c==1 else (None,w) if c==-1 else (None,None))
            if i!=j and dot(simple[i],simple[j])==0:
                for op in (raise_at,lower_at):
                    a=op(w,j); ij=op(a,i) if a is not None else None
                    a=op(w,i); ji=op(a,j) if a is not None else None
                    assert ij==ji
            if i!=j and dot(simple[i],simple[j])==-1:
                for op in (raise_at,lower_at):
                    a=op(w,i)
                    a=op(a,j) if a is not None else None
                    a=op(a,i) if a is not None else None
                    assert a is None
def pin_sign(a,i,j):
    v=[]
    for w in weights:
        x=raise_at(w,i)
        if x is not None:
            k,s=a[wi[w]]; l,t=a[wi[x]]
            assert raise_at(weights[k],j)==weights[l]
            v.append(s*t)
    assert len(set(v))==1
    return v[0]
for i in range(7):
    assert power(n[i],2)==h(simple[i])
    for j in range(i):
        if dot(simple[i],simple[j])==0: assert mul(n[i],n[j])==mul(n[j],n[i])
        else: assert word([i+1,j+1,i+1])==word([j+1,i+1,j+1])

d4=[1,2,3,4]
centre=[add(unit(5),unit(4,-2)),add(unit(5),unit(4)),add(unit(7),unit(6,-2))]
assert all(dot(x,y)==(2 if i==j else 0) for i,x in enumerate(centre) for j,y in enumerate(centre))
assert all(dot(x,simple[i])==0 for x in centre for i in d4)
counts=[]
for k in range(1,4):
    zero={r for r in roots if all(dot(r,c)==0 for c in centre[:k])}
    mod={r for r in roots if all(dot(r,c)%3==0 for c in centre[:k])}
    assert zero==mod
    counts.append(len(zero))
assert counts==[60,26,24]
mu=tuple(sum(x) for x in zip(*centre))
mu_values=Counter(dot(r,mu) for r in roots)
assert mu_values=={0:72,2:27,-2:27}

# Transport the simple root alpha7 to each perpendicular A1 root.
trans={simple[6]:()}; queue=deque([simple[6]])
while queue:
    r=queue.popleft()
    for i,s in enumerate(simple):
        z=reflect(r,s)
        if z not in trans: trans[z]=(i+1,)+trans[r]; queue.append(z)
assert len(trans)==126
flips=[]
for r in centre:
    a=word(trans[r]); b=mul(mul(a,n[6]),inv(a))
    assert power(b,2)==h(r)
    assert all(pin_sign(b,i,i)==1 for i in d4)
    flips.append(b)
q1=word([1,3,4,5,2,4,3,1]); q2=word([6,5,4,3,2,4,5,6])
assert power(q1,2)==ID and power(q2,2)==ID
sigma=mul(q1,q2)
assert sigma!=ID and power(sigma,3)==ID
for q in (q1,q2):
    images=[]
    for i in d4:
        found=[]
        for j in d4:
            try:
                if pin_sign(q,i,j)==1: found.append(j)
            except AssertionError: pass
        assert len(found)==1
        images.extend(found)
    assert sorted(images)==d4
gens=flips+[q1,q2]
carrier={ID}; queue=deque([ID])
while queue:
    a=queue.popleft()
    for b in gens:
        z=mul(a,b)
        if z not in carrier: carrier.add(z); queue.append(z)
assert len(carrier)==384
unsigned={tuple(j for j,s in x) for x in carrier}
assert len(unsigned)==48
toral={x for x in carrier if tuple(j for j,s in x)==tuple(range(56))}
assert len(toral)==8

# Coxeter half-powers give the actual longest Tits elements (reduced words).
long_e7=power(word([1,2,3,4,5,6,7]),9)
long_d4=power(word([2,3,4,5]),3)
assert all(weights[j]==scale(-1,w) for w,(j,s) in zip(weights,long_e7))
m=mul(long_e7,inv(long_d4))
z=tuple((i,-1) for i in range(56))
assert power(m,2)==z
assert all(pin_sign(m,i,i)==1 for i in d4)
assert all(mul(m,x)==mul(x,m) for x in carrier)
assert all(any(dot(r,s)%3 for s in simple) for r in roots)

out={"model":"E8 Euclidean slice orthogonal to e7+e8; doubled coordinates",
     "root_count":126,"weight_count":56,"centralizer_root_counts":counts,
     "centre_roots_doubled":centre,"fixed_line_doubled":mu,
     "fixed_line_root_pairings":dict(mu_values),"full_carrier_order":len(carrier),
     "toral_2kernel":len(toral),"relative_order":len(unsigned),
     "q1_q2_relations":True,"opposition_square":"minus identity on the faithful 56-space",
     "Chevalley_Serre_integral_relations":True,
     "opposition_centralizes_carrier":True,"opposition_fixes_D4_pinning":True,
     "scope":"uniform integral relations, all q and all valuations"}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
