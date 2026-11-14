"""Exact E7 toral action and central-twist transport, independent node 41.

Uses the Bourbaki Cartan matrix directly.  No finite-field q is sampled.
The geometric coweight points are not asserted to be rational SC elements;
their actual rational representatives live in the connected-centre cover.
"""
from collections import deque, Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
import json
import sympy as S

C = 2*S.eye(7)
for i,j in [(0,2),(1,3),(2,3),(3,4),(4,5),(5,6)]:
    C[i,j] = C[j,i] = -1
assert C.det() == 2
Ci = C.inv()
refl=[]
for i in range(7):
    m=S.eye(7)
    for j in range(7): m[i,j] -= C[i,j]
    refl.append(m)

def tup(v): return tuple(int(x) for x in v)
roots=set(tup(S.eye(7)[:,i]) for i in range(7))
todo=deque(roots)
while todo:
    v=S.Matrix(todo.popleft())
    for r in refl:
        w=tup(r*v)
        if w not in roots: roots.add(w); todo.append(w)
assert len(roots)==126

def rank3(mat):
    a=[[int(x)%3 for x in row] for row in mat.tolist()]
    m,n=len(a),len(a[0]); k=0
    for j in range(n):
        pivot=next((i for i in range(k,m) if a[i][j]),None)
        if pivot is None: continue
        a[k],a[pivot]=a[pivot],a[k]
        inv=1 if a[k][j]==1 else 2
        a[k]=[(inv*x)%3 for x in a[k]]
        for i in range(m):
            if i!=k:
                q=a[i][j]
                a[i]=[(x-q*y)%3 for x,y in zip(a[i],a[k])]
        k+=1
        if k==m: break
    return k

def perm(n,cyc):
    p=list(range(n))
    for seq in cyc:
        for i,j in zip(seq,seq[1:]+seq[:1]): p[i]=j
    m=S.zeros(n)
    for i,j in enumerate(p): m[j,i]=1
    return m

def augmentation(m):
    n=m.rows
    base=S.Matrix.hstack(*(S.eye(n)[:,i]-S.eye(n)[:,n-1] for i in range(n-1)))
    image=m*base
    assert all(sum(image[:,j])==0 for j in range(n-1))
    return image[:n-1,:]

action={}
models={
    'row5_A7':[augmentation(perm(8,[(0,1,2)])),augmentation(perm(8,[(3,4,5)]))],
    'row12_A3_A3_A1':[
        S.diag(augmentation(perm(4,[(0,1,2)])),S.eye(3),S.eye(1)),
        S.diag(S.eye(3),augmentation(perm(4,[(0,1,2)])),S.eye(1))],
    'row3_D6_A1':[perm(7,[(0,1,2)]),perm(7,[(3,4,5)])],
    'row6_C4':[perm(4,[(0,1,2)])],
    'row14_A3_A1':[S.diag(augmentation(perm(4,[(0,1,2)])),S.eye(1))]
}
for label,gens in models.items():
    n=gens[0].rows; identity=S.eye(n); values=[]
    for powers in product(range(3),repeat=len(gens)):
        g=identity
        for h,k in zip(gens,powers): g=g*h**k
        if g==identity: continue
        values.append(rank3((g-identity)**2))
    assert all(x>0 for x in values)
    action[label]={'dimension':n,'sylow_image_order':3**len(gens),
                   'square_ranks':values,'rank_counts':dict(Counter(values))}
negative=augmentation(perm(3,[(0,1,2)]))
assert rank3((negative-S.eye(2))**2)==0
action['negative_control_A2_augmentation']={'dimension':2,'square_rank':0}

def components(rs):
    unseen=set(rs); out=[]
    while unseen:
        x=unseen.pop(); found={x}; q=[x]
        while q:
            y=S.Matrix(q.pop())
            neighbors={z for z in unseen if (y.T*C*S.Matrix(z))[0]!=0}
            unseen-=neighbors; found|=neighbors; q.extend(neighbors)
        out.append(found)
    return sorted(out,key=lambda a:(-len(a),min(a)))

twists={}
for label,node,den,expected in [('A7',1,2,56),('A3_A3_A1',3,4,26),('D6_A1',0,2,62)]:
    point=Ci[:,node]/den; center=Ci[:,6]
    modulus=S.ilcm(*[x.q for x in point],*[x.q for x in center])
    start=tuple(int(x*modulus)%int(modulus) for x in point)
    target=tuple(int(x*modulus)%int(modulus) for x in point+center)
    previous={start:None}; q=deque([start])
    while q:
        v=q.popleft()
        for i,r in enumerate(refl):
            w=tuple(int(x)%int(modulus) for x in r*S.Matrix(v))
            if w not in previous:
                previous[w]=(v,i); q.append(w)
    rs={r for r in roots if (S.Matrix(r).T*C*point)[0].q==1}
    assert len(rs)==expected
    comps=components(rs)
    entry={'point_node':node+1,'point_denominator':den,'modulus':int(modulus),
           'weyl_orbit_size':len(previous),'root_component_sizes':[len(x) for x in comps],
           'central_twist_conjugate':target in previous}
    if target in previous:
        back=[]; v=target
        while previous[v] is not None:
            v,i=previous[v]; back.append(i)
        path=list(reversed(back)); w=S.eye(7)
        for i in path: w=refl[i]*w
        assert all(x.q==1 for x in w*point-point-center)
        images=[]
        for comp in comps:
            img={tup(w*S.Matrix(r)) for r in comp}
            images.append(next(i for i,c in enumerate(comps) if c==img))
        entry.update({'simple_reflections_in_application_order':[i+1 for i in path],
                      'component_permutation':images,
                      'weyl_matrix':[[int(x) for x in w.row(i)] for i in range(7)]})
    twists[label]=entry
assert twists['A7']['central_twist_conjugate']
assert twists['A3_A3_A1']['central_twist_conjugate']
assert twists['A3_A3_A1']['component_permutation']==[1,0,2]
assert not twists['D6_A1']['central_twist_conjugate']
result={'Cartan':[[int(x) for x in C.row(i)] for i in range(7)],
        'roots':len(roots),'action_certificates':action,'central_twist_transport':twists,
        'scope':'Integral root datum, valid for both matching Ennola signs and every a>=1; rational lifts supplied by KM25 connected-centre source rows.'}
path=Path(__file__).with_suffix('.json')
path.write_text(json.dumps(result,indent=2)+'\n')
for label,value in action.items(): print(label,value)
for label,value in twists.items():
    print(label,{k:v for k,v in value.items() if k!='weyl_matrix'})
print('All exact node-41 certificates passed.')
