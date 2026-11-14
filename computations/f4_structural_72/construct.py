"""Independent integral F4 certificate for whole 3-adic source families.

Only the Cartan matrix is entered. No prior root list, Weyl matrix, or
certificate is imported. Small valuations below a proved coefficient bound
are exhaustive congruence classes of all q, not sample finite fields.
"""
from pathlib import Path
from itertools import product
import json

n=4
I=tuple(tuple(int(i==j) for j in range(n)) for i in range(n))
C=((2,-1,0,0),(-1,2,-1,0),(0,-2,2,-1),(0,0,-1,2))
G=tuple(tuple((2,2,1,1)[i]*C[i][j] for j in range(n)) for i in range(n))
def mm(A,B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(n))
                       for j in range(n)) for i in range(n))
def mv(A,v):
    return tuple(sum(A[i][k]*v[k] for k in range(n)) for i in range(n))
def add(A,B):
    return tuple(tuple(A[i][j]+B[i][j] for j in range(n)) for i in range(n))
def scale(c,A):
    return tuple(tuple(c*x for x in row) for row in A)
def minus(A,B):return add(A,scale(-1,B))
def dot(a,b):return sum(a[i]*mv(G,b)[i] for i in range(n))
def refl(r):
    row=mv(G,r)
    norm=dot(r,r)
    assert all((2*x)%norm==0 for x in row)
    coeff=tuple(2*x//norm for x in row)
    return tuple(tuple(I[i][j]-r[i]*coeff[j] for j in range(n)) for i in range(n))
simple=tuple(tuple(int(i==j) for i in range(n)) for j in range(n))
S=tuple(refl(r) for r in simple)
roots=set(simple); queue=list(simple)
while queue:
    r=queue.pop()
    for h in S:
        x=mv(h,r)
        if x not in roots:roots.add(x);queue.append(x)
roots=tuple(sorted(roots))
assert len(roots)==48
W={I};queue=[I]
while queue:
    w=queue.pop()
    for h in S:
        x=mm(h,w)
        if x not in W:W.add(x);queue.append(x)
W=tuple(sorted(W));assert len(W)==1152
def zero_mod(A,m):return all(x%m==0 for row in A for x in row)
def nonzero(A):return any(x for row in A for x in row)
def maxabs(A):return max(abs(x) for row in A for x in row)

report={'cartan':C,'twice_gram':G,'roots':roots,'weyl_order':len(W),
        'pure':{},'rank_one':[],'mixed':[]}
report['pure']['trivial_roots_mod3']=[r for r in roots if all(x%3==0 for x in r)]
report['pure']['pointwise_weyl_mod3']=[w for w in W if zero_mod(minus(w,I),3)]
assert not report['pure']['trivial_roots_mod3']
assert report['pure']['pointwise_weyl_mod3']==[I]
order3rank2=[w for w in W if w!=I and mm(mm(w,w),w)==I
             and sum(w[i][i] for i in range(n))==1]
assert all(not zero_mod(mm(minus(w,I),minus(w,I)),3) for w in order3rank2)
models=[mm(S[0],S[1]),mm(S[2],S[3])]
classes=[[w for w in order3rank2 if any(mm(h,v)==mm(w,h) for h in W)]
         for v in models]
assert set(classes[0]).isdisjoint(classes[1])
assert set(classes[0]+classes[1])==set(order3rank2)
report['pure']['order3_rank2_class_sizes']=[len(c) for c in classes]
report['pure']['all_rank2_order3_actions_noncentral_commutator_mod3']=True
for i in (0,3):
    trivial=[r for r in roots if all(r[j]%3==0 for j in range(n) if j!=i)]
    pointwise=[w for w in W if all((w[j][k]-I[j][k])%3==0
               for j in range(n) if j!=i for k in range(n))]
    assert set(trivial)=={simple[i],tuple(-x for x in simple[i])}
    assert set(pointwise)=={I,S[i]}
    # v fixes alpha; its relative action on ker alpha has class > 2 at p=3.
    candidates=[w for w in W if w!=I and mm(mm(w,w),w)==I and mv(w,simple[i])==simple[i]]
    assert candidates and all(any(mm(minus(w,I),minus(w,I))[j][k]%3
            for j in range(n) if j!=i for k in range(n)) for w in candidates)
    v=next(w for w in candidates if any(mm(minus(w,I),minus(w,I))[j][k]%3
            for j in range(n) if j!=i for k in range(n)))
    delta=mm(minus(v,I),minus(v,I))
    report['rank_one'].append({'simple_root':i+1,'trivial_roots_mod3':trivial,
      'pointwise_weyl_mod3':pointwise,'order3_actions_checked':len(candidates),
      'all_order3_actions_noncentral_commutator_mod3':True,'v':v,'delta':delta})

for name,w in [('long_A2',mm(S[0],S[1])),('short_A2',mm(S[2],S[3]))]:
    w2=mm(w,w);assert mm(w2,w)==I
    N=add(add(I,w),w2);K=add(scale(2,w2),w)
    perp=[r for r in roots if mv(w,r)==r]
    assert len(perp)==6
    b,c=next((b,c) for b in perp for c in perp
             if dot(b,b)==dot(c,c) and 2*dot(b,c)==-dot(b,b))
    v=mm(refl(b),refl(c));assert mm(v,w)==mm(w,v) and mm(mm(v,v),v)==I
    delta=mm(minus(v,I),minus(v,I))
    n_delta=mm(N,delta)
    bounds={'root':max(abs(x) for r in roots for x in mv(N,r)),
      'weyl':max(maxabs(mm(N,minus(h,I))) for h in W)}
    cut=1
    while 3**cut<=max(bounds.values()):cut+=1
    small=[]
    for a,u in product(range(1,cut),(1,2)):
        R=add(N,scale(3**a*u,K));mod=3**(a+1)
        trivial=[r for r in roots if all(x%mod==0 for x in mv(R,r))]
        pointwise=[h for h in W if zero_mod(mm(R,minus(h,I)),mod)]
        commutator=not zero_mod(mm(R,delta),mod)
        assert not trivial and pointwise==[I] and commutator
        small.append({'a':a,'u_mod3':u,'trivial_roots':trivial,
          'pointwise_weyl':pointwise,'commutator_noncentral':commutator})
    large_roots=[r for r in roots if not any(mv(N,r)) and
                  all(x%3==0 for x in mv(K,r))]
    large_pointwise=[h for h in W if not nonzero(mm(N,minus(h,I))) and
                     zero_mod(mm(K,minus(h,I)),3)]
    assert not large_roots and large_pointwise==[I]
    # This entry has exact 3-adic valuation 1, so no 3^a correction, a>=2,
    # can make N*delta divisible by 3^a.
    witness=next((i,j,n_delta[i][j]) for i in range(n) for j in range(n)
                if n_delta[i][j]%3==0 and n_delta[i][j]%9!=0)
    report['mixed'].append({'name':name,'w':w,'N':N,'K':K,
      'perpendicular_simple_pair':[b,c],'v':v,'N_delta':n_delta,
      'noncentral_witness':witness,'bounds':bounds,'cut':cut,
      'small_valuations':small,'large_trivial_roots':large_roots,
      'large_pointwise_weyl':large_pointwise})

out=Path(__file__).with_name('certificate.json')
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'root_count':len(roots),'weyl_order':len(W),
 'rank_one':[{'root':c['simple_root'],'centralizer_roots':len(c['trivial_roots_mod3']),
  'pointwise_W':len(c['pointwise_weyl_mod3'])} for c in report['rank_one']],
 'mixed':[{'name':c['name'],'bounds':c['bounds'],'cut':c['cut'],
 'all_trivial_root_sets_empty':True,'all_pointwise_W_trivial':True,
 'commutator_noncentral_all_a':True,'witness':c['noncentral_witness']}
 for c in report['mixed']]},indent=2))
