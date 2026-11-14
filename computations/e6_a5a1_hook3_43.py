"""Formal hook valuations for all type-A factors used in the 3-height inverse."""
from pathlib import Path
from math import factorial
import json

def parts(n,b=None):
    if not n:
        yield ()
    else:
        for a in range(min(n,b or n),0,-1):
            for q in parts(n-a,a):yield (a,)+q

def vp(n):
    a=0
    while n%3==0 and n:
        n//=3;a+=1
    return a

def hooks(p):
    return [p[i]-j+sum(x>j for x in p[i+1:]) for i in range(len(p)) for j in range(p[i])]

data={}
for n in range(1,7):
    rows=[]
    for p in parts(n):
        hs=hooks(p)
        den=1
        for h in hs:den*=h
        deg=factorial(n)//den
        match=vp(factorial(n))-sum(vp(h) for h in hs)
        even=[h for h in hs if h%2==0]
        # Nonmatching 3-valuation is a*linear+constant, a=v3(q+1).
        linear=n//2-len(even)
        constant=vp(factorial(n//2))-sum(vp(h//2) for h in even)
        assert match==vp(deg) and linear>=0
        rows.append({'partition':p,'symmetric_degree':deg,'matching_v3':match,
                     'nonmatching_v3':[linear,constant]})
    data[str(n)]=rows
assert len([r for r in data['6'] if r['matching_v3']==0])==9
assert len([r for r in data['6'] if r['nonmatching_v3']==[0,0]])==6
assert next(r for r in data['6'] if r['partition']==(3,2,1))['nonmatching_v3']==[3,1]
# Every nontrivial 3^k rational orbit in the matching six-space has length 3.
# Its maximal centralizer is one length-3 orbit plus one scalar eigenspace of dim 3.
# v3(GL6)=(6a+2); v3(maximal such centralizer)=(a+1)+(3a+1).
data['nontoral_matching_minimal_index_v3']=[2,0]
data['cubic_self_twist_toral_minimal_index_v3']=2
data['scope']='Identities in the variable a>=1; no finite-q evaluation.'
p=Path(__file__).with_suffix('.json');p.write_text(json.dumps(data,indent=2)+'\n')
print('ALL_PARTITIONS_1_TO_6|MATCHING_A5_PRIME3=9|NONMATCHING_A5_PRIME3=6|321=3a+1|PASSED')
