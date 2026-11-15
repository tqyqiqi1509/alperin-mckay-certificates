"""All 127 proper standard E7 Levis from the Cartan matrix alone.

This verifies the finite lattice assertions in e7_whole_positive_install_76.md.
It reads no predecessor arrays and is independent of character/block claims.
"""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

C = 2*s.eye(7)
for i,j in [(0,2),(2,3),(3,4),(4,5),(5,6),(1,3)]:
    C[i,j]=C[j,i]=-1
assert C.det()==2

def invariants(A):
    D=smith_normal_form(A,domain=s.ZZ)
    return [abs(int(D[i,i])) for i in range(min(D.shape)) if abs(D[i,i])>1]

def root_type(I):
    todo=set(I); labels=[]
    while todo:
        stack=[min(todo)]; comp=set()
        while stack:
            i=stack.pop()
            if i in comp: continue
            comp.add(i); todo.discard(i)
            stack += [j for j in I if C[i,j]==-1 and j not in comp]
        adj={i:[j for j in comp if C[i,j]==-1] for i in comp}
        branch=[i for i in comp if len(adj[i])==3]
        if not branch: label='A'+str(len(comp))
        else:
            assert len(branch)==1
            b=branch[0]; arms=[]
            for v in adj[b]:
                old=b; length=1
                while len(adj[v])==2:
                    w=next(w for w in adj[v] if w!=old)
                    old,v=v,w; length+=1
                arms.append(length)
            arms.sort()
            if arms[:2]==[1,1]: label='D'+str(len(comp))
            else:
                assert arms==[1,2,2]
                label='E6'
        labels.append(label)
    return '+'.join(sorted(labels)) if labels else 'torus'

# The unique nonzero mod-two null vector controls every saturation defect.
k=s.Matrix([0,1,0,0,1,0,1])
assert all(int(x)%2==0 for x in C*k)
assert sum(all(int(y)%2==0 for y in C*s.Matrix(x))
           for x in product(range(2),repeat=7))==2
support={1,4,6}
rows=[]
for r in range(7):
    for I in combinations(range(7),r):
        if not I:
            rows.append(dict(nodes=[],type='torus',centre_components=[],multiplication_kernel=[]))
            continue
        comp=invariants(C[:,I])
        assert comp==([2] if support.issubset(I) else [])
        A=C.extract(I,I)
        if comp:
            v=C[:,I]*k.extract(I,[0])/2
            assert all(x.q==1 for x in v)
            A=A.row_join(v.extract(I,[0]))
        ker=invariants(A)
        assert s.prod(comp)*s.prod(ker)==C.extract(I,I).det()
        rows.append(dict(nodes=[i+1 for i in I],type=root_type(I),
                         centre_components=comp,multiplication_kernel=ker))
assert len(rows)==127
assert sum(bool(x['centre_components']) for x in rows)==15

# Independent root closure under the seven integral reflections.
roots={tuple(s.eye(7)[:,i]) for i in range(7)}
queue=list(roots)
while queue:
    x=s.Matrix(queue.pop())
    for i in range(7):
        y=x.copy(); y[i]-=(C*x)[i]
        t=tuple(y)
        if t not in roots: roots.add(t);queue.append(t)
assert len(roots)==126
assert all((s.Matrix(x).T*C*s.Matrix(x))[0]==2 for x in roots)

counts=Counter((r['type'],tuple(r['centre_components']),tuple(r['multiplication_kernel'])) for r in rows)
out=dict(cartan=[list(map(int,C.row(i))) for i in range(7)],root_count=len(roots),
         proper_levi_count=len(rows),disconnected_schematic_centres=15,
         saturation_support=[2,5,7],
         summary=[dict(type=t,centre_components=list(c),multiplication_kernel=list(k),standard_subsets=n)
                  for (t,c,k),n in sorted(counts.items())],proper_levis=rows)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['proper_levis','cartan']},indent=2))
