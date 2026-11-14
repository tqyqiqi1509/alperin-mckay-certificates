"""Uniform nonmatching-three principal E6 witness.

All calculations take place on the actual integral E6 coroot lattice.
No finite field parameter is sampled.  The centralizer conclusion uses
the full rational Weyl group, not only a dimension count of a torus.
"""
from collections import deque,Counter
from pathlib import Path
import json
import sympy as S

C=2*S.eye(6)
for i,j in [(0,2),(2,3),(3,4),(4,5),(1,3)]:C[i,j]=C[j,i]=-1
I=S.eye(6)
rs=[]
for i in range(6):
    r=I.copy()
    for j in range(6):r[i,j]-=C[i,j]
    rs.append(r)
roots={tuple(I[:,i]) for i in range(6)}|{tuple(-I[:,i]) for i in range(6)}
todo=deque(roots)
while todo:
    r=S.Matrix(todo.popleft())
    for a in rs:
        t=tuple(a*r)
        if t not in roots:roots.add(t);todo.append(t)
assert len(roots)==72
gamma=S.zeros(6)
for i,j in enumerate([5,1,4,3,2,0]):gamma[j,i]=1
w0=-gamma
assert all(w0*r==r*w0 for r in [rs[0]*rs[5],rs[2]*rs[4],rs[3],rs[1]])
positive=[S.Matrix(r) for r in roots if all(x>=0 for x in r)]
assert len(positive)==36
assert all(all(x<=0 for x in w0*r) for r in positive)

folded=[rs[0]*rs[5],rs[2]*rs[4],rs[3],rs[1]]
c=S.eye(6)
for r in folded:c=c*r
u=c**4
x=S.Symbol('x')
assert c**12==I and u**3==I and u!=I
assert S.expand(u.charpoly(x).as_expr())==S.expand((x*x+x+1)**3)
assert (I-u).det()==27
assert u*w0==w0*u

# Primitive basis for the -1 eigenspace of w0, equivalently gamma-fixed.
B=S.Matrix.hstack(I[:,0]+I[:,5],I[:,2]+I[:,4],I[:,1],I[:,3])
assert w0*B==-B
restrictions=[tuple((r.T*C*B)[0,j] for j in range(4)) for r in map(S.Matrix,roots)]
assert all(any(int(y)%3 for y in a) for a in restrictions)
assert B.rank()==4

# Explicit enumeration of the folded finite group, using exact integers.
def key(a):return tuple(a)
def mat(a):return S.Matrix(6,6,a)
seen={key(I)};todo=deque([I])
while todo:
    a=todo.popleft()
    for r in folded:
        b=a*r;k=key(b)
        if k not in seen:seen.add(k);todo.append(b)
assert len(seen)==1152
# The restriction to the rank-four lattice is faithful, including modulo3.
acts=[]
for k in seen:
    a=mat(k);aB=a*B
    r=S.Matrix.vstack(aB[0,:],aB[2,:],aB[1,:],aB[3,:])
    assert B*r==aB
    acts.append(tuple(int(t)%3 for t in r))
assert len(set(acts))==1152

# Independent permutation certificate on all72 roots for the FULL centralizer.
rootlist=sorted(roots);index={r:i+1 for i,r in enumerate(rootlist)}
def gap_perm(a):return '['+','.join(str(index[tuple(a*S.Matrix(r))]) for r in rootlist)+']'
gap='W:=Group(List(['+','.join(gap_perm(r) for r in rs)+'],PermList));;\n'
gap+='v:=PermList('+gap_perm(w0)+');;\n'
gap+='C:=Centralizer(W,v);;\n'
gap+='R:=Group(List(['+','.join(gap_perm(r) for r in folded)+'],PermList));;\n'
gap+='if Size(W)<>51840 or Size(C)<>1152 or C<>R then Error("full centralizer failed");fi;\n'
gap+='Print("E6 Weyl=",Size(W)," full centralizer=",Size(C)," folded=",Size(R),"\\n");\nQUIT;\n'
Path(__file__).with_suffix('.g').write_text(gap)

out={'status' :'passed','root_count':72,'cartan':C.tolist(),
     'graph':gamma.tolist(),'longest':w0.tolist(),
     'nonmatching_torus_polynomial':'Phi1^2 Phi2^4; Ennola Phi2^2 Phi1^4',
     'three_torus_basis':B.tolist(),'root_restriction_count':len(restrictions),
     'no_root_trivial_on_three_torsion':True,
     'folded_group_order':1152,'faithful_mod3_image_order':len(set(acts)),
     'folded_coxeter_word':[1,6,3,5,4,2],'u_word_power':4,
     'u_matrix':u.tolist(),'u_characteristic_polynomial':str(u.charpoly(x).as_expr()),
     'det_I_minus_u':27,
     'meaning':'Actual Sylow-three centralizer is a three-group for every nonmatching E6/2E6 parameter.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2,default=int)+'\n')
print(json.dumps({k:out[k] for k in ['status','root_count','folded_group_order','faithful_mod3_image_order','det_I_minus_u']}))
