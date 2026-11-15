from itertools import product,combinations
from pathlib import Path
import json
import sympy as s
from sympy.combinatorics import Permutation as P,PermutationGroup as PG
C=2*s.eye(6)
for i,j in [(0,2),(2,3),(3,4),(4,5),(1,3)]:C[i,j]=C[j,i]=-1
simple=[tuple(int(v) for v in s.eye(6)[:,i]) for i in range(6)]
def dot(a,b):return int((s.Matrix(a).T*C*s.Matrix(b))[0])
def ref(a,b):return tuple(b[i]-dot(a,b)*a[i] for i in range(6))
roots=set(simple);todo=list(simple)
while todo:
 b=todo.pop()
 for a in simple:
  z=ref(a,b)
  if z not in roots:roots.add(z);todo.append(z)
roots=sorted(roots);idx={r:i for i,r in enumerate(roots)};n=len(roots);ID=P(list(range(2*n)))
def eps(a,b):return (sum(a[i]*b[i] for i in range(6))+sum(a[i]*int(C[i,j])*b[j] for i in range(6) for j in range(i)))%2
def lift(a):
 im=[]
 for b in roots:
  v=dot(a,b);r=ref(a,b)
  sign=eps(a if v==-1 else tuple(-x for x in a),b) if abs(v)==1 else 0
  im.append(idx[r]+n*sign)
 return P(im+[(j+n)%(2*n) for j in im])
def tor(v):
 im=[i+n*(dot(r,v)%2) for i,r in enumerate(roots)]
 return P(im+[(j+n)%(2*n) for j in im])
N=[lift(a) for a in simple]
assert all(N[i]**2==tor(simple[i]) for i in range(6))
for i,j in combinations(range(6),2):assert N[i]*N[j]==N[j]*N[i] if C[i,j]==0 else N[i]*N[j]*N[i]==N[j]*N[i]*N[j]
# Pinned graph automorphism: root-vector signs determined recursively by brackets.
perm=[5,1,4,3,2,0]
def graph(r):return tuple(r[perm[i]] for i in range(6))
gsign={a:0 for a in simple}
for r in sorted([r for r in roots if min(r)>=0],key=sum):
 if r in gsign:continue
 a=next(a for a in simple if tuple(r[i]-a[i] for i in range(6)) in gsign)
 b=tuple(r[i]-a[i] for i in range(6))
 gsign[r]=(gsign[a]+gsign[b]+eps(a,b)+eps(graph(a),graph(b)))%2
for r in list(gsign):gsign[tuple(-v for v in r)]=gsign[r]
im=[idx[graph(r)]+n*gsign[r] for r in roots];B=P(im+[(j+n)%(2*n) for j in im])
assert B**2==ID
assert all(B**-1*N[i]*B==N[perm[i]] for i in range(6))
H=[(0,0,1,1,1,1),(0,1,0,1,1,1),(1,0,1,1,1,0),(1,1,1,1,0,0)]
K=[(0,0,0,1,0,0),(0,1,1,1,1,0)]
base=H+K
T=s.Matrix.hstack(*map(s.Matrix,base));assert abs(T.det())==3
assert abs((T.T*C).det())==9
Araw=N[4]*N[2]*N[1]
def correct(op):
 found=[]
 for v in product([0,1],repeat=6):
  x=op*tor(v)
  if all(x(idx[r])<n for r in base):found.append((v,x))
 assert len(found)==1
 return found[0]
av,A=correct(Araw);bv,B=correct(B)
assert A**2==B**2==ID and A*B==B*A
assert all(A(idx[r])%n in [idx[x] for x in base] for r in base)
assert all(B(idx[r])%n in [idx[x] for x in base] for r in base)
# longest SL3 representative for perpendicular A2 (orientation fixed by lifts)
J=lift(K[0])*lift(K[1])*lift(K[0]);assert J**2==ID
Tau=J**-1*A
assert Tau**2==ID
# Signed permutation transpositions are -P12 and -P23. Correct root lift by
# coroot 2-torsion so it induces the expected diagonal sign change on third basis.
# Search SL3's four torus2 corrections and retain two involutions centralizing Tau.
choices=[]
for k in K:
 cs=[]
 for u,v in product([0,1],repeat=2):
  z=tuple(u*K[0][i]+v*K[1][i] for i in range(6));r=lift(k)*tor(z)
  if r**2==ID and r*Tau==Tau*r:cs.append(r)
 choices.append(cs)
rs=next((x,y) for x in choices[0] for y in choices[1] if (x*y)**3==ID and x!=y)
Q=PG([*rs,Tau]);assert Q.order()==12
assert all(B*x==x*B for x in [*rs,Tau])
frames=[ID,Tau,Tau*B,B]
assert all(all(f*x==x*f for x in [*rs,Tau]) for f in frames)
# Source component group action via root lattice quotient: torsion character
# is omega1 modulo span(source roots). Check a acts identity, graph negates.
W=s.eye(6)
for i in [4,2,1]:
 Ri=s.eye(6);Ri[i,:]=Ri[i,:]-C[i,:];W=Ri*W
G=s.zeros(6)
for i,j in enumerate(perm):G[i,j]=1
# fundamental-weight coordinate lattice, source roots columns C H
HS=C*s.Matrix.hstack(*map(s.Matrix,H));o=s.eye(6)[:,0]
def in_source(v):
 try:sol=HS.gauss_jordan_solve(v)[0]
 except ValueError:return False
 return all(x.q==1 for x in sol)
# W defined in root coordinates; action on weights C*W*C^-1
assert in_source(C*W*C.inv()*o-o)
assert in_source(C*G*C.inv()*o+o)
out={'roots':n,'root_index':abs(int(T.det())),'character_determinant':abs(int((T.T*C).det())),'a_correction':av,'b_correction':bv,'A_order':int(A.order()),'B_order':int(B.order()),'tau_order':int(Tau.order()),'Q12_order':int(Q.order()),'four_frames_centralize_Q12':True,'component_action_a_plus_b_minus':True,'simple_tits_squares_braid_graph':True,'scope':'Signed root-vector faithful adjoint action; all tested discrepancies lie in 2-kernel, which injects through E6 central kernel C3. Thus these relation checks lift to simply connected Tits/graph carriers.'}
Path('outputs/e6_carrier_round3_checks.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
