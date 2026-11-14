"""Uniform E8 source model for Table8 E7A1 e=4 at ell=5.
Independent enumeration uses all D4 signed permutations and pinning automorphisms.
No finite group of a sampled field size is used.
"""
from itertools import product, permutations, combinations
from collections import Counter
from pathlib import Path
import json
import sympy as S

I4=S.eye(4)
def key(m): return tuple(m)
def M(k): return S.Matrix(4,4,k)
def act(m,v): return tuple(m*S.Matrix(v))
def rootline(v):
 v=tuple(v)
 return min(v,tuple(-x for x in v))
roots4=set()
for i,j in combinations(range(4),2):
 for a,b in product((-1,1),repeat=2):
  r=[0]*4;r[i]=a;r[j]=b;roots4.add(tuple(r))
W0={}
for p in permutations(range(4)):
 for sig in product((-1,1),repeat=4):
  if S.prod(sig)!=1: continue
  a=S.zeros(4)
  for j in range(4): a[p[j],j]=sig[j]
  W0[key(a)]=a
assert len(W0)==192
# Bourbaki D4: central root is position1, outer roots positions0,2,3.
B=S.Matrix([[1,0,0,0],[-1,1,0,0],[0,-1,1,1],[0,0,-1,1]])
outer=[]
for p in permutations((0,2,3)):
 perm=(p[0],1,p[1],p[2]);P=S.zeros(4)
 for j in range(4):P[perm[j],j]=1
 g=B*P*B.inv()
 assert {act(g,r) for r in roots4}==roots4
 outer.append(g)
Aut={}
for u in W0.values():
 for g in outer:
  a=u*g
  assert key(a) not in Aut
  Aut[key(a)]=(a,g)
assert len(Aut)==1152
J=S.Matrix([[0,-1],[1,0]])
w=S.diag(J,J)
assert key(w) in W0 and w*w==-I4
ennola_h=S.diag(1,-1,1,-1)
assert key(ennola_h) in W0 and ennola_h*w*ennola_h.inv()==-w
C=[(a,g) for a,g in Aut.values() if a*w==w*a]
C0=[a for a in W0.values() if a*w==w*a]
assert len(C)==96 and len(C0)==16
assert len({key(g) for a,g in C})==6
# Full E8 roots in doubled Euclidean coordinates.
roots8=set()
for i,j in combinations(range(8),2):
 for a,b in product((-2,2),repeat=2):
  r=[0]*8;r[i]=a;r[j]=b;roots8.add(tuple(r))
for r in product((-1,1),repeat=8):
 if sum(x<0 for x in r)%2==0:roots8.add(r)
assert len(roots8)==240
alpha=S.Matrix([1,-1,0,0])
svec=alpha.col_join(S.zeros(4,1))
Cs=[r for r in roots8 if (S.Matrix(1,8,r)*svec)[0]%4==0]
assert len(Cs)==128
basisA1=[S.Matrix(x) for x in [(1,-1,0,0),(1,1,0,0),(0,0,1,-1),(0,0,1,1)]]
assert {rootline(r) for r in roots4 if (S.Matrix(1,4,r)*alpha)[0]%2==0}=={rootline(r) for r in basisA1}
# For each actual centralizer element correct its source diagram part to fix alpha.
# Any two corrections differ by W(A1^4), acting trivially on the four unipotent bits.
records=[]
for a,g in C:
 v=next(u for u in W0.values() if u*g*alpha==alpha)
 d=v*g
 lift=S.diag(d,a)
 assert {tuple(lift*S.Matrix(r)) for r in roots8}==roots8
 assert lift*S.diag(I4,w)==S.diag(I4,w)*lift
 assert lift*svec==svec
 perm=tuple(next(i for i,b in enumerate(basisA1) if rootline(d*x)==rootline(b)) for x in basisA1)
 assert perm[0]==0
 records.append((a,perm,lift))
assert len({p for a,p,l in records})==6
assert Counter(p for a,p,l in records)==Counter({(0,)+p:16 for p in permutations((1,2,3))})
# All four A1 factors are individually fixed by F=w F_q.
assert all(S.diag(I4,w)*x.col_join(S.zeros(4,1))==x.col_join(S.zeros(4,1)) for x in basisA1)
# Bit orbits and actual stabilizer sizes.
def perm_bits(p,b):
 out=[0]*4
 for j in range(4):out[p[j]]=b[j]
 return tuple(out)
bits=set(product((0,1),repeat=4));orbits=[]
while bits:
 b=min(bits);orb={perm_bits(p,b) for a,p,l in records};stab=[a for a,p,l in records if perm_bits(p,b)==b]
 bits-=orb
 orbits.append({'representative':b,'orbit':sorted(orb),'size':len(orb),'relative_stabilizer_order':len(stab)})
assert Counter((len(x['orbit']),x['relative_stabilizer_order']) for x in orbits)==Counter({(1,96):4,(3,32):4})
# Principal5-torsion of the center is the F-fixed kernel on the FULL D4 lattice.
wY=B.inv()*w*B
assert all(x.q==1 for x in wY)
mod5=[]
for qres in (2,3):
 tor=[S.Matrix(z) for z in product(range(5),repeat=4) if all(x%5==0 for x in (qres*wY-I4)*S.Matrix(z))]
 assert len(tor)==25
 vanished=[]
 for r in roots8:
  functional=S.Matrix(1,4,r[4:])*B/2
  assert all(x.q==1 for x in functional)
  if all((functional*z)[0]%5==0 for z in tor): vanished.append(r)
 assert len(vanished)==24 and all(not any(r[4:]) for r in vanished)
 mod5.append({'q_mod5':qres,'fixed_torsion_size':25,'centralizing_roots':24})
out={'model':'E8 Euclidean integral lattice, source first D4, complementary D4 Frobenius',
 'source_D4_roots':24,'ambient_E8_roots':240,'involution_cocharacter':list(svec),
 'involution_centralizer_roots':128,'source_centralizer_A1_roots':[list(b) for b in basisA1],
 'center_lattice_basis':B.tolist(),'Frobenius_matrix':w.tolist(),
 'Frobenius_characteristic_polynomial':str(w.charpoly().as_expr()),
 'source_Frobenius_identity':True,'Ennola_conjugator':ennola_h.tolist(),'relative_group_order':96,'relative_D4_kernel_order':16,
 'relative_action_on_four_A1_factors':sorted({p for a,p,l in records}),
 'A1_bit_orbits':orbits,'center_five_torsion':mod5,
 'scope':'All odd q with ord_5(q)=4, all v5(q^2+1); no finite field sample'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2,default=int)+'\n')
print(json.dumps({'passed':True,'relative_order':96,'kernel':16,'bit_orbits':[(x['size'],x['relative_stabilizer_order']) for x in orbits],'centralizers_mod5':mod5}))
