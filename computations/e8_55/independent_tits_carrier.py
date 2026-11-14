"""Independent signed-root carrier in node55's recovered E8 simple basis.

All Tits relations and root-sum brackets are checked. The generated GAP4
file checks the full finite pinning quotient and transports the actual
order-six cubic source to it. No parent Tits arrays are imported.
"""
from pathlib import Path
from itertools import product
import json
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

OUT=Path(__file__).parent
data=json.loads((OUT/'independent_geometry.json').read_text())
C=S.Matrix(data['Cartan']);E=S.eye(8)
sim=[]
for i in range(8):
 m=E.copy();m[i,:]-=C[i,:];sim.append(m)
roots={tuple(int(z) for z in E[:,i]) for i in range(8)};todo=list(roots)
while todo:
 r=todo.pop()
 for i in range(8):
  z=list(r);z[i]-=sum(int(C[i,j])*r[j] for j in range(8));z=tuple(z)
  if z not in roots:roots.add(z);todo.append(z)
roots=sorted(roots);ri={r:i for i,r in enumerate(roots)}
assert len(roots)==240
def eps(a,b):
 return (-1)**(sum(a[i]*b[i] for i in range(8))+sum(a[i]*int(C[i,j])*b[j] for i in range(8) for j in range(i)))
ID=tuple((i,1) for i in range(240))
def mul(a,b):return tuple((a[j][0],sg*a[j][1]) for j,sg in b)
def inv(a):
 ans=[None]*240
 for i,(j,sg) in enumerate(a):ans[j]=(i,sg)
 return tuple(ans)
def power(a,n):
 ans=ID
 for _ in range(n):ans=mul(ans,a)
 return ans
def tor(bits):
 return tuple((i,(-1)**sum(r[j]*int(C[j,k])*bits[k] for j in range(8) for k in range(8))) for i,r in enumerate(roots))
ns=[]
for i in range(8):
 al=tuple(int(z) for z in E[:,i]);act=[]
 for r in roots:
  m=sum(int(C[i,j])*r[j] for j in range(8));z=list(r);z[i]-=m
  act.append((ri[tuple(z)],eps(al,r) if abs(m)==1 else 1))
 ns.append(tuple(act))
def word(w):
 ans=ID
 for i in w:ans=mul(ans,ns[i])
 return ans
def reduced(mat):
 m=mat.copy();steps=[]
 while m!=E:
  i=next(i for i in range(8) if all(z<=0 for z in m[:,i]))
  m=m*sim[i];steps.append(i);assert len(steps)<=120
 return list(reversed(steps))
def root_tits(beta):
 b=S.Matrix(beta);seq=[]
 if all(z<=0 for z in b):b=-b
 while sum(b)>1:
  i=next(i for i in range(8) if (C*b)[i]>0 and b!=E[:,i])
  b=sim[i]*b;seq.append(i)
 i=next(i for i in range(8) if b==E[:,i]);w=word(seq)
 return mul(mul(w,ns[i]),inv(w)),[j+1 for j in seq],i+1
def perm(a):
 return [j+1 if sg==1 else j+241 for j,sg in a]+[j+241 if sg==1 else j+1 for j,sg in a]
def wp(mat):return [ri[tuple(mat*S.Matrix(r))]+1 for r in roots]
def wpart(a):return [j+1 for j,sg in a]
def base(rs):
 pos=[r for r in rs if all(z>=0 for z in r)];ps=set(map(tuple,pos))
 return [r for r in pos if not any(tuple(r-t) in ps for t in pos)]

def addv(a,b):
 z=dict(a)
 for k,v in b.items():z[k]=z.get(k,0)+v
 return {k:v for k,v in z.items() if v}
def scalev(a,c):return {k:c*v for k,v in a.items() if c*v}
def adroot(al,v):
 z={}
 for k,c in v.items():
  if k<240:
   be=roots[k]
   if tuple(-x for x in al)==be:
    z=addv(z,{240+i:-c*al[i] for i in range(8) if al[i]})
   else:
    ab=tuple(x+y for x,y in zip(al,be))
    if ab in ri:z=addv(z,{ri[ab]:c*eps(al,be)})
  else:
   val=-c*sum(int(C[k-240,j])*al[j] for j in range(8))
   if val:z=addv(z,{ri[al]:val})
 return z
def expad(al,v):
 d=adroot(al,v);dd=adroot(al,d)
 assert not adroot(al,dd)
 return addv(addv(v,d),scalev(dd,S.Rational(1,2)))

brackets=0;opposite_brackets=0;adjoint_products=0
for i in range(8):
 al=tuple(int(z) for z in E[:,i]);neg=tuple(-z for z in al)
 # Explicit integral Chevalley unipotent product on all 248 basis vectors.
 # The cocycle convention is [e_alpha,e_-alpha]=-h_alpha.
 for j in range(248):
  v=expad(al,expad(neg,expad(al,{j:1})))
  if j<240:
   k,sg=ns[i][j];assert v=={k:sg}
  else:
   want={240+k:int(sim[i][k,j-240]) for k in range(8) if sim[i][k,j-240]}
   assert v==want
  assert all(S.Rational(x).q==1 for x in v.values());adjoint_products+=1
 bits=[0]*8;bits[i]=1;assert power(ns[i],2)==tor(bits)
 for j in range(8):
  if i!=j and C[i,j]==0:assert word([i,j])==word([j,i])
  if C[i,j]==-1:assert word([i,j,i])==word([j,i,j])
 for r in roots:
  ir,sr=ns[i][ri[r]];rr=roots[ir]
  # The Cartan action is the corresponding root reflection. Together
  # with this signed opposite-root identity this checks [e_r,e_-r]=h_r.
  it,st=ns[i][ri[tuple(-x for x in r)]]
  assert roots[it]==tuple(-x for x in rr) and sr*st==1
  assert tuple(sim[i]*S.Matrix(r))==rr
  opposite_brackets+=1
  for t in roots:
   rt=tuple(x+y for x,y in zip(r,t))
   if rt not in ri:continue
   it,st=ns[i][ri[t]];tt=roots[it]
   irt,srt=ns[i][ri[rt]]
   assert tuple(x+y for x,y in zip(rr,tt))==roots[irt]
   assert eps(r,t)*srt==sr*st*eps(rr,tt);brackets+=1

hi,hj=next((i,j) for i in range(8) for j in range(i+1,8) if C[i,j]==-1)
H=[E[:,hi],E[:,hj]]
perp=[S.Matrix(r) for r in roots if all((S.Matrix(r).T*C*h)[0]==0 for h in H)]
K=base(perp);assert len(perp)==72 and len(K)==6
KB=S.Matrix.hstack(*K);assert abs(S.Matrix.hstack(*(H+K)).det())==3
Ksmith=[abs(int(z)) for z in smith_normal_form(KB,domain=ZZ).diagonal() if z]
assert Ksmith==[1]*6
nk=[];kw=[]
for b in K:
 n,seq,si=root_tits(b);nk.append(n);kw.append({'conjugator':seq,'simple':si})
 for h in H:assert n[ri[tuple(h)]]==(ri[tuple(h)],1)
n0=word(reduced(-E));nH=word([hi,hj,hi]);a=mul(n0,inv(nH))
assert power(n0,2)==ID and power(nH,2)==ID and power(a,2)==ID
assert a[ri[tuple(H[0])]]==(ri[tuple(H[1])],1)
assert a[ri[tuple(H[1])]]==(ri[tuple(H[0])],1)
assert all(mul(a,n)==mul(n,a) for n in nk)
c=ID
for n in nk:c=mul(c,n)
u=power(c,4);assert power(u,3)==ID and u!=ID
assert power(mul(a,u),6)==ID and power(mul(a,u),3)==a
um=S.Matrix.hstack(*[S.Matrix(roots[u[ri[tuple(E[:,i])]][0]]) for i in range(8)])
ux=KB.gauss_jordan_solve(um*KB)[0]
assert S.expand(ux.charpoly().as_expr())==S.expand((S.Symbol('lambda')**2+S.Symbol('lambda')+1)**3)
assert abs((ux-S.eye(6)).det())==27
out={'H_simple_indices':[hi+1,hj+1],'K_base':[list(r) for r in K],
 'X_cocharacter_smith':Ksmith,
 'joint_lattice_index':3,'K_tits_words':kw,'longest_word':[i+1 for i in reduced(-E)],
 'u_X_matrix':ux.tolist(),'u_minus_identity_det':int((ux-S.eye(6)).det()),
 'bracket_checks':brackets,'opposite_root_Cartan_checks':opposite_brackets,
 'explicit_248_Chevalley_products':adjoint_products,
 'a_order':2,'u_order':3,'au_order':6,
 'a_positive_H_graph':True,'a_centralizes_X_Tits':True}
(OUT/'independent_tits_carrier.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
print(out)
g=['nx:=List('+str([perm(n) for n in nk])+',PermList);; V:=Group(nx,());;',
 'AA:=Group(List(nx,x->x*x),());;',
 'TT:=Group(List('+str([perm(power(n,2)) for n in ns])+',PermList),());;',
 'u:=PermList('+str(perm(u))+');;a:=PermList('+str(perm(a))+');;',
 'Q0:=Centralizer(V,u);;Q:=Group(Concatenation(GeneratorsOfGroup(Q0),[a]),());;',
 'Print("CARRIER|",Size(V),"|",Size(AA),"|",Size(TT),"|",Size(Q0),"|",Size(Q),"|",Size(Intersection(Q,TT)),"\\n");',
 'Print("RETURNS|",Order(u),"|",Order(a),"|",Order(a*u),"|",ForAll(GeneratorsOfGroup(Q),x->x*u=u*x and x*a=a*x),"\\n");',
 'W:=Group(List('+str([wpart(n) for n in ns])+',PermList));;',
 'WX:=Group(List('+str([wpart(n) for n in nk])+',PermList));;',
 'up:=PermList('+str(wpart(u))+');;',
 'Print("WEYL|",Size(WX),"|",Size(Centralizer(WX,up)),"\\n");']
src=data['sources']['six_cubic'];lbs=[S.Matrix(r) for r in src['L_base']]
vsrc=S.Matrix(src['F']);assert all(vsrc*r==r for r in lbs)
g+=['lbase:='+str([ri[tuple(r)]+1 for r in lbs])+';;hbase:='+str([ri[tuple(r)]+1 for r in H])+';;',
 'tr:=RepresentativeAction(W,lbase,hbase,OnTuples);;',
 'vs:=PermList('+str(wp(vsrc))+');;',
 'vv:=vs^tr;; zz:=RepresentativeAction(WX,vv,up);;',
 'Print("ACTUAL_SOURCE_TRANSPORT|",tr<>fail,"|",vv in WX,"|",zz<>fail,"\\n");',
 'P0:=Stabilizer(W,hbase,OnTuples);;P1:=Stabilizer(W,Set(hbase),OnSets);;',
 'ap:=PermList('+str(wpart(a))+');;',
 'Print("FULL_PINNING_QUOTIENT|",Size(P0),"|",Size(P1),"|",Size(Centralizer(P1,up)),"|",Size(Centralizer(P1,ap*up)),"\\n");',
 'Print("TRANSPORT_ROOT_PERM|",List([1..240],i->i^tr),"\\n");',
 'Print("AUXILIARY_CONJUGATOR|",List([1..240],i->i^zz),"\\n");']
(OUT/'independent_tits_carrier.g').write_text('\n'.join(g)+'\n')
