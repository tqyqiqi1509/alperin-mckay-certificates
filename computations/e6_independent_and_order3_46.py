"""Independent root-closure and formal valuation checks for node46.
No node43 data or code is imported. New order3 source geometry is included.
"""
from pathlib import Path
from collections import deque
from itertools import product
import json, math
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
base=Path(__file__).parent
C=2*S.eye(6)
for i,j in [(0,2),(2,3),(3,4),(4,5),(1,3)]:C[i,j]=C[j,i]=-1
I=S.eye(6);ref=lambda r:I-r*r.T*C;pair=lambda r,t:(r.T*C*t)[0]
gens=[ref(I[:,i]) for i in range(6)]
seen={tuple(I[:,0])};queue=deque(seen)
while queue:
 r=S.Matrix(queue.popleft())
 for w in gens:
  z=tuple(w*r)
  if z not in seen:seen.add(z);queue.append(z)
roots=[S.Matrix(r) for r in sorted(seen)]
assert len(roots)==72
gamma=S.zeros(6)
for i,j in enumerate([5,1,4,3,2,0]):gamma[j,i]=1
w0=-gamma
theta=max(roots,key=lambda r:sum(r));assert list(theta)==[1,2,2,3,2,1]
def rankmod(M,p):
 a=[[int(t)%p for t in M.row(i)] for i in range(M.rows)];row=0
 for col in range(M.cols):
  k=next((i for i in range(row,len(a)) if a[i][col]),None)
  if k is None:continue
  a[row],a[k]=a[k],a[row];inv=pow(a[row][col],-1,p)
  a[row]=[(inv*t)%p for t in a[row]]
  for j in range(len(a)):
   if j!=row:
    t=a[j][col];a[j]=[(x-t*y)%p for x,y in zip(a[j],a[row])]
  row+=1
 return row
def simples(rs):
 pos=[r for r in rs if all(t>=0 for t in r)];ps=set(map(tuple,pos))
 return [r for r in pos if not any(tuple(r-t) in ps for t in pos)]
def comps(bs):
 todo=set(range(len(bs)));ans=[]
 while todo:
  co={todo.pop()}
  while True:
   ns={j for i in co for j in todo if pair(bs[i],bs[j])}
   if not ns:break
   todo-=ns;co|=ns
  ans.append(sorted(co))
 return ans
z=3*C.inv()[:,0]
assert all(t.q==1 for t in z) and rankmod(z,3)==1
u=gens[0]*gens[2];v=gens[4]*gens[5]
assert u**3==v**3==I and u*v==v*u
comm=S.Matrix.hstack(u-I,v-I)
assert rankmod(comm,3)==rankmod(comm.row_join(z),3)
bad3=[]
for a,b in product(range(3),repeat=2):
 if not (a or b):continue
 m=(u**a*v**b-I)**2
 assert rankmod(m,3)==1 and rankmod(m.row_join(z),3)==2
 bad3.append({'powers':[a,b],'square_rank':1,'with_central_vector':2})
# A5 model reconstruction gives a separate interpretation of these ranks:
# in sum-zero F3^6, the two 3-cycles yield constant vectors on each
# triple; their product square has image (1,1,1,-1,-1,-1), never Z.
E=S.eye(6);P1=S.zeros(6);P2=S.zeros(6)
for i in range(6):P1[[1,2,0,3,4,5][i],i]=1;P2[[0,1,2,4,5,3][i],i]=1
Aug=S.Matrix.hstack(*[E[:,i]-E[:,5] for i in range(5)])
for a,b in product(range(3),repeat=2):
 if a or b:
  M=(P1**a*P2**b-I)**2*Aug
  assert rankmod(M,3)==1 and rankmod(M.row_join(S.ones(6,1)),3)==2

# Formal centralizer valuations for every orbit shape on six dimensions.
def val3(n):
 v=0
 while n and n%3==0:n//=3;v+=1
 return v
def parts(n,m=None):
 if n==0:yield ();return
 if m is None:m=n
 for k in range(min(n,m),0,-1):
  for rest in parts(n-k,k):yield (k,)+rest
valuation=[]
for m in [1,2]:
 rem=6-3*m
 for tail in parts(rem):
  # GL_m(q^3) plus scalar rational eigenspaces of multiplicities tail.
  ca=m+sum(tail);cb=m+val3(math.factorial(m))+sum(val3(math.factorial(k)) for k in tail)
  loss_a=6-ca;loss_b=2-cb
  assert loss_a+loss_b>=2
  valuation.append({'orbit3_multiplicity':m,'rational_multiplicities':tail,
    'index_v3_coeff_a':loss_a,'index_v3_constant':loss_b})
assert 6-2==4 # Two distinct length-three orbits, index=4a.

# New order3 point and its actual A2^3 subsystem.
s=C.inv()[:,3]/3
assert all(t.q==1 for t in 3*s) and any(t.q!=1 for t in s)
cr=[r for r in roots if pair(r,s).q==1];cb=simples(cr)
assert len(cr)==18 and sorted(map(len,comps(cb)))==[2,2,2]
assert all(t.q==1 for t in w0*s+s)
Ftor=w0;Freg=gens[1]*w0
F6=(gens[0]*gens[2])*(gens[1]*ref(theta))*w0
cases={}
for name,F,d in [('e1_A1cube_regular',Freg,1),('e1_D4_phi21',Ftor,1),
                  ('e2_toral',Ftor,2),('e6_candidate',F6,6)]:
 poly={1:F-I,2:F+I,6:F*F-F+I}[d];V=poly.nullspace()
 lr=[r for r in roots if all(pair(r,t)==0 for t in V)];lb=simples(lr)
 kr=[r for r in lr if pair(r,s).q==1];kb=simples(kr)
 LB=S.Matrix.hstack(*lb) if lb else S.zeros(6,0)
 smith=[abs(int(t)) for t in smith_normal_form(LB,domain=ZZ).diagonal() if t] if lb else []
 derived=rankmod(LB,3)==rankmod(LB.row_join(3*s),3)
 jr=[r for r in cr if all(pair(r,k)==0 for k in kb)]
 cases[name]={'F':F.tolist(),'d':d,'F_characteristic':str(F.charpoly().as_expr().factor()),
  'active_torus_dimension':len(V),'active_torus_basis':[list(t) for t in V],
  'L_roots':len(lr),'L_simple':[list(r) for r in lb],'L_component_ranks':list(map(len,comps(lb))),
  'K_roots':len(kr),'K_simple':[list(r) for r in kb], 'K_component_ranks':list(map(len,comps(kb))),
  'L_center_smith':smith,'s_in_dual_L_derived':derived,'J_roots':len(jr),
  'J_components':list(map(len,comps(simples(jr))))}
# Root restriction certificates hold for every relevant prime, not samples.
for name,dat in cases.items():
 F=S.Matrix(dat['F']);d=dat['d'];lr={tuple(r) for r in roots if all(pair(r,S.Matrix(t))==0 for t in dat['active_torus_basis'])}
 if d in [1,2]:
  V=[S.Matrix(t) for t in dat['active_torus_basis']]
  rows=[]
  for r in roots:
   if tuple(r) in lr:continue
   vals=[pair(r,t) for t in V];den=S.ilcm(*[x.q for x in vals])
   h=math.gcd(*[int(x*den) for x in vals]);assert h and all(p in [2,3] for p in S.factorint(abs(h*int(den))))
   rows.append(abs(h))
  dat['root_restriction_gcd_values']=sorted(set(rows))
 else:
  rt=S.sqrt(-3);zeta=(1+rt)/2
  V=(F-zeta*I).nullspace();assert len(V)==2
  V=[t.applyfunc(S.simplify) for t in V]
  norms=[]
  for r in roots:
   vals=[S.expand(S.simplify(pair(r,t))) for t in V]
   coeff=[]
   for v in vals:
    a=S.simplify(S.re(v));b=S.simplify(S.im(v)/S.sqrt(3))
    assert a.is_Rational and b.is_Rational
    coeff.append((a-b,2*b))
   den=1
   for a,b in coeff:den=S.ilcm(den,a.q,b.q)
   ns=[int(den**2*(a*a+a*b+b*b)) for a,b in coeff]
   h=math.gcd(*ns)
   assert h and all(p in [2,3] for p in S.factorint(abs(h*int(den))))
   norms.append(h)
  dat['Eisenstein_gcd_norm_values']=sorted(set(norms))
  dat['Eisenstein_eigenbasis']=[[str(t) for t in v] for v in V]
# Every A2 factor sees the nontrivial order3 kernel of its SC product.
CB=S.Matrix.hstack(*cb);assert abs(CB.det())==3
kernel=CB.inv()*I[:,3]
assert all(any(t.q==3 for t in [kernel[i] for i in co]) for co in comps(cb))
out={'root_closure':72,'Cartan':C.tolist(),'z':list(z),'bad3_square_checks':bad3,
 'formal_nontoral_valuation':valuation,'order3_point':list(s),'Cs_simple':[list(r) for r in cb],
 'order3_cases':cases,'SC_A2cube_kernel_coordinates':list(kernel)}
base.joinpath('e6_independent_and_order3_46.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
print(json.dumps({k:{x:y[x] for x in ['F_characteristic','L_roots','K_roots','s_in_dual_L_derived']} for k,y in cases.items()}))
# Independent finite Weyl images, using root closure rather than imported lists.
keys=[tuple(r) for r in roots];index={r:i+1 for i,r in enumerate(keys)}
def gp(w):return 'PermList('+str([index[tuple(w*r)] for r in roots])+')'
g=[]
g+=['W:=Group(['+','.join(gp(w) for w in gens)+']);;']
g+=['R:=Group(['+','.join(gp(ref(r)) for r in cb)+']);;']
g+=['Print("ROOT_CLOSURE|W=",Size(W),"|R=",Size(R),"\\n");']
# The radical of the E6 root form modulo3 is the natural central line.
# Independently build the five-dimensional quotient and certify faithfulness.
B=S.Matrix.hstack(z,*[I[:,j] for j in range(1,6)])
def fp3(x):
 return (int(x.p)*pow(int(x.q),-1,3))%3
qgens=[]
for w in gens:
 m=B.inv()*w*B
 assert all(fp3(m[j,0])==(1 if j==0 else 0) for j in range(6))
 qgens.append([[fp3(m[i,j]) for j in range(1,6)] for i in range(1,6)])
g+=['Q:=Group(['+','.join(str(m)+'*One(GF(3))' for m in qgens)+']);;']
g+=['Print("QUOTIENT_MOD3|dimension=5|image=",Size(Q),"\\n");']
g+=['if Size(Q)<>51840 then Error("Nonfaithful quotient"); fi;']
for name,dat in cases.items():
 F=S.Matrix(dat['F']);bs=[S.Matrix(r) for r in dat['L_simple']]
 g+=['f:='+gp(F)+';; Cw:=Centralizer(W,f);; Cr:=Centralizer(R,f);;']
 g+=['L:=Group(['+','.join(gp(ref(r)) for r in bs)+'],());;']
 g+=['Print("CASE|'+name+'|CW=",Size(Cw),"|CR=",Size(Cr),"|WL=",Size(L),"|inter=",Size(Intersection(L,Cr)),"\\n");']
 g+=['K:=Intersection(L,Cr);; Print("RELATIVE|'+name+'|",StructureDescription(Cr/K),"\\n");']
g+=['QUIT;']
base.joinpath('e6_independent_and_order3_46.g').write_text('\n'.join(g)+'\n')
