"""Independent E7 integral source geometry, not a finite-q sample.

Generate roots first, recover a simple system from a generic height,
then discover isolated vertices by their actual centralizing roots.
No previous source JSON or source-root list is imported.
"""
from pathlib import Path
from itertools import combinations, product, permutations
from functools import reduce
from math import gcd, lcm
import json
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

I8=S.eye(8)
roots=[]
for i,j in combinations(range(8),2):
 for a,b in product((-1,1),repeat=2):
  r=S.zeros(8,1);r[i]=a;r[j]=b
  if r[6]+r[7]==0:roots.append(r)
for signs in product((-1,1),repeat=8):
 if sum(x<0 for x in signs)%2==0 and signs[6]+signs[7]==0:
  roots.append(S.Matrix(signs)/2)
assert len(roots)==126
height=S.Matrix([1,4,16,64,256,1024,4096,-4096])
dot=lambda a,b:(a.T*b)[0]
pos=[r for r in roots if dot(height,r)>0]
assert len(pos)==63
def simples(rs):
 pp=[r for r in rs if dot(height,r)>0]; keys=set(map(tuple,pp))
 return [r for r in pp if not any(tuple(r-t) in keys for t in pp)]
def components(bs):
 rem=set(range(len(bs)));out=[]
 while rem:
  co={rem.pop()}
  while True:
   more={j for j in rem if any(dot(bs[i],bs[j]) for i in co)}
   if not more:break
   co|=more;rem-=more
  out.append(sorted(co))
 return out
base=simples(roots);assert len(base)==7
B=S.Matrix.hstack(*base);C=B.T*B;Ci=C.inv(); left=Ci*B.T
assert C.det()==2
coords=[left*r for r in roots]
assert all(x.q==1 for r in coords for x in r)
assert sorted(map(len,components(base)))==[7]
I=S.eye(7)
ref=lambda r:I8-r*r.T
orderQ=lambda v:lcm(*[x.q for x in left*v])
orderP=lambda v:lcm(*[dot(r,v).q for r in base])
highest=max(pos,key=lambda r:sum(left*r));marks=left*highest
vertices=[]
for j,m in enumerate(marks):
 p=B*Ci[:,j]/m
 for k in range(1,int(m)+1):
  s=k*p;cr=[r for r in roots if dot(r,s).q==1]
  vertices.append((len(cr),orderQ(s),s))
point3=next(s for n,o,s in vertices if n==36 and o==3)
point2=next(s for n,o,s in vertices if n==62 and o==2)
# The two SC central translates have the same root centralizer. Select
# the trace +8 anchor on the independently generated minuscule orbit.
mj=next(j for j,m in enumerate(marks) if m==1)
minuscule=B*Ci[:,mj]
weights={tuple(minuscule)};todo=[minuscule]
while todo:
 w=todo.pop()
 for r in base:
  v=w-dot(w,r)*r
  if tuple(v) not in weights:weights.add(tuple(v));todo.append(v)
assert len(weights)==56
trace=lambda s:sum(1 if dot(S.Matrix(w),s).q==1 else -1 for w in weights)
assert trace(point2) in (-8,8)
if trace(point2)==-8:point2+=minuscule
assert trace(point2)==8 and orderQ(point2)==2

def chain(bs,inds):
 adj={i:[j for j in inds if i!=j and dot(bs[i],bs[j])!=0] for i in inds}
 ends=[i for i in inds if len(adj[i])==1];assert len(ends)==2
 seq=[min(ends)]
 while len(seq)<len(inds):seq.append(next(j for j in adj[seq[-1]] if j not in seq))
 return [bs[i] for i in seq]

def mod_span(M,b,p):
 a=[[int(x)%p for x in M.row(i)]+[int(b[i])%p] for i in range(M.rows)]
 rr=0;piv=[]
 for j in range(M.cols):
  k=next((k for k in range(rr,M.rows) if a[k][j]),None)
  if k is None:continue
  a[k],a[rr]=a[rr],a[k];iv=pow(a[rr][j],-1,p)
  a[rr]=[(x*iv)%p for x in a[rr]]
  for k in range(M.rows):
   if k!=rr:
    aa=a[k][j];a[k]=[(x-aa*y)%p for x,y in zip(a[k],a[rr])]
  piv.append(j);rr+=1
 ok=not any(all(x==0 for x in row[:-1]) and row[-1] for row in a)
 ans=[0]*M.cols
 if ok:
  for k,j in enumerate(piv):ans[j]=a[k][-1]
 return ok,ans

cases={}
def examine(name,s,F,d,modulus):
 Fc=left*F*B;assert all(x.q==1 for x in Fc)
 poly={1:Fc-I,2:Fc+I,3:Fc**2+Fc+I,4:Fc**2+I,6:Fc**2-Fc+I}[d]
 active=poly.nullspace();LR=[r for r in roots if all(dot(r,B*v)==0 for v in active)]
 LB=simples(LR);LC=left*S.Matrix.hstack(*LB) if LB else S.zeros(7,0)
 CR=[r for r in roots if dot(r,s).q==1]
 KR=[r for r in LR if dot(r,s).q==1];KB=simples(KR)
 JR=[r for r in CR if all(dot(r,k)==0 for k in KB)];JB=simples(JR)
 smith=[abs(int(t)) for t in smith_normal_form(LC,domain=ZZ).diagonal() if t] if LB else []
 inderived,witness=mod_span(LC,modulus*left*s,modulus) if LB else (False,[])
 norms=[]
 for rc in coords:
  if d in (1,2):
   row=rc.T*C*(I+Fc if d==1 else I-Fc)
   n=reduce(gcd,[abs(int(t)) for t in row],0)
  else:
   # Kill both inactive signs first; harmless extra unit at d=3.
   aa=rc.T*C*(Fc**3-Fc);bb=rc.T*C*(I-Fc**2)
   sign={3:-1,6:1,4:0}[d]
   vals=[int(a*a+sign*a*b+b*b) for a,b in zip(aa,bb)]
   n=reduce(gcd,vals,0)
  norms.append(n)
  assert (n==0)==(B*rc in LR)
  assert all(p in (2,3) for p in S.factorint(n)) if n else True
 # Correct F to a source positive system; retain its diagram, not just order.
 Fb=F;LP={tuple(r) for r in LR if dot(height,r)>0}
 while any(tuple(Fb*r) not in LP for r in LB):
  r=next(r for r in LB if tuple(Fb*r) not in LP);Fb=Fb*ref(r)
 diagram=[next(j for j,t in enumerate(LB) if Fb*r==t) for r in LB]
 cases[name]={'s_Euclidean':list(s),'F_on_Q':Fc.tolist(),'d':d,
 'active_rank':len(active),'L_root_count':len(LR),'L_ranks':list(map(len,components(LB))),
 'L_basis_Euclidean':[list(r) for r in LB],'L_basis_Q':LC.tolist(),
 'L_center_smith':smith,'derived_membership_modulus':modulus,'s_in_derived':inderived,
 'derived_membership_witness':witness,'K_root_count':len(KR),
 'K_ranks':list(map(len,components(KB))),'K_basis_Euclidean':[list(r) for r in KB],
 'J_identity_root_count':len(JR),'J_identity_ranks':list(map(len,components(JB))),
 'J_basis_Euclidean':[list(r) for r in JB],'source_F_diagram':diagram,
 'K_root_smith_in_SC_character_P':([abs(int(x)) for x in
   smith_normal_form(C*left*S.Matrix.hstack(*KB),domain=ZZ).diagonal() if x] if KB else []),
 'restriction_norm_gcds':sorted(set(norms))}
 print(name,cases[name]['L_ranks'],cases[name]['K_ranks'],cases[name]['J_identity_ranks'],
       'center',smith,'derived',inderived,'norms',sorted(set(norms)))

cr3=[r for r in roots if dot(r,point3).q==1];cb3=simples(cr3);co3=components(cb3)
a5=chain(cb3,next(c for c in co3 if len(c)==5))
a2=chain(cb3,next(c for c in co3 if len(c)==2))
w2=ref(a5[0])*ref(a5[2])*ref(a5[4])
w3=ref(a5[0])*ref(a5[1])*ref(a5[3])*ref(a5[4])*ref(a2[0])*ref(a2[1])
for eps in (1,-1):
 examine('A5A2_A1cube_'+str(eps),point3,eps*w2*ref(a2[0]),2 if eps==1 else 1,3)
 examine('A5A2_D4_'+str(eps),point3,eps*w2,2 if eps==1 else 1,3)
 examine('A5A2_cubic_'+str(eps),point3,eps*w3,3 if eps==1 else 6,3)

cr2=[r for r in roots if dot(r,point2).q==1];cb2=simples(cr2);co2=components(cb2)
dc=next(c for c in co2 if len(c)==6);D=S.Matrix.hstack(*[cb2[i] for i in dc])
I6=S.eye(6);std=S.Matrix.hstack(*[I6[:,i]-I6[:,i+1] for i in range(5)],I6[:,4]+I6[:,5])
per=next(p for p in permutations(range(6)) if std[:,p].T*std[:,p]==D.T*D)
T=D*std[:,per].inv();assert T.T*T==I6
P3=S.zeros(6);P4=S.eye(6)
for cyc in ((0,1,2),(3,4,5)):
 for i,j in zip(cyc,cyc[1:]+cyc[:1]):P3[j,i]=1
for i,j in ((0,1),(2,3)):
 P4[i,i]=P4[j,j]=0;P4[j,i]=1;P4[i,j]=-1
for d,P in ((3,P3),(6,-P3),(4,P4)):
 F=T*P*T.T+(I8-T*T.T)
 examine('D6A1_d'+str(d),point2,F,d,2)

# The old D4 marked source, independently recovered from its active
# three-dimensional center; the two d-marks remain distinct below.
PD4=S.diag(-1,-1,-1,-1,1,1)
FD4=T*PD4*T.T+(I8-T*T.T)
for eps in (1,-1):
 examine('D6A1_D4_'+str(eps),point2,eps*FD4,1 if eps==1 else 2,2)

# Isolated centralizer central products in the SC lattice: the finite
# gluing, not merely the abstract type or order, is recorded.
central_products={}
for name,cb,co in (('A5A2',cb3,co3),('D6A1',cb2,co2)):
 ordered=sum(co,[]);joint_cp=left*S.Matrix.hstack(*[cb[i] for i in ordered])
 index=abs(int(joint_cp.det()));ji_cp=joint_cp.inv()
 residues={tuple(x-S.floor(x) for x in ji_cp*S.Matrix(v))
           for v in product(range(index),repeat=7)}
 cut=len(co[0]);projections=[{g[:cut] for g in residues},{g[cut:] for g in residues}]
 assert len(residues)==index
 assert all(len(g)==index for g in projections)
 central_products[name]={'lattice_index':index,'component_ranks':list(map(len,co)),
  'glue':sorted(residues),'glue_projection_sizes':list(map(len,projections)),
  'factor_primitive_smith':[[abs(int(x)) for x in
    smith_normal_form(left*S.Matrix.hstack(*[cb[i] for i in c]),domain=ZZ).diagonal() if x]
    for c in co]}
 assert all(all(x==1 for x in ss) for ss in central_products[name]['factor_primitive_smith'])

# Actual central gluing for the quartic A1-cube, including its component.
dat=cases['D6A1_d4'];KB=[S.Matrix(r) for r in dat['K_basis_Euclidean']]
perp=[r for r in roots if all(dot(r,k)==0 for k in KB)];PB=simples(perp)
joint=left*S.Matrix.hstack(*KB,*PB);assert joint.shape==(7,7)
idx=abs(int(joint.det()));assert idx==4
Jinv=joint.inv();frac=lambda x:tuple(t-S.floor(t) for t in x)
glue={frac(Jinv*S.Matrix(v)) for v in product(range(2),repeat=7)}
assert len(glue)==4
Kglue={g[:3] for g in glue};Dglue={g[3:] for g in glue}
assert len(Kglue)==len(Dglue)==4
z=next(Ci[:,i] for i in range(7) if any(t.q==2 for t in Ci[:,i]))
zj=Jinv*z
assert frac(zj[:3,0]) not in Kglue
gluing={'joint_lattice_index':idx,'glue':sorted(glue),'K_projection_size':len(Kglue),
 'D4_projection_size':len(Dglue),'global_center_Q':list(z),'global_center_joint':list(zj),
 'global_center_outside_D4_component':True,
 'D4_primitive_smith':[abs(int(t)) for t in smith_normal_form(left*S.Matrix.hstack(*PB),domain=ZZ).diagonal() if t]}

# The reverse centralizer of the D4 factor has no extra component:
# its full D4 center is exhausted by the order-four gluing, leaving
# the embedded simply connected A1 cube. This is the full J for d1/2.
old=cases['D6A1_D4_1'];oldK=[S.Matrix(r) for r in old['K_basis_Euclidean']]
oldP=[r for r in roots if all(dot(r,k)==0 for k in oldK)];oldPB=simples(oldP)
oldjoint=left*S.Matrix.hstack(*oldK,*oldPB);oldji=oldjoint.inv()
oldglue={frac(oldji*S.Matrix(v)) for v in product(range(2),repeat=7)}
assert abs(oldjoint.det())==4 and len(oldglue)==4
assert len({g[:4] for g in oldglue})==4 and len({g[4:] for g in oldglue})==4
assert len(oldP)==6 and all(dot(r,point2).q==1 for r in oldP)
gluing['old_D4_reverse_full_J']={'joint_lattice_index':4,
 'D4_center_projection_size':4,'A1cube_projection_size':4,
 'full_J':'SL2^3 connected','global_center_in_A1cube':True}

# E7 adjoint quasi-isolated order six: independent standard-alcove model.
Std=S.Matrix.hstack(S.Matrix([1,-1,-1,-1,-1,-1,-1,1])/2,
 I8[:,0]+I8[:,1],*[I8[:,i]-I8[:,i-1] for i in range(1,6)])
CS=Std.T*Std;fund=Std*CS.inv();s6=(fund[:,2]+fund[:,4])/6
cr6=[r for r in roots if dot(r,s6).q==1];kb6=simples(cr6)
assert len(cr6)==18 and sorted(map(len,components(kb6)))==[2,2,2]
center=(S.Matrix.vstack(*[r.T for r in kb6])).col_join(S.Matrix([[0,0,0,0,0,0,1,1]])).nullspace()
assert len(center)==1
lr6=[r for r in roots if dot(r,center[0])==0];lb6=simples(lr6)
assert len(lr6)==72 and len(lb6)==6
L6=S.Matrix.hstack(*lb6);C6=L6.T*L6
s60=L6*C6.inv()*L6.T*s6;c60=C6.inv()*L6.T*s60
assert lcm(*[(C6*c60)[i].q for i in range(6)])==3
zc=center[0]/2
assert all(dot(r,zc).q==1 for r in base)
assert orderQ(zc)==2 and orderP(zc)==1
assert s6-s60==center[0]/12
order6={'point_Euclidean':list(s6),'point_in_independent_Q':list(left*s6),
 'adjoint_order':orderP(s6),'SC_lift_order':orderQ(s6),'centralizer_roots':18,
 'Cs_component_ranks':list(map(len,components(kb6))),
 'minimal_Levi_roots':len(lr6),'minimal_Levi_basis':[list(r) for r in lb6],
 'minimal_Levi_center_direction':list(center[0]),
 'global_E7_center_representative':list(zc),
 'global_E7_center_lies_in_connected_Levi_center':True,
 'central_torus_projection':list(s6-s60),
 'derived_projection_Euclidean':list(s60),'derived_projection_E6_Q':list(c60),
 'derived_projection_in_E6adj_order':lcm(*[t.q for t in C6*c60]),
 'derived_projection_in_E6sc_order':lcm(*[t.q for t in c60]),
 'derived_centralizer_roots':len([r for r in lr6 if dot(r,s60).q==1]),
 'primal_derived_embedding_smith':[abs(int(t)) for t in smith_normal_form(left*L6,domain=ZZ).diagonal() if t]}
assert order6['adjoint_order']==6
print('ORDER6',order6)
# Uniform nonabelian toral defect at ell=5, every valuation.  The
# calculation is on the seven-dimensional integral torus lattice.
w5a=ref(a5[0])*ref(a5[1])*ref(a5[2])*ref(a5[3])
P5=S.zeros(6);P5[5,5]=1
for i in range(5):P5[(i+1)%5,i]=1
w5d=T*P5*T.T+(I8-T*T.T)
def rankmod52(M,p):
 a=[[int(x)%p for x in M.row(i)] for i in range(M.rows)];rr=0
 for j in range(M.cols):
  k=next((k for k in range(rr,M.rows) if a[k][j]),None)
  if k is None:continue
  a[k],a[rr]=a[rr],a[k];iv=pow(a[rr][j],-1,p)
  a[rr]=[(x*iv)%p for x in a[rr]]
  for k in range(M.rows):
   if k!=rr:
    z=a[k][j];a[k]=[(x-z*y)%p for x,y in zip(a[k],a[rr])]
  rr+=1
 return rr
fives={}
for name,F in [('A5A2',w5a),('D6A1',w5d)]:
 M=left*F*B;assert M**5==I
 ranks=[rankmod52((M**k-I)**2,5) for k in range(1,5)]
 assert ranks==[3]*4
 fives[name]={'matrix_on_Q':M.tolist(),'ranks_square_g_minus_one_mod5':ranks}
out={'E7_simple_system_independently_recovered':[list(r) for r in base],
 'Cartan':C.tolist(),'highest_marks':list(marks),'cases':cases,
 'isolated_central_product_gluings':central_products,
 'toral_five_cycle_characteristic_test':fives,
 'quartic_full_J_gluing':gluing,'order6_E6T_carrier':order6}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
