"""Fresh Cartan/lattice reconstruction of the six matching E6 source geometries.

No arrays from node50 are imported. Frobenius classes are integral operators,
not chosen finite-field examples. The character arguments live in the note.
"""
from pathlib import Path
from itertools import product
from math import gcd
from functools import reduce
import json
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

outdir=Path(__file__).parent
C=2*S.eye(6)
for i,j in ((0,2),(2,3),(3,4),(4,5),(1,3)):C[i,j]=C[j,i]=-1
I=S.eye(6)
ref=lambda r:I-r*r.T*C
pair=lambda r,t:(r.T*C*t)[0]
key=lambda r:tuple(r)
roots={key(I[:,j]) for j in range(6)}
pending=list(roots)
simple=[ref(I[:,j]) for j in range(6)]
while pending:
 r=S.Matrix(pending.pop())
 for w in simple:
  u=key(w*r)
  if u not in roots:roots.add(u);pending.append(u)
roots=[S.Matrix(r) for r in sorted(roots)]
assert len(roots)==72 and C.det()==3
def based(rs):
 pos=[r for r in rs if all(v>=0 for v in r)]
 keys=set(map(key,pos))
 return [r for r in pos if not any(key(r-t) in keys for t in pos)]
def components(bs):
 todo=set(range(len(bs)));ret=[]
 while todo:
  cur={todo.pop()}
  while True:
   new={j for i in cur for j in todo if pair(bs[i],bs[j])}
   if not new:break
   cur|=new;todo-=new
  ret.append(sorted(cur))
 return ret
def smith(B):
 if not B.cols:return []
 return [abs(int(t)) for t in smith_normal_form(B,domain=ZZ).diagonal() if t]
def rankmod(M,p):
 rows=[[int(v)%p for v in M.row(i)] for i in range(M.rows)]
 rank=0
 for col in range(M.cols):
  piv=next((i for i in range(rank,M.rows) if rows[i][col]),None)
  if piv is None:continue
  rows[rank],rows[piv]=rows[piv],rows[rank]
  inv=pow(rows[rank][col],-1,p);rows[rank]=[(v*inv)%p for v in rows[rank]]
  for i in range(M.rows):
   if i!=rank:
    v=rows[i][col];rows[i]=[(a-v*b)%p for a,b in zip(rows[i],rows[rank])]
  rank+=1
 return rank

s=C.inv()[:,3]/3;z=C.inv()[:,0]
cr=[r for r in roots if pair(r,s).q==1];cb=based(cr);co=components(cb)
assert len(cr)==18 and sorted(map(len,co))==[2,2,2]
chains=[[cb[j] for j in x] for x in co]
Bs=S.Matrix.hstack(*sum(chains,[]))
cycs=[]
for flips in product((0,1),repeat=3):
 target=[]
 for j in range(3):target+=chains[(j+1)%3][::(-1 if flips[j] else 1)]
 w=S.Matrix.hstack(*target)*Bs.inv()
 if w**3==I and all(v.q==1 for v in w) and all(key(w*r) in set(map(key,roots)) for r in roots):
  shifts=[k for k in (1,2) if all(v.q==1 for v in w*s-s-k*z)]
  if shifts:cycs.append((w,shifts[0]))
assert cycs
c,shift=cycs[0]
theta=max((r for r in roots if all(v>=0 for v in r)),key=lambda r:sum(r))
w3=(simple[0]*simple[2])*(simple[4]*simple[5])*(simple[1]*ref(theta))
assert w3**3==I and c**3==I
cases=[('split_d1',I,1),('split_d2_regular',simple[0]*simple[4]*simple[1],2),
 ('split_d2_marked',simple[0]*simple[4],2),('split_d3',w3,3),
 ('cyclic_d1',c,1),('cyclic_d3',c,3)]
ri={key(r):i+1 for i,r in enumerate(roots)}
gp=lambda m:'PermList('+str([ri[key(m*r)] for r in roots])+')'
refperms=[gp(ref(r)) for r in roots]
gap=['refs:=['+','.join(refperms)+'];;',
 'W:=Group(List('+str([ri[key(I[:,j])] for j in range(6)])+',i->refs[i]));;',
 'WC:=Group(List('+str([ri[key(r)] for r in cb])+',i->refs[i]));;',
 'c:='+gp(c)+';; Wbar:=ClosureGroup(WC,c);;',
 'Print("AMBIENT|",Size(W),"|",Size(WC),"|",Size(Wbar),"\\n");']
result={'cartan':C.tolist(),'s':list(s),'z':list(z),'centralizer_base':Bs.tolist(),
        'centralizer_lattice_index':abs(Bs.det()),'cyclic_F':c.tolist(),'cyclic_shift':shift,'cases':{}}
for name,v,d in cases:
 if d==1: M=v-I
 elif d==2:M=v+I
 else:M=v*v+v+I
 tor=M.nullspace();lr=[r for r in roots if all(pair(r,t)==0 for t in tor)];lb=based(lr)
 LB=S.Matrix.hstack(*lb) if lb else S.zeros(6,0)
 kr=[r for r in lr if pair(r,s).q==1];kb=based(kr)
 per=[r for r in roots if all(pair(r,t)==0 for t in lb)] if lb else []
 pb=based(per);PB=S.Matrix.hstack(*pb) if pb else S.zeros(6,0)
 jrs=[r for r in cr if all(pair(r,t)==0 for t in kb)];jb=based(jrs)
 if d==3:
  P0=I-v;P1=v*v-v
  norms=[]
  for r in roots:
   a=r.T*C*P0;b=r.T*C*P1
   ns=[int(a[j]**2-a[j]*b[j]+b[j]**2) for j in range(6)]
   norms.append(reduce(gcd,ns))
 else:
  P=I+v+v*v if v**3==I and v!=I else (I-v if d==2 else I+v)
  norms=[reduce(gcd,[abs(int(t)) for t in r.T*C*P]) for r in roots]
 inderived=rankmod(LB,3)==rankmod(LB.row_join(3*s),3) if lb else False
 rec={'F':v.tolist(),'effective_d':d,'active_rank':len(tor),'L_roots':len(lr),
  'L_components':list(map(len,components(lb))),'L_base':LB.tolist(),
  'natural_center_smith':smith(C*LB),'root_lattice_smith':smith(LB),
  's_in_derived_mod_Q':inderived,'K_roots':len(kr),'K_base':[list(r) for r in kb],
  'J_roots':len(jrs),'J_components':list(map(len,components(jb))),
  'restriction_norm_gcds':sorted(set(norms)),
  'perpendicular_roots':len(per),'perpendicular_base':PB.tolist()}
 if lb:
  joint=LB.row_join(PB)
  rec['joint_lattice_index']=abs(joint.det()) if joint.cols==6 else None
  rec['joint_rank']=joint.rank()
  rec['perp_primitive_smith']=smith(PB)
  rec['global_z_in_perp']=rankmod(PB,3)==rankmod(PB.row_join(3*z),3)
 result['cases'][name]=rec
 li=[ri[key(r)] for r in lr];lbi=[ri[key(r)] for r in lb]
 gap+=['f:='+gp(v)+';;','lroots:='+str(li)+';;',
  'WL:=Group(List('+str(lbi)+',i->refs[i]),());;',
  'NW:=Stabilizer(W,lroots,OnSets);;',
  'if ForAll(GeneratorsOfGroup(NW),x->x^f*x^-1 in WL) then FW:=NW; else fix:=Filtered(Elements(NW),x->x^f*x^-1 in WL);FW:=Group(fix);fi;',
  'Q:=FactorGroup(FW,WL);;',
  'Print("'+name+'|FULL|",Size(Q),"\\n");']
 if lb:
  gap+=['Print("ACTION|'+name+'|FULL|",List(GeneratorsOfGroup(FW),g->List('+str([ri[key(I[:,j])] for j in range(6)])+',i->i^g)),"\\n");']
 for label,G in [('REG','WC'),('NAT','Wbar')]:
  gap+=['NC:=Stabilizer('+G+',lroots,OnSets);;K:=Intersection(NC,WL);;',
   'if ForAll(GeneratorsOfGroup(NC),x->x^f*x^-1 in K) then FC:=NC; else fix:=Filtered(Elements(NC),x->x^f*x^-1 in K);FC:=Group(fix);fi;Q:=FactorGroup(FC,K);;',
   'Print("'+name+'|'+label+'|",Size(Q),"|",Length(Irr(Q)),"|",StructureDescription(Q),"\\n");']
  if lb:
   gap+=['Print("ACTION|'+name+'|'+label+'|",List(GeneratorsOfGroup(FC),g->List('+str([ri[key(I[:,j])] for j in range(6)])+',i->i^g)),"\\n");']
 print(name,{k:rec[k] for k in ('L_roots','L_components','natural_center_smith','K_roots','J_roots','restriction_norm_gcds','s_in_derived_mod_Q')})
(outdir/'independent_geometry.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
(outdir/'independent_geometry.g').write_text('\n'.join(gap)+'\nQUIT;\n')
