"""Node55 independent E8 order-five/order-six geometry.

Build all Euclidean roots and recover a new simple system from a generic
height. No root/source/normalizer arrays from the root manuscripts are read.
"""
from pathlib import Path
from itertools import combinations, product
from functools import reduce
from math import gcd
import json
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

OUT=Path(__file__).parent
E=S.eye(8); eu=[]
for i,j in combinations(range(8),2):
 for a,b in product((-1,1),repeat=2):
  r=S.zeros(8,1);r[i]=a;r[j]=b;eu.append(r)
for signs in product((-1,1),repeat=8):
 if sum(x<0 for x in signs)%2==0:eu.append(S.Matrix(signs)/2)
assert len(eu)==240
height=S.Matrix([3**i for i in range(8)])
pos=[r for r in eu if (r.T*height)[0]>0];ps=set(map(tuple,pos))
eb=[r for r in pos if not any(tuple(r-t) in ps for t in pos)]
assert len(eb)==8
B=S.Matrix.hstack(*eb);Bi=B.inv();C=B.T*B;I=S.eye(8)
assert C.det()==1
roots=[Bi*r for r in eu];rootset=set(map(tuple,roots))
dot=lambda a,b:(a.T*C*b)[0]
ref=lambda r:I-r*r.T*C
def base(rs):
 p=[r for r in rs if all(z>=0 for z in r)];sp=set(map(tuple,p))
 return [r for r in p if not any(tuple(r-t) in sp for t in p)]
def comps(bs):
 todo=set(range(len(bs)));out=[]
 while todo:
  cc={todo.pop()}
  while True:
   new={j for i in cc for j in todo if dot(bs[i],bs[j])}
   if not new:break
   cc|=new;todo-=new
  out.append(sorted(cc))
 return out
def chain(bs,cc):
 if len(cc)==1:return [bs[cc[0]]]
 ad={i:[j for j in cc if i!=j and dot(bs[i],bs[j])] for i in cc}
 seq=[min(i for i in cc if len(ad[i])==1)]
 while len(seq)<len(cc):seq.append(next(j for j in ad[seq[-1]] if j not in seq))
 return [bs[i] for i in seq]
def mul(ms):
 w=I
 for m in ms:w=w*m
 return w
def longest(bs):
 w=I
 while any(all(z>=0 for z in w*b) for b in bs):
  b=next(b for b in bs if all(z>=0 for z in w*b));w=w*ref(b)
 return w
def based(w,rs,bs):
 p={tuple(r) for r in rs if all(z>=0 for z in r)}
 while any(tuple(w*b) not in p for b in bs):
  b=next(b for b in bs if tuple(w*b) not in p);w=w*ref(b)
 return w
def snf(bs):
 if not bs:return []
 return [abs(int(z)) for z in smith_normal_form(S.Matrix.hstack(*bs),domain=ZZ).diagonal() if z]

highest=max(roots,key=lambda r:sum(r))
points={}
for n in (5,6):
 js=[i for i in range(8) if highest[i]==n];assert len(js)==1
 ss=C.inv()[:,js[0]]/n
 cs=[r for r in roots if dot(r,ss).q==1];bs=base(cs)
 sizes=sorted(map(len,comps(bs)))
 assert sizes==([4,4] if n==5 else [1,2,5])
 points[n]=(ss,cs,bs)

settings=[]
s5,cr5,cb5=points[5]; a5,b5=[chain(cb5,cc) for cc in comps(cb5)]
CB=S.Matrix.hstack(*(a5+b5));assert abs(CB.det())==5
R=S.zeros(4)
for i in range(4):R[i,3-i]=1
gam=CB*S.BlockMatrix([[S.zeros(4),S.eye(4)],[R,S.zeros(4)]]).as_explicit()*CB.inv()
assert gam**4==I and all(x.q==1 for x in gam)
assert all(tuple(gam*r) in rootset for r in roots)
res={k:next(q for q in range(1,5) if all(z.q==1 for z in q*gam**k*s5-s5)) for k in range(4)}
w0a=longest(a5);w0b=longest(b5);assert gam**2==-w0a*w0b
v3=mul([ref(r) for r in a5[:2]+b5[:2]])
v4=mul([ref(r) for r in a5[:3]+b5[:3]])
v5=mul([ref(r) for r in a5+b5])
settings += [('five_split_matching',5,I,1,[],1),('five_unitary_matching',5,-I,2,[],4)]
for sign,d,qres in [(1,1,4),(-1,2,1)]:
 for label,core in [('regular',[]),('one_core',a5[1:3]),('two_core',a5[1:3]+b5[1:3])]:
  settings.append(('five_nonmatching_'+str(d)+'_'+label,5,sign*gam**2,d,core,qres))
for label,v,d in [('cubic',v3,3),('quartic',v4,4),('quintic',v5,5)]:
 settings.append(('five_split_'+label,5,v,d,[],1))
 settings.append(('five_unitary_'+label,5,-v,6 if d==3 else 10 if d==5 else 4,[],4))
for k in (1,3):
 for d in (1,2):
  settings.append(('five_exchange_'+str(k)+'_d'+str(d),5,gam**k,d,[],res[k]))
  settings.append(('five_exchange_'+str(k)+'_d'+str(d)+'_core',5,gam**k,d,a5[1:3]+b5[1:3],res[k]))
 vv=w0a*gam
 assert vv**2==-I
 settings.append(('five_exchange_'+str(k)+'_quartic',5,vv if k==1 else vv.inv(),4,[],res[k]))

s6,cr6,cb6=points[6];cc6=comps(cb6)
a6=chain(cb6,next(cc for cc in cc6 if len(cc)==5))
b6=chain(cb6,next(cc for cc in cc6 if len(cc)==2))
c6=chain(cb6,next(cc for cc in cc6 if len(cc)==1))
assert abs(S.Matrix.hstack(*(a6+b6+c6)).det())==6
sig6=-longest(cb6)
assert sig6**2==I
settings += [('six_split_matching',6,I,1,[],1),('six_unitary_matching',6,-I,2,[],5)]
for sign,d,qres in [(1,1,5),(-1,2,1)]:
 for label,core in [('regular',[]),('A2_core',b6),('A5_core',a6),('A5A2_core',a6+b6)]:
  settings.append(('six_nonmatching_'+str(d)+'_'+label,6,sign*sig6,d,core,qres))
v6=mul([ref(r) for r in a6[:2]+a6[3:5]+b6])
assert v6**3==I
settings += [('six_cubic',6,v6,3,[],1),('six_sextic',6,-v6,6,[],5)]

ri={tuple(r):i+1 for i,r in enumerate(roots)}
gapperm=lambda v:'PermList('+str([ri[tuple(v*r)] for r in roots])+')'
gap=['refs:=['+','.join(gapperm(ref(r)) for r in roots)+'];;']
out={'Cartan':C.tolist(),'Euclidean_simple_roots':[list(r) for r in eb],
 's5':list(s5),'s6':list(s6),'gamma5':gam.tolist(),'gamma5_residues':res,'sources':{}}
x=S.Symbol('x')
for name,n,v,d,core,qres in settings:
 ss,cr,cb=points[n]
 assert v.T*C*v==C and all(z.q==1 for z in v)
 assert all(z.q==1 for z in qres*v*ss-ss)
 cp=S.Poly(S.cyclotomic_poly(d,x),x);mat=S.zeros(8)
 for (e,),c in cp.terms():mat+=c*v**e
 con=mat.col_join(S.Matrix.vstack(*[r.T*C for r in core])) if core else mat
 az=con.nullspace();AZ=S.Matrix.hstack(*az) if az else S.zeros(8,0)
 lr=[r for r in roots if all(dot(r,t)==0 for t in az)];lb=base(lr)
 kr=[r for r in lr if dot(r,ss).q==1];kb=base(kr)
 jr=[r for r in cr if all(dot(r,k)==0 for k in kb)];jb=base(jr)
 vl=based(v,lr,lb);fg=[next(i for i,t in enumerate(lb) if vl*r==t) for r in lb]
 vj=based(v,jr,jb);jfg=[next(i for i,t in enumerate(jb) if vj*r==t) for r in jb]
 poly=str(AZ.gauss_jordan_solve(v*AZ)[0].charpoly().as_expr().factor()) if az else '1'
 row={'n':n,'d':d,'q_residue':qres,'F':v.tolist(),'active_dimension':len(az),
 'active_basis':[list(t) for t in az],'active_polynomial':poly,
 'L_roots':len(lr),'L_base':[list(r) for r in lb],'L_components':comps(lb),
 'L_smith':snf(lb),'L_graph':fg,'K_roots':len(kr),'K_base':[list(r) for r in kb],
 'J_roots':len(jr),'J_base':[list(r) for r in jb],'J_graph':jfg}
 assert all(z==1 for z in row['L_smith'])
 if d in (3,4,6) and name not in ('five_exchange_1_quartic','five_exchange_3_quartic'):
  if d==3:aa=v*v-v;bb=I-v
  elif d==6:aa=v*v+v;bb=-v-I
  else:aa=v**3-v;bb=I-v*v
  norms=[]
  for r in roots:
   av=r.T*C*aa;bv=r.T*C*bb
   vals=[abs(int(a*a-a*b+b*b)) if d==3 else abs(int(a*a+a*b+b*b)) if d==6 else abs(int(a*a+b*b)) for a,b in zip(av,bv)]
   norms.append(reduce(gcd,vals,0))
  row['norm_gcds']=sorted(set(norms))
  lrset=set(map(tuple,lr))
  assert all((nrm==0)==(tuple(r) in lrset) for r,nrm in zip(roots,norms))
  row['norm_zero_iff_source_root']=True
 out['sources'][name]=row
 print(name,'L',len(lr),'K',len(kr),'J',len(jr),'F_L',fg,'Z',poly,'norms',row.get('norm_gcds'))
 if name in ('five_split_cubic','five_split_quartic','six_cubic','six_sextic'):
  cbi=[ri[tuple(r)] for r in cb];li=sorted(ri[tuple(r)] for r in lr);kbi=[ri[tuple(r)] for r in kb]
  gap+=['WC:=Group(List('+str(cbi)+',i->refs[i]),());;',
   'NC:=Stabilizer(WC,'+str(li)+',OnSets);;',
   'WK:=Group(List('+str(kbi)+',i->refs[i]),());;',
   'vp:='+gapperm(v)+';; NFF:=Group(Concatenation(GeneratorsOfGroup(NC),[vp]),());;',
   'ph:=NaturalHomomorphismByNormalSubgroup(NFF,WK);;Q:=Image(ph,NC);;',
   'R:=Centralizer(Q,Image(ph,vp));;',
   'Print("'+name+'|",Size(WC),"|",Size(NC),"|",Size(WK),"|",Size(Q),"|",Size(R),"|",StructureDescription(R),"\\n");']

(OUT/'independent_geometry.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
(OUT/'normalizers.g').write_text('\n'.join(gap)+'\n')
