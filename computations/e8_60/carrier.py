"""Fresh integral E8 Chevalley carrier for the actual order-three A2 source."""
from pathlib import Path
from itertools import combinations,product
from functools import reduce
import json
import sympy as S

DIR=Path(__file__).resolve().parent
geo=json.loads((DIR/'geometry.json').read_text())
B=S.Matrix(geo['lattice_basis']);H=S.Matrix(geo['gram'])
E=S.Matrix.hstack(*(S.Matrix(r) for r in geo['simple_basis']));C=E.T*H*E
EI=E.inv();R=[]
for i,j in combinations(range(8),2):
  for a,b in product((-1,1),repeat=2):
    z=[0]*8;z[i]=a;z[j]=b;R.append(EI*B.inv()*S.Matrix(z))
for ss in product((-1,1),repeat=8):
  if sum(a<0 for a in ss)%2==0:R.append(EI*B.inv()*S.Matrix(ss)/2)
R=sorted(set(tuple(map(int,r)) for r in R));idx={r:i for i,r in enumerate(R)}
dot=lambda a,b:sum(a[i]*int(C[i,j])*b[j] for i in range(8) for j in range(8))
neg=lambda a:tuple(-x for x in a)
plus=lambda a,b:tuple(x+y for x,y in zip(a,b))
def epsilon(a,b):
  z=sum(a[i]*b[i] for i in range(8))+sum(a[i]*int(C[i,j])*b[j] for i in range(8) for j in range(i))
  return (-1)**(z%2)
units=[tuple(int(i==j) for i in range(8)) for j in range(8)]
ident=list(range(480))
def compose(a,b):return [a[b[i]] for i in range(len(a))]
def inverse(a):
  b=[0]*len(a)
  for i,j in enumerate(a):b[j]=i
  return b
def power(a,n):
  b=list(range(len(a)))
  while n:
    if n&1:b=compose(b,a)
    a=compose(a,a);n//=2
  return b
def tperm(a):
  out=[]
  for b in R:
    k=dot(a,b)
    if b==a or b==neg(a):target=neg(b);sgn=1
    elif k==0:target=b;sgn=1
    elif k==-1:target=plus(a,b);sgn=epsilon(a,b)
    elif k==1:target=plus(neg(a),b);sgn=epsilon(neg(a),b)
    else:raise ValueError((a,b,k))
    out.append(idx[target]+(240 if sgn==-1 else 0))
  return out+[v+240 if v<240 else v-240 for v in out]
gen=[tperm(a) for a in units]
for i,a in enumerate(units):
  square=compose(gen[i],gen[i])
  expected=[k+(240 if dot(a,b)%2 else 0) for k,b in enumerate(R)]
  expected+= [v+240 if v<240 else v-240 for v in expected]
  assert square==expected
  for j in range(i):
    left,right=compose(gen[i],gen[j]),compose(gen[j],gen[i])
    if C[i,j]==0:assert left==right
    else:assert compose(left,gen[i])==compose(right,gen[j])

# Direct sparse exp(ad e_a) products on all 248 basis vectors.
def addscaled(out,v,c):
  for k,x in v.items():
    out[k]=out.get(k,0)+c*x
    if out[k]==0:del out[k]
def ad(a,v):
  out={}
  for k,c in v.items():
    if k>=240:
      z=-dot(a,units[k-240]);addscaled(out,{idx[a]:1},c*z)
    else:
      b=R[k]
      if b==neg(a):addscaled(out,{240+j:-a[j] for j in range(8) if a[j]},c)
      elif plus(a,b) in idx:addscaled(out,{idx[plus(a,b)]:epsilon(a,b)},c)
  return out
def exponential(a,v):
  out=dict(v);term=dict(v)
  for n in range(1,5):
    term=ad(a,term)
    if not term:return out
    addscaled(out,term,S.Rational(1,S.factorial(n)))
  raise ValueError('nilpotence exceeded')
for i,a in enumerate(units):
  for k in range(248):
    v=exponential(a,exponential(neg(a),exponential(a,{k:1})))
    assert all(S.sympify(z).q==1 for z in v.values())
    if k<240:
      target=gen[i][k];expected={target%240:1 if target<240 else -1}
    else:
      j=k-240;expected={240+j:1};addscaled(expected,{240+i:1},-int(C[i,j]))
    assert v==expected,(i,k,v,expected)
brackets=0
for i in range(8):
  g=gen[i]
  for aa,a in enumerate(R):
    ga=R[g[aa]%240];sa=1 if g[aa]<240 else -1
    for bb,b in enumerate(R):
      if plus(a,b) not in idx:continue
      gb=R[g[bb]%240];sb=1 if g[bb]<240 else -1
      ab=idx[plus(a,b)];gab=R[g[ab]%240];sab=1 if g[ab]<240 else -1
      assert plus(ga,gb)==gab and sa*sb*epsilon(ga,gb)==sab*epsilon(a,b)
      brackets+=1

source=next(z for z in geo['settings'] if z['family']=='A8' and z['eps']==1 and z['d']==3)
hb=[tuple(map(int,EI*S.Matrix(a))) for a in source['L_base']]
assert len(hb)==2 and dot(hb[0],hb[1])==-1
perp=[r for r in R if all(dot(r,h)==0 for h in hb)]
ps=[r for r in perp if all(x>=0 for x in r)]
pset=set(ps)
pb=[r for r in ps if not any(tuple(a-b for a,b in zip(r,t)) in pset for t in ps)]
pb=sorted(pb,reverse=True);assert len(pb)==6
pg=[tperm(r) for r in pb]
for g in pg:
  for h in hb:assert g[idx[h]]==idx[h]
c=reduce(compose,pg,ident);u=power(c,4)
assert power(u,3)==ident and u!=ident
def longest(bb):
  perm=ident;word=[];gens=[tperm(b) for b in bb]
  for _ in range(240):
    j=next((j for j,b in enumerate(bb) if all(v>=0 for v in R[perm[idx[b]]%240])),None)
    if j is None:return perm,word
    perm=compose(perm,gens[j]);word.append(j)
  raise ValueError('longest did not terminate')
n0,w0word=longest(units);nh,whword=longest(hb)
a=compose(n0,inverse(nh));assert power(a,2)==ident
assert all(compose(a,g)==compose(g,a) for g in pg)
assert all(a[idx[hb[j]]]==idx[hb[1-j]] for j in range(2))
assert power(compose(a,u),6)==ident and power(compose(a,u),3)==a

# Actual based F: multiply by source reflections until the H base is positive.
v=EI*S.Matrix(source['F'])*E
for _ in range(8):
  j=next((j for j,h in enumerate(hb) if any(z<0 for z in v*S.Matrix(h))),None)
  if j is None:break
  h=S.Matrix(hb[j]);v=v*(S.eye(8)-h*(h.T*C))
assert all(v*S.Matrix(h)==S.Matrix(h) for h in hb)
vf=[idx[tuple(map(int,v*S.Matrix(r)))]+1 for r in R]
shadow=lambda g:[g[j]%240+1 for j in range(240)]
gapperm=lambda a:'PermList('+str([j+1 for j in a])+')'
pgap=lambda a:'PermList('+str(a)+')'
gg=['tits60:=['+','.join(gapperm(g) for g in gen)+'];;',
 'perp60:=['+','.join(gapperm(g) for g in pg)+'];;',
 'vx60:=Group(perp60,());;',
 'u60:='+gapperm(u)+';;','a60:='+gapperm(a)+';;',
 'q060:=Centralizer(vx60,u60);;','q60:=Group(Concatenation(GeneratorsOfGroup(q060),[a60]));;',
 'tor260:=Group(List(tits60,g->g^2),());;',
 'Print("CARRIER60|",Size(vx60),"|",Order(u60),"|",Order(a60),"|",Size(q060),"|",Size(q60),"|",Size(Intersection(q60,tor260)),"\\n");',
 'Assert(0,Size(q060)=648 and Size(q60)=1296 and Size(Intersection(q60,tor260))=1);',
 'we60:=Group(['+','.join(pgap(shadow(g)) for g in gen)+'],());;',
 'wx60:=Group(['+','.join(pgap(shadow(g)) for g in pg)+'],());;',
 'wf60:='+pgap(vf)+';;','wu60:='+pgap(shadow(u))+';;',
 'Assert(0,wf60 in wx60);;','conj60:=RepresentativeAction(wx60,wf60,wu60);;',
 'Assert(0,conj60<>fail);;',
 'hbase60:='+str([idx[h]+1 for h in hb])+';;',
 'pin60:=Stabilizer(we60,Set(hbase60),OnSets);;',
 'point60:=Stabilizer(we60,hbase60,OnTuples);;',
 'Print("GEOMETRIC60|",Size(wx60),"|",Size(point60),"|",Size(pin60),"|",Size(Centralizer(pin60,wu60)),"\\n");',
 'Print("TRANSPORTER60|",conj60,"\\n");',
 'QUIT;']
(DIR/'carrier.g').write_text('\n'.join(gg)+'\n')
out={'direct_exp_columns':1984,'bracket_identities':brackets,'source_H_base':hb,'perpendicular_base':pb,
 'E8_longest_word_length':len(w0word),'H_longest_word_length':len(whword),'actual_based_F':v.tolist(),
 'u_order':3,'a_order':2,'au_order':6,'positive_graph':True,'centralizes_E6_Tits':True,
 'u':u,'a':a,'E8_generators':gen,'E6_generators':pg}
(DIR/'carrier.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
print('Fresh carrier complete:',1984,'integral Chevalley columns;',brackets,'bracket identities; actual cubic F recorded.',flush=True)
