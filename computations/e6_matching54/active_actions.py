"""Reconstruct exact ordinary coefficient actions from actual E6 Weyl lifts."""
from pathlib import Path
import ast,io,contextlib,runpy,json
from itertools import product
import sympy as S

p=Path(__file__).parent
C=2*S.eye(6)
for i,j in ((0,2),(2,3),(3,4),(4,5),(1,3)):C[i,j]=C[j,i]=-1
I=S.eye(6)
rr={tuple(I[:,j]) for j in range(6)};todo=list(rr)
refs=[I-I[:,j]*I[:,j].T*C for j in range(6)]
while todo:
 r=S.Matrix(todo.pop())
 for w in refs:
  t=tuple(w*r)
  if t not in rr:rr.add(t);todo.append(t)
roots=[S.Matrix(r) for r in sorted(rr)];assert len(roots)==72
cases=json.loads((p/'independent_geometry.json').read_text())['cases']
print('ROOTS|72',flush=True)
lines=(p/'independent_geometry_gap_final.txt').read_text().splitlines()
actions={};i=0
while i<len(lines):
 if lines[i].startswith('ACTION|'):
  _,name,label,tail=lines[i].split('|',3)
  buf=tail;balance=buf.count('[')-buf.count(']')
  while '[' not in buf or balance:
   i+=1;buf+='\n'+lines[i];balance=buf.count('[')-buf.count(']')
  rows=ast.literal_eval(buf)
  actions[name,label]=[S.Matrix.hstack(*[roots[j-1] for j in row]) for row in rows]
 i+=1

def normalize(M):return M.applyfunc(S.simplify)
def key(M):return tuple(M)
def closure(gs):
 n=gs[0].rows;one=S.eye(n);seen={key(one):one};todo=[one]
 gens=list({key(g):g for g in gs}.values())
 while todo:
  g=todo.pop()
  for h in gens:
   m=normalize(g*h);k=key(m)
   if k not in seen:seen[k]=m;todo.append(m)
 return list(seen.values())
def act(U,m):
 rows=list(U.T.rref()[1]);cols=list(range(U.cols))
 return normalize(U.extract(rows,cols).inv()*(m*U).extract(rows,cols))

out={}
for name in ['split_d2_regular','split_d2_marked']:
 F=S.Matrix(cases[name]['F']);U=S.Matrix.hstack(*(F+I).nullspace())
 gs=[act(U,m) for m in actions[name,'REG']]
 eig=[]
 for signs in product((-1,1),repeat=len(gs)):
  ns=S.Matrix.vstack(*[g-e*S.eye(U.cols) for g,e in zip(gs,signs)]).nullspace()
  if len(ns)==1:eig.append(ns[0])
 E=S.Matrix.hstack(*eig);assert E.cols==U.cols and E.det()!=0
 Z=U*E;norm=[(Z[:,j].T*C*Z[:,j])[0] for j in range(Z.cols)]
 # Normalize by equal root-form lengths. The needed ratios are rational squares.
 fac=[S.sqrt(norm[0]/x) for x in norm];assert all(x.is_Rational for x in fac)
 Z=Z*S.diag(*fac)
 rec={'basis':Z.tolist(),'norm':norm[0],'groups':{}}
 for lab in ['REG','NAT','FULL']:
  ms=closure([act(Z,m) for m in actions[name,lab]])
  mon=all(all(sum(v!=0 for v in m.row(j))==1 and sum(v!=0 for v in m.col(j))==1
     for j in range(m.rows)) and all(v in (-1,0,1) for v in m) for m in ms)
  assert mon
  perms={tuple(next(i for i in range(m.rows) if m[i,j]) for j in range(m.cols)) for m in ms}
  rec['groups'][lab]={'order':len(ms),'permutation_image_order':len(perms),'signed_monomial':mon}
 out[name]=rec;print(name,rec['groups'],flush=True)

for name in ['cyclic_d1','cyclic_d3']:
 F=S.Matrix(cases[name]['F'])
 U=S.Matrix.hstack(*((F-I) if name=='cyclic_d1' else (F*F+F+I)).nullspace())
 assert U.cols==(2 if name=='cyclic_d1' else 4)
 rec={'basis':U.tolist(),'groups':{}}
 reps={}
 for lab in ['REG','NAT','FULL']:
  ms=closure([act(U,m) for m in actions[name,lab]])
  reps[lab]=ms
  rec['groups'][lab]={'order':len(ms),'traces':sorted({str(S.simplify(S.trace(m))) for m in ms})}
 if name=='cyclic_d1':
  assert key(-S.eye(2)) in {key(m) for m in reps['FULL']}
  assert {key(m) for m in reps['FULL']}=={key(e*m) for m in reps['REG'] for e in (-1,1)}
  rec['ambient_is_S3_times_inversion']=True
 else:
  cf=act(U,F)
  assert cf*cf+cf+S.eye(4)==S.zeros(4)
  assert all(cf*m==m*cf for m in reps['REG'])
  expected={key(cf**j*m) for m in reps['REG'] for j in range(3)}
  assert {key(m) for m in reps['NAT']}==expected
  assert sorted([S.trace(m) for m in reps['REG']])==[-2,-2,0,0,0,4]
  rec['natural_is_scalar_C3_times_standard_S3']=True
  rec['cyclotomic_scalar']=cf.tolist()
 out[name]=rec;print(name,rec['groups'],flush=True)
(p/'active_actions.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
print('DONE|all actual actions verified')
