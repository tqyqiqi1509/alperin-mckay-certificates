"""Independent exact gluing and type-A core exhaustion, no finite-q samples."""
from pathlib import Path
from itertools import product
import json
import sympy as S
O=Path(__file__).parent
D=json.loads((O/'independent_geometry.json').read_text())
C=S.Matrix(D['Cartan']);E=S.eye(8)
roots={tuple(int(z) for z in E[:,i]) for i in range(8)};todo=list(roots)
while todo:
 r=todo.pop()
 for i in range(8):
  x=list(r);x[i]-=sum(int(C[i,j])*r[j] for j in range(8));x=tuple(x)
  if x not in roots:roots.add(x);todo.append(x)
rs=[S.Matrix(r) for r in sorted(roots)]
def base(rr):
 pos=[r for r in rr if all(z>=0 for z in r)];ps=set(map(tuple,pos))
 return [r for r in pos if not any(tuple(r-t) in ps for t in pos)]
def chains(bs):
 rem=set(range(len(bs)));out=[]
 while rem:
  cc={rem.pop()}
  while True:
   new={j for i in cc for j in rem if (bs[i].T*C*bs[j])[0]}
   if not new:break
   cc|=new;rem-=new
  adj={i:[j for j in cc if j!=i and (bs[i].T*C*bs[j])[0]] for i in cc}
  sq=[min(i for i in cc if len(adj[i])<=1)]
  while len(sq)<len(cc):sq.append(next(j for j in adj[sq[-1]] if j not in sq))
  out.append([bs[i] for i in sq])
 return sorted(out,key=len,reverse=True)
answer={}
for n in [5,6]:
 ss=S.Matrix(D['s'+str(n)])
 cs=[r for r in rs if (r.T*C*ss)[0].q==1];cc=chains(base(cs))
 B=S.Matrix.hstack(*sum(cc,[]));Bi=B.inv()
 mods=[len(x)+1 for x in cc]
 offsets=[sum(map(len,cc[:i])) for i in range(len(cc))]
 gens=[]
 for j in range(8):
  col=Bi*E[:,j]
  row=tuple(int(mod*col[off])%mod for mod,off in zip(mods,offsets))
  assert all((mod*col[off]).q==1 for mod,off in zip(mods,offsets))
  gens.append(row)
 subgroup={(0,)*len(mods)};changed=True
 while changed:
  new={tuple((x+y)%m for x,y,m in zip(a,b,mods)) for a in subgroup for b in gens}
  changed=not new.issubset(subgroup);subgroup|=new
 assert len(subgroup)==n==abs(B.det())
 full=[g for g in subgroup if all(len({(i*x)%m for i in range(n)})==m for x,m in zip(g,mods))]
 assert full
 answer['order'+str(n)]={'component_centers':mods,'root_lattice_index':n,
  'ambient_generators':gens,'actual_discriminant_glue':sorted(subgroup),
  'full_projection_generators':full}
def partitions(n,cap=None):
 if n==0:yield ();return
 for k in range(min(n,cap or n),0,-1):
  for tail in partitions(n-k,k):yield (k,)+tail
def core(part,e):
 m=12 # fixed number of abacus beads, divisible by each e used below except5
 b={part[i]+m-i-1 if i<len(part) else m-i-1 for i in range(m)}
 while True:
  movable=[x for x in b if x>=e and x-e not in b]
  if not movable:break
  x=min(movable);b.remove(x);b.add(x-e)
 bb=sorted(b,reverse=True)
 return tuple(v for i,x in enumerate(bb) if (v:=x-(m-i-1))>0)
# The abacus is independent of the bead count; cross-check with hook divisibility.
def hooks(part):
 return [part[i]-j+sum(part[k]>j for k in range(i+1,len(part))) for i in range(len(part)) for j in range(part[i])]
cores={}
for n in [2,3,5,6]:
 for e in range(1,7):
  out={}
  for lam in partitions(n):
   co=core(lam,e);assert not any(h%e==0 for h in hooks(co))
   wt=(n-sum(co))//e;assert sum(co)+e*wt==n
   out.setdefault(str(co),{'core':co,'weight':wt,'partitions':[]})['partitions'].append(lam)
  cores[str((n,e))]=list(out.values())
answer['all_type_A_core_packets']=cores
# Positive-rank, noncyclic packets: product of partitions' core choices.
def packets(sizes,e):
 ans=[]
 for cc in product(*[cores[str((n,e))] for n in sizes]):
  wt=sum(c['weight'] for c in cc)
  if wt>=2:ans.append({'cores':[c['core'] for c in cc],'total_weight':wt})
 return ans
answer['five_noncyclic_effective_packets']={str(e):packets([5,5],e) for e in [2,3,4,5]}
answer['six_noncyclic_effective_packets']={str(e):packets([6,3,2],e) for e in [2,3,4,5,6]}
(O/'gluing_and_cores.json').write_text(json.dumps(answer,indent=2)+'\n')
for k in ['order5','order6','five_noncyclic_effective_packets','six_noncyclic_effective_packets']:
 print(k,answer[k])
