from itertools import permutations
from collections import Counter
from pathlib import Path
import json,math
z=(1,2,0,4,5,3)
def cyc(p):
 todo=set(range(len(p)));v=[]
 while todo:
  x=next(iter(todo));j=x;k=0
  while j in todo:todo.remove(j);j=p[j];k+=1
  v.append(k)
 return tuple(sorted(v))
def vp(n):
 k=0
 while n%3==0:n//=3;k+=1
 return k
forms=Counter(cyc(p) for p in permutations(range(6)) if all(p[z[i]]==z[p[i]] for i in range(6)))
rows=[]
for ds,count in sorted(forms.items()):
 # product GL1(Q^d), valuations sum(a+v3(d)); GL6 valuation6a+2
 delta=(6-len(ds),2-sum(vp(d) for d in ds))
 assert delta[0]+delta[1]>=2
 rows.append({'multiplicities':[1]*6,'F_cycles':ds,'multiplicity':count,'index_valuation_a_coefficient':delta[0],'constant':delta[1]})
for ds in [(1,1,1),(3,)]:
 # product GL2(Q^d), each 2a+2v3(d)
 delta=(6-2*len(ds),2-2*sum(vp(d) for d in ds));assert sum(delta)>=2
 rows.append({'multiplicities':[2]*3,'F_cycles':ds,'index_valuation_a_coefficient':delta[0],'constant':delta[1]})
Path('outputs/e6_cubic_index_round3_checks.json').write_text(json.dumps({'rows':rows,'all_a_at_least_one_index_v3_at_least_2':True,'matching_Q': 'Q=epsilon*q congruent1 mod3; sign included in rational GL/GU factor'},indent=2));print(rows)
