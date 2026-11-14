"""Compact 36-setting certificate, classifying actual Dynkin graphs."""
from pathlib import Path
import json
import sympy as S
P=Path(__file__).parent
D=json.loads((P/'independent_geometry.json').read_text());C=S.Matrix(D['Cartan'])
def types(raw):
 bs=[S.Matrix(r) for r in raw];n=len(bs);rem=set(range(n));answer=[]
 adj={i:[j for j in range(n) if i!=j and (bs[i].T*C*bs[j])[0]] for i in range(n)}
 assert all((b.T*C*b)[0]==2 for b in bs)
 assert all((bs[i].T*C*bs[j])[0]==-1 for i in adj for j in adj[i])
 while rem:
  comp={rem.pop()};front=list(comp)
  while front:
   i=front.pop()
   for j in adj[i]:
    if j not in comp:comp.add(j);rem.remove(j);front.append(j)
  branches=[i for i in comp if len(adj[i])==3]
  if not branches:
   assert all(len(adj[i])<=2 for i in comp);answer.append('A'+str(len(comp)))
  else:
   assert len(branches)==1;c=branches[0];arms=[]
   for j in adj[c]:
    prev=c;now=j;le=1
    while len(adj[now])==2:
     nex=next(k for k in adj[now] if k!=prev);prev,now=now,nex;le+=1
    assert len(adj[now])==1;arms.append(le)
   arms.sort()
   if arms[:2]==[1,1]:answer.append('D'+str(len(comp)))
   else:
    assert arms in [[1,2,2],[1,2,3],[1,2,4]]
    answer.append('E'+str(len(comp)))
 return sorted(answer)
rows=[]
for name,r in D['sources'].items():
 rows.append({'setting':name,'d':r['d'],'q_residue':r['q_residue'],
  'L':types(r['L_base']),'K':types(r['K_base']),'J_derived':types(r['J_base']),
  'actual_L_F_diagram':r['L_graph'],'actual_J_F_diagram':r['J_graph'],
  'source_center_polynomial':r['active_polynomial'],
  'source_center_connected':all(z==1 for z in r['L_smith']),
  'norm_zero_iff_source_root':r.get('norm_zero_iff_source_root')})
assert len(rows)==36
core=json.loads((P/'gluing_and_cores.json').read_text())
carrier=json.loads((P/'independent_tits_carrier.json').read_text())
out={'settings':rows,'number_settings':36,
 'order5_actual_glue':core['order5']['actual_discriminant_glue'],
 'order6_actual_glue':core['order6']['actual_discriminant_glue'],
 'five_noncyclic_core_exhaustion':core['five_noncyclic_effective_packets'],
 'six_noncyclic_core_exhaustion':core['six_noncyclic_effective_packets'],
 'carrier':{k:carrier[k] for k in ['joint_lattice_index','X_cocharacter_smith','bracket_checks','opposite_root_Cartan_checks','explicit_248_Chevalley_products','a_order','u_order','au_order']}}
(P/'source_exhaustion.json').write_text(json.dumps(out,indent=2)+'\n')
for r in rows:
 print(r['setting'],'| L','+'.join(r['L']) or 'T','| K','+'.join(r['K']) or '1',
 '| Jder','+'.join(r['J_derived']) or '1','| Z',r['source_center_polynomial'])
print('CARRIER',out['carrier'])
