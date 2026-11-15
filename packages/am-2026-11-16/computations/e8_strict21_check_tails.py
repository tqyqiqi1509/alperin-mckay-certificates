"""Exact generic defects of every cyclic/zero candidate, in actual E8 forms."""
from pathlib import Path
import ast,json,re
import sympy as sp
from e8_strict21_parse_actual import cyc,phi_valuation,mat_eval

HERE=Path(__file__).parent
data=json.loads((HERE/'e8_strict21_actual_sources.json').read_text())
roots=data['roots'];cartan=sp.Matrix(data['simple_coroots'])
raw=(HERE/'e8_strict21_cyclic_tails.out').read_text()
assert 'DONE_TAILS21' in raw
raw=re.sub(r'(?m)^(?:warning:.*|#.*)$','',raw)
starts=list(re.finditer(r'(?:^|\n|gap> )(ZERO_TAIL|CYCLIC_TAIL)\|',raw))
zero=[];cyclic=[];cache={};primes=(3,5,7,11,13,17,19,23,29,31,41,61)
for i,m in enumerate(starts):
 end=starts[i+1].start() if i+1<len(starts) else raw.index('DONE_TAILS21',m.end())
 rec=raw[m.start(1):end].strip().split('|')
 common={'form_index':int(rec[1]),'s_order':int(rec[2]),'q_residue':int(rec[3]),'e':int(rec[4])}
 if rec[0]=='ZERO_TAIL':
  row={**common,'mark':rec[5],'kappa_degree':cyc(rec[6]),'M_order':cyc(rec[7]),'checks':[]}
  for ell in primes:
   if row['s_order']%ell==0 or (ell-1)%row['e']:continue
   a=phi_valuation(row['M_order'],ell,row['e']);b=phi_valuation(row['kappa_degree'],ell,row['e'])
   diff=(a[0]-b[0],a[1]-b[1])
   row['checks'].append({'ell':ell,'defect':diff,'valid_central_defect_source':diff==(0,0)})
  zero.append(row)
 else:
  row={**common,'M_type':rec[5],'mark':rec[6],'W_order':int(rec[7]),
       'K_indices':ast.literal_eval(rec[8]),'F_matrix':ast.literal_eval(rec[9]),
       'M_connected_centre':cyc(rec[10]),'kappa_degree':cyc(rec[11]),'M_order':cyc(rec[12]),'checks':[]}
  key=(tuple(row['K_indices']),str(row['F_matrix']),row['e'])
  if key not in cache:
   f=sp.Matrix(row['F_matrix']);k=sp.Matrix([roots[j-1] for j in row['K_indices']]) if row['K_indices'] else sp.zeros(0,8)
   cond=mat_eval(sp.cyclotomic_poly(row['e']),f).T.col_join(k*cartan)
   active=sp.Matrix.hstack(*cond.nullspace())
   assert active.cols==sp.totient(row['e'])
   lr=[r for r in roots if sp.Matrix([r])*cartan*active==sp.zeros(1,active.cols)]
   lm=sp.Matrix(lr) if lr else sp.zeros(0,8)
   z=sp.Matrix.hstack(*(lm*cartan).nullspace())
   h=z.gauss_jordan_solve(f.T*z)[0]
   x=sp.Symbol('x');p=h.charpoly(x).as_expr()
   ph={}
   for j in range(1,31):
    count=0;fac=sp.cyclotomic_poly(j,x)
    while sp.degree(p,x)>0 and sp.rem(p,fac,x)==0:
     p=sp.div(p,fac,x)[0];count+=1
    if count:ph[j]=count
   assert p==1
   cache[key]={'ambient_Levi_rank':lm.rank(),'L_centre_dimension':z.cols,
               'L_connected_centre':{'coefficient':'1','phi':ph,'display':str(h.charpoly(x).as_expr())}}
  row.update(cache[key])
  for ell in primes:
   if row['s_order']%ell==0 or (ell-1)%row['e']:continue
   a=phi_valuation(row['M_order'],ell,row['e']);b=phi_valuation(row['kappa_degree'],ell,row['e'])
   diff=(a[0]-b[0],a[1]-b[1]);mc=phi_valuation(row['M_connected_centre'],ell,row['e'])
   lc=phi_valuation(row['L_connected_centre'],ell,row['e'])
   valid=diff==mc
   if valid:
    assert lc==mc==(1,0),(row['M_type'],row['mark'],ell,diff,mc,lc)
    assert row['W_order']%ell
   row['checks'].append({'ell':ell,'source_defect':diff,'M_centre_defect':mc,'L_centre_defect':lc,
                        'valid_central_defect_source':valid})
  cyclic.append(row)
out={'method':'all actual rational source cosets, generic degree and centre valuations; no finite q',
     'cyclic_candidate_records':len(cyclic),'zero_candidate_records':len(zero),'cyclic_geometries':len(cache),
     'cyclic_sources':cyclic,'zero_sources':zero}
(HERE/'e8_strict21_cyclic_tails.json').write_text(json.dumps(out,indent=2,default=int))
print(json.dumps({k:v for k,v in out.items() if k not in ('cyclic_sources','zero_sources')},indent=2))
