"""Interpret all actual order-six rational pinning normalizers from GAP."""
import independent_e7 as m
from pathlib import Path
import ast,json
import sympy as S

text=Path(__file__).with_name('normalizers.txt').read_text().splitlines()
rows=[];pending=None
for line in text:
 if line.startswith('RATIONAL6|'):
  _,q,data=line.split('|',2);pending=[int(q),data]
 elif pending:pending[1]+=line
 if pending and pending[1].count('[')==pending[1].count(']'):
  rows.append((pending[0],S.Matrix(ast.literal_eval(pending[1])).T));pending=None
assert len(rows)==8 and pending is None
LC=m.C6.inv()*m.L6.T
Cb=m.components(m.kb6)
Lpos={tuple(r) for r in m.lr6 if m.dot(m.height,r)>0}
z6=next(m.C6.inv()[:,j] for j in range(6) if any(x.q==3 for x in m.C6.inv()[:,j]))
z6E=m.L6*z6
out=[]
frac=lambda x:tuple(t-S.floor(t) for t in x)
for q,M in rows:
 F=m.B*M*m.left+(m.I8-m.B*m.left)
 assert all(m.dot(r,q*F*m.s6-m.s6).q==1 for r in m.base)
 sig=m.dot(m.center[0],F*m.center[0])/m.dot(m.center[0],m.center[0])
 assert sig in (1,-1)
 graphC=[next(j for j,t in enumerate(m.kb6) if F*r==t) for r in m.kb6]
 factor=[next(j for j,c in enumerate(Cb) if set(graphC[i] for i in co)==set(c)) for co in Cb]
 Fb=F
 while any(tuple(Fb*r) not in Lpos for r in m.lb6):
  r=next(r for r in m.lb6 if tuple(Fb*r) not in Lpos);Fb=Fb*m.ref(r)
 graphL=[next(j for j,t in enumerate(m.lb6) if Fb*r==t) for r in m.lb6]
 eps=1 if graphL==list(range(6)) else -1
 assert eps==sig
 fix=[]
 for k in range(3):
  delta=LC*(q*F*(m.s60+k*z6E)-(m.s60+k*z6E))
  if all(t.q==1 for t in delta):fix.append(k)
 obstruction=frac(LC*(q*F*m.s60-m.s60))
 natural_obstruction=frac(m.left*(q*F*m.s6-m.s6))
 # Rational types of C_s from the actual based root permutation.
 used=set();factor_types=[]
 for i in range(3):
  if i in used:continue
  cyc=[i];used.add(i);j=factor[i]
  while j!=i:cyc.append(j);used.add(j);j=factor[j]
  power=len(cyc);Fp=F**power
  nodes=Cb[i]
  innergraph=all(Fp*m.kb6[t]==m.kb6[t] for t in nodes)
  factor_types.append(('A2' if innergraph else '2A2')+'(q^'+str(power)+')')
 assert sorted(len(c) for c in [Cb])==[3]
 out.append({'q_mod12':q,'F_matrix_Q':M.tolist(),'center_sign':int(sig),
 'E6_form_sign':eps,'Cs_component_permutation':factor,'Cs_root_permutation':graphC,
 'Cs_actual_types':factor_types,'E6_source_center_natural_size':3 if (q-eps)%3==0 else 1,
 'E6_SC_rational_lift_choices':fix,'E6_SC_rational_obstruction':list(obstruction),
 'natural_E7_central2_faithful':any(natural_obstruction),
 'E7_SC_obstruction_in_Q':list(natural_obstruction)})
 print(q, 'E6 sign',eps,'Cs',factor_types,'SC choices',fix,'E7 faithful',any(natural_obstruction))
 assert fix
 assert not any(len(set([i,factor[i],factor[factor[i]]]))==3 for i in range(3))
Path(__file__).with_suffix('.json').write_text(json.dumps({'rows':out,'no_A2_q3_branch':True,
 'all_derived_E6_parameters_admit_SC_rational_lifts':True},indent=2,default=str)+'\n')
