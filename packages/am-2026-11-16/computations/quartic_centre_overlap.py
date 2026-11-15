import sys,json,pickle
from pathlib import Path
base=Path(__file__).resolve().parent;sys.path.insert(0,str(base))
from e8_quartic19_phase import *
out=[]
for row in ROWS:
 d=pickle.loads((base/('e8_quartic19_'+row['name']+'.pickle')).read_bytes())
 if 'signed' in row:
  suffix='_corrected_phase.pickle' if row['name']=='A7A1_external_A1' else '_phase.pickle'
  c=pickle.loads((base/('e8_quartic19_'+row['name']+suffix)).read_bytes())
  P=PermutationGroup([Permutation(g) for g in c['Pgens']]); zs=[Permutation(x) for y,x,v in c['zlist']]
 else:
  P=PermutationGroup(list(d['E0'].args));zs=[torus2(coeff_y(v)) for v in combos(d['ker2'],2)]
 W=PermutationGroup([shadow(g) for g in P.generators]);Z=PermutationGroup(zs)
 assert all(P.contains(z) and z.order()<=2 and shadow(z)==I240 for z in zs)
 assert P.order()==W.order()*Z.order()
 out.append({'row':row['name'],'P_order':int(P.order()),'W_order':int(W.order()),'kernel_order':int(Z.order()),'kernel_equals_listed_central_two_torsion':True})
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
