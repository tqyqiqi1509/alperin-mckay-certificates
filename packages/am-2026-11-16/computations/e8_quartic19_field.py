from e8_quartic19_phase import *
out=[]
for row in ROWS:
 d=pickle.loads((OUT/('e8_quartic19_'+row['name']+'.pickle')).read_bytes());E=PermutationGroup(list(d['E0'].args));a,_=longest(SIMPLE);b,_=longest(row['base']);h=a*b**-1
 if 'signed' not in row:
  P=E;n=h;inv=ID;hp=h
  item={'row':row['name'],'regular_signed_carrier_primefield_fixed':True,'h_central':all(h*g==g*h for g in P.generators)}
 else:
  fn='e8_quartic19_'+row['name']+('_corrected_phase.pickle' if row['name']=='A7A1_external_A1' else '_phase.pickle');c=pickle.loads((OUT/fn).read_bytes());P=PermutationGroup([Permutation(g) for g in c['Pgens']]);n=Permutation(c['n']);hp=Permutation(c['h']);inv=IOTA
  p4gens=[phase_lift(g) for g in E.generators]+[torus4(v) for v in d['ker4']]
  item={'row':row['name'],'h_equals_iota_on_P4':all(hp**-1*g*hp==IOTA**-1*g*IOTA for g in p4gens),'iota_fixes_signed_returns':inv*n==n*inv and inv*hp==hp*inv}
 item['n_central_on_P']=all(n*g==g*n for g in P.generators)
 item['q1_alpha_f_matches_return']=all((n**-1)*g*n==g for g in P.generators)
 u=n*hp
 item['q3_alpha_f_matches_return']=all((u**-1)*g*u==inv*g*inv for g in P.generators)
 item['return_orders']=[int(n.order()),int(u.order())]
 assert all(v for k,v in item.items() if isinstance(v,bool))
 print(json.dumps(item),flush=True);out.append(item)
(OUT/'e8_quartic19_field_certificate.json').write_text(json.dumps(out,indent=2))
