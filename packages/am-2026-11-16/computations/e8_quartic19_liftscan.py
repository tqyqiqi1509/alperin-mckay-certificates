from e8_quartic19_phase import *
row=ROWS[-1]
d=pickle.loads((OUT/('e8_quartic19_'+row['name']+'.pickle')).read_bytes())
E=PermutationGroup(list(d['E0'].args)); n0=d['n']; a,_=longest(SIMPLE); b,_=longest(row['base']); h=a*b**-1; hp=phase_lift(h)
P4=PermutationGroup([phase_lift(g) for g in E.generators]+[torus4(v) for v in d['ker4']])
zall=[]
for v in combos(d['ker2'],2):
 z=torus2(coeff_y(v))
 if z*n0==n0*z:zall.append((v,z,2*zeta(v,row)))
results=[]
for no,v in enumerate(combos(d['ker2'],2)):
 n=torus2(coeff_y(v))*n0; np=phase_lift(n)
 P=P4.centralizer(np);Q=PermutationGroup(P.generators+[IOTA])
 chk=extension_check(Q,[(phase_lift(z),val) for vv,z,val in zall],[(np,2),(hp,2)])
 result={'index':no,'torus_coefficients':v,'n_order':int(n.order()),'nh_order':int((n*h).order()),'P_order':int(P.order()),'Q_order':int(Q.order()),'signed_centralizer_image':int(PermutationGroup([shadow(g) for g in E.centralizer(n).generators]).order()),'extension':chk}
 print(json.dumps(result),flush=True);results.append(result)
 if chk['assignment_consistent'] and chk['kills_derived_intersection']:
  cache={'row':row,'Pgens':[g.array_form for g in P.generators],'Qgens':[g.array_form for g in Q.generators],'n':np.array_form,'h':hp.array_form,'zlist':[(vv,phase_lift(z).array_form,val) for vv,z,val in zall],'correcting_torus_coefficients':v}
  (OUT/'e8_quartic19_A7A1_external_A1_corrected_phase.pickle').write_bytes(pickle.dumps(cache))
  print('SUCCESS',no,flush=True)
  break
(OUT/'e8_quartic19_external_lift_scan.json').write_text(json.dumps(results,indent=2))
