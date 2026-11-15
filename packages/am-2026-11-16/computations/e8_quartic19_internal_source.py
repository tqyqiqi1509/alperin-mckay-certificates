"""Actual internal A3 and classical Springer symbol; source centralizer kernels."""
from e8_quartic19_phase import *
from sympy import Symbol, factor, eye
row=ROWS[3];core=row['base'][:3];roots,pos=source_data(core)
cr,cb,C=get_C('D');Lroots,Lpos=source_data(row['base']);w=signed_weyl(row['signed']);cw=based_projection(w,core,pos);ids=[INDEX[a] for a in core]
A=Matrix([[dot(a,b)//4 for b in SIMPLE] for a in core]);V=Matrix.hstack(*A.nullspace());M=Matrix.hstack(*[Matrix(action_vector(w,b)) for b in SIMPLE]);M=BM.T.inv()*M;W=V.gauss_jordan_solve(M*V)[0];q=Symbol('q')
# XSP(2,0,4,true): .31 has symbol [[0,2],[1,5]].
# It is distinguished because its interleaving is already increasing.
symbol=[[0,2],[1,5]];combined=sorted(sum(symbol,[]));orbit=sorted([2*(v-i)+1 for i,v in enumerate(combined)],reverse=True)
ldeg=(q*q-1)*(q**4-1)**2*(q**6-1)*(q*q+1)**2
cdeg=(q*q-1)*(q**3+1)*(q**4-1)*(q+1)*(q*q+1)**2
phi22=q*q*(q*q+1)
out={'L_C_intersection_roots':len(set(Lroots)&{INDEX[r] for r in cr}),'core_root_count':len(roots),'core_based_frobenius_action':[ids.index(cw(i)) for i in ids],'centralizer_torus_polynomial':str(factor((q*W-eye(W.rows)).det())),'ordinary_springer_symbol':symbol,'symbol_orbit':orbit,'phi22_degree':str(phi22),'Jordan_source_degree':str(factor(ldeg/cdeg*phi22)),'opposite_degree':str(factor((ldeg/cdeg*phi22).subs(q,-q)))}
(OUT/'e8_quartic19_internal_source.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
kernels=[]
for row in ROWS:
 lr,lp=source_data(row['base']);cr,cb,wc=get_C(row['kind']);inter=[r for r in cr if INDEX[r] in lr];group=PermutationGroup([shadow(lift(r)) for r in inter] or [I240]);kernels.append({'row':row['name'],'source_semisimple_centralizer_root_count':len(inter),'source_semisimple_centralizer_Weyl_order':int(group.order())})
(OUT/'e8_quartic19_source_kernel.json').write_text(json.dumps(kernels,indent=2))
