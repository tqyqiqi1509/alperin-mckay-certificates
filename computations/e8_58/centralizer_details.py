"""Intrinsic cocharacter lattices of the full perpendicular centralizers."""
from pathlib import Path
from itertools import product
import json
import sympy as S
from sympy.matrices.normalforms import smith_normal_form

out=Path(__file__).parent
j=json.loads((out/'independent_geometry.json').read_text());Q=S.Matrix(j['lattice_basis'])
records={}
for name,d in j['cases'].items():
    if not d['K']:continue
    P=S.Matrix.hstack(*[S.Matrix(r) for r in d['perpendicular_base']])
    D=S.Matrix.hstack(*[S.Matrix(r) for r in d['J0_base']])
    Pinv=(P.T*P).inv()*P.T
    yy=Pinv*D
    assert all(t.q==1 for t in yy)
    zbas=(D.T*P).nullspace()
    assert len(zbas)<=1
    zs=[]
    for z in zbas:
        z=z*S.ilcm(*[v.q for v in z]);z=z/S.igcd(*[int(v) for v in z]);zs.append(z)
    joint=S.Matrix.hstack(yy,*zs)
    assert joint.rows==joint.cols
    smith=smith_normal_form(joint.applyfunc(int),domain=S.ZZ)
    rec={'derived_rank':D.cols,'perpendicular_rank':P.cols,
         'product_kernel_Smith':[abs(int(smith[i,i])) for i in range(smith.rows)],
         'product_kernel_order':abs(int(joint.det())),
         'central_cocharacters_intrinsic':[z.tolist() for z in zs],
         'central_cocharacters_E8':[(P*z).tolist() for z in zs]}
    K=S.Matrix.hstack(*[S.Matrix(r) for r in d['K_base']]);Z=S.Matrix.hstack(*K.T.nullspace());F=S.Matrix(d['F'])
    rec['C_Ls_center_poly']=str(S.factor(((Z.T*Z).inv()*Z.T*F*Z).charpoly().as_expr()))
    if name=='E7A1_quartic':
        alpha=sum([K[:,i] for i in range(K.cols)],S.zeros(8,1))/2
        assert all(v.q==1 for v in Q.inv()*alpha)
        assert (alpha.T*S.Matrix(d['s']))[0]==S.Rational(1,2)
        rec['nonidentity_component_character']=alpha.tolist()
        rec['component_value_on_s']='1/2'
    records[name]=rec
(out/'centralizer_details.json').write_text(json.dumps(records,indent=2,default=str)+'\n')
for name,r in records.items():print(name,'KERNEL',r['product_kernel_Smith'],'Z',r['central_cocharacters_E8'],'CLsCENTER',r['C_Ls_center_poly'])
assert records['D5A3_matching_D4']['product_kernel_order']==2
assert records['A7A1_cubic']['product_kernel_order']==3
assert records['A7A1_nonmatching_D6']['product_kernel_order']==2
