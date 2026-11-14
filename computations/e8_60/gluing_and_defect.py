"""Intrinsic full-J gluing and the A8 seven-defect normalizer argument."""
from pathlib import Path
from itertools import product
import json
import sympy as S
from sympy.matrices.normalforms import hermite_normal_form,smith_normal_form
DIR=Path(__file__).resolve().parent
geo=json.loads((DIR/'geometry.json').read_text());C=S.Matrix(geo['gram'])
def integer_kernel(A):
    zz=A.nullspace()
    if not zz:return S.zeros(A.cols,0)
    N=S.Matrix.hstack(*zz);den=int(S.ilcm(*[z.q for z in N]))
    vv=[]
    for y in product(range(den),repeat=N.cols):
        v=N*S.Matrix(y)
        if any(v) and all(z.q==1 for z in v):vv.append(v)
    vv.extend(N[:,j]*den for j in range(N.cols))
    return hermite_normal_form(S.Matrix.hstack(*vv).applyfunc(int))
rows=[]
for i,r in enumerate(geo['settings']):
    if r['K_type']=='1':continue
    P=S.Matrix.hstack(*(S.Matrix(z) for z in r['perpendicular_base']))
    J=S.Matrix.hstack(*(S.Matrix(z) for z in r['J_base']))
    jc=P.gauss_jordan_solve(J)[0];assert all(z.q==1 for z in jc)
    pairing=(J.T*C*P).applyfunc(int)
    cz=integer_kernel(pairing)
    joint=S.Matrix.hstack(jc,cz);index=abs(int(joint.det()))
    sd=smith_normal_form(pairing,domain=S.ZZ)
    center_component=[abs(int(sd[j,j])) for j in range(J.cols)]
    rows.append({'row':i,'family':r['family'],'eps':r['eps'],'d':r['d'],
                 'P':r['perpendicular_type'],'J':r['J_type'],
                 'derived_times_connected_center_index':index,
                 'center_character_smith':center_component})
    print(i,r['family'],r['eps'],r['d'],'fullJ',r['J_type'],'product index',index,'center components',center_component)

# A8 weight/root indices are 9/3, hence units at seven.  This is its
# actual seven-cycle on the augmentation lattice with basis e_i-e_9.
pi=S.zeros(8)
for j in range(8):pi[(j+1)%7 if j<7 else j,j]=1
def rankmod(M,p):
    a=[[int(v)%p for v in M.row(i)] for i in range(M.rows)];r=0
    for j in range(M.cols):
        k=next((k for k in range(r,M.rows) if a[k][j]),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r];inv=pow(a[r][j],-1,p);a[r]=[(v*inv)%p for v in a[r]]
        for k in range(M.rows):
            if k!=r:
                c=a[k][j];a[k]=[(v-c*w)%p for v,w in zip(a[k],a[r])]
        r+=1
    return r
ranks=[rankmod((pi**j-S.eye(8))**2,7) for j in range(1,7)]
assert ranks==[5]*6
out={'full_J_gluing':rows,'seven_cycle_square_ranks':ranks,
     'deduction':'Every abelian normal subgroup B of D has (pi-1)^2 A=0 for each image pi of B; nonzero ranks force B<=A. Thus A is characteristic, for every valuation and without a split extension assumption.'}
(DIR/'gluing_and_defect.json').write_text(json.dumps(out,indent=2)+'\n')
print('A8 ell7 nonidentity square ranks:',ranks)
