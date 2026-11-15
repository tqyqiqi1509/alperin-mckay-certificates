from pathlib import Path
import sys,pickle,json
from sympy import Matrix,Symbol,factor,eye
folder=Path(__file__).resolve().parent;sys.path.insert(0,str(folder))
from e8_quartic19 import ROOTS,INDEX,ROWS,SIMPLE,BM,source_data,get_C,dot,signed_weyl,based_projection
row=ROWS[3];base=row['base'];core=base[:3];lr,lp=source_data(base);kr,kp=source_data(core);cr,_,_=get_C('D')
d=pickle.loads((folder/'e8_quartic19_D5A3_internal_A3.pickle').read_bytes())
phase=pickle.loads((folder/'e8_quartic19_D5A3_internal_A3_phase.pickle').read_bytes())
ids=[INDEX[r] for r in base]
actions=[]
for arr in phase['Qgens']:
    action=[ids.index(arr[i]%240) for i in ids]
    coefficients=[arr[i]//240 for i in ids]
    assert coefficients==[0]*4,coefficients
    assert action in [[0,1,2,3],[1,0,2,3]],action
    actions.append(action)
assert [1,0,2,3] in actions
assert set(kr)==set(lr)&{INDEX[r] for r in cr}
w=d['w'];assert [w(i) for i in ids]==ids
raw=signed_weyl(row['signed']);c=based_projection(raw,core,kp)
assert [ids.index(c(i)) for i in ids[:3]]==[1,0,2]
def act(p,y):
    coords=BM.T.inv()*Matrix(y)
    return sum((coords[j]*Matrix(ROOTS[p(INDEX[SIMPLE[j]])]) for j in range(8)),Matrix.zeros(8,1))
q=Symbol('q')
def torus_poly(roots,p):
    eq=Matrix([[dot(a,b)//4 for b in SIMPLE] for a in roots])
    V=Matrix.hstack(*eq.nullspace());M=Matrix.hstack(*[act(p,a) for a in SIMPLE]);M=BM.T.inv()*M
    R=V.gauss_jordan_solve(M*V)[0]
    return factor(R.charpoly(q).as_expr())
zL=torus_poly(base,w);zC=torus_poly(core,raw)
assert zL==(q*q+1)**2
assert zC==(q+1)*(q*q+1)**2
# Hook formula for partition22 in GL4, followed by Ennola substitution.
hooks=[3,2,2,1];degree=factor(q**2*__import__('functools').reduce(lambda a,b:a*b,[q**i-1 for i in range(1,5)],1)/__import__('functools').reduce(lambda a,b:a*b,[q**h-1 for h in hooks],1))
assert degree==q*q*(q*q+1)
Lprime=(q**2-1)*(q**4-1)**2*(q**6-1)*zL
Cprime=(q**2-1)*(q**3+1)*(q**4-1)*zC
source=factor(Lprime/Cprime*degree)
assert source==q*q*(q-1)**2*(q*q+q+1)*(q*q+1)**2
rec={'root_gram':[[dot(a,b)//4 for b in base] for a in base],'L_root_count':len(lr),'centralizer_root_count':len(kr),'Q_based_actions':actions,'Q_simple_root_coefficients':'all +1','L_frobenius_based_action':[0,1,2,3],'centralizer_based_action':[1,0,2],'L_center_polynomial':str(zL),'centralizer_center_polynomial':str(zC),'hook22_degree':str(degree),'Jordan_source_degree':str(source),'opposite_degree':str(factor(source.subs(q,-q)))}
Path(__file__).with_suffix('.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec))
