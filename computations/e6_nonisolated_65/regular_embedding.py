"""Actual connected-centre E6 Levi in E7, for both Frobenius signs."""
from itertools import product
from pathlib import Path
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

C=2*s.eye(7)
for i,j in [(0,2),(2,3),(3,4),(4,5),(5,6),(1,3)]:
    C[i,j]=C[j,i]=-1
B=s.eye(7)[:,:6]
C6=B.T*C*B
g=s.zeros(6)
for j,i in enumerate([5,1,4,3,2,0]):
    g[i,j]=1
w0H=s.eye(7)+B*(-g-s.eye(6))*C6.inv()*B.T*C
w=-w0H
assert all(z.q==1 for z in w)
assert w.T*C*w==C and w*w==s.eye(7)
assert w*B==B*g
centre=(B.T*C).nullspace()[0]
assert w*centre==-centre
SNF=smith_normal_form(C*B,domain=ZZ)
assert list(SNF.diagonal())==[1]*6
roots=[]
for v in product(range(3),range(3),range(4),range(5),range(4),range(3),range(2)):
    a=s.Matrix(v)
    if (a.T*C*a)[0]==2:
        roots.extend([a,-a])
assert len(roots)==126
rs={tuple(a) for a in roots}
assert {tuple(w*a) for a in roots}==rs
out={'E7_cartan':C.tolist(),'E7_root_count':126,
     'Levi_centre_smith':[1]*6,'centre_connected':True,
     'graph_normalizer_matrix':w.tolist(),'normalizer_order':2,
     'E6_action':'pinned graph on root system',
     'central_torus_action':-1,
     'conclusion':'An actual F-stable E6.T1 Levi of E7_sc supplies the connected-centre regular group for both epsilon signs; KM15 Theorem A(e) applies directly.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
print(json.dumps(out,indent=2,default=str))
