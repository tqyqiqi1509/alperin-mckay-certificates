"""Exact rational equivalent of e67_source_regular_phases.wl.

Computes the Cartan determinant, Gram matrix and displayed phases only.
"""
from pathlib import Path
import json
import sympy as s

C=2*s.eye(6)
for i,j in [(0,2),(2,3),(3,4),(4,5),(1,3)]:
    C[i,j]=C[j,i]=-1
H=s.Matrix([[0,0,1,1,1,1],[0,1,0,1,1,1],
            [1,0,1,1,1,0],[1,1,1,1,0,0]])
K=s.Matrix([[0,0,0,1,0,0],[0,1,1,1,1,0]])
positive=H[0,:].col_join(H[3,:]).col_join(H[0,:]+H[3,:])
positive=positive.col_join(H[1,:]).col_join(H[2,:]).col_join(H[1,:]+H[2,:])
unit=s.eye(6)
parameters=[C.inv()*unit[:,3]/3,
            C.inv()*(unit[:,0]+unit[:,5])/3,
            C.inv()*(unit[:,1]+unit[:,2]+unit[:,4])/6]
single=K.col_join(K[0,:]+K[1,:])
phase=lambda rows,v: [x % 1 for x in rows*C*v]
all_phases=[phase(positive,v) for v in parameters]
single_three=phase(single,parameters[0])
single_six=phase(single,parameters[2])
quadratic=phase(positive,3*parameters[2])
assert C.det()==3
assert H*C*H.T==s.Matrix([[2,0,0,-1],[0,2,-1,0],[0,-1,2,0],[-1,0,0,2]])
assert all_phases==[[s.Rational(1,3),s.Rational(1,3),s.Rational(2,3)]*2]*3
assert single_three==[s.Rational(1,3),s.Rational(1,3),s.Rational(2,3)]
assert single_six==[0,s.Rational(1,2),s.Rational(1,2)]
assert quadratic==[0]*6
result={"CartanDeterminant":int(C.det()),
        "TwoA2PositiveGram":[[int(x) for x in row] for row in (H*C*H.T).tolist()],
        "TwoA2RegularPhases":[list(map(str,row)) for row in all_phases],
        "SingleA2IsolatedThreePhases":list(map(str,single_three)),
        "SingleA2OrderSixPhases":list(map(str,single_six)),
        "QuadraticWholeLeviPhases":list(map(str,quadratic))}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
