import sympy as s,json,itertools,hashlib
from pathlib import Path
O=Path(__file__).parent
I=s.eye(8);zero=s.zeros(8)
def E(i,j):
 a=s.zeros(8);a[i,j]=1;return a
def M(i,j):return E(i,j)-E(j+4,i+4)
def P(i,j):return E(i,j+4)-E(j,i+4)
roots=[]
for i,j in itertools.combinations(range(4),2): roots.extend([M(i,j),M(j,i),P(i,j),P(i,j).T])
hs=[E(i,i)-E(i+4,i+4) for i in range(4)]
basis=roots+hs
x=M(0,1)+M(1,2)+P(1,2);xd=-x.T;H=s.diag(4,2,0,0,-4,-2,0,0)
flat=lambda a:s.Matrix(list(a))
def grade(a):
 b=H*a-a*H
 return next(k for k in range(-8,9) if b==k*a)
grades=[grade(a) for a in basis]
g0=[a for a,k in zip(basis,grades) if k==0];g2=[M(0,1),M(1,2),P(1,2),M(1,3),P(1,3)]
B=s.Matrix.hstack(*map(flat,g2))
A=s.Matrix.hstack(*[B.gauss_jordan_solve(flat(a*x-x*a))[0] for a in g0])
minor=next((list(c),int(A[:,c].det())) for c in itertools.combinations(range(A.cols),5) if abs(A[:,c].det())==2)
T=s.Matrix([[(a*b).trace() for b in basis] for a in basis])
C=s.Matrix.hstack(*[flat(a*x-x*a) for a in basis]);centraldim=28-C.rank()
chain=[I[:,4]]
for _ in range(4):chain.append(x*chain[-1])
chain += [I[:,2]-I[:,6],I[:,3],I[:,7]]
g=I[:,[0,1,2,7,4,5,6,3]];u=(I+x)*(I-x).inv()
coeff=[int((xd*b).trace()) for b in g2]
assert all(k%2==0 for k in grades)
assert coeff==[-2,-2,-2,0,0] and abs(T.det())==2**28
assert centraldim==8 and A.rank()==5
assert g*x*g==x and g*xd*g==xd and g*H*g==H
assert (u-I)*(u+I).inv()==x and g*u*g==u
assert abs(s.Matrix.hstack(*chain).det())==8 and x**5==zero and x**4!=zero
r={'grading_dimensions':{str(k):grades.count(k) for k in sorted(set(grades))},'centralizer_lie_dimension_Q':centraldim,'bracket_matrix':A.tolist(),'minor_columns_and_value':minor,'trace_determinant':int(T.det()),'chain_determinant':int(s.Matrix.hstack(*chain).det()),'dual_functional':coeff,'graph_permutation':[g2.index(g*b*g)+1 for b in g2],'graph_and_cayley_checks':True,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(O/'cubic_functional_exact.json').write_text(json.dumps(r,indent=2,default=int))
print(json.dumps(r,default=int))

