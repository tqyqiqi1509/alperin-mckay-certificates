from pathlib import Path
import numpy as np, sympy as s,json
C=2*np.eye(6,dtype=np.int64)
for i,j in [(0,2),(2,3),(3,4),(4,5),(1,3)]:C[i,j]=C[j,i]=-1
wt={(1,0,0,0,0,0)};todo=list(wt)
while todo:
 v=todo.pop()
 for i in range(6):
  z=tuple(np.array(v)-v[i]*C[:,i])
  if z not in wt:wt.add(z);todo.append(z)
wt=sorted(wt);ix={v:i for i,v in enumerate(wt)};d=len(wt);I=np.eye(d,dtype=np.int64)
E=[];N=[]
for i in range(6):
 e=np.zeros((d,d),dtype=np.int64)
 for j,v in enumerate(wt):
  if v[i]==-1:e[ix[tuple(np.array(v)+C[:,i])],j]=1
 E.append(e);N.append((I+e)@(I-e.T)@(I+e))
count=0
for i in range(6):
 for j in range(6):
  assert np.array_equal(E[i]@E[j].T-E[j].T@E[i],np.diag([v[i] for v in wt]) if i==j else np.zeros_like(I));count+=1
  if i!=j:
   if C[i,j]==0:assert np.array_equal(E[i]@E[j],E[j]@E[i])
   else:assert not np.any(E[i]@E[j]@E[i])
   count+=1

def word(w):
 x=I.copy()
 for i in w:x=x@N[i-1]
 return x
q1=word([1,3,4,2,5,4,3,1]);q2=word([6,5,4,2,3,4,5,6])
assert np.array_equal(q1@q1,I) and np.array_equal(q2@q2,I) and np.array_equal(np.linalg.matrix_power(q1@q2,3),I)
pins=[]
for q in [q1,q2]:
 pins.append({str(i+1):[j+1 for j in [1,2,3,4] if np.array_equal(q@E[i]@q.T,E[j])] for i in [1,2,3,4]})
assert pins==[{'2':[5],'3':[3],'4':[4],'5':[2]},{'2':[3],'3':[2],'4':[4],'5':[5]}]
# Central-character lattice and actual torus maps from the Weyl words.
S=[np.eye(6,dtype=np.int64)-np.outer(np.eye(6,dtype=np.int64)[:,i],C[i]) for i in range(6)]
def ww(w):
 x=np.eye(6,dtype=np.int64)
 for i in w:x=x@S[i-1]
 return s.Matrix(x)
b=s.Matrix([[2,0],[1,1],[2,1],[2,2],[1,2],[0,2]])
m1=b.gauss_jordan_solve(ww([1,3,4,2,5,4,3,1])*b)[0];m2=b.gauss_jordan_solve(ww([6,5,4,2,3,4,5,6])*b)[0]
G=s.eye(6)[[5,1,4,3,2,0],:];mg=b.gauss_jordan_solve(G*b)[0]
assert (m1*m2*m1).inv()*mg==-s.eye(2)
out={'weights':d,'Chevalley_Serre_checks':count,'q1q2_genuine_S3':True,'positive_D4_pins':pins,'torus_q1':m1.tolist(),'torus_q2':m2.tolist(),'corrected_graph_torus_minus_identity':True}
Path('outputs/e6_D4_carrier_round3_checks.json').write_text(json.dumps(out,indent=2,default=int));print(out)
