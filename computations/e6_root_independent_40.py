"""Independent integral E6/F4 certificate, no sampled finite-field parameters.

Roots are found by the positive norm-two equation in the highest-root box,
not by copying the root orbit routine of node37.  GAP receives these roots
and separately certifies the full Weyl centralizer and its mod-three image.
"""
from itertools import product
from pathlib import Path
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

here=Path(__file__).resolve().parent
C=s.Matrix([[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],
            [0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
highest=(1,2,2,3,2,1)
pos=[]
for a in product(*(range(n+1) for n in highest)):
    norm=2*sum(x*x for x in a)-2*(a[0]*a[2]+a[1]*a[3]+a[2]*a[3]+a[3]*a[4]+a[4]*a[5])
    if norm==2:pos.append(a)
assert len(pos)==36
roots=sorted(pos+[tuple(-x for x in a) for a in pos])
index={a:i+1 for i,a in enumerate(roots)}
Id=s.eye(6)
ref=[]
for i in range(6):
    r=Id-Id[:,i]*C[i,:]
    assert r*r==Id and r.T*C*r==C
    assert all(tuple(r*s.Matrix(a)) in index for a in roots)
    ref.append(r)
gamma=s.zeros(6)
for i,j in enumerate((5,1,4,3,2,0)):gamma[j,i]=1
longest=-gamma
assert all(all(x<=0 for x in longest*s.Matrix(a)) for a in pos)
B=s.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1],[0,1,0,0],[1,0,0,0]])
assert gamma*B==B
assert list(smith_normal_form(B,domain=s.ZZ).diagonal())==[1,1,1,1]
restrictions=[tuple((s.Matrix(a).T*C*B)) for a in roots]
assert all(any(int(x)%3 for x in r) for r in restrictions)
fold=[ref[0]*ref[5],ref[2]*ref[4],ref[3],ref[1]]
coef=[s.Matrix.vstack((r*B)[0,:],(r*B)[2,:],(r*B)[1,:],(r*B)[3,:]) for r in fold]
assert all(B*a==r*B for a,r in zip(coef,fold))
c=fold[0]*fold[1]*fold[2]*fold[3]
u=c**4
assert c**12==Id and u**3==Id and u!=Id
assert u*longest==longest*u
char=u.charpoly().as_expr()
x=u.charpoly().gen
assert s.expand(char-(x*x+x+1)**3)==0
snf=[abs(int(v)) for v in smith_normal_form(Id-u,domain=s.ZZ).diagonal()]
assert snf==[1,1,1,3,3,3]

def perm(r):return 'PermList('+str([index[tuple(r*s.Matrix(a))] for a in roots])+')'
def mat(r):return str([[int(v) for v in row] for row in r.tolist()])
gap='W:=Group(['+','.join(perm(r) for r in ref)+']);;\n'
gap+='v:='+perm(longest)+';;\n'
gap+='C:=Centralizer(W,v);;\n'
gap+='R:=Group(['+','.join(perm(r) for r in fold)+']);;\n'
gap+='R3:=Group(List(['+','.join(mat(r) for r in coef)+'],r->r*One(GF(3))));;\n'
gap+='if Size(W)<>51840 or Size(C)<>1152 or C<>R or Size(R3)<>1152 then Error("independent E6 certificate");fi;\n'
gap+='Print("INDEPENDENT_E6|W=",Size(W),"|centralizer=",Size(C),"|mod3=",Size(R3),"|PASSED\\n");\nQUIT;\n'
(here/'e6_root_independent_40.g').write_text(gap)
out={'positive_roots_from_norm_equation':len(pos),'total_roots':len(roots),
     'primitive_fixed_lattice_rank':4,'all_roots_nonzero_mod3':True,
     'u_characteristic_polynomial':str(char),'smith_I_minus_u':snf,
     'u_matrix':u.tolist(),'folded_generators_mod3':[r.tolist() for r in coef],
     'conclusion':'T^u is mu_3^3, stronger than only det(1-u)=27; full Weyl assertions are separate GAP checks.'}
(here/'e6_root_independent_40.json').write_text(json.dumps(out,indent=2,default=int)+'\n')
print(json.dumps({k:out[k] for k in ['positive_roots_from_norm_equation','total_roots','primitive_fixed_lattice_rank','all_roots_nonzero_mod3','smith_I_minus_u']}))
