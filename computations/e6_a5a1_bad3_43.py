"""Uniform E6 A5A1 three-block witnesses, reconstructed at node43.

Integral root datum only; no finite fields q or valuations are sampled.
The GF(3) calculation certifies the characteristic-torus assertion for
every valuation via the order-three subgroup of the actual torus.
"""
from pathlib import Path
from itertools import product
import json
import sympy as sp

outdir=Path(__file__).resolve().parent
I=sp.eye(6)
C=2*I
for i,j in [(0,2),(2,3),(3,4),(4,5),(1,3)]:
    C[i,j]=C[j,i]=-1
theta=sp.Matrix([1,2,2,3,2,1])
roots=[]
for a in product(*[range(int(x)+1) for x in theta]):
    x=sp.Matrix(a)
    if (x.T*C*x)[0]==2:
        roots.extend([x,-x])
roots=sorted(roots,key=lambda x:tuple(x))
assert len(roots)==72
idx={tuple(x):i+1 for i,x in enumerate(roots)}
ref=[I-I[:,i]*(I[:,i].T*C) for i in range(6)]
gamma=sp.zeros(6)
for i,j in enumerate([5,1,4,3,2,0]):gamma[j,i]=1
w0=-gamma
assert C*theta==I[:,1] and gamma*theta==theta
rt=I-theta*(theta.T*C)
Rgens=[ref[i] for i in [0,2,3,4,5]]+[rt]
Rbasis=sp.Matrix.hstack(*[I[:,i] for i in [0,2,3,4,5]],theta)
assert abs(Rbasis.det())==2
involution_roots=[x for x in roots if int((x.T*C*theta)[0])%2==0]
assert len(involution_roots)==32

def rank3(M):
    a=[[int(x)%3 for x in r] for r in M.tolist()]
    n,m=len(a),len(a[0]) if a else 0
    k=0
    for j in range(m):
        rows=[r for r in range(k,n) if a[r][j]]
        if not rows:continue
        r=rows[0];a[k],a[r]=a[r],a[k]
        inv=pow(a[k][j],-1,3);a[k]=[(x*inv)%3 for x in a[k]]
        for r in range(n):
            if r!=k:
                c=a[r][j];a[r]=[(a[r][l]-c*a[k][l])%3 for l in range(m)]
        k+=1
        if k==n:break
    return k

z=next(sp.Matrix(a) for a in product(range(3),repeat=6)
       if any(a) and all(int(x)%3==0 for x in C*sp.Matrix(a)))
assert rank3(C)==5
assert all(all(int(x)%3==0 for x in (r-I)*z) for r in ref)
assert all(int(x)%3==0 for x in gamma*z+z)
u=ref[0]*ref[2]
v=ref[4]*ref[5]
assert u**3==I and v**3==I and u*v==v*u
assert w0*u*w0==v**-1
comm=sp.Matrix.hstack(u-I,v-I)
assert rank3(comm)==rank3(sp.Matrix.hstack(comm,z))
square_checks=[]
for i,j in product(range(3),repeat=2):
    if i==j==0:continue
    g=u**i*v**j
    sq=(g-I)**2
    rank=rank3(sq)
    rz=rank3(sp.Matrix.hstack(z,sq))
    assert rz>1
    square_checks.append({'powers':[i,j],'square_rank':rank,'rank_with_center':rz})

# A basis with central first column; quotient matrices are integral mod 3.
cols=[z]
for i in range(6):
    if rank3(sp.Matrix.hstack(*cols,I[:,i]))>len(cols):cols.append(I[:,i])
    if len(cols)==6:break
B=sp.Matrix.hstack(*cols)
assert int(B.det())%3
den=int(B.det())
Binv=(pow(den%3,-1,3)*B.adjugate()).applyfunc(lambda x:int(x)%3)
quot=[]
for r in ref:
    m=(Binv*r*B).applyfunc(lambda x:int(x)%3)
    assert list(m[:,0])==[1,0,0,0,0,0]
    quot.append(m[1:,1:])

E=sp.Matrix.hstack(I[:,0]+I[:,5],I[:,2]+I[:,4],I[:,1],I[:,3])
assert gamma*E==E
assert all(any(int(t)%3 for t in (x.T*C*E)) for x in roots)
un=u*v**-1
assert un*w0==w0*un
unE=E.gauss_jordan_solve(un*E)[0]
assert (unE-sp.eye(4))**3==sp.zeros(4) or unE**3==sp.eye(4)
assert rank3((unE-sp.eye(4))**2)>0

def mat(M):return str([[int(x) for x in row] for row in M.tolist()])
def perm(M):return 'PermList('+str([idx[tuple(M*x)] for x in roots])+')'
g=[]
g.append('wp:=['+','.join(perm(r) for r in ref)+'];; W:=Group(wp);;')
g.append('R:=Group(['+','.join(perm(r) for r in Rgens)+']);;')
g.append('w0:='+perm(w0)+';; R2:=Centralizer(R,w0);; WF:=Centralizer(W,w0);;')
g.append('qgens:=List(['+','.join(mat(r.T) for r in quot)+'],x->x*One(GF(3)));;')
g.append('Wbar:=Group(qgens);;')
g.append('if Size(W)<>51840 or Size(R)<>1440 or Size(R2)<>96 or Size(WF)<>1152 or Size(Wbar)<>51840 then Error("root or quotient action");fi;')
g.append('Print("ROOT_AND_CENTRAL_QUOTIENT|W=",Size(W),"|R=",Size(R),"|R2=",Size(R2),"|WF=",Size(WF),"|Wbar=",Size(Wbar),"|PASSED\\n");')
g.append('wmQ:=['+','.join(mat(r.T) for r in ref)+'];; rho:=GroupHomomorphismByImages(W,Group(wmQ),wp,wmQ);;')
g.append('vm:=Image(rho,w0);; pe:=(IdentityMat(6)-vm)/2;;')
g.append('cnt:=0;; for cl in ConjugacyClassesSubgroups(R2) do')
g.append(' H:=Representative(cl); ph:=Sum(Elements(H),h->Image(rho,h))/Size(H)*pe;')
g.append(' closure:=Filtered(Elements(R2),x->ph*Image(rho,x)=ph);')
g.append(' if Length(closure)=Size(H) then')
g.append('  N:=Normalizer(R2,H); for eta in Irr(H) do J:=InertiaSubgroup(N,eta); if not ForAny(Irr(J),x->RestrictedClassFunction(x,H)=eta) then Error("ordinary coefficient");fi; od;')
g.append('  cnt:=cnt+1; Print("ORDINARY_PARABOLIC|",StructureDescription(H),"|N/H=",StructureDescription(N/H),"|PASSED\\n");')
g.append(' fi; od; Print("ORDINARY_PARABOLIC_CLASSES|",cnt,"|PASSED\\n");')
g.append('QUIT;')
(outdir/'e6_a5a1_bad3_43.g').write_text('\n'.join(g)+'\n')
result={'root_count':72,'involution_root_count':32,'A5A1_root_index':2,
        'center_vector_mod3':list(z),'center_in_sylow_commutator':True,
        'matching_square_checks':square_checks,
        'full_W_quotient_dimension':5,'nonmatching_square_rank':rank3((unE-sp.eye(4))**2),
        'nonmatching_root_restrictions_nonzero_mod3':True,
        'all_q_and_valuation_scope':True}
(outdir/'e6_a5a1_bad3_43.json').write_text(json.dumps(result,indent=2,default=int)+'\n')
print(json.dumps(result,default=int))
