"""Reconstruct the marked D4 and toral bad-three geometry in fresh A7 coordinates."""
from pathlib import Path
import contextlib, io, json
from itertools import product
from functools import reduce
from math import gcd
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
with contextlib.redirect_stdout(io.StringIO()):
    import independent_carrier as c
OUT=Path(__file__).resolve().parent
roots=sorted(c.roots); ri={r:i for i,r in enumerate(roots)}
coords={r:s.Matrix(c.simple_coords[r]) for r in roots}
C=c.C; ci=C.inv();
fixed=[r for r in roots if coords[r][0]%2==0]
assert len(fixed)==62
phases=[]
for w in c.weights:
    fc=s.Matrix([c.pairing(w,a) for a in c.alpha])
    a=(fc.T*ci[:,0])[0]/2; b=(fc.T*ci[:,6])[0]
    assert (2*a).q==1 and (2*b).q==1
    # Parity evaluation stays integral even for negative exponents.
    phases.append((1-2*(int(2*a)%2),1-2*(int(2*(a+b))%2)))
tr=[sum(z[i] for z in phases) for i in [0,1]]
assert tr==[8,-8]
H=[1,2,3,4]
hroots=[r for r in roots if all(coords[r][i]==0 for i in [0,5,6])]
assert len(hroots)==24 and all(r in fixed for r in hroots)
ker3=[r for r in roots if all(coords[r][i]%3==0 for i in [0,5,6])]
assert set(ker3)==set(hroots)
perp=[r for r in roots if all(sum(x*y for x,y in zip(r,c.alpha[i]))==0 for i in H)]
assert len(perp)==6 and all(r in fixed for r in perp)
kp=[r for r in perp if all(z>=0 for z in coords[r])]
assert len(kp)==3
assert all(sum(x*y for x,y in zip(kp[i],kp[j]))==0 for i,j in product(range(3),repeat=2) if i!=j)

def snf(M):
    D=smith_normal_form(M,domain=ZZ)
    return [abs(int(D[i,i])) for i in range(min(D.shape)) if D[i,i]]
QH=s.Matrix.hstack(*(coords[c.alpha[i]] for i in H))
QK=s.Matrix.hstack(*(coords[r] for r in kp))
J=s.Matrix.hstack(QH,QK)
assert abs(J.det())==4
assert snf(QK)==[1,1,1]
# Characters alpha_i have fundamental-weight coordinates C[:,i].
assert snf(C*QH)==[1,1,1,1]
# Finite gluing Q(E7)/(Q(D4)+Q(A1^3)); project to the D4 weight lattice.
columns=[tuple(z%1 for z in J.inv()[:,i]) for i in range(7)]
glue={(s.Rational(0),)*7};todo=list(glue)
while todo:
    v=todo.pop()
    for w in columns:
        z=tuple((a+b)%1 for a,b in zip(v,w))
        if z not in glue:glue.add(z);todo.append(z)
assert len(glue)==4
hp={v[:4] for v in glue};kpgl={v[4:] for v in glue}
assert len(hp)==len(kpgl)==4

def perm_for_root(r):
    ans=[]
    for x in roots:
        k=sum(a*b for a,b in zip(x,r))//4
        y=tuple(a-k*b for a,b in zip(x,r))
        ans.append(ri[y]+1)
    return ans
def gapgroup(rs):return 'Group('+','.join('PermList('+str(perm_for_root(r))+')' for r in rs)+')'
code='W:='+gapgroup(c.alpha)+';;\n'
code+='C:='+gapgroup(fixed)+';;\n'
code+='H:='+gapgroup([c.alpha[i] for i in H])+';;\n'
code+='P:='+gapgroup(kp)+';;\n'
code+='basis:='+str([ri[c.alpha[i]]+1 for i in H])+';;\n'
code+='H2:='+gapgroup([c.alpha[i] for i in [0,2]])+';;\n'
def rootperm_of_signed_weight_operator(a):
    ans=[]
    for r in roots:
        w=next(w for w in c.weights if c.pairing(w,r)==-1)
        t=c.move(w,r)
        wa=c.weights[abs(a[c.wi[w]])-1]
        ta=c.weights[abs(a[c.wi[t]])-1]
        image=tuple((x-y)//2 for x,y in zip(ta,wa))
        assert image in ri
        ans.append(ri[image]+1)
    return ans
code+='v3:=PermList('+str(rootperm_of_signed_weight_operator(c.v3))+');;\n'
code+='v6:=PermList('+str(rootperm_of_signed_weight_operator(c.v6))+');;\n'
code+=r'''
Pw:=Stabilizer(W,basis,OnTuples);;
Nc:=Normalizer(C,H);; Nw:=Normalizer(W,H);;
Print("W,C,H,P orders=",List([W,C,H,P],Size),"\n");
Print("D4 pointwise stabilizer order=",Size(Pw)," equals perpendicular=",Pw=P,"\n");
Print("Marked D4 normalizer quotients in C,W=",[Index(Nc,H),Index(Nw,H)],"\n");
if Pw<>P or Index(Nc,H)<>16 or Index(Nw,H)<>48 then Error("Geometry mismatch"); fi;
N2:=Normalizer(W,H2);; hom:=NaturalHomomorphismByNormalSubgroup(N2,H2);;
Q2:=Image(hom);;
Print("A2 geometric normalizer quotient=",Size(Q2),"\n");
for v in [v3,v6] do
  Fcentral56:=Centralizer(Q2,Image(hom,v));;
  Print("A2 rational fixed quotient=",Size(Fcentral56)," structure=",StructureDescription(Fcentral56),"\n");
  if Size(Fcentral56)<>36 then Error("Wrong A2 rational quotient"); fi;
od;
QUIT;
'''
(OUT/'independent_bad3.g').write_text(code)
record={'centralizer_roots':len(fixed),'minuscule_traces_s_zs':tr,
        'D4_source_roots':len(hroots),'D4_mod3_kernel_exact':True,
        'perpendicular_roots':len(perp),'perpendicular_bases_scaled2':kp,
        'D4_source_center_smith_in_P':snf(C*QH),
        'A1cube_coroot_primitive_smith':snf(QK),'joint_lattice_index':abs(int(J.det())),
        'glue_D4_projection_order':len(hp),'glue_A1cube_projection_order':len(kpgl),
        'glue_fraction_coordinates':[[str(z) for z in v] for v in sorted(glue)]}
(OUT/'independent_bad3.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
