"""Uniform A2-source certificate for the order-two E6 A5A1 series.

The two non-toral source characters have degrees Phi3,q*Phi3 (Ennola:
Phi6,q*Phi6).  This script checks the root interface, every nonzero
root restriction at all ell>=7, the actual relative-Weyl quotient,
and all ordinary coefficient inertia extensions.
"""
from pathlib import Path
from functools import reduce
from math import gcd
import runpy,json
import sympy as s
here=Path(__file__).resolve().parent
d=runpy.run_path(str(here/'e6_root_independent_40.py'))
C,Id,ref,roots,index,gamma,longest=[d[k] for k in ['C','Id','ref','roots','index','gamma','longest']]
theta=s.Matrix([1,2,2,3,2,1]); v=ref[0]*ref[2]*ref[4]*ref[5]
assert v**3==Id and v!=Id
assert longest*v*longest==v**2
P=(2*Id-v-v**2)/3
assert P*P==P and P.rank()==4
fixed=[a for a in roots if s.Matrix(a).T*C*P==s.zeros(1,6)]
simpleL=[Id[:,1],theta-Id[:,1]]
assert len(fixed)==6
assert s.Matrix.hstack(*simpleL).T*C*s.Matrix.hstack(*simpleL)==s.Matrix([[2,-1],[-1,2]])
assert theta==sum(simpleL,s.zeros(6,1))
assert all(v*a==a for a in simpleL)
# Three times the zeta3-projector is (I-v)+zeta3*(v^2-v).
# If a root vanished after any reduction at ell, every coordinate norm
# a^2-a*b+b^2 would be divisible by ell.  Their gcds have only prime3.
norm_gcds=[]
for a in roots:
    ra=s.Matrix(a).T*C*(Id-v); rb=s.Matrix(a).T*C*(v**2-v)
    vals=[int(aa*aa-aa*bb+bb*bb) for aa,bb in zip(ra,rb)]
    g=reduce(gcd,vals,0)
    if a in fixed:assert g==0
    else:
        assert g>0 and set(s.factorint(g))<={3}
        norm_gcds.append(g)
Lgens=[Id-a*a.T*C for a in simpleL]
Rgens=[ref[i] for i in [0,2,3,4,5]]+[Id-theta*theta.T*C]
def perm(r):return 'PermList('+str([index[tuple(r*s.Matrix(a))] for a in roots])+')'
def mat(r):return str([[int(vv) for vv in row] for row in r.tolist()])
gap='wp:=['+','.join(perm(r) for r in ref)+'];; W:=Group(wp);;\n'
gap+='R:=Group(['+','.join(perm(r) for r in Rgens)+']);;\n'
gap+='L:=Group(['+','.join(perm(r) for r in Lgens)+']);;\n'
gap+='v:='+perm(v)+';; CW:=Centralizer(W,v);; CR:=Centralizer(R,v);;\n'
gap+='wm:=['+','.join(mat(r.T) for r in ref)+'];;\n'
gap+='rho:=GroupHomomorphismByImages(W,Group(wm),wp,wm);;\n'
gap+='vm:=Image(rho,v);; pe:=(2*IdentityMat(6)-vm-vm^2)/3;;\n'
gap+='ker:=Group(Filtered(Elements(CW),g->pe*Image(rho,g)=pe));;\n'
gap+='if Size(CW)<>108 or Size(CR)<>36 or Size(L)<>6 or not IsNormal(CW,L) or ker<>L or Size(Intersection(CR,L))<>2 or ClosureGroup(CR,L)<>CW then Error("A2 source quotient");fi;\n'
gap+='qmap:=NaturalHomomorphismByNormalSubgroup(CW,L);; Q:=Image(qmap);;\n'
gap+='if Size(Q)<>18 or IdGroup(Q)<>[18,3] then Error("monomial relative group");fi;\n'
gap+='Print("E6_A5A1_E3|C_W=",Size(CW),"|C_R=",Size(CR),"|W_L=",Size(L),"|Q=",StructureDescription(Q),"|PASSED\\n");\n'
gap+='classes:=ConjugacyClassesSubgroups(Q);; count:=0;;\n'
gap+='for cl in classes do\n'
gap+=' H:=Representative(cl); HH:=PreImage(qmap,H); ph:=Sum(Elements(HH),h->Image(rho,h))/Size(HH)*pe;\n'
gap+=' closure:=Filtered(Elements(CW),g->ph*Image(rho,g)=ph);\n'
gap+=' if Length(closure)=6*Size(H) then\n'
gap+='  N:=Normalizer(Q,H);\n'
gap+='  for eta in Irr(H) do\n'
gap+='   J:=InertiaSubgroup(N,eta);\n'
gap+='   if not ForAny(Irr(J),x->RestrictedClassFunction(x,H)=eta) then Error("ordinary inertia extension");fi;\n'
gap+='  od; count:=count+1;\n'
gap+='  Print("E3_PARABOLIC|H=",StructureDescription(H),"|N/H=",StructureDescription(N/H),"|complexfixedrank=",RankMat(ph)/2,"|ALL_ORDINARY_EXTEND\\n");\n'
gap+=' fi;\nod;\n'
gap+='Print("E3_ALL_COEFFICIENTS|classes=",count,"|PASSED\\n");\nQUIT;\n'
(here/'e6_a5a1_e3_40.g').write_text(gap)
out={'v_matrix':v.tolist(),'v_charpoly':str(v.charpoly().as_expr()),
     'fixed_root_count':6,'simple_A2_roots':[list(a) for a in simpleL],
     'theta_in_source_derived_lattice':True,'root_restriction_norm_gcds':sorted(set(norm_gcds)),
     'no_non_A2_root_vanishes_at_any_ell_ge_7':True,
     'w0_conjugates_v_to_inverse':True,
     'scope':'split e3, twisted e6, odd q, good ell>=7'}
(here/'e6_a5a1_e3_40.json').write_text(json.dumps(out,indent=2,default=int)+'\n')
print(json.dumps({k:out[k] for k in ['fixed_root_count','root_restriction_norm_gcds','no_non_A2_root_vanishes_at_any_ell_ge_7']}))
