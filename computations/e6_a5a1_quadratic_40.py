"""Exact order-two A5A1 root and local coefficient interface in E6.

Uses the independently generated norm-two roots from e6_root_independent_40.
All calculations are on one integral root datum, with no finite-q sampling.
"""
from pathlib import Path
import runpy,json
import sympy as s
here=Path(__file__).resolve().parent
d=runpy.run_path(str(here/'e6_root_independent_40.py'))
C,Id,ref,roots,index,gamma,longest,B=[d[k] for k in ['C','Id','ref','roots','index','gamma','longest','B']]
theta=s.Matrix([1,2,2,3,2,1])
assert C*theta==Id[:,1]
assert gamma*theta==theta
fixed=[a for a in roots if int((s.Matrix(a).T*C*theta)[0])%2==0]
assert len(fixed)==32
perp=[a for a in roots if int((s.Matrix(a).T*C*theta)[0])==0]
assert len(perp)==30
assert set(fixed)==set(perp)|{tuple(theta),tuple(-theta)}
alpha=[Id[:,i] for i in [0,2,3,4,5]]+[-theta]
gram=s.Matrix.hstack(*alpha).T*C*s.Matrix.hstack(*alpha)
assert gram[5,5]==2 and gram[:5,5]==s.zeros(5,1)
assert abs(s.Matrix.hstack(*alpha).det())==2
rt=Id-theta*theta.T*C
Rgens=[ref[i] for i in [0,2,3,4,5]]+[rt]

def perm(r):return 'PermList('+str([index[tuple(r*s.Matrix(a))] for a in roots])+')'
def mat(r):return str([[int(v) for v in row] for row in r.tolist()])
gap='W:=Group(['+','.join(perm(r) for r in ref)+']);;\n'
gap+='R:=Group(['+','.join(perm(r) for r in Rgens)+']);;\n'
gap+='v:='+perm(longest)+';;\n'
gap+='R2:=Centralizer(R,v);;\n'
gap+='wp:=GeneratorsOfGroup(W);;\n'
gap+='wm:=List(['+','.join(mat(r.T) for r in ref)+'],r->r*One(GF(2)));;\n'
gap+='Wmod2:=Group(wm);; hom:=GroupHomomorphismByImages(W,Wmod2,wp,wm);;\n'
gap+='s0:='+str(list(theta))+'*One(GF(2));;\n'
gap+='St:=Stabilizer(Wmod2,s0,OnPoints);;\n'
gap+='if Size(W)<>51840 or Size(R)<>1440 or not v in R or Normalizer(W,R)<>R or Size(R2)<>96 or Size(Wmod2)<>51840 or Image(hom,R)<>St then Error("quadratic root anchor");fi;\n'
gap+='Print("QUADRATIC_ROOT|R=",Size(R),"|R2=",Size(R2),"|normalizer=",Size(Normalizer(W,R)),"|s_stabilizer=",Size(St),"|PASSED\\n");\n'
# Independent exact real-rank-four fixed-space closure and extension test.
gap+='wmQ:=['+','.join(mat(r.T) for r in ref)+'];;\n'
gap+='rho:=GroupHomomorphismByImages(W,Group(wmQ),wp,wmQ);;\n'
gap+='vm:=Image(rho,v);; pe:=(IdentityMat(6)-vm)/2;;\n'
gap+='classes:=ConjugacyClassesSubgroups(R2);; count:=0;;\n'
gap+='for cl in classes do\n'
gap+=' H:=Representative(cl); ph:=Sum(Elements(H),h->Image(rho,h))/Size(H)*pe;\n'
gap+=' closure:=Filtered(Elements(R2),g->ph*Image(rho,g)=ph);\n'
gap+=' if Length(closure)=Size(H) then\n'
gap+='  N:=Normalizer(R2,H);\n'
gap+='  for eta in Irr(H) do\n'
gap+='   J:=InertiaSubgroup(N,eta);\n'
gap+='   if not ForAny(Irr(J),x->RestrictedClassFunction(x,H)=eta) then Error("ordinary inertia extension");fi;\n'
gap+='  od; count:=count+1;\n'
gap+='  Print("A5A1_FOLDED_PARABOLIC|H=",StructureDescription(H),"|N/H=",StructureDescription(N/H),"|fixedrank=",RankMat(ph),"|ALL_ORDINARY_EXTEND\\n");\n'
gap+=' fi;\nod;\n'
gap+='Print("A5A1_FOLDED_ALL|classes=",count,"|PASSED\\n");\nQUIT;\n'
(here/'e6_a5a1_quadratic_40.g').write_text(gap)
out={'root_count_fixed_by_s':32,'orthogonal_A5_root_count':30,
     'theta_is_fundamental_coweight_2':True,'gamma_fixes_s':True,
     'root_lattice_index_R_in_E6':2,'cartan_R':gram.tolist(),
     'scope':'Order-two A5A1 pseudo-Levi, odd q; good ell e=1 or 2 family interface.'}
(here/'e6_a5a1_quadratic_40.json').write_text(json.dumps(out,indent=2,default=int)+'\n')
print(json.dumps({k:out[k] for k in ['root_count_fixed_by_s','orthogonal_A5_root_count','root_lattice_index_R_in_E6']}))
