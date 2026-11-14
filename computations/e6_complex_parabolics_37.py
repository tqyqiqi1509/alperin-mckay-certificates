"""Exact ordinary coefficient calculation in E6 generic Sylow-torus groups.

This is a finite root-datum calculation, not a finite-q experiment. It
determines every pointwise stabilizer in the characteristic-zero regular
eigenspaces, its full normalizer, and actual ordinary-character extensions.
"""
from pathlib import Path
from collections import deque
from collections import Counter
import json
import sympy as S
C=2*S.eye(6)
for i,j in [(0,2),(2,3),(3,4),(4,5),(1,3)]: C[i,j]=C[j,i]=-1
I=S.eye(6);rs=[]
for i in range(6):
    r=I.copy()
    for j in range(6):r[i,j]-=C[i,j]
    rs.append(r)
roots={tuple(I[:,i]) for i in range(6)}|{tuple(-I[:,i]) for i in range(6)}
todo=deque(roots)
while todo:
    r=S.Matrix(todo.popleft())
    for a in rs:
        t=tuple(a*r)
        if t not in roots:roots.add(t);todo.append(t)
assert len(roots)==72
rootlist=sorted(roots);index={r:i+1 for i,r in enumerate(rootlist)}
def gp(a):return '['+','.join(str(index[tuple(a*S.Matrix(r))]) for r in rootlist)+']'
def gm(a):return '['+','.join('['+','.join(str(a[i,j]) for j in range(6))+']' for i in range(6))+']'
c=rs[0]*rs[5]*rs[2]*rs[4]*rs[3]*rs[1]
vs=[c**4,c**3,c**2]
for v in vs: assert v.det()==1
data={'status':'passed','root_count':72,'generic_tori':[]}
for d,v in zip([3,4,6],vs):
    pe=I if d==3 else ((I-v**2)/2 if d==4 else (I-v**3)/2)
    assert pe*pe==pe
    assert not [a for a in roots if S.Matrix(a).T*C*pe==S.zeros(1,6)]
    data['generic_tori'].append({'e':d,'weyl_matrix':v.tolist(),
       'characteristic_polynomial':str(v.charpoly().as_expr()),
       'cyclotomic_projector':pe.tolist(),
       'no_root_vanishes_on_regular_eigenspace':True})
B=S.Matrix([[2,0],[1,1],[2,1],[2,2],[1,2],[0,2]])
source=[1,2,3,4]
assert C[source,:]*B==S.zeros(4,2)
assert not [a for a in roots if all((S.Matrix(a).T*C)[0,j]==0 for j in source)]
restr=Counter(tuple(S.Matrix(a).T*C*B) for a in roots)
assert restr==Counter({(0,0):24,(2,-1):8,(-2,1):8,(-1,2):8,
                      (1,-2):8,(1,1):8,(-1,-1):8})
data['d4_centralizer_has_no_roots']=True
data['d4_center_basis']=B.tolist()
data['d4_center_root_restrictions']={str(k):v for k,v in restr.items()}
Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2,default=str)+'\n')
gap='''# Root permutations and exact rational row matrices.
'''
gap+='wp:=List(['+','.join(gp(r) for r in rs)+'],PermList);;\n'
gap+='wm:=['+','.join(gm(r.T) for r in rs)+'];;\n'
gap+='W:=Group(wp);; WM:=Group(wm);; rho:=GroupHomomorphismByImages(W,WM,wp,wm);;\n'
gap+='if Size(W)<>51840 or Size(Kernel(rho))<>1 then Error("root representation");fi;\n'
gap+='vs:=List(['+','.join(gp(v) for v in vs)+'],PermList);;\n'
gap+='''
ds:=[3,4,6];;
for j in [1..3] do
  v:=vs[j]; G:=Centralizer(W,v); vm:=Image(rho,v); one:=IdentityMat(6);
  if j=1 then pe:=one;
  elif j=2 then pe:=(one-vm^2)/2;
  else pe:=(one-vm^3)/2;fi;
  classes:=ConjugacyClassesSubgroups(G);
  Print("GROUP|e=",ds[j],"|order=",Size(G),"|subgroupclasses=",Length(classes),"\\n");
  for cl in classes do
    H:=Representative(cl);
    ph:=Sum(Elements(H),h->Image(rho,h))/Size(H)*pe;
    rank:=RankMat(ph);
    if rank>0 then
      closure:=Filtered(Elements(G),g->ph*Image(rho,g)=ph);
      if Length(closure)=Size(H) then
        N:=Normalizer(G,H); ctH:=CharacterTable(H); ih:=Irr(ctH);
        exts:=[];
        for eta in ih do
          J:=InertiaSubgroup(N,eta);
          ij:=Irr(J);
          ext:=Filtered(ij,x->RestrictedClassFunction(x,H)=eta);
          if Length(ext)=0 then Error("ordinary extension failed");fi;
          Add(exts,[eta[1],Size(J),Length(ext)]);
        od;
        Print("PARABOLIC|order=",Size(H),"|id=",IdGroup(H),"|fixedrank=",rank/2,
          "|normalizer=",Size(N),"|quotient=",StructureDescription(N/H),
          "|ordinaryextensions=",exts,"\\n");
      fi;
    fi;
  od;
od;
Print("ALL ORDINARY PARABOLIC COEFFICIENTS EXTEND: PASSED\\n");
QUIT;
'''
Path(__file__).with_suffix('.g').write_text(gap)
print('wrote',Path(__file__).with_suffix('.g'))
