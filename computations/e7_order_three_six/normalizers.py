"""Emit independent GAP4 root-permutation tests from node52's own roots."""
import independent_e7 as m
from pathlib import Path
import sympy as S

ri={tuple(r):i+1 for i,r in enumerate(m.roots)}
rp=lambda A:'PermList('+str([ri[tuple(A*r)] for r in m.roots])+')'
refs=[rp(m.ref(r)) for r in m.roots]
base=[ri[tuple(r)] for r in m.base]
g=['refs52:=['+','.join(refs)+'];;',
   'W52:=Group(List('+str(base)+',i->refs52[i]));;',
   'Print("WE7|",Size(W52),"\\n");',
   'Based52:=function(x,bb,pos) local b; while ForAny(bb,b->not b^x in pos) do b:=First(bb,b->not b^x in pos);x:=refs52[b]*x;od;return x;end;;']
for name,dat in m.cases.items():
 s=S.Matrix(dat['s_Euclidean']);Croots=[r for r in m.roots if m.dot(r,s).q==1]
 cb=m.simples(Croots);kc=[S.Matrix(r) for r in dat['K_basis_Euclidean']]
 lc=[S.Matrix(r) for r in dat['L_basis_Euclidean']]
 Kroots=[r for r in Croots if all(m.dot(r,m.B*z)==0 for z in ({1:S.Matrix(dat['F_on_Q'])-m.I,2:S.Matrix(dat['F_on_Q'])+m.I,3:S.Matrix(dat['F_on_Q'])**2+S.Matrix(dat['F_on_Q'])+m.I,4:S.Matrix(dat['F_on_Q'])**2+m.I,6:S.Matrix(dat['F_on_Q'])**2-S.Matrix(dat['F_on_Q'])+m.I}[dat['d']]).nullspace())]
 # Use the complete source root set recovered from the simple roots.
 Lspan=S.Matrix.hstack(*lc) if lc else S.zeros(8,0)
 Lroots=[r for r in m.roots if Lspan.row_join(r).rank()==len(lc)]
 Fc=S.Matrix(dat['F_on_Q']);F=m.B*Fc*m.left+(m.I8-m.B*m.left)
 kgens='['+','.join('refs52['+str(ri[tuple(r)])+']' for r in kc)+']'
 lgens='['+','.join('refs52['+str(ri[tuple(r)])+']' for r in lc)+']'
 # GAP Group([]) requires an identity; use Group([()]) for empty source.
 if not kc:kgens='[()]'
 if not lc:lgens='[()]'
 g += ['Print("CASE|'+name+'\\n");',
       'WC:=Group(List('+str([ri[tuple(r)] for r in cb])+',i->refs52[i]));;',
       'WK:=Group('+kgens+');;WL:=Group('+lgens+');;',
       'F52:='+rp(F)+';;',
       'NC:=Normalizer(WC,WK);;NN:=Group(Concatenation(GeneratorsOfGroup(NC),[F52]));;',
       'pi:=NaturalHomomorphismByNormalSubgroup(NN,WK);;',
       'Q:=Intersection(Image(pi,NC),Centralizer(Image(pi),Image(pi,F52)));;',
       'Print("SOURCE_REL|",Size(WC),"|",Size(WK),"|",Size(Q),"|",StructureDescription(Q),"\\n");',
       'NL:=Normalizer(W52,WL);;NNL:=Group(Concatenation(GeneratorsOfGroup(NL),[F52]));;',
       'pil:=NaturalHomomorphismByNormalSubgroup(NNL,WL);;',
       'QL:=Intersection(Image(pil,NL),Centralizer(Image(pil),Image(pil,F52)));;',
       'Print("AMBIENT_REL|",Size(QL),"|",StructureDescription(QL),"\\n");']
 if name=='D6A1_d4':
  pairs=sorted(sorted([ri[tuple(r)],ri[tuple(-r)]]) for r in kc)
  g+=['Kpairs:='+str(pairs)+';;QQ:=PreImage(pi,Q);;',
      'act:=Action(QQ,Kpairs,OnSets);;',
      'Print("K_FACTOR_ACTION|",Size(act),"|",GeneratorsOfGroup(act),"\\n");',
      'mixed:=Filtered(Elements(Q),x->OnSets(Kpairs[2],PreImagesRepresentative(pi,x))=Kpairs[2]);;',
      'mixed:=Subgroup(Q,mixed);;',
      'Print("MIXED_MARK_KERNEL|",Size(mixed),"|",StructureDescription(mixed),"\\n");']

# Quasi-isolated order6 rational normalizer and exact based representatives.
cb=m.kb6;cr=m.cr6;cbidx=[ri[tuple(r)] for r in cb]
pos=[ri[tuple(r)] for r in cr if m.dot(m.height,r)>0]
labels=[int(12*m.dot(r,m.s6)) for r in m.roots]
rootcoords=[[int(x) for x in m.left*r] for r in m.roots]
g+=['CR6:='+str(sorted(ri[tuple(r)] for r in cr))+';;',
    'NR6:=Stabilizer(W52,CR6,OnSets);;',
    'WC6:=Group(List('+str(cbidx)+',i->refs52[i]));;',
    'VV6:=Group(List(GeneratorsOfGroup(NR6),x->Based52(x,'+str(cbidx)+','+str(pos)+')));;',
    'Print("ORDER6_NORMALIZER|",Size(NR6),"|",Size(WC6),"|",Size(VV6),"\\n");',
    'labels6:='+str(labels)+';; rootcoords52:='+str(rootcoords)+';; base52:='+str(base)+';;',
    'for q in [1,5,7,11] do for v in Elements(VV6) do',
    'if ForAll([1..126],i->(labels6[i^v]-q*labels6[i]) mod 12=0) then',
    'Print("RATIONAL6|",q,"|",List(base52,i->rootcoords52[i^v]),"\\n");',
    'fi;od;od;',
    'QUIT;']
Path(__file__).with_suffix('.g').write_text('\n'.join(g)+'\n')
