LoadPackage("ctbllib");;LoadPackage("perfgrp");;
SetInfoLevel(InfoWarning,0);;
gg:=PerfectGroup(IsPermGroup,40320,4);;aut:=AutomorphismGroup(gg);;
inn:=InnerAutomorphismsAutomorphismGroup(aut);;s2:=SylowSubgroup(aut,2);;
comps:=Filtered(ConjugacyClassesSubgroups(s2),c->Size(Representative(c))=4 and Size(Intersection(Representative(c),inn))=1);;
aa:=Representative(comps[1]);;yy:=SemidirectProduct(aa,gg);;
PrimeExpected63:=function(p) if p=5 then return 4;else return 5;fi;end;;
for prime in [5,7] do
 Print("START_PRIME|",prime,"\n");
  eg:=Embedding(yy,2);;
  g0:=Image(eg,gg);; pp:=SylowSubgroup(gg,prime);;p0:=Image(eg,pp);;
Print("STEP|embeddedSylow|",prime,"\n");
  m0:=Image(eg,Normalizer(gg,pp));;hh:=m0;;ea:=Embedding(yy,1);;
  for autgen in GeneratorsOfGroup(aa) do
    pa:=Image(autgen,pp);;conj:=RepresentativeAction(gg,pa,pp,function(a,x)return a^x;end);;
    Assert(0,conj<>fail);;hhgen:=Image(ea,autgen)*Image(eg,conj);;
    Assert(0,p0^hhgen=p0);;hh:=ClosureGroup(hh,hhgen);;
  od;
Print("STEP|normalizers|",prime,"\n");
  Assert(0,Size(hh)=4*Size(m0) and Centralizer(yy,g0)=Centre(g0));
  Assert(0,Centralizer(g0,p0)=ClosureGroup(p0,Centre(g0)));
  tg:=CharacterTable(g0);; tm:=CharacterTable(m0);;Print("STEP|lowerTables|",prime,"\n");
  ty:=CharacterTable(yy);; th:=CharacterTable(hh);;Print("STEP|upperTables|",prime,"\n");
  fg:=FusionConjugacyClasses(tg,ty);;
  fm:=FusionConjugacyClasses(tm,th);;
  Assert(0,fg<>fail and fm<>fail);
  zl:=First(Elements(Centre(g0)),x->Order(x)=2);;
  cg:=ConjugacyClasses(g0);; cm:=ConjugacyClasses(m0);;
  zg:=PositionProperty(cg,c->zl in c);;
  zm:=PositionProperty(cm,c->zl in c);;
  ig:=Irr(tg);; im:=Irr(tm);;
  gx:=Filtered([1..Length(ig)],j->ig[j][1] mod prime<>0 and ig[j][zg]=-ig[j][1]);;
  mx:=Filtered([1..Length(im)],j->im[j][1] mod prime<>0 and im[j][zm]=-im[j][1]);;
  Assert(0,Length(gx)=Length(mx) and Length(gx)=PrimeExpected63(prime));
  bg3:=PrimeBlocks(tg,prime);;bm3:=PrimeBlocks(tm,prime);;
  Assert(0,Length(Set(gx,j->bg3.block[j]))=1 and Length(Set(mx,j->bm3.block[j]))=1);
  qh:=NaturalHomomorphismByNormalSubgroup(hh,m0);;
  qo:=Image(qh);; qe:=Elements(qo);;
  Assert(0,StructureDescription(qo)="C2 x C2");
  lifts:=List(qe,q->PreImagesRepresentative(qh,q));;
  RowPerm63:=function(ct,cls,h)
    local cp,j,vals;
    cp:=List(cls,c->PositionProperty(cls,e->Representative(c)^h in e));
    Assert(0,not fail in cp);
    vals:=List(Irr(ct),ValuesOfClassFunction);
    return List(vals,v->Position(vals,v{cp}));
  end;;
  gaperms:=List(lifts,h->RowPerm63(tg,cg,h));;
  maperms:=List(lifts,h->RowPerm63(tm,cm,h));;
  gorbs:=Orbits(Group(List(gaperms,PermList)),gx);;
  morbs:=Orbits(Group(List(maperms,PermList)),mx);;
  printDescriptor63:=function(tag,ct,up,f,perms,orbs)
    local o,j,stab,mults,r,ext;
    for o in orbs do
      j:=o[1];
      stab:=Filtered([1..Length(qe)],q->perms[q][j]=j);
      ext:="cyclic-inertia";
      if Length(stab)=4 then
        mults:=List(Irr(up),chi->ScalarProduct(Irr(ct)[j],ClassFunction(ct,ValuesOfClassFunction(chi){f})));
        mults:=Filtered(mults,x->x>0);
        Assert(0,mults=[1,1,1,1] or mults=[2]);
        ext:=mults;
      fi;
      Print(tag,"|orbit=",o,"|degrees=",List(o,i->Irr(ct)[i][1]),"|stabilizer=",stab,"|full_inertia_multiplicities=",ext,"\n");
    od;
  end;;
  Print("LOCAL|prime=",prime,"|M=",Size(m0),"|H=",Size(hh),"|Mstructure=",StructureDescription(m0),"|quotient=",StructureDescription(qo),"\n");
  printDescriptor63("GLOBAL_PATTERN",tg,ty,fg,gaperms,gorbs);
  printDescriptor63("LOCAL_PATTERN",tm,th,fm,maperms,morbs);
  Descriptor63:=function(ct,up,f,perms,o)
    local j,stab,mults;
    j:=o[1];stab:=Filtered([1..Length(qe)],q->perms[q][j]=j);
    if Length(stab)<4 then return [stab,0];fi;
    mults:=List(Irr(up),chi->ScalarProduct(Irr(ct)[j],ClassFunction(ct,ValuesOfClassFunction(chi){f})));
    return [stab,Maximum(mults)];
  end;;
  omega:=[];;used:=[];;
  for go in gorbs do
    descr:=Descriptor63(tg,ty,fg,gaperms,go);;
    pos:=First([1..Length(morbs)],i->not i in used and Descriptor63(tm,th,fm,maperms,morbs[i])=descr);;
    Assert(0,pos<>fail);Add(used,pos);lo:=morbs[pos];;
    for q in [1..Length(qe)] do
      x:=gaperms[q][go[1]];;y:=maperms[q][lo[1]];;
      if IsBound(omega[x]) then Assert(0,omega[x]=y);else omega[x]:=y;fi;
    od;
  od;
  Assert(0,Set(omega{gx})=Set(mx));
  Assert(0,ForAll([1..Length(qe)],q->ForAll(gx,x->omega[gaperms[q][x]]=maperms[q][omega[x]])));
  Print("EXPLICIT_BIJECTION|",List(gx,x->[x,omega[x],ig[x][1],im[omega[x]][1]]),"\n");
  Print("Q_ELEMENTS|",qe,"\nGLOBAL_ACTION|",List(gaperms,r->r{gx}),"\nLOCAL_ACTION|",List(maperms,r->r{mx}),"\n");
od;
Print("FRESH_L34_ALL_ODD_PROJECTIVE_COMPLETE\n");
QUIT;
