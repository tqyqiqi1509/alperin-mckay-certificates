LoadPackage("ctbllib");;LoadPackage("perfgrp");;
SetInfoLevel(InfoWarning,0);;
Rows4C:=function(lo,hi)
 local fu,cp,j,inds,vals;
 fu:=FusionConjugacyClasses(lo,hi);Assert(0,fu<>fail);cp:=[];
 for j in [1..Length(fu)] do
  inds:=Positions(fu,fu[j]);Assert(0,Length(inds)<=2);
  if Length(inds)=1 then cp[j]:=j;else cp[j]:=First(inds,x->x<>j);fi;
 od;
 vals:=List(Irr(lo),ValuesOfClassFunction);
 return List(vals,v->Position(vals,v{cp}));
end;;
Run4C:=function(k)
 local g,z,zc,aut,inn,ac,s2,ao,a,y,eg,g0,p,m,h,tg,tm,ty,th,rg,rm,
  zg,zm,bg,bm,c,gi,mi,go,mo,omega,used,o,pos,r,x,matched;
 g:=PerfectGroup(IsPermGroup,80640,k);z:=Centre(g);
 Assert(0,IsCyclic(z) and Size(z)=4);zc:=First(Elements(z),x->Order(x)=4);
 aut:=AutomorphismGroup(g);inn:=InnerAutomorphismsAutomorphismGroup(aut);
 ac:=Stabilizer(aut,zc,function(x,a)return Image(a,x);end);
 Assert(0,Size(ac)=2*Size(inn));
 s2:=SylowSubgroup(ac,2);ao:=First(Elements(s2),a->Order(a)=2 and not a in inn);
 Assert(0,ao<>fail);a:=Group(ao);y:=SemidirectProduct(a,g);eg:=Embedding(y,2);
 g0:=Image(eg,g);zc:=Image(eg,zc);p:=SylowSubgroup(g0,3);
 m:=Normalizer(g0,p);h:=Normalizer(y,p);
 Assert(0,Size(m)=288 and Size(h)=576 and Centralizer(y,g0)=Centre(g0));
 Assert(0,Centralizer(g0,p)=ClosureGroup(p,Centre(g0)));
 tg:=CharacterTable(g0);tm:=CharacterTable(m);ty:=CharacterTable(y);th:=CharacterTable(h);
 matched:=Filtered(["4_1.L3(4)","4_2.L3(4)"],nm->TransformingPermutationsCharacterTables(tg,CharacterTable(nm))<>fail);
 Assert(0,Length(matched)=1);
 rg:=Rows4C(tg,ty);rm:=Rows4C(tm,th);
 zg:=PositionProperty(ConjugacyClasses(g0),cl->zc in cl);zm:=PositionProperty(ConjugacyClasses(m),cl->zc in cl);
 bg:=PrimeBlocks(tg,3);bm:=PrimeBlocks(tm,3);
 Print("FOUR_BASE|perf_id=",k,"|table=",matched[1],"|Aut=",Size(aut),"|C_Aut_Z=",Size(ac),"|N=",Size(m),"|H=",Size(h),"|NmodP=",StructureDescription(m/p),"\n");
 for c in [E(4),-E(4)] do
  gi:=Filtered([1..Length(Irr(tg))],j->Irr(tg)[j][1] mod 3<>0 and Irr(tg)[j][zg]=c*Irr(tg)[j][1]);
  mi:=Filtered([1..Length(Irr(tm))],j->Irr(tm)[j][1] mod 3<>0 and Irr(tm)[j][zm]=c*Irr(tm)[j][1]);
  Assert(0,Length(gi)=6 and Length(mi)=6 and Length(Set(gi,j->bg.block[j]))=1 and Length(Set(mi,j->bm.block[j]))=1);
  go:=Orbits(Group(PermList(rg)),gi);mo:=Orbits(Group(PermList(rm)),mi);
  omega:=[];used:=[];
  for o in go do
   pos:=First([1..Length(mo)],j->not j in used and Length(mo[j])=Length(o));Assert(0,pos<>fail);Add(used,pos);
   omega[o[1]]:=mo[pos][1];
   if Length(o)=2 then omega[rg[o[1]]]:=rm[mo[pos][1]];fi;
  od;
  Assert(0,Set(omega{gi})=Set(mi) and ForAll(gi,j->omega[rg[j]]=rm[omega[j]]));
  Print("FOUR_COLOUR|perf_id=",k,"|colour=",c,"|global_orbits=",Collected(List(go,Length)),"|local_orbits=",Collected(List(mo,Length)),
  "|bijection=",List(gi,j->[j,omega[j],Irr(tg)[j][1],Irr(tm)[omega[j]][1]]),"\n");
 od;
end;;
Run4C(3);
Run4C(4);
Print("FRESH_L34_FOURFOLD_ACTUAL_COMPLETE\n");
QUIT;
