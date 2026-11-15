LoadPackage("ctbllib");;
SetInfoLevel(InfoWarning,0);;

# A table row permutation from an ACTUAL index-two fusion.  Every class
# fibre has size at most two, so the full outer action is determined.
OuterRows63 := function(lo,hi)
  local f,n,cperm,bucket,j,vals,rperm;
  f:=FusionConjugacyClasses(lo,hi);
  n:=NrConjugacyClasses(lo); cperm:=[];
  for j in [1..n] do
    bucket:=Positions(f,f[j]);
    Assert(0,Length(bucket)<=2);
    if Length(bucket)=1 then cperm[j]:=j;
    else cperm[j]:=First(bucket,x->x<>j); fi;
  od;
  vals:=List(Irr(lo),ValuesOfClassFunction);
  rperm:=List(vals,x->Position(vals,x{cperm}));
  Assert(0,not fail in rperm);
  Assert(0,ForAll([1..n],j->rperm[rperm[j]]=j));
  return rperm;
end;;

ff:=GF(5);;
qi:=ImmutableMatrix(ff,[[0,1],[-1,0]]*One(ff));;
qj:=ImmutableMatrix(ff,[[2,0],[0,-2]]*One(ff));;
qu:=Group(qi,qj);;
Assert(0,Size(qu)=8 and StructureDescription(qu)="Q8");
im:=Normalizer(GL(2,5),qu);;
Assert(0,Size(im)=96);
im0s:=Filtered(NormalSubgroups(im),h->Size(h)=48);;
Assert(0,Length(im0s)=1);;
im0:=im0s[1];;
vv:=Cartesian(Elements(ff),Elements(ff));;
pts:=Cartesian([0,1,2],vv);;
tr1:=PermList(List(pts,x->Position(pts,[x[1],x[2]+[1,0]*One(ff)])));;
tr2:=PermList(List(pts,x->Position(pts,[x[1],x[2]+[0,1]*One(ff)])));;
zc:=PermList(List(pts,x->Position(pts,[(x[1]+1) mod 3,x[2]])));;
EmbedI63:=function(m)
  local sign;
  if m in im0 then sign:=1; else sign:=-1; fi;
  return PermList(List(pts,x->Position(pts,[(sign*x[1]) mod 3,x[2]*m])));
end;;
ii:=Group(List(GeneratorsOfGroup(im),EmbedI63));;
dd:=Group(tr1,tr2);;
nn:=Group(Concatenation([tr1,tr2,zc],GeneratorsOfGroup(ii)));;
sig:=PermList(List(pts,x->Position(pts,[(-x[1]) mod 3,x[2]])));;
yy:=ClosureGroup(nn,sig);;
Assert(0,Size(nn)=7200 and Size(yy)=14400 and IsNormal(yy,nn));
Assert(0,Centralizer(nn,dd)=ClosureGroup(dd,zc));
Print("FIVE|I=",StructureDescription(im),"|Iab=",AbelianInvariants(im),
 "|I0=",StructureDescription(im0),"|N/D=",[Size(nn/dd),AbelianInvariants(nn/dd)],
 "|sigma_order=",Order(sig),"|sigma_pointfix_DI=",ForAll(GeneratorsOfGroup(ClosureGroup(dd,ii)),g->g^sig=g),"\n");
ct:=CharacterTable(nn);; yt:=CharacterTable(yy);;
rp:=OuterRows63(ct,yt);;
bl:=PrimeBlocks(ct,5);;
for b in Set(bl.block) do
  ids:=Positions(bl.block,b);;
  orb:=Orbits(Group(PermList(rp)),ids);;
  Print("FIVE_BLOCK|id=",b,"|def=",bl.defect[b],"|k=",Length(ids),
   "|outer_orbits=",Collected(List(orb,Length)),"|rows=",ids,"|row_images=",rp{ids},"\n");
od;

# The central-three nonabelian row: an actual Heisenberg group,
# with the order-four symplectic automorphism, rather than a Q8 model.
hp:=Cartesian([0,1,2],[0,1,2],[0,1,2]);;
HProd63:=function(x,y)
 return [(x[1]+y[1]) mod 3,(x[2]+y[2]) mod 3,
 (x[3]+y[3]+2*(x[1]*y[2]-x[2]*y[1])) mod 3];
end;;
ht:=List([[1,0,0],[0,1,0],[0,0,1]],y->PermList(List(hp,x->Position(hp,HProd63(x,y)))));;
hh:=Group(ht);;
Assert(0,Exponent(hh)=3 and Size(DerivedSubgroup(hh))=3);
zz:=Centre(hh);;
u4:=PermList(List(hp,x->Position(hp,[(-x[2]) mod 3,x[1],x[3]])));;
Assert(0,Order(u4)=4 and ForAll(GeneratorsOfGroup(zz),z->z^u4=z));
ext:=Group(Concatenation(ht,[u4]));;
Assert(0,Size(ext)=108);
degrees:=List(Irr(ext),x->x[1]);;
Print("EXTRASPECIAL_C4|D=",[Size(hh),Size(Centre(hh))],"|Dexp=",Exponent(hh),
 "|center=",Size(zz),"|degrees=",Collected(degrees),
 "|k0=",Number(degrees,x->x mod 3<>0),"|k=",Length(degrees),"\n");
Assert(0,Collected(degrees)=[[1,4],[3,8],[4,2]]);

# Small defining-characteristic witness for the actual maximal-three
# centralizer used inside 2E6(2): Malle's O7(3) embedding contains N(P).
oo:=Omega(7,3);;
pp:=SylowSubgroup(oo,3);;
cc:=Centralizer(oo,pp);;
Print("O7_THREE|order=",Size(oo),"|P=",Size(pp),"|C=",Size(cc),
 "|ZP=",Size(Centre(pp)),"|actual_C_equals_ZP=",cc=Centre(pp),"\n");
Assert(0,cc=Centre(pp));
Print("FRESH_LOCAL_MATRIX_CARRIERS_COMPLETE\n");
QUIT;
