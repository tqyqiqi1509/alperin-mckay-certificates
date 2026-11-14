LoadPackage("ctbllib");;
SetInfoLevel(InfoWarning,0);;
f5:=GF(5);; f3:=GF(3);;
a5:=ImmutableMatrix(f5,[[0,1],[-1,0]]*One(f5));;
b5:=ImmutableMatrix(f5,[[2,0],[0,-2]]*One(f5));;
q8:=Group(a5,b5);; ei:=Normalizer(GL(2,5),q8);;
o2:=PCore(ei,2);; Assert(0,Size(o2)=16);;
quo:=NaturalHomomorphismByNormalSubgroup(ei,o2);;
t3:=ImmutableMatrix(f3,[[1,0],[1,1]]*One(f3));;
r3:=ImmutableMatrix(f3,[[1,0],[0,-1]]*One(f3));;
s3:=Group(t3,r3);; Assert(0,Size(s3)=6);;
iso:=IsomorphismGroups(Image(quo),s3);; Assert(0,iso<>fail);;
action:=CompositionMapping(iso,quo);;
vec5:=Cartesian(Elements(f5),Elements(f5));;
vec3:=Cartesian(Elements(f3),Elements(f3));;
points:=Cartesian(vec3,vec5);;
Perm6:=function(c,m)
 return PermList(List(points,x->Position(points,[x[1]*c,x[2]*m])));
end;;
iip:=Group(List(GeneratorsOfGroup(ei),m->Perm6(Image(action,m),m)));;
Trans6:=function(a,b)
 return PermList(List(points,x->Position(points,[x[1]+a,x[2]+b])));
end;;
z3:=Trans6([1,0]*One(f3),[0,0]*One(f5));;
w3:=Trans6([0,1]*One(f3),[0,0]*One(f5));;
d1:=Trans6([0,0]*One(f3),[1,0]*One(f5));;
d2:=Trans6([0,0]*One(f3),[0,1]*One(f5));;
p5:=Group(d1,d2);;
ng3:=Group(Concatenation([z3,w3,d1,d2],GeneratorsOfGroup(iip)));;
Assert(0,Size(ng3)=21600 and Size(Centre(ng3))=3 and z3 in DerivedSubgroup(ng3));
Assert(0,Centralizer(ng3,p5)=Group(z3,w3,d1,d2));
Assert(0,Size(ng3/p5)=864); # Identifier not retested: smallgrp absent.
Assert(0,AbelianInvariants(ng3/p5)=[4]);
p2i:=SylowSubgroup(ei,2);;
Assert(0,Size(p2i)=32 and IsomorphismGroups(p2i,WreathProduct(CyclicGroup(4),SymmetricGroup(2)))<>fail);
sig3:=Perm6(-IdentityMat(2,f3),IdentityMat(2,f5));;
yg3:=ClosureGroup(ng3,sig3);;
Assert(0,Size(yg3)=43200 and z3^sig3=z3^-1);
Assert(0,ForAll(GeneratorsOfGroup(ClosureGroup(p5,iip)),g->g^sig3=g));
ct3:=CharacterTable(ng3);; yt3:=CharacterTable(yg3);;
fu:=FusionConjugacyClasses(ct3,yt3);;
cp:=List([1..Length(fu)],j->j);;
for j in [1..Length(fu)] do
 ids:=Positions(fu,fu[j]);;
 Assert(0,Length(ids)<=2);
 if Length(ids)=2 then cp[j]:=First(ids,x->x<>j); fi;
od;
vals:=List(Irr(ct3),ValuesOfClassFunction);;
rp:=List(vals,x->Position(vals,x{cp}));;
Assert(0,not fail in rp);
zc:=PositionProperty(ConjugacyClasses(ct3),c->z3 in c);;
Assert(0,zc<>fail);
bl:=PrimeBlocks(ct3,5);;
reg5:=Filtered([1..Length(vals)],j->OrdersClassRepresentatives(ct3)[j] mod 5<>0);;
Print("SIX_FIVE|Norder=",Size(ng3),"|NmodD=",[Size(ng3/p5),AbelianInvariants(ng3/p5)],"|NmodDab=",AbelianInvariants(ng3/p5),"|Ialpha=",StructureDescription(p2i),"|Zsize=",Size(Centre(ng3)),
 "|C_D_order=",Size(Centralizer(ng3,p5)),"|stem=",z3 in DerivedSubgroup(ng3),
 "|O2ab=",AbelianInvariants(Kernel(action)),"|outer_inverts_Z=",z3^sig3=z3^-1,"\n");
for b in Set(bl.block) do
 ids:=Positions(bl.block,b);;
 cols:=Set(ids,j->Irr(ct3)[j][zc]/Irr(ct3)[j][1]);;
 imp:=Set(rp{ids},j->bl.block[j]);;
 Print("SIX_FIVE_BLOCK|id=",b,"|def=",bl.defect[b],"|k=",Length(ids),
 "|l=",RankMat(List(ids,j->vals[j]{reg5})),"|colour=",cols,"|sigma_image_block=",imp,"|degree_multiset=",Collected(List(ids,j->Irr(ct3)[j][1])),"\n");
od;
Print("FRESH_SIXFOLD_FIVE_NORMALIZER_COMPLETE\n");
QUIT;
