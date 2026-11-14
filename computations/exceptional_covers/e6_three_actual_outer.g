LoadPackage("ctbllib");;
SetInfoLevel(InfoWarning,0);;
yy:=SO(7,3);; gg:=DerivedSubgroup(yy);;
Print("ORTHOGONAL_INITIAL|SO=",Size(yy),"|derived=",Size(gg),"|perfect=",IsPerfectGroup(gg),"\n");
Assert(0,IsSubgroup(yy,gg) and Size(yy)=2*Size(gg));
s2:=SylowSubgroup(yy,2);;si:=First(Elements(s2),x->Order(x)=2 and not x in gg);;
Assert(0,si<>fail);;Print("OUTER_INVOLUTION|order=",Order(si),"|outside_derived=",not si in gg,"\n");
pp:=SylowSubgroup(gg,3);;
mm:=Normalizer(gg,pp);; hh:=Normalizer(yy,pp);;
Assert(0,IsNormal(hh,mm) and Size(hh)=2*Size(mm));
pd:=DerivedSubgroup(pp);;
qh:=NaturalHomomorphismByNormalSubgroup(hh,pd);;
qg:=Image(qh,mm);;qy:=Image(qh);;
Print("E6_THREE_ACTUAL|Omega=",Size(gg),"|SO=",Size(yy),"|P=",Size(pp),"|Pderived=",Size(pd),
 "|N=",Size(mm),"|H=",Size(hh),"|NmodPderived=",Size(qg),"|HmodPderived=",Size(qy),"|C=",Size(Centralizer(gg,pp)),"\n");
ct:=CharacterTable(qg);;yt:=CharacterTable(qy);;fu:=FusionConjugacyClasses(ct,yt);;
cp:=[];;
for j in [1..Length(fu)] do
 ids:=Positions(fu,fu[j]);;Assert(0,Length(ids)<=2);
 if Length(ids)=1 then cp[j]:=j;else cp[j]:=First(ids,x->x<>j);fi;
od;
vals:=List(Irr(ct),ValuesOfClassFunction);;
rp:=List(vals,v->Position(vals,v{cp}));;
Assert(0,not fail in rp and ForAll(vals,x->x[1] mod 3<>0));
Print("E6_THREE_PPRIME|k0=",Length(vals),"|degrees=",Collected(List(vals,x->x[1])),
 "|outer_orbits=",Collected(List(Orbits(Group(PermList(rp)),[1..Length(vals)]),Length)),"|images=",rp,"\n");
Assert(0,Length(vals)=18 and Number([1..Length(vals)],j->rp[j]=j)=12);
Read("computations/exceptional_covers/actual_c2_pairing.g");;
ActualC2Pair63("2.2E6(2)",16,ct,[1..Length(vals)],rp);
Print("SO_POINTWISE_SYLOW|C=",Size(Centralizer(yy,pp)),"|contained_inner=",IsSubgroup(gg,Centralizer(yy,pp)),"\n");
Assert(0,IsSubgroup(gg,Centralizer(yy,pp)));
Print("FRESH_E6_THREE_ACTUAL_OUTER_COMPLETE\n");
QUIT;
