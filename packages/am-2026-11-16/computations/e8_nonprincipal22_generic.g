RequirePackage("chevie");;
Run22:=function() local w22,uc22,deg22,names22,h22,k22,d22,e22,pp22,pair22,ss22,rr22,nn22,wD22,ppD22,idsD22,a22,p22;
SizeScreen([1000,1000]);
w22:=CoxeterGroup("E",8);
uc22:=UnipotentCharacters(w22);
deg22:=CycPolUnipotentDegrees(w22);
names22:=CharNames(uc22);
Print("CHEVIE_COMMIT|64365f30757f5374c86511b377d7d98eb15c1e02\n");
for h22 in uc22.harishChandra do
 Print("HC|",h22.cuspidalName,"|",h22.charNumbers,"\n");
od;
Print("ENNOLA|",List([1..Length(deg22)],i->i^Ennola(w22)),"\n");
for k22 in [1..Length(deg22)] do
 d22:=deg22[k22];
 Print("DEG|",k22,"|",names22[k22],"|",d22.coeff,"|",List([1..30],j->Valuation(d22,j)),"\n");
od;
for e22 in [1,2,3,4,5,6,7,8,9,10,12,14,15,18,20,24,30] do
 pp22:=CuspidalPairs(Spets(w22),e22);
 Print("E_BEGIN|",e22,"|",Length(pp22),"\n");
 for pair22 in pp22 do
  if SemisimpleRank(pair22[1])<8 and SemisimpleRank(pair22[1])>0 then
   ss22:=SeriesNC(Spets(w22),pair22[1],pair22[2],e22);
   rr22:=RelativeGroup(ss22);
   nn22:=SeriesOps.CharNumbers(ss22);
   Print("SER|",e22,"|",ReflectionName(pair22[1]),"|",CharNames(UnipotentCharacters(pair22[1]))[pair22[2]],"|",Size(rr22),"|",nn22,"\n");
  fi;
 od;
 Print("E_END|",e22,"\n");
od;
wD22:=CoxeterGroup("D",4);
ppD22:=Filtered(CuspidalPairs(Spets(wD22),4),x->SemisimpleRank(x[1])=4);
idsD22:=List(ppD22,x->x[2]);
Print("D4_MARKS|",idsD22,"|",CharNames(UnipotentCharacters(wD22)){idsD22},"\n");
for a22 in [(1,2),(1,2,4)] do
 p22:=PermutationOnUnipotents(wD22,Spets(wD22,a22).phi);
 Print("D4_GRAPH|",a22,"|",p22,"|",List(idsD22,i->i^p22),"\n");
od;
Print("DONE22\n");
end;;
Run22();
quit;
