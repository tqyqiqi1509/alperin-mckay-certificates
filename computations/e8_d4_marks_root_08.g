RequirePackage("chevie");
WD4R:=CoxeterGroup("D",4);
UD4R:=UnipotentCharacters(WD4R);
namesD4R:=CharNames(UD4R);
pairsD4R:=Filtered(CuspidalPairs(Spets(WD4R),4),x->SemisimpleRank(x[1])=4);
idsD4R:=List(pairsD4R,x->x[2]);
Print("CUSPIDAL4|",idsD4R,"|",namesD4R{idsD4R},"\n");
for autD4R in [(1,2),(1,2,4)] do
 phiD4R:=Spets(WD4R,autD4R).phi;
 puD4R:=PermutationOnUnipotents(WD4R,phiD4R);
 Print("GRAPH|",autD4R,"|ACTION|",puD4R,"|MARK_IMAGES|",List(idsD4R,i->namesD4R[i^puD4R]),"\n");
od;
