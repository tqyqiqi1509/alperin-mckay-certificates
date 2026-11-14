RequirePackage("chevie");
for spec75 in [["A",8],["E",6,"A",2]] do
  W75:=ApplyFunc(CoxeterGroup,spec75);
  for eps75 in [1,-1] do
    C75:=Spets(W75,eps75*IdentityMat(8));
    for e75 in [1,2,4] do
      Print("SOURCE75|",spec75,"|",eps75,"|",e75,"\n");
      for pair75 in CuspidalPairs(C75,e75) do
        L75:=pair75[1]; U75:=UnipotentCharacters(L75);
        Print("MARK75|",L75,"|",CharNames(U75)[pair75[2]],"\n");
      od;
    od;
  od;
od;
WD475:=CoxeterGroup("D",4);
UD475:=UnipotentCharacters(WD475);
namesD475:=CharNames(UD475);
autD475:=Spets(WD475,(1,2,4)).phi;
permD475:=PermutationOnUnipotents(WD475,autD475);
autD275:=Spets(WD475,(1,2)).phi;
permD275:=PermutationOnUnipotents(WD475,autD275);
for markD475 in ["D4","1.21"] do
  posD475:=Position(namesD475,markD475);
  Print("D4MARK75|",markD475,"|",posD475,"|triality_image|",
        namesD475[posD475^permD475],"|graph2_image|",namesD475[posD475^permD275],"\n");
od;
