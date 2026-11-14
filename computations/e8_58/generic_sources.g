RequirePackage("chevie");
for shape58 in [["D",8],["E",7,"A",1],["D",5,"A",3],["A",7,"A",1]] do
  W58:=ApplyFunc(CoxeterGroup,shape58);
  Print("CENTRALIZER|",shape58,"\n");
  for d58 in [1,2,3,4,5,6,7,8,9,10,12,14,18] do
    cp58:=CuspidalPairs(Spets(W58),d58);
    Print("D|",d58,"|",Length(cp58),"\n");
    for pa58 in cp58 do
      Print("SOURCE|",pa58[1],"|",CharNames(UnipotentCharacters(pa58[1]))[pa58[2]],"\n");
    od;
  od;
od;
DD58:=CoxeterGroup("D",4);
UD58:=UnipotentCharacters(DD58);
Print("D4_NAMES|",CharNames(UD58),"\n");
Print("D4_DEGREES|",CycPolUnipotentDegrees(DD58),"\n");
