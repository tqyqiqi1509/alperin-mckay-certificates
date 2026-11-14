RequirePackage("chevie");
for spec60 in [["A",8],["E",6,"A",2]] do
  grp60:=ApplyFunc(CoxeterGroup,spec60);
  Print("GROUP60|",spec60,"\n");
  for eps60 in [1,-1] do
    cos60:=Spets(grp60,eps60*IdentityMat(8));
    for d60 in [1,2,3,4,5,6,7,8,9,10,12,14,18] do
      pairs60:=CuspidalPairs(cos60,d60);
      Print("CATALOGUE60|",eps60,"|",d60,"|",Length(pairs60),"\n");
      for pair60 in pairs60 do
        lev60:=pair60[1];
        Print("PAIR60|",eps60,"|",d60,"|",lev60,"|",CharNames(UnipotentCharacters(lev60))[pair60[2]],"\n");
        if SemisimpleRank(lev60)<8 then
          Print("K60|",lev60.reflectionGroup.simpleRoots,"\n");
          Print("F60|",lev60.F0Mat,"\n");
        fi;
        Print("END60\n");
      od;
    od;
  od;
od;
quit;
