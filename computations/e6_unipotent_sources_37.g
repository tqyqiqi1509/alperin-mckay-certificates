RequirePackage("chevie");
W:=CoxeterGroup("E",6);
for eps in [1,-1] do
  WF:=Spets(W,eps*IdentityMat(6));
  Print("FORM|",eps,"\n");
  for d in [1,2,3,4,5,6,8,9,10,12,18] do
    ps:=CuspidalPairs(WF,d);
    Print("D|",d,"|ALLCOUNT|",Length(ps),"\n");
    for pair in ps do
      if SemisimpleRank(pair[1])<6 then
        Print("SOURCE|",d,"|",pair[1],"|",CharNames(UnipotentCharacters(pair[1]))[pair[2]],"\n");
      fi;
    od;
  od;
od;
