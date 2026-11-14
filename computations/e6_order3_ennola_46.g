RequirePackage("chevie");
W:=CoxeterGroup("E",6);
R:=ReflectionSubgroup(W,[1,3,5,6,2,36]);
f:=LongestCoxeterElement(W);
RF:=Spets(R,-MatXPerm(W,f));
Print("ENNOLA_SPETS|",RF,"\n");
for d in [1,2,3,4,6] do
  ps:=CuspidalPairs(RF,d);
  Print("D|",d,"|COUNT|",Length(ps),"\n");
  for pair in ps do
    Print("SOURCE|",d,"|",pair[1],"|",CharNames(UnipotentCharacters(pair[1]))[pair[2]],"\n");
  od;
od;
