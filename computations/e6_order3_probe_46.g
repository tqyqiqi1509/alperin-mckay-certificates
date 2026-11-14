RequirePackage("chevie");
W:=CoxeterGroup("E",6);
Print("HIGHEST|",W.roots[36],"\n");
R:=ReflectionSubgroup(W,[1,3,5,6,2,36]);
Print("SUB|",R,"\n");
f:=LongestCoxeterElement(W);
Print("F|",f,"\n");
RF:=Spets(R,f);
Print("SPETS|",RF,"\n");
for d in [1,2,3,4,5,6,8,9,10,12,18] do
  ps:=CuspidalPairs(RF,d);
  Print("D|",d,"|COUNT|",Length(ps),"\n");
  for pair in ps do Print("SOURCE|",d,"|",pair[1],"|",CharNames(UnipotentCharacters(pair[1]))[pair[2]],"\n");od;
od;
