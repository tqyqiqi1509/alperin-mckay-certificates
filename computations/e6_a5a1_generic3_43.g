RequirePackage("chevie");
for eps in [1,-1] do
  R:=Spets(CoxeterGroup("A",5,"A",1),eps*IdentityMat(6));
  Print("A5A1_FORM|",eps,"\n");
  for d in [1,2] do
    ps:=CuspidalPairs(R,d);
    for pair in ps do
      Print("SOURCE|",d,"|",pair[1],"|",CharNames(UnipotentCharacters(pair[1]))[pair[2]],"\n");
    od;
  od;
od;
for n in [2,3,6] do
  W:=CoxeterGroup("A",n-1);
  Print("GENERIC_A|",n,"\n");
  uc:=UnipotentCharacters(W);
  ns:=CharNames(uc);
  ds:=CycPolUnipotentDegrees(W);
  for i in [1..Length(ns)] do Print("DEGREE|",ns[i],"|",ds[i],"\n");od;
od;
quit;
