RequirePackage("chevie");
# Formal root data and all relevant cyclotomic orders; no q is chosen.
for shape52 in [["A",5,"A",2],["D",6,"A",1]] do
  W52:=ApplyFunc(CoxeterGroup,shape52);
  if shape52[1]="A" then signs52:=[1,-1]; else signs52:=[1];fi;
  for sign52 in signs52 do
    WF52:=Spets(W52,sign52*IdentityMat(7));
    Print("CENTRALIZER|",shape52,"|SIGN|",sign52,"\n");
    for order52 in [1,2,3,4,5,6,7,8,9,10,12,14,18] do
      sources52:=CuspidalPairs(WF52,order52);
      Print("D|",order52,"|NUMBER_OF_RAW_MARKED_PAIRS|",Length(sources52),"\n");
      for source52 in sources52 do
        names52:=CharNames(UnipotentCharacters(source52[1]));
        Print("ACTUAL_SOURCE|",source52[1],"|ORDINARY_MARK|",names52[source52[2]],"\n");
      od;
    od;
  od;
od;
