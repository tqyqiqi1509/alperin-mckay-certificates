RequirePackage("chevie");;
d:=false;; w:=false;;
# Fresh direct use of the CHEVIE generic degree objects, without a printed
# cyclotomic factor parser or previously stored valuation certificates.
for rank in [6,8] do
  w:=CoxeterGroup("E",rank);
  uc:=UnipotentCharacters(w);
  deg:=CycPolUnipotentDegrees(w);
  names:=CharNames(uc);
  Print("TYPE|E",rank,"\n");
  for hc in uc.harishChandra do
    Print("HC|",hc.cuspidalName,"|",hc.charNumbers,"\n");
  od;
  if rank=8 then Print("ENNOLA|",List([1..Length(deg)],i->i^Ennola(w)),"\n"); fi;
  for k in [1..Length(deg)] do
    d:=deg[k];
    Print("DEG|",k,"|",names[k],"|",d.coeff,"|",
      List([1..30],j->Valuation(d,j)),"\n");
  od;
od;
quit;
