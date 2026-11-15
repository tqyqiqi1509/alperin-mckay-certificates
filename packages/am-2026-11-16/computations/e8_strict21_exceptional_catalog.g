RequirePackage("chevie");;
RunStrict21:=function()
 local n,w,wf,uc,deg,names,h,k,d,e,pp,pair,ss,rr,nn,go,ad,tw;
 SizeScreen([3000,3000]);
 Print("CHEVIE_COMMIT|64365f30757f5374c86511b377d7d98eb15c1e02\n");
 for n in [6,7] do
  w:=CoxeterGroup("E",n);
  for tw in [1,-1] do
   wf:=Spets(w,tw*IdentityMat(n));
   Print("GROUP_BEGIN|",ReflectionName(wf),"\n");
   uc:=UnipotentCharacters(wf);
   deg:=CycPolUnipotentDegrees(wf);
   names:=CharNames(uc);
   Print("COUNT|",Length(deg),"\n");
   for k in [1..Length(deg)] do
    d:=deg[k];
    Print("DEG|",k,"|",names[k],"|",d.coeff,"|",List([1..30],j->Valuation(d,j)),"\n");
   od;
   for e in [1..30] do
    pp:=CuspidalPairs(wf,e);
    Print("E_BEGIN|",e,"|",Length(pp),"\n");
    for pair in pp do
     ss:=SeriesNC(wf,pair[1],pair[2],e);
     rr:=RelativeGroup(ss);
     nn:=SeriesOps.CharNumbers(ss);
     ad:=Number(ReflectionDegrees(pair[1]),x->x[1]=1 and x[2]=E(e));
     go:=CycPol(GenericOrder(pair[1],X(Cyclotomics)));
     Print("SER|",e,"|",ReflectionName(pair[1]),"|",SemisimpleRank(pair[1]),"|",
      CharNames(UnipotentCharacters(pair[1]))[pair[2]],"|",Size(rr),"|",nn,"|",go,"|",ad,"\n");
    od;
    Print("E_END|",e,"\n");
   od;
   Print("GROUP_END|",ReflectionName(wf),"\n");
  od;
 od;
 Print("DONE_STRICT21\n");
end;;
RunStrict21();
quit;
