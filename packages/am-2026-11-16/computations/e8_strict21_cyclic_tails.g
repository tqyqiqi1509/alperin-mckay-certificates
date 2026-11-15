RequirePackage("chevie");;
Tails21:=function()
 local w,ss,s,m,idx,c,cf,qres,a,e,cp,pp,pair,ad,k,f,cen,poly,degree,order,rel;
 SizeScreen([3000,3000]);w:=CoxeterGroup("E",8);ss:=QuasiIsolatedRepresentatives(w);
 for idx in [1..Length(ss)] do
  s:=ss[idx];m:=Lcm(List(s.v,Denominator));
  if m>1 then
   c:=Centralizer(w,s).group;
   for qres in PrimeResidues(m) do
    a:=ToAlcove(s^qres);cf:=Spets(c,a[1]);cp:=CycPol(GenericOrder(cf,X(Cyclotomics)));
    for e in [1..20] do
     if Valuation(cp,e)>0 then
      pp:=CuspidalPairs(cf,e);
      for pair in pp do
       ad:=Number(ReflectionDegrees(pair[1]),x->x[1]=1 and x[2]=E(e));
       if ad<=1 then
        degree:=CycPolUnipotentDegrees(pair[1])[pair[2]];
        order:=CycPol(GenericOrder(pair[1],X(Cyclotomics)));
        if ad=0 then
         Print("ZERO_TAIL|",idx,"|",m,"|",qres,"|",e,"|",CharNames(UnipotentCharacters(pair[1]))[pair[2]],"|",degree,"|",order,"\n");
        else
         k:=Group(pair[1]);f:=pair[1].phi;cen:=GetRelativeAction(w,k,f);
         poly:=CycPol(DeterminantMat(X(Cyclotomics)*IdentityMat(Length(cen))-cen));
         rel:=RelativeGroup(SeriesNC(cf,pair[1],pair[2],e));
         Print("CYCLIC_TAIL|",idx,"|",m,"|",qres,"|",e,"|",ReflectionName(pair[1]),"|",CharNames(UnipotentCharacters(pair[1]))[pair[2]],"|",Size(rel),"|",InclusionGens(k),"|",MatXPerm(w,f),"|",poly,"|",degree,"|",order,"\n");
        fi;
       fi;
      od;
     fi;
    od;
   od;
  fi;
 od;
 Print("DONE_TAILS21\n");
end;;
Tails21();
quit;
