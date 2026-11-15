RequirePackage("chevie");;
ActualStrict21:=function()
 local w,ss,s,m,res,qres,a,c,cf,cp,e,pp,pair,k,kr,f,fm,kb,kfm,jr,j,jf,jgo,kpoly,
       ser,rel,ad,cent,centpoly,rc,idx,rs,nm,georel,kdegree,morder;
 SizeScreen([5000,5000]);
 w:=CoxeterGroup("E",8);
 ss:=QuasiIsolatedRepresentatives(w);
 Print("AM21_ACTUAL_ROOT_DATA|",w.simpleRoots,"|",w.simpleCoroots,"|",w.roots,"\n");
 Print("ISOLATED_COUNT|",Length(ss),"\n");
 for idx in [1..Length(ss)] do
  s:=ss[idx]; m:=Lcm(List(s.v,Denominator));
  if m>1 then
   c:=Centralizer(w,s).group;
   Print("ISOLATED|",idx,"|",m,"|",s.v,"|",ReflectionName(c),"|",InclusionGens(c),"\n");
   res:=PrimeResidues(m);
   for qres in res do
    a:=ToAlcove(s^qres);
    if a[2]<>s then Error("Power moves isolated vertex");fi;
    cf:=Spets(c,a[1]);
    Print("FORM_BEGIN|",idx,"|",qres,"|",ReflectionName(cf),"|",MatXPerm(w,cf.phi),"\n");
    cp:=CycPol(GenericOrder(cf,X(Cyclotomics)));
    Print("FORM_ORDER|",cp,"\n");
    for e in [1..20] do
     if Valuation(cp,e)>0 then
      pp:=CuspidalPairs(cf,e);
      Print("ACTUAL_E|",e,"|",Length(pp),"\n");
      for pair in pp do
       ad:=Number(ReflectionDegrees(pair[1]),x->x[1]=1 and x[2]=E(e));
       if ad>=2 then
        k:=Group(pair[1]); kr:=InclusionGens(k);
        f:=pair[1].phi; fm:=MatXPerm(w,f);
        kb:=w.roots{kr};
        if Length(kb)=0 then kpoly:=CycPol(1);
        else
         kfm:=List(kb,x->SolutionMat(kb,x*fm));
         kpoly:=CycPol(DeterminantMat(X(Cyclotomics)*IdentityMat(Length(kb))-kfm));
        fi;
        jr:=Filtered(c.rootInclusion{[1..c.N]},i->ForAll(kb,b->w.roots[i]*w.cartan*b=0));
        j:=ReflectionSubgroup(w,jr); jf:=Spets(j,f);
        jgo:=CycPol(GenericOrder(jf,X(Cyclotomics)))/kpoly;
        cent:=GetRelativeAction(w,k,f);
        centpoly:=CycPol(DeterminantMat(X(Cyclotomics)*IdentityMat(Length(cent))-cent));
        ser:=SeriesNC(cf,pair[1],pair[2],e);
        rel:=RelativeGroup(ser);
        nm:=CharNames(UnipotentCharacters(pair[1]))[pair[2]];
        georel:=Size(Normalizer(c,k))/Size(k);
        kdegree:=CycPolUnipotentDegrees(pair[1])[pair[2]];
        morder:=CycPol(GenericOrder(pair[1],X(Cyclotomics)));
        Print("ACTUAL_SOURCE|",e,"|",ad,"|",ReflectionName(pair[1]),"|",nm,"|",Size(rel),"|",kr,"|",fm,"|",centpoly,"|",ReflectionName(jf),"|",jgo,"|",SeriesOps.CharNumbers(ser),"|",georel,"|",kdegree,"|",morder,"\n");
       fi;
      od;
     fi;
    od;
    Print("FORM_END|",idx,"|",qres,"\n");
   od;
  fi;
 od;
 Print("DONE_ACTUAL_STRICT21\n");
end;;
ActualStrict21();
quit;
