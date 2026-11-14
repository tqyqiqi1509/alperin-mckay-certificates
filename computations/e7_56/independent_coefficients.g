A:=Group(PermList([22, 27, 28, 46, 48, 49, 7, 44, 45, 4, 47, 5, 6, 54, 15, 16, 17, 55, 19, 20, 21, 34, 56, 24, 25, 26, 39, 43, 14, 18, 31, 32, 33, 1, 23, 36, 37, 38, 2, 40, 41, 42, 3, 51, 52, 10, 53, 12, 13, 50, 8, 9, 11, 29, 30, 35]),PermList([2, 3, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 23, 24, 25, 26, 14, 15, 16, 17, 27, 18, 19, 20, 21, 28, 22, 35, 29, 36, 37, 38, 39, 30, 40, 41, 42, 43, 31, 32, 33, 34, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 56, 54, 55]),PermList([56, 35, 23, 40, 41, 42, 50, 24, 25, 36, 26, 37, 38, 43, 51, 52, 53, 28, 44, 45, 47, 55, 3, 8, 9, 11, 30, 18, 39, 27, 46, 48, 49, 54, 2, 10, 12, 13, 29, 4, 5, 6, 14, 19, 20, 31, 21, 32, 33, 7, 15, 16, 17, 34, 22, 1]),PermList([56, 55, 54, 51, 52, 53, 50, 46, 48, 44, 49, 45, 47, 43, 40, 41, 42, 39, 36, 37, 38, 35, 34, 31, 32, 33, 30, 29, 28, 27, 24, 25, 26, 23, 22, 19, 20, 21, 18, 15, 16, 17, 14, 10, 12, 8, 13, 9, 11, 7, 4, 5, 6, 3, 2, 1]));;

Print("Actual quotient order=",Size(A)," structure=",StructureDescription(A),"\n");
if Size(A)<>36 then Error("Wrong actual quotient"); fi;
subclasses:=ConjugacyClassesSubgroups(A);;
count:=0;; pairs:=0;; failed:=[];; witnesses:=[];;
for sc in subclasses do
  S:=Representative(sc);;
  for R in NormalSubgroups(S) do
    if Index(S,R) in [1,2] then
      pairs:=pairs+1;;
      K:=Intersection(Normalizer(A,S),Normalizer(A,R));;
      ir:=Irr(R);; is:=Irr(S);;
      for j in [1..Length(ir)] do
        x:=ir[j];;
        stab:=Filtered(Elements(K),g->x^g=x);;
        cand:=Filtered([1..Length(is)],i->ScalarProduct(RestrictedClassFunction(is[i],R),x)>0);;
        wit:=Filtered(cand,i->ForAll(stab,g->is[i]^g=is[i]));;
        count:=count+1;;
        if Length(wit)=0 then
          Add(failed,[Size(S),Size(R),j]);;
        else
          Add(witnesses,[Size(S),Size(R),j,wit[1],Size(K),Length(stab)]);;
        fi;
      od;
    fi;
  od;
od;
Print("Subgroup classes=",Length(subclasses)," pairs=",pairs," ordinary entries=",count,"\n");
Print("Failed=",failed,"\n");
Print("Witnesses (|S|,|R|,Rlabel,Slabel,|K|,|K_eta|)=",witnesses,"\n");
if Length(failed)>0 then Error("Ordinary coefficient condition failed"); fi;
QUIT;
