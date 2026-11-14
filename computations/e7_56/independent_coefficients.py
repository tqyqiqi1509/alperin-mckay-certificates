"""Build the ordinary coefficient check from node56's fresh natural carrier."""
from pathlib import Path
import contextlib
import io
with contextlib.redirect_stdout(io.StringIO()):
    import independent_carrier as c
p=Path(__file__).resolve().parent
gens=[tuple(abs(j) for j in g) for g in c.gens]
head="A:=Group("+",".join("PermList("+str(list(g))+")" for g in gens)+");;\n"
body=r'''
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
'''
(p/'independent_coefficients.g').write_text(head+body)
print(p/'independent_coefficients.g')
