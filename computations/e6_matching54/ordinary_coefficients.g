# Exact finite coefficient groups, independent of q and ell.
B:=Group((1,2),(3,4),(5,6));;
c:=(1,3,5)(2,4,6);; h:=(1,3)(2,4);;
S:=ClosureGroup(B,c);; K:=ClosureGroup(S,h);;
Print("A1_AMBIENT|",Size(B),"|",Size(S),"|",Size(K),"\n");
for a in Irr(B) do
 found:=[];
 for nu in Irr(K) do
  eta:=RestrictedClassFunction(nu,S);
  if ScalarProduct(eta,eta)=1 and ScalarProduct(RestrictedClassFunction(eta,B),a)>0 then
   Add(found,nu[1]);
  fi;
 od;
 Assert(0,Length(found)>0);
 Print("A1_ZERO_FIBRE|",ValuesOfClassFunction(a),"|EXTENDING_DEGREES|",Set(found),"\n");
od;
# The nonzero diagonal line has coefficient pair 1 < C3. Its trivial
# coefficient extends to the full actual normalizer of this pair.
S3diag:=Group(c);; N3:=Normalizer(K,S3diag);;
Assert(0,RestrictedClassFunction(TrivialCharacter(N3),S3diag)=TrivialCharacter(S3diag));
Print("A1_NONZERO_DIAGONAL|",Size(N3),"|TRIVIAL_EXTENSION|true\n");

CharStab:=function(U,R,a)
 local cls,vals,els;
 cls:=ConjugacyClasses(R);vals:=ValuesOfClassFunction(a);
 els:=Filtered(Elements(U),g->ForAll([1..Length(cls)],i->
   vals[i]=vals[PositionProperty(cls,z->Representative(cls[i])^g in z)]));
 return Group(els);
end;;
for r in [0..2] do
 R:=Group(List([1..r],i->(2*i-1,2*i)),());
 N:=Normalizer(K,R);
 for a in Irr(R) do
  U:=CharStab(N,R,a);
  good:=Filtered(Irr(U),nu->nu[1]=a[1] and RestrictedClassFunction(nu,R)=a);
  Assert(0,Length(good)>0);
 od;
 Print("A1_COORDINATE|",r,"|",Size(R),"|",Size(N),"|ALL_ORDINARY|true\n");
od;

# Marked A3: the actual rank-two signed reflection group and coordinate
# stabilizers; the central inversion alone is not a stabilizer of nonzero t.
B2:=Group((1,2),(3,4));; K2:=ClosureGroup(B2,(1,3)(2,4));;
for R in [Group(()),Group((1,2)),B2] do
 N:=Normalizer(K2,R);
 for a in Irr(R) do
  U:=CharStab(N,R,a);
  good:=Filtered(Irr(U),nu->nu[1]=a[1] and RestrictedClassFunction(nu,R)=a);
  Assert(0,Length(good)>0);
 od;
 Print("MARKED_A3|",Size(R),"|",Size(N),"|ALL_ORDINARY|true\n");
od;

# Faithful cubic d3: the zero fibre S3 < C3 x S3, and the nonzero
# eigenline 1 < C3. Central-C3-trivial tensor extensions work for all
# three ordinary S3 characters, including its degree-two character.
U3:=Group((1,2),(1,2,3));; Z3:=Group((4,5,6));; P3:=ClosureGroup(U3,GeneratorsOfGroup(Z3));;
for a in Irr(U3) do
 good:=Filtered(Irr(P3),nu->nu[1]=a[1] and RestrictedClassFunction(nu,U3)=a and
  RestrictedClassFunction(nu,Z3)=a[1]*TrivialCharacter(Z3));
 Assert(0,Length(good)=1);
 Print("FAITHFUL_D3_ZERO|",a[1],"|CANONICAL_EXTENSION|true\n");
od;
Print("DONE|all ordinary finite coefficient checks passed\n");
QUIT;
