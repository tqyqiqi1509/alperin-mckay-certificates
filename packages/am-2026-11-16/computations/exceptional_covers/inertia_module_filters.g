LoadPackage("ctbllib");;
SetInfoLevel(InfoWarning,0);;
g3:=GL(2,3);;
v3:=Filtered(Cartesian(Elements(GF(3)),Elements(GF(3))),v->v<>[Zero(GF(3)),Zero(GF(3))]);;
classes3:=Filtered(ConjugacyClassesSubgroups(g3),c->Size(Representative(c)) mod 3<>0);;
for cl in classes3 do
  ee:=Representative(cl);;
  orb:=Orbits(ee,v3,OnRight);;
  delta:=Sum(orb,o->Size(Stabilizer(ee,o[1],OnRight)));;
  invlin:=Filtered(Cartesian(Elements(GF(3)),Elements(GF(3))),w->ForAll(GeneratorsOfGroup(ee),m->w*m=w));;
  Print("ORDER9_INERTIA|E=",StructureDescription(ee),"|order=",Size(ee),
  "|orbits=",Collected(List(orb,Length)),"|k_minus_l=",delta,
  "|fixed_vectors=",Length(invlin),"|det_orders=",Set(Elements(ee),x->Order(DeterminantMat(x))),"\n");
od;
Print("FRESH_INERTIA_MODULE_FILTERS_COMPLETE\n");
QUIT;
