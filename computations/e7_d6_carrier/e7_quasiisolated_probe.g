RequirePackage("chevie");
W:=CoxeterGroup("E",7);
ss:=QuasiIsolatedRepresentatives(W);
for s in ss do
  C:=Centralizer(W,s);
  o:=Lcm(List(s.v,Denominator));
  Print("SS|",s,"|ORDER|",o,"|C|",C,"\n");
  if o=4 and C.group.semisimpleRank=6 then
    Print("TARGET|",s,"|COMPONENT|",C.phis,"\n");
    for a in C.phis do
      Print("COMPONENT_SPETS|",Spets(C.group,a),"\n");
    od;
    Print("IDENTITY_SPETS|",Spets(C.group),"\n");
    w0:=LongestCoxeterElement(W);
    Print("INVERSE_POWER|",s^w0=s^-1,"\n");
    Print("Q3_BASE_SPETS|",Spets(C.group,w0),"\n");
    for a in C.phis do
      Print("Q3_COMPONENT_SPETS|",Spets(C.group,a*w0),"\n");
    od;
  fi;
od;
