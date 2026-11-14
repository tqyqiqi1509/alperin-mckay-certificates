RequirePackage("chevie");;
w:=CoxeterGroup("E",8);;
ss:=Series(w,4);;
for s in ss do
 Print("SERIES|",ReflectionName(s.levi),"|",s.cuspidal,"|",CharNames(UnipotentCharacters(s.levi))[s.cuspidal],"|",s.charNumbers,"\n");
od;
quit;
