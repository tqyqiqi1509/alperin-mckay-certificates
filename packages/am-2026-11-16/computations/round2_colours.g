LoadPackage("ctbllib");;
SetInfoLevel(InfoWarning,0);;
ColourBlock60:=function(nm,p,ids)
 local t,v,bb,z,cp,b,rows,h0,cols,c,subset;
 t:=CharacterTable(nm);; v:=List(Irr(t),ValuesOfClassFunction);;
 bb:=PrimeBlocksOp(t,p);; z:=ClassPositionsOfCentre(t);;
 cp:=List(v,r->List(z,j->r[j]/r[1]));;
 for b in ids do
  rows:=Positions(bb.block,b);; h0:=Filtered(rows,j->bb.height[j]=0);;
  cols:=Set(cp{rows});;
  Print("PER_COLOUR_BLOCK|",nm,"|p=",p,"|block=",b,
        "|defect=",bb.defect[b],"|k=",Length(rows),"|k0=",Length(h0),"\n");
  for c in cols do
   subset:=Filtered(rows,j->cp[j]=c);;
   Print("COLOUR|",c,"|heights=",Collected(bb.height{subset}),"|rows=",subset,"\n");
  od;
  if nm="6.2E6(2)" and p=3 and b=17 then
   Assert(0,Length(cols)=3 and Length(rows)=27 and Length(h0)=27);
   Assert(0,ForAll(cols,c->Number(rows,j->cp[j]=c)=9));
  fi;
  if nm="6.2E6(2)" and p=3 and b=18 then
   Assert(0,Length(cols)=3 and Length(rows)=14 and Length(h0)=6);
   Assert(0,Length(Set(cp{h0}))=1);
  fi;
  if p=5 and nm="6.2E6(2)" then
   Assert(0,Length(cols)=1 and Length(rows)=20 and Length(h0)=20 and bb.defect[b]=2);
  fi;
  if p=7 and nm="6.2E6(2)" then
   Assert(0,Length(cols)=1 and Length(rows)=27 and Length(h0)=27 and bb.defect[b]=2);
  fi;
 od;
end;;
ColourBlock60("6.U6(2)",3,[7]);;
ColourBlock60("6.2E6(2)",3,[17,18]);;
ColourBlock60("6.2E6(2)",5,[257,258]);;
ColourBlock60("6.2E6(2)",7,[231,232]);;
Print("C60_PER_COLOUR_ANCHORS_COMPLETE\n");
QUIT;
