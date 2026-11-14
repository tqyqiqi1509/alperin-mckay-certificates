LoadPackage("ctbllib");;
SetInfoLevel(InfoWarning,0);;
QuotientBlocks63:=function(uppername,lowername)
 local up,lo,fu,uv,lv,map,ub,lb,b,ids,upb;
 up:=CharacterTable(uppername);lo:=CharacterTable(lowername);fu:=GetFusionMap(up,lo);
 if fu=fail then Print("NO_STORED_QUOTIENT|",uppername,"->",lowername,"\n");return;fi;
 Assert(0,Size(up)=3*Size(lo));uv:=List(Irr(up),ValuesOfClassFunction);lv:=List(Irr(lo),ValuesOfClassFunction);
 map:=List(lv,v->Position(uv,v{fu}));Assert(0,not fail in map and Length(Set(map))=Length(map));
 ub:=PrimeBlocks(up,3);lb:=PrimeBlocks(lo,3);
 for b in Set(lb.block) do
  ids:=Positions(lb.block,b);upb:=Set(map{ids},j->ub.block[j]);Assert(0,Length(upb)=1);
  Assert(0,ub.defect[upb[1]]=1+lb.defect[b]);
  if lb.defect[b]>0 then Print("ACTUAL_CENTRAL_BLOCK|",uppername,"->",lowername,"|lower_block=",b,"|upper_block=",upb[1],"|lower_def=",lb.defect[b],"|upper_def=",ub.defect[upb[1]],"|inflated_rows=",map{ids},"\n");fi;
 od;
end;;
QuotientBlocks63("6.L3(4)","2.L3(4)");
QuotientBlocks63("12_1.L3(4)","4_1.L3(4)");
QuotientBlocks63("12_1.L3(4)","4_2.L3(4)");
QuotientBlocks63("12_2.L3(4)","4_1.L3(4)");
QuotientBlocks63("12_2.L3(4)","4_2.L3(4)");
QuotientBlocks63("6.U6(2)","2.U6(2)");
QuotientBlocks63("6.2E6(2)","2.2E6(2)");
Print("FRESH_CENTRAL_BLOCK_QUOTIENTS_COMPLETE\n");
QUIT;
