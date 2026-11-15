LoadPackage("ctbllib");;
SetInfoLevel(InfoWarning,0);;
families60:=["2.L3(4)","4_1.L3(4)","4_2.L3(4)","6.L3(4)",
 "12_1.L3(4)","12_2.L3(4)","2.U6(2)","6.U6(2)",
 "2.O8+(2)","2.2E6(2)","6.2E6(2)","2.2E6(2).2"];;
for name60 in families60 do
 table60:=CharacterTable(name60);;
 Print("TABLE_INPUT|",name60,"|order_factors=",Collected(FactorsInt(Size(table60))),
       "|centre_orders=",OrdersClassRepresentatives(table60){ClassPositionsOfCentre(table60)},"\n");
 for prime60 in [3,5,7] do
  wascached60:=IsBound(ComputedPrimeBlockss(table60)[prime60]);;
  fresh60:=PrimeBlocksOp(table60,prime60);;
  usual60:=PrimeBlocks(table60,prime60);;
  Assert(0,fresh60.block=usual60.block and fresh60.defect=usual60.defect and fresh60.height=usual60.height);
  Print("UNCACHED_BLOCK_CHECK|",name60,"|p=",prime60,"|before_cached=",wascached60,
        "|blocks=",Length(fresh60.defect),"|defects=",Collected(fresh60.defect),"\n");
 od;
od;

pairs60:=[["6.L3(4)","2.L3(4)"],["12_1.L3(4)","4_1.L3(4)"],
 ["12_2.L3(4)","4_2.L3(4)"],["6.U6(2)","2.U6(2)"],
 ["6.2E6(2)","2.2E6(2)"]];;
for pair60 in pairs60 do
 up60:=CharacterTable(pair60[1]);; lo60:=CharacterTable(pair60[2]);;
 fusion60:=GetFusionMap(up60,lo60);; Assert(0,fusion60<>fail);
 kernel60:=Positions(fusion60,1);;
 Assert(0,Length(kernel60)=3 and ForAll(kernel60,j->SizesConjugacyClasses(up60)[j]=1));
 Assert(0,SortedList(OrdersClassRepresentatives(up60){kernel60})=[1,3,3]);
 uv60:=List(Irr(up60),ValuesOfClassFunction);; lv60:=List(Irr(lo60),ValuesOfClassFunction);;
 pull60:=List(lv60,v->Position(uv60,v{fusion60}));;
 kills60:=Filtered([1..Length(uv60)],j->ForAll(kernel60,k->uv60[j][k]=uv60[j][1]));;
 Assert(0,not fail in pull60 and Length(Set(pull60))=Length(pull60) and Set(pull60)=Set(kills60));
 ub60:=PrimeBlocksOp(up60,3);; lb60:=PrimeBlocksOp(lo60,3);;
 Assert(0,ForAll([1..Length(pull60)],j->ub60.height[pull60[j]]=lb60.height[j]));
 for block60 in Set(lb60.block) do
  rows60:=Positions(lb60.block,block60);; upperblocks60:=Set(pull60{rows60},j->ub60.block[j]);;
  Assert(0,Length(upperblocks60)=1 and ub60.defect[upperblocks60[1]]=lb60.defect[block60]+1);
 od;
 Print("EXACT_CENTRAL_KERNEL|",pair60,"|kernel_classes=",kernel60,
       "|all_killed_rows=",Length(kills60),"|all_pullbacks=",Length(pull60),
       "|height_preservation=true|all_dominated_blocks=true\n");
 if pair60[1]="6.2E6(2)" then
  h0up60:=Filtered(Positions(ub60.block,18),j->ub60.height[j]=0);;
  h0lo60:=Filtered(Positions(lb60.block,18),j->lb60.height[j]=0);;
  Assert(0,Set(h0up60)=Set(pull60{h0lo60}) and Length(h0up60)=6);
  Assert(0,Collected(ub60.height{Positions(ub60.block,18)})=[[0,6],[1,8]]);
  Print("C4_LIFT_EXACT_HEIGHTS|",Collected(ub60.height{Positions(ub60.block,18)}),
        "|inflated_h0=",h0up60,"\n");
 fi;
od;

low60:=CharacterTable("2.2E6(2)");; high60:=CharacterTable("2.2E6(2).2");;
f60:=GetFusionMap(low60,high60);; Assert(0,f60<>fail);;
lowblocks60:=PrimeBlocksOp(low60,3);; highblocks60:=PrimeBlocksOp(high60,3);;
highvals60:=List(Irr(high60),ValuesOfClassFunction);;
reg60:=Filtered([1..NrConjugacyClasses(high60)],j->OrdersClassRepresentatives(high60)[j] mod 3<>0);;
for index60 in [20,21,22] do
 rows60:=Positions(highblocks60.block,index60);;
 brauer60:=RankMat(List(rows60,j->highvals60[j]{reg60}));;
 covered60:=[];
 for j60 in rows60 do
  restricted60:=ClassFunction(low60,highvals60[j60]{f60});;
  for k60 in [1..NrConjugacyClasses(low60)] do
   if ScalarProduct(Irr(low60)[k60],restricted60)>0 then AddSet(covered60,lowblocks60.block[k60]);fi;
  od;
 od;
 Assert(0,highblocks60.defect[index60]=2 and Length(rows60)=9 and brauer60=5);
 Assert(0,ForAll(rows60,j->highblocks60.height[j]=0));
 if index60=22 then Assert(0,covered60=[18]);else Assert(0,covered60=[17]);fi;
 Print("ACTUAL_UPPER_BLOCK_ANCHOR|",index60,"|covered=",covered60,
       "|defect=2|k=9|k0=9|l=",brauer60,"\n");
od;
Print("C60_INDEPENDENT_BLOCK_KERNEL_COMPLETE\n");
QUIT;
