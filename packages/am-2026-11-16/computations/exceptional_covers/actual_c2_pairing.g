# Exact orbit pairing at a specified global block and an actual local table.
# rp is computed from the actual local subgroup inclusion, not guessed.
ActualC2Pair63:=function(nm,blockno,ct,ids,rp)
 local gt,ut,fu,cp,j,ps,gr,gv,gb,gi,go,mo,omega,used,o,k;
 gt:=CharacterTable(nm);ut:=CharacterTable(Concatenation(nm,".2"));
 fu:=GetFusionMap(gt,ut);Assert(0,fu<>fail);cp:=[];
 for j in [1..Length(fu)] do
  ps:=Positions(fu,fu[j]);Assert(0,Length(ps)<=2);
  if Length(ps)=1 then cp[j]:=j;else cp[j]:=First(ps,x->x<>j);fi;
 od;
 gv:=List(Irr(gt),ValuesOfClassFunction);gr:=List(gv,v->Position(gv,v{cp}));
 gb:=PrimeBlocks(gt,3);gi:=Filtered(Positions(gb.block,blockno),j->gv[j][1] mod 3<>0);
 go:=Orbits(Group(PermList(gr)),gi);mo:=Orbits(Group(PermList(rp)),ids);
 omega:=[];used:=[];
 for o in go do
  k:=First([1..Length(mo)],j->not j in used and Length(mo[j])=Length(o));
  Assert(0,k<>fail);Add(used,k);omega[o[1]]:=mo[k][1];
  if Length(o)=2 then omega[gr[o[1]]]:=rp[mo[k][1]];fi;
 od;
 Assert(0,Set(omega{gi})=Set(ids) and ForAll(gi,j->omega[gr[j]]=rp[omega[j]]));
 Print("ACTUAL_C2_BIJECTION|group=",nm,"|block=",blockno,"|rows=",List(gi,j->[j,omega[j],gv[j][1],Irr(ct)[omega[j]][1]]),"\n");
end;;
