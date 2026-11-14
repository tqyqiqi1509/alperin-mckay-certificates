"""Fresh integral Chevalley/Tits pinning for the difficult order-four sources.

The source geometries are node58's independent HNF/signed-cycle output.
The signed adjoint representation is reconstructed from a root-lattice
cocycle in a freshly recovered E8 simple basis, not read from Geck data.
"""
from pathlib import Path
from itertools import product, combinations
import json
import sympy as S

OUT=Path(__file__).parent;E=S.eye(8)
data=json.loads((OUT/'independent_geometry.json').read_text())
eu=[]
for i,j in combinations(range(8),2):
    for a,b in product((-1,1),repeat=2):
        r=S.zeros(8,1);r[i]=a;r[j]=b;eu.append(r)
for v in product((-1,1),repeat=8):
    if sum(t<0 for t in v)%2==0:eu.append(S.Matrix(v)/2)
h=S.Matrix([3**i for i in range(8)]);pos=[r for r in eu if (r.T*h)[0]>0];ps=set(map(tuple,pos))
bs=[r for r in pos if not any(tuple(r-t) in ps for t in pos)]
B=S.Matrix.hstack(*bs);Bi=B.inv();C=B.T*B
rr=[tuple(int(t) for t in Bi*r) for r in eu];ri={r:i for i,r in enumerate(rr)}
def eps(a,b):return (-1)**(sum(a[i]*b[i] for i in range(8))+sum(a[i]*int(C[i,j])*b[j] for i in range(8) for j in range(i)))
ns=[];ws=[];checks=0
def add(u,v):
    z=dict(u)
    for k,t in v.items():z[k]=z.get(k,0)+t
    return {k:t for k,t in z.items() if t}
def scale(u,t):return {k:t*v for k,v in u.items() if t*v}
def ad(a,v):
    z={}
    for k,c in v.items():
        if k<240:
            b=rr[k]
            if tuple(-t for t in a)==b:z=add(z,{240+i:-c*a[i] for i in range(8) if a[i]})
            else:
                ab=tuple(t+u for t,u in zip(a,b))
                if ab in ri:z=add(z,{ri[ab]:c*eps(a,b)})
        else:
            m=-c*sum(int(C[k-240,j])*a[j] for j in range(8))
            if m:z=add(z,{ri[a]:m})
    return z
def expad(a,v):
    d=ad(a,v);dd=ad(a,d);assert not ad(a,dd)
    return add(add(v,d),scale(dd,S.Rational(1,2)))
for i in range(8):
    a=tuple(int(t) for t in E[:,i]);act=[];w=E.copy();w[i,:]-=C[i,:];ws.append(w)
    for r in rr:
        m=sum(int(C[i,j])*r[j] for j in range(8));z=list(r);z[i]-=m
        act.append((ri[tuple(z)],eps(a,r) if abs(m)==1 else 1))
    ns.append(act)
    for j in range(248):
        got=expad(a,expad(tuple(-t for t in a),expad(a,{j:1})))
        want={act[j][0]:act[j][1]} if j<240 else {240+k:int(w[k,j-240]) for k in range(8) if w[k,j-240]}
        assert got==want and all(S.Rational(t).q==1 for t in got.values());checks+=1
def signed_mul(a,b):return [(a[j][0],sg*a[j][1]) for j,sg in b]
ident=[(i,1) for i in range(240)]
bracket_checks=0
for i in range(8):
    assert signed_mul(ns[i],ns[i])==[(j,int((-1)**sum(int(C[i,k])*r[k] for k in range(8)))) for j,r in enumerate(rr)]
    for j in range(i):
        if C[i,j]==0:assert signed_mul(ns[i],ns[j])==signed_mul(ns[j],ns[i])
        else:assert signed_mul(ns[i],signed_mul(ns[j],ns[i]))==signed_mul(ns[j],signed_mul(ns[i],ns[j]))
    for r in rr:
        a,sg=ns[i][ri[r]]
        for t in rr:
            rt=tuple(u+v for u,v in zip(r,t))
            if rt not in ri:continue
            b,sh=ns[i][ri[t]];c,sk=ns[i][ri[rt]]
            assert tuple(u+v for u,v in zip(rr[a],rr[b]))==rr[c]
            assert eps(r,t)*sk==sg*sh*eps(rr[a],rr[b]);bracket_checks+=1
def perm(a):return [j+1+(240 if sg<0 else 0) for j,sg in a]+[j+1+(240 if sg>0 else 0) for j,sg in a]
def wperm(w):return [ri[tuple(w*S.Matrix(r))]+1 for r in rr]
def tor(y,m=2):
    shifts=[int((S.Matrix(r).T*C*y)[0])%m for r in rr]
    return [a+1+240*((j+shifts[a])%m) for j in range(m) for a in range(240)]
g=['Reset(GlobalMersenneTwister,580001);;Reset(GlobalRandomSource,580001);;',
   'nn58:=List('+str([perm(n) for n in ns])+',PermList);;',
   'ww58:=List('+str([wperm(w) for w in ws])+',PermList);;',
   'T58:=Group(nn58);;W58:=Group(ww58);;',
   'phi58:=GroupHomomorphismByImages(T58,W58,nn58,ww58);;Assert(0,phi58<>fail);',
   'Print("CHEVALLEY_PRODUCTS|'+str(checks)+'|TITS_KERNEL|",Size(Kernel(phi58)),"\\n");',
   'ref58:=List('+str([wperm(E-S.Matrix(r)*S.Matrix(r).T*C) for r in rr])+',PermList);;',
   'Based58:=function(w,lb,ps) local b;while ForAny(lb,b->not b^w in ps) do b:=First(lb,b->not b^w in ps);w:=ref58[b]*w;od;return w;end;;',
   'Lift4:=function(x) local a,b,j,k,pp;pp:=[];for j in [0..3] do for a in [1..240] do b:=a^x;k:=0;if b>240 then b:=b-240;k:=2;fi;Add(pp,b+240*((j+k) mod 4));od;od;return PermList(pp);end;;',
   'phase58:=PermList(Concatenation(List([0..3],j->List([1..240],a->a+240*((-j) mod 4)))));;',
   'ChooseLine58:=function(G,Z,vals) local ab,map,ch,rho;map:=NaturalHomomorphismByNormalSubgroup(G,DerivedSubgroup(G));ab:=Image(map);ch:=First(Irr(ab),z->ForAll([1..Length(Z)],i->Image(map,Z[i])^z=vals[i]));Assert(0,ch<>fail);rho:=List(GeneratorsOfGroup(G),z->Image(map,z)^ch);return rec(values:=rho,map:=map,character:=ch);end;;']
summary={}
keys=['D5A3_nonmatching_regular','A7A1_nonmatching_regular','D5A3_quartic31','D5A3_quartic33','A7A1_quartic36']
for key in keys:
    d=data['cases'][key];s=Bi*S.Matrix(d['s']);w=Bi*S.Matrix(d['F'])*B
    lb=[Bi*S.Matrix(r) for r in d['L_base']]
    # Use the root span directly, without an output-root list.
    LM=S.Matrix.hstack(*lb);rank=LM.cols
    lr=[S.Matrix(r) for r in rr if S.Matrix.hstack(LM,S.Matrix(r)).rank()==rank]
    cr=[S.Matrix(r) for r in rr if (S.Matrix(r).T*C*s)[0].q==1]
    cpos=[r for r in cr if all(t>=0 for t in r)];cp=set(map(tuple,cpos))
    cb=[r for r in cpos if not any(tuple(r-t) in cp for t in cpos)]
    lbi=[ri[tuple(r)]+1 for r in lb];lp=[ri[tuple(r)]+1 for r in lr if all(t>=0 for t in r)]
    lw=E
    while any(all(t>=0 for t in lw*b) for b in lb):
        b=next(b for b in lb if all(t>=0 for t in lw*b));lw=lw*(E-b*b.T*C)
    hhmat=-lw
    g += ['Print("CASE|'+key+'\\n");',
          'WC:=Group(List('+str([ri[tuple(r)]+1 for r in cb])+',i->ref58[i]));;',
          'NC:=Stabilizer(WC,'+str(sorted(ri[tuple(r)]+1 for r in lr))+',OnSets);;',
          'lb:='+str(lbi)+';;ps:='+str(lp)+';;',
          'V:=Group(List(GeneratorsOfGroup(NC),v->Based58(v,lb,ps)));;',
          'E0:=Stabilizer(PreImage(phi58,V),Set(lb),OnSets);;pe:=RestrictedMapping(phi58,E0);;',
          'aa:=Based58(PermList('+str(wperm(w))+'),lb,ps);;',
          ('nn:=PreImagesRepresentative(RestrictedMapping(pe,Centre(E0)),aa);;Assert(0,nn<>fail);' if d['d']==2 else 'nn:=PreImagesRepresentative(pe,aa);;'),
          'hh:=PreImagesRepresentative(RestrictedMapping(pe,Centre(E0)),PermList('+str(wperm(hhmat))+'));;Assert(0,hh<>fail);',
          'CV:=Centralizer(V,aa);;CE:=Centralizer(E0,nn);;',
          'Print("BASE|",Size(V),"|",Size(E0),"|",Size(CV),"|",Size(CE),"|",Size(Image(pe,CE)),"\\n");']
    zrat=[];zsg=[];period=[]
    for yy in product(range(2),repeat=8):
        y=S.Matrix(yy)
        if any(int(t)%2 for t in LM.T*C*y) or any(int(t)%2 for t in (w-E)*y):continue
        zz=(s.T*C*(w-E)*y/2)[0];assert (2*zz).q==1
        sign=int((-1)**int(2*zz));zrat.append(tor(y));zsg.append(sign)
        for sg,q in [(1,1),(1,5),(-1,3),(-1,7)]:
            z=(s.T*C*(sg*q*w-E)*y/2)[0];assert (2*z).q==1 and (-1)**int(2*z)==sign
    summary[key]={'rational_two_kernel':len(zrat),'zeta_values':zsg,'source_rank':rank}
    g += ['zz:=List('+str(zrat)+',PermList);;sv:='+str(zsg)+';;']
    if d['d']==2:
        g+=['Assert(0,Size(Image(pe,CE))=Size(CV));',
            'rho:=ChooseLine58(CE,zz,sv);;',
            'Print("CENTRAL_LINE|",AbelianInvariants(CE),"|",rho.values,"\\n");',
            'Print("ACTUAL_RETURN|",Order(nn),"|",Image(rho.map,nn)^rho.character,"|",Order(nn*hh),"|",Image(rho.map,nn*hh)^rho.character,"\\n");',
            'Print("PINNING_ACTION|",Size(Action(CE,Set(lb))),"\\n");']
    else:
        M=LM.T*C
        piv=next(cc for cc in combinations(range(8),rank) if int(M[:,cc].det())%2)
        free=[i for i in range(8) if i not in piv];mi=M[:,piv].inv_mod(4);z4=[]
        for j in free:
            y=S.zeros(8,1);y[j]=1;v=-mi*M[:,j]
            for a,k in enumerate(piv):y[k]=int(v[a])%4
            assert all(int(t)%4==0 for t in M*y);z4.append(tor(y,4))
        g += ['zg4:=List('+str(z4)+',PermList);;eg4:=List(GeneratorsOfGroup(E0),Lift4);;',
              'P4:=Group(Concatenation(zg4,eg4));;',
              'p4:=GroupHomomorphismByImages(P4,V,Concatenation(zg4,eg4),Concatenation(List(zg4,z->()),List(GeneratorsOfGroup(E0),z->Image(pe,z))));;',
              'P:=Centralizer(P4,Lift4(nn));;QQ:=Group(Concatenation(GeneratorsOfGroup(P),[phase58]));;',
              'Assert(0,Size(Image(p4,P))=Size(CV));',
              'rho:=ChooseLine58(QQ,List(zz,Lift4),sv);;',
              'Print("FULL_PIN|",Size(P4),"|",Size(P),"|",Size(Image(p4,P)),"|",Size(QQ),"|",AbelianInvariants(QQ),"|",rho.values,"\\n");',
              'Assert(0,ForAll(GeneratorsOfGroup(P),v->v^phase58=v^Lift4(hh)));',
              'Print("ACTUAL_RETURN|",Order(nn),"|",Image(rho.map,Lift4(nn))^rho.character,"|",Order(nn*hh),"|",Image(rho.map,Lift4(nn*hh))^rho.character,"\\n");',
              'Print("PINNING_ACTION|",Size(Action(P,Set(lb))),"|",List(GeneratorsOfGroup(P),v->List(lb,b->Position(lb,b^v))),"\\n");']
g.append('QUIT;')
(OUT/'independent_pinning.g').write_text('\n'.join(g)+'\n')
(OUT/'independent_pinning.json').write_text(json.dumps({'Chevalley_product_checks':checks,'root_bracket_checks':bracket_checks,'Tits_braid_and_square_checks':True,'basis':B.tolist(),'sources':summary},indent=2,default=str)+'\n')
print('fresh E8 Chevalley carrier',checks,'products;',bracket_checks,'brackets; prepared',len(keys),'whole-center source pinning groups')
