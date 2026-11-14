"""Independent order-two/four E8 root data, node 58.

Input consists solely of the Euclidean E8 roots, four explicit torus
points, and signed-cycle/Coxeter constructions.  No predecessor result
array or Frobenius matrix is loaded.  An HNF basis of the complete E8
lattice is used instead of the predecessor Bourbaki basis.
"""
from itertools import combinations, product
from functools import reduce
from math import gcd
from pathlib import Path
import json
import sympy as S
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ

I=S.eye(8); x=S.Symbol('x'); OUT=Path(__file__).parent
rr=[]
for i,j in combinations(range(8),2):
    for a,b in product((-1,1),repeat=2):
        v=S.zeros(8,1);v[i]=a;v[j]=b;rr.append(v)
for v in product((-1,1),repeat=8):
    if sum(t<0 for t in v)%2==0:rr.append(S.Matrix(v)/2)
assert len(rr)==240
Q=hermite_normal_form(S.Matrix.hstack(*[2*r for r in rr]).applyfunc(S.Integer))/2
assert abs(Q.det())==1
Qi=Q.inv(); rid={tuple(r):i+1 for i,r in enumerate(rr)}
generic=S.Matrix([3**j for j in range(8)])
dot=lambda u,v:(u.T*v)[0]
rf=lambda r:I-r*r.T
def cols(a):return S.Matrix.hstack(*a) if a else S.zeros(8,0)
def pos(r):return dot(generic,r)>0
def base(rs):
    pp=[r for r in rs if pos(r)];ss=set(map(tuple,pp))
    return [r for r in pp if not any(tuple(r-t) in ss for t in pp)]
def smith(m):
    return [abs(int(t)) for t in smith_normal_form(m,domain=ZZ).diagonal() if t] if m.cols else []
def components(bs):
    todo=set(range(len(bs)));out=[]
    while todo:
        cc={todo.pop()}
        while True:
            added={j for i in cc for j in todo if dot(bs[i],bs[j])}
            if not added:break
            cc|=added;todo-=added
        out.append(sorted(cc))
    return out
def diagram(bs):
    types=[]
    for cc in components(bs):
        n=len(cc);adj={i:[j for j in cc if i!=j and dot(bs[i],bs[j])] for i in cc}
        branch=[i for i in cc if len(adj[i])==3]
        if not branch:types.append('A'+str(n));continue
        assert len(branch)==1
        b=branch[0];arms=[]
        for start in adj[b]:
            prev=b;cur=start;narm=1
            while len(adj[cur])==2:
                nxt=next(t for t in adj[cur] if t!=prev);prev,cur=cur,nxt;narm+=1
            arms.append(narm)
        arms.sort()
        if arms[:2]==[1,1]:types.append('D'+str(n))
        else:
            assert arms in ([1,2,2],[1,2,3],[1,2,4]);types.append('E'+str(n))
    return sorted(types)
def based(w,bs):
    while any(not pos(w*b) for b in bs):
        b=next(b for b in bs if not pos(w*b));w=w*rf(b)
    assert all(w*b in bs for b in bs)
    return w
def graph(w,bs):
    if not bs:return []
    v=based(w,bs);return [bs.index(v*b) for b in bs]
def prodref(bs):
    w=I
    for b in bs:w=w*rf(b)
    return w
def negcycle(seq):
    w=I.copy()
    for i,j in zip(seq,seq[1:]):w[:,i]=I[:,j]
    w[:,seq[-1]]=-I[:,seq[0]]
    return w
def cycle(seq):
    w=I.copy()
    for i,j in zip(seq,seq[1:]+seq[:1]):w[:,i]=I[:,j]
    return w
def dbase(seq):
    return [I[:,i]-I[:,j] for i,j in zip(seq,seq[1:])]+[I[:,seq[-2]]+I[:,seq[-1]]]
def chain(bs,cc):
    adj={i:[j for j in cc if j!=i and dot(bs[i],bs[j])] for i in cc}
    order=[min(i for i in cc if len(adj[i])==1)]
    while len(order)<len(cc):order.append(next(i for i in adj[order[-1]] if i not in order))
    return [bs[i] for i in order]
def subpoly(w,bs):
    if not bs:return '1'
    M=cols(bs);A=M.gauss_jordan_solve(w*M)[0];assert M*A==w*M
    return str(S.factor(A.charpoly(x).as_expr()))
def active(w,d):
    cp=S.Poly(S.cyclotomic_poly(d,x),x);A=S.zeros(8)
    for (e,),c in cp.terms():A+=c*w**e
    return A.nullspace()
def perm(w):return 'PermList('+str([rid[tuple(w*r)] for r in rr])+')'

theta=I[:,6]+I[:,7]
points={'D8':I[:,7], 'E7A1':theta/2,
        'D5A3':S.Matrix([0,0,0,0,0,1,1,1])/2,
        'A7A1':S.Matrix([-1,1,1,1,1,1,1,7])/8}
settings=[]
def add(f,n,w,d,z=None):settings.append((f,n,w,d,z))

# Involution D8: signed cycles give all active cyclotomic packets.
add('D8','toral1',I,1);add('D8','toral2',-I,2)
for d,sg in [(1,1),(2,-1)]:add('D8','D4_'+str(d),sg*I,d,[I[:,j] for j in range(4,8)])
w3=cycle([0,1,2])*cycle([3,4,5])
add('D8','cubic',w3,3);add('D8','sextic',-w3,6)
add('D8','quartic_toral',negcycle([0,1])*negcycle([2,3])*negcycle([4,5])*negcycle([6,7]),4)
add('D8','quartic_D4',negcycle([0,1])*negcycle([2,3]),4)
add('D8','octic_toral',negcycle([0,1,2,3])*negcycle([4,5,6,7]),8)

# Involution E7A1: E7 Coxeter cube packets and an orthogonal D3^2.
e7=[r for r in rr if dot(r,theta)==0];e7b=base(e7)
c7=prodref(e7b);assert c7**18==I
add('E7A1','toral1',I,1);add('E7A1','toral2',-I,2)
kD=dbase([0,1,2,3]);zD=cols(kD).T.nullspace()
for d,sg in [(1,1),(2,-1)]:add('E7A1','D4_'+str(d),sg*I,d,zD)
# An E6 Levi: delete an endpoint from the actual based E7 diagram,
# selecting the first deletion with E6 graph, independently of node labels.
e6b=next([r for i,r in enumerate(e7b) if i!=j] for j in range(7)
         if diagram([r for i,r in enumerate(e7b) if i!=j])==['E6'])
for d,sg in [(1,1),(2,-1)]:add('E7A1','E6_'+str(d),sg*I,d,cols(e6b).T.nullspace())
add('E7A1','cubic',c7**6,3);add('E7A1','sextic',-(c7**6),6)
add('E7A1','quartic',prodref(dbase([0,1,2]))*prodref(dbase([3,4,5])),4)

# Order-four D5A3 in coordinate D5 x D3, and A7A1 from its graph.
for f in ['D5A3','A7A1']:
    s=points[f];cr=[r for r in rr if dot(r,s).q==1];cb=base(cr)
    w0=-based(-I,cb)
    add(f,'matching_toral',I,1);add(f,'nonmatching_regular',w0,2)
    if f=='D5A3':
        db=dbase([0,1,2,3,4]);ab=chain(cb,next(c for c in components(cb) if len(c)==3))
        add(f,'matching_D4',I,1,cols(dbase([0,1,2,3])).T.nullspace())
        add(f,'cubic',cycle([0,1,2])*prodref(ab[:2]),3)
        add(f,'quartic31',prodref(db)**2*prodref(ab),4)
        add(f,'quartic32',prodref(db)**2,4)
        add(f,'quartic33',negcycle([0,1])*negcycle([2])*prodref(ab),4)
    else:
        ab=chain(cb,next(c for c in components(cb) if len(c)==7));a1=cb[next(c[0] for c in components(cb) if len(c)==1)]
        add(f,'nonmatching_D6',rf(ab[0])*rf(a1),2)
        add(f,'cubic',prodref(ab[:2])*prodref(ab[3:5]),3)
        add(f,'quartic36',prodref(ab[:3])*prodref(ab[4:]),4)

SB=cols(base(rr));SBi=SB.inv()
g=['refs58:=['+','.join(perm(rf(r)) for r in rr)+'];;',
   'coords58:='+str([[int(t) for t in SBi*r] for r in rr])+';;',
   'simp58:='+str([rid[tuple(SB[:,i])] for i in range(8)])+';;',
   'Based58:=function(w,bs,ps) local r;while ForAny(bs,r->not r^w in ps) do r:=First(bs,r->not r^w in ps);w:=refs58[r]*w;od;return w;end;;']
out={'lattice_basis':Q.tolist(),'cases':{}}
for f,n,w,d,z in settings:
    s=points[f];cr=[r for r in rr if dot(r,s).q==1];cb=base(cr)
    assert w.T*w==I and all(t.q==1 for t in Qi*w*Q)
    assert all(t.q==1 for t in Qi*(w*s-s))
    z=active(w,d) if z is None else z
    lr=[r for r in rr if all(dot(r,a)==0 for a in z)];lb=base(lr)
    kr=[r for r in lr if dot(r,s).q==1];kb=base(kr)
    pr=[r for r in rr if all(dot(r,k)==0 for k in kb)];pb=base(pr)
    jr=[r for r in pr if dot(r,s).q==1];jb=base(jr)
    zfull=cols(lb).T.nullspace();assert len(zfull)==len(z)
    sn=smith(Qi*cols(lb));assert sn==[1]*len(lb)
    rec={'family':f,'name':n,'d':d,'s':list(s),'F':w.tolist(),
         'Cs_type':diagram(cb),'Cs_root_count':len(cr),'Cs_index':abs(int((Qi*cols(cb)).det())),
         'L':diagram(lb),'L_root_count':len(lr),'L_base':[list(r) for r in lb],
         'K':diagram(kb),'K_root_count':len(kr),'K_base':[list(r) for r in kb],
         'L_center_Smith':sn,'K_kernel_Smith':smith(Qi*cols(kb)),
         'perpendicular':diagram(pb),'perpendicular_base':[list(r) for r in pb],
         'perpendicular_Smith':smith(Qi*cols(pb)),'J0':diagram(jb),
         'J0_base':[list(r) for r in jb],'J0_Smith':smith(Qi*cols(jb)),
         'full_center_poly':subpoly(w,zfull),'L_F_graph':graph(w,lb),
         'K_F_graph':graph(w,kb),'perp_F_graph':graph(w,pb),'J0_F_graph':graph(w,jb),
         'J0_center_poly':subpoly(w,cols(kb+jb).T.nullspace()),
         'L_Ennola_graph':graph(-w,lb),'J0_Ennola_graph':graph(-w,jb)}
    CB=Qi*cols(cb);cinv=CB.inv()
    glue=next((cinv*I[:,j] for j in range(8)
               if int(S.ilcm(1,*[t.q for t in cinv*I[:,j]]))==rec['Cs_index']),None)
    if glue is not None:
        rec['Cs_glue_projection_orders']=[int(S.ilcm(1,*[glue[j].q for j in cc])) for cc in components(cb)]
    order=2 if f in ['D8','E7A1'] else 4
    target=tuple(int(t)%order for t in order*Qi*s)
    sm=Qi*cols(lb)
    wit=next((list(v) for v in product(range(order),repeat=len(lb))
              if tuple(int(t)%order for t in sm*S.Matrix(v))==target),None) if lb else None
    rec['s_in_Lderived_numerator']=wit
    # Full complement components: test the involution in its identity torus.
    if f=='E7A1' and n=='quartic':
        pm=Qi*cols(pb)
        rec['s_in_perp_derived']=any(tuple(int(t)%2 for t in pm*S.Matrix(v))==target for v in product(range(2),repeat=len(pb)))
        assert not rec['s_in_perp_derived']
        assert rec['K_kernel_Smith']==[1,1,1,2]
    # A denominator-cleared full-center projector, independently obtained
    # by orthogonal linear algebra, detects all relevant odd roots.
    Z=cols(zfull);proj=Z*(Z.T*Z).inv()*Z.T if Z.cols else S.zeros(8)
    den=int(S.ilcm(1,*[t.q for t in Qi*proj*Q]));pp=den*proj
    norms={}
    for r in rr:
        vals=[int(t) for t in r.T*pp*Q]
        g0=reduce(gcd,map(abs,vals),0);norms[g0]=norms.get(g0,0)+1
    rec['projector_denominator']=den;rec['restriction_integer_gcds']=norms
    if d in [3,4,6,8]:
        qq=S.Symbol('q');tt=S.Symbol('t');pol=S.cyclotomic_poly(d,tt)
        quotient=S.div(pol,tt-qq,tt)[0];P=S.zeros(8)
        for (power,),coeff in S.Poly(quotient,tt).terms():P+=coeff*w**power
        pmat=P*pp*Q;pn=S.cyclotomic_poly(d,qq);ng={};norm_cache={}
        for r in rr:
            vals=[S.Poly(t,qq).as_expr() for t in r.T*pmat]
            nv=[]
            for t in vals:
                if t not in norm_cache:norm_cache[t]=abs(int(S.resultant(pn,t,qq)))
                nv.append(norm_cache[t])
            v=reduce(gcd,nv,0);ng[v]=ng.get(v,0)+1
            assert (v==0)==(r in lr)
            if v:assert set(S.factorint(v))<=({3} if d in [3,6] else {2})
        rec['exact_cyclotomic_root_norm_gcds']=ng
    key=f+'_'+n;out['cases'][key]=rec
    print(key,'L',rec['L'],'K',rec['K'],'P',rec['perpendicular'],'J0',rec['J0'],
          'KSmith',rec['K_kernel_Smith'],'poly',rec['full_center_poly'],flush=True)
    li=sorted(rid[tuple(r)] for r in lr);lbi=[rid[tuple(r)] for r in lb]
    kbi=[rid[tuple(r)] for r in kb];cbi=[rid[tuple(r)] for r in cb]
    lp=[rid[tuple(r)] for r in lr if pos(r)]
    g += ['Print("CASE|'+key+'\\n");',
          'WC:=Group(List('+str(cbi)+',i->refs58[i]),());;',
          'NC:=Stabilizer(WC,'+str(li)+',OnSets);;',
          'WK:=Group(List('+str(kbi)+',i->refs58[i]),());;',
          'pq:=NaturalHomomorphismByNormalSubgroup(NC,WK);;V:=Image(pq);;',
          'ww:='+perm(w)+';;',
          'QQ:=Centralizer(V,Image(pq,ww));;',
          'Print("REL|",Size(WC),"|",Size(NC),"|",Size(WK),"|",Size(V),"|",Size(QQ),"\\n");']
    if d in [1,2] and ((f in ['D8','E7A1'] and not lb) or f in ['D5A3','A7A1']):
        mm=SBi*pp*SB;assert all(t.q==1 for t in mm)
        g += ['project58:='+str([[int(t) for t in row] for row in mm.tolist()])+';;',
              'for ell58 in [3,5,7] do if Size(QQ) mod ell58=0 then',
              'P58:=SylowSubgroup(QQ,ell58);;ff58:=GF(ell58);;rank58:=[];;',
              'for z58 in Elements(P58) do if z58<>One(P58) then',
              'zp58:=PreImagesRepresentative(pq,z58);;',
              'mat58:=TransposedMat(List(simp58,b->List(coords58[b^zp58],t->t*One(ff58))));;',
              'Add(rank58,RankMat((mat58-IdentityMat(8,ff58))^2*(project58*One(ff58))));',
              'fi;od;',
              'Assert(0,ForAll(rank58,v->v>0));',
              'Print("SQUARE_RANKS|",ell58,"|",Size(P58),"|",Collected(rank58),"\\n");',
              'fi;od;']
    if len(kb)>0 and diagram(kb)==['A1']*len(kb):
        kp=[rid[tuple(r)] for r in kr if pos(r)]
        g += ['ng:=GeneratorsOfGroup(NC);; kb:='+str(kbi)+';;kp:='+str(kp)+';;',
              'imgs:=List(ng,z->PermList(List(kb,r->Position(kb,r^Based58(z,kb,kp)))));;',
              'mh:=GroupHomomorphismByImages(NC,Group(imgs),ng,imgs);;Assert(0,mh<>fail);',
              'pre:=PreImage(pq,QQ);;',
              'Print("MARK|",Size(Image(mh,pre)),"|",Size(Intersection(pre,Kernel(mh)))/Size(WK),"\\n");']
g.append('QUIT;')
(OUT/'independent_geometry.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
(OUT/'independent_geometry.g').write_text('\n'.join(g)+'\n')
