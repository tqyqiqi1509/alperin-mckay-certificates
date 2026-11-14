"""Independent order-three E8 root geometry.  No predecessor arrays imported."""
from pathlib import Path
from itertools import combinations, product
from functools import reduce
from math import gcd
import ast, json, re
import sympy as S
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

DIR=Path(__file__).resolve().parent
eu=[]
for i,j in combinations(range(8),2):
    for a,b in product((-1,1),repeat=2):
        r=[0]*8;r[i]=a;r[j]=b;eu.append(S.Matrix(r))
for signs in product((-1,1),repeat=8):
    if sum(a<0 for a in signs)%2==0:eu.append(S.Matrix(signs)/2)
assert len(eu)==240
B=hermite_normal_form(S.Matrix.hstack(*(2*r for r in eu)).applyfunc(int))/2
assert B.det()==1
C=B.T*B;I=S.eye(8);BI=B.inv()
roots=sorted((BI*r for r in eu),key=lambda r:tuple(r))
assert all(v.q==1 for r in roots for v in r)
RSET={tuple(r) for r in roots};INDEX={tuple(r):i+1 for i,r in enumerate(roots)}
height=S.Matrix([[5**i for i in range(8)]])*B
dot=lambda r,t:(r.T*C*t)[0]
def simple(rr):
    positive=[r for r in rr if (height*r)[0]>0]
    ps={tuple(r) for r in positive}
    return sorted([r for r in positive if not any(tuple(r-a) in ps for a in positive)],key=lambda r:tuple(r),reverse=True)
def comp(bb):
    remaining=set(range(len(bb)));out=[]
    while remaining:
        todo=[max(remaining)];cc=[];remaining.remove(todo[0])
        while todo:
            k=todo.pop();cc.append(k)
            new=[j for j in remaining if dot(bb[k],bb[j])!=0]
            todo+=new;remaining.difference_update(new)
        out.append(sorted(cc))
    return sorted(out,key=lambda cc:(-len(cc),cc))
def graph_type(bb):
    out=[]
    for cc in comp(bb):
        adj={i:[j for j in cc if j!=i and dot(bb[i],bb[j])!=0] for i in cc}
        branch=[i for i in cc if len(adj[i])==3]
        if not branch:out.append('A'+str(len(cc)));continue
        assert len(branch)==1
        m=branch[0];arms=[]
        for x in adj[m]:
            past=m;length=1
            while len(adj[x])>1:
                nxt=next(y for y in adj[x] if y!=past);past,x=x,nxt;length+=1
            arms.append(length)
        arms=sorted(arms)
        out.append(('D' if arms[:2]==[1,1] else 'E')+str(len(cc)))
    return '+'.join(out) if out else '1'
def smith(bb):
    if not bb:return []
    z=smith_normal_form(S.Matrix.hstack(*bb),domain=S.ZZ)
    return [abs(int(z[i,i])) for i in range(len(bb))]
def reflection(r):return I-r*(r.T*C)
def based(v,rr,bb):
    ps={tuple(r) for r in rr if (height*r)[0]>0};out=v
    for _ in range(2000):
        neg=next((r for r in bb if tuple(out*r) not in ps),None)
        if neg is None:return out
        out=out*reflection(neg)
    raise ValueError('based action did not terminate')
def restricted_charpoly(v,bb):
    V=S.Matrix.hstack(*bb) if bb else S.zeros(8,0)
    return str(V.gauss_jordan_solve(v*V)[0].charpoly().as_expr().factor()) if V.cols else '1'
def nullorth(bb):return S.Matrix.vstack(*(r.T*C for r in bb)).nullspace() if bb else list(I.columnspace())
def chain(bb,cc):
    adj={i:[j for j in cc if i!=j and dot(bb[i],bb[j])] for i in cc}
    seq=[max(i for i in cc if len(adj[i])==1)]
    while len(seq)<len(cc):seq.append(next(i for i in adj[seq[-1]] if i not in seq))
    return [bb[i] for i in seq]

EB=simple(roots);EM=S.Matrix.hstack(*EB);EC=EM.T*C*EM
assert len(EB)==8 and EC.det()==1
highest=max((EM.inv()*r for r in roots),key=lambda r:sum(r))
families={}
for j in range(8):
    if highest[j]!=3:continue
    s=EM*EC.inv()[:,j]/3
    cr=[r for r in roots if dot(r,s).q==1]
    cb=simple(cr);ct=graph_type(cb)
    if ct=='A8':ordered=chain(cb,comp(cb)[0]);family='A8'
    else:
        assert sorted(ct.split('+'))==['A2','E6'];family='E6A2'
        cc=comp(cb);ec=next(x for x in cc if len(x)==6);ac=next(x for x in cc if len(x)==2)
        adj={i:[j for j in ec if i!=j and dot(cb[i],cb[j])] for i in ec}
        mid=next(i for i in ec if len(adj[i])==3)
        short=next(i for i in adj[mid] if len(adj[i])==1)
        near=sorted((i for i in adj[mid] if i!=short),reverse=True)
        far=[next(j for j in adj[i] if j!=mid) for i in near]
        ordered=[cb[i] for i in [far[0],short,near[0],mid,near[1],far[1]]]+chain(cb,ac)
    CB=S.Matrix.hstack(*ordered)
    assert abs(CB.det())==3
    families[family]=(s,cr,ordered,CB)
assert len(families)==2

raw=(DIR/'generic_sources.txt').read_text();x=S.Symbol('x');settings={};tails=[]
for chunk in raw.split('GROUP60|')[1:]:
    fam='A8' if chunk.lstrip().startswith('[ "A"') else 'E6A2'
    ss,cr,cb,CB=families[fam]
    for m in re.finditer(r'PAIR60\|(-?\d+)\|(\d+)\|([^\n]*)\nK60\|(.*?)\nF60\|(.*?)\nEND60',chunk,re.S):
        eps,d=int(m[1]),int(m[2]);desc,mark=m[3].rsplit('|',1)
        kb0=[CB*S.Matrix(r) for r in ast.literal_eval(m[4])]
        v=CB*S.Matrix(ast.literal_eval(m[5])).T*CB.inv()
        assert v.T*C*v==C and all(z.q==1 for z in v)
        assert all(z.q==1 for z in eps*v*ss-ss)
        poly=S.Poly(S.cyclotomic_poly(d,x),x);eq=S.zeros(8)
        for (j,),coef in poly.terms():eq+=coef*v**j
        if kb0:eq=eq.col_join(S.Matrix.vstack(*(r.T*C for r in kb0)))
        active=eq.nullspace();mult=len(active)//int(S.totient(d))
        if mult<2:
            tails.append({'family':fam,'q_mod3':eps%3,'d':d,'source':desc,'mark':mark,'active_multiplicity':mult})
            continue
        lr=[r for r in roots if all(dot(r,a)==0 for a in active)]
        kr=[r for r in lr if dot(r,ss).q==1];lb=simple(lr);kb=simple(kr)
        assert len(kb)==len(kb0)
        assert S.Matrix.hstack(*(kb+kb0)).rank()==len(kb) if kb else not kb0
        key=(fam,eps,d,tuple(INDEX[tuple(r)] for r in lr),tuple(v))
        if key in settings:
            settings[key]['marks'].append(mark);continue
        pr=[r for r in roots if all(dot(r,k)==0 for k in kb)];pb=simple(pr)
        jr=[r for r in cr if all(dot(r,k)==0 for k in kb)];jb=simple(jr)
        lf=based(v,lr,lb);jf=based(v,jr,jb)
        lgraph=[lb.index(lf*r) for r in lb];jgraph=[jb.index(jf*r) for r in jb]
        V=S.Matrix.hstack(*active);P=V*(V.T*C*V).inv()*V.T*C
        den=S.ilcm(*[z.q for z in list(P)+list(P*v)])
        U=den*P*v;W=-den*P;norms=[]
        if d in (3,4,6):
            for r in roots:
                rr=r.T*C*U;tt=r.T*C*W
                if d==3:ns=[int(a*a-a*b+b*b) for a,b in zip(rr,tt)]
                elif d==6:ns=[int(a*a+a*b+b*b) for a,b in zip(rr,tt)]
                else:ns=[int(a*a+b*b) for a,b in zip(rr,tt)]
                g=reduce(gcd,ns,0);norms.append(g)
                assert (g==0)==(r in lr)
                assert all(k in (2,3,5) for k in S.factorint(g)) if g else True
        else:
            assert d in (1,2)
            for r in roots:
                vals=[abs(int(z)) for z in r.T*C*(den*P)]
                g=reduce(gcd,vals,0);norms.append(g)
                assert (g==0)==(r in lr)
                assert all(k in (2,3,5) for k in S.factorint(g)) if g else True
        witness=None
        for coeff in product(range(3),repeat=len(lb)):
            z=ss-sum((a*r for a,r in zip(coeff,lb)),S.zeros(8,1))/3
            if all(t.q==1 for t in z):witness=list(coeff);break
        item={'family':fam,'q_mod3':eps%3,'eps':eps,'d':d,'source':desc,'marks':[mark],
          's':list(ss),'F':v.tolist(),'active_multiplicity':mult,'active_basis':[list(a) for a in active],
          'L_type':graph_type(lb),'K_type':graph_type(kb),'perpendicular_type':graph_type(pb),'J_type':graph_type(jb),
          'L_roots':[INDEX[tuple(r)] for r in lr],'L_base':[list(r) for r in lb],'K_base':[list(r) for r in kb],
          'perpendicular_base':[list(r) for r in pb],'J_base':[list(r) for r in jb],
          'L_smith':smith(lb),'K_smith':smith(kb),'perpendicular_smith':smith(pb),
          'L_graph':lgraph,'J_graph':jgraph,'L_center':restricted_charpoly(v,nullorth(lb)),
          'J_center':restricted_charpoly(v,nullorth(kb+jb)),
          'derived_witness_mod3':witness,'projector_den':int(den),'restriction_gcds':sorted(set(norms))}
        settings[key]=item
        print(fam,eps,d,desc,'L=',item['L_type'],'K=',item['K_type'],'J=',item['J_type'],'center=',item['L_center'],'F=',lgraph,'norms=',sorted(set(norms)),flush=True)

items=list(settings.values());assert len(items)==30,len(items)
for item in items:
    assert all(t==1 for t in item['L_smith'])
    assert all(t==1 for t in item['K_smith'])
    assert all(t==1 for t in item['perpendicular_smith'])
output={'lattice_basis':B.tolist(),'gram':C.tolist(),'simple_basis':[list(r) for r in EB],
        'highest_root':list(highest),'families':{k:{'s':list(v[0]),'base':[list(r) for r in v[2]]} for k,v in families.items()},
        'settings':items,'tails':tails}
(DIR/'geometry.json').write_text(json.dumps(output,indent=2,default=str)+'\n')
def perm(mat):return 'PermList('+str([INDEX[tuple(mat*r)] for r in roots])+')'
gg=['refs60:=['+','.join(perm(reflection(r)) for r in roots)+'];;']
for idx,item in enumerate(items):
    s,cr,cb,CB=families[item['family']];v=S.Matrix(item['F']);kb=[S.Matrix(r) for r in item['K_base']]
    cbi=[INDEX[tuple(r)] for r in cb];kbi=[INDEX[tuple(r)] for r in kb]
    gg+=['wc60:=Group(List('+str(cbi)+',j->refs60[j]),());;',
       'nn60:=Stabilizer(wc60,'+str(item['L_roots'])+',OnSets);;',
       'wk60:=Group(List('+str(kbi)+',j->refs60[j]),());;',
       'f60:='+perm(v)+';;',
       'ext60:=Group(Concatenation(GeneratorsOfGroup(nn60),[f60]));;',
       'map60:=NaturalHomomorphismByNormalSubgroup(ext60,wk60);;',
       'rel60:=Intersection(Image(map60,nn60),Centralizer(Image(map60),Image(map60,f60)));;',
       'Print("ROW60|'+str(idx)+'|",Size(wc60),"|",Size(nn60),"|",Size(wk60),"|",Size(rel60),"|",StructureDescription(rel60),"\\n");']
gg+=['QUIT;']
(DIR/'normalizers.g').write_text('\n'.join(gg)+'\n')
print('COMPLETE',len(items),'settings;',len(tails),'proper cyclic/zero marked tails',flush=True)
