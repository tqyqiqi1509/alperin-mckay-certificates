"""Fresh integral geometry for the whole isolated order-three E8 family at five.

Only root lattices and finite Weyl groups occur. No finite-q sample is used.
"""
from itertools import combinations, product, permutations
from pathlib import Path
import json
import sympy as S
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

HERE = Path(__file__).resolve().parent
I = S.eye(8)
roots = []
for a,b in combinations(range(8),2):
    for u,v in product((-1,1),repeat=2):
        x=S.zeros(8,1); x[a]=u; x[b]=v; roots.append(x)
for signs in product((-1,1),repeat=8):
    if signs.count(-1)%2==0: roots.append(S.Matrix(signs)/2)
assert len(roots)==240
simple=[S.Matrix([1,-1,-1,-1,-1,-1,-1,1])/2,I[:,0]+I[:,1]]
simple += [I[:,i]-I[:,i-1] for i in range(1,7)]
B=S.Matrix.hstack(*simple); Bi=B.inv(); gram=B.T*B
assert gram.det()==1
key=lambda v:tuple(v)
rootkeys={key(r) for r in roots}
def reflection(r): return I-r*r.T
def positive(rr): return [r for r in rr if sum(Bi*r)>0]
def base(rr):
    pos=positive(rr); keys={key(r) for r in pos}
    return [r for r in pos if not any(key(r-t) in keys for t in pos)]
def components(bb):
    bb=S.Matrix.hstack(*bb); cc=bb.T*bb; left=set(range(bb.cols)); ans=[]
    while left:
        comp={left.pop()}; old=set()
        while old!=comp:
            old=comp.copy(); new={j for i in comp for j in left if cc[i,j]!=0}
            comp |= new; left -= new
        ans.append(sorted(comp))
    return ans
def root_type(rr):
    if not rr:return '1'
    bb=base(rr); cc=components(bb); ans=[]
    for ids in cc:
        n=len(ids); g=S.Matrix.hstack(*(bb[i] for i in ids)); d=int((g.T*g).det())
        if d==n+1: t=f'A{n}'
        elif d==4: t=f'D{n}'
        elif (n,d) in ((6,3),(7,2),(8,1)):t=f'E{n}'
        else:raise AssertionError((n,d))
        ans.append(t)
    return '+'.join(sorted(ans))
def rootspan(rr):
    if not rr:return []
    n=S.Matrix.hstack(*rr).T.nullspace()
    return [r for r in roots if all((r.T*v)[0]==0 for v in n)]
def perpendicular(rr):
    return [r for r in roots if all((r.T*k)[0]==0 for k in rr)]
def centralizer(s):return [r for r in roots if (r.T*s)[0].q==1]
def longest(rr):
    bb=base(rr); pos={key(r) for r in positive(rr)}; w=I; count=0
    while True:
        poss=[r for r in bb if key(w*r) in pos]
        if not poss:break
        w=w*reflection(poss[0]);count+=1
    assert count==len(pos)
    return w
def diag_action(f,rr):
    if not rr:return []
    bb=base(rr); keys={key(r) for r in positive(rr)}; w=f; count=0
    while True:
        bad=[r for r in bb if key(w*r) not in keys]
        if not bad:break
        w=w*reflection(bad[0]);count+=1;assert count<=len(keys)
    return [next(i for i,z in enumerate(bb) if w*r==z) for r in bb]
def smith_nonzero(m):
    d=smith_normal_form(m.applyfunc(int),domain=S.ZZ)
    return [abs(int(d[i,i])) for i in range(min(d.shape)) if d[i,i]]
def saturated_kernel(m):
    ns=m.nullspace()
    if not ns:return S.zeros(m.cols,0)
    nn=S.Matrix.hstack(*ns); den=int(S.ilcm(*[v.q for v in nn]))
    vv=[nn[:,j]*den for j in range(nn.cols)]
    for a in product(range(den),repeat=nn.cols):
        v=nn*S.Matrix(a)
        if any(v) and all(x.q==1 for x in v):vv.append(v)
    return hermite_normal_form(S.Matrix.hstack(*vv).applyfunc(int))
def path(bb,ids,n):
    g=S.Matrix.hstack(*bb); c=g.T*g
    for js in permutations(ids,n):
        if all(c[js[i],js[j]]==(-1 if abs(i-j)==1 else 0)
               for i in range(n) for j in range(i+1,n)):
            return list(js)
    raise AssertionError('No path')

out={'scope':'Integral all-q, all five-adic valuations; no finite-q samples',
     'simple_basis':B.tolist(),'gram':gram.tolist(),'centralizers':[], 'five_cycles':[]}
models={}
for family,node in [('A8',1),('E6A2',6)]:
    s=B*gram.inv()[:,node]/3; cs=centralizer(s); bb=base(cs)
    assert root_type(cs)==('A8' if family=='A8' else 'A2+E6')
    f=-longest(cs)
    assert all(x.q==1 for x in Bi*(f*s+s))
    models[family]=(s,cs,bb,f)

sa,ca,ba,fa=models['A8']
chain=path(ba,range(8),8)
se,ce,be,fe=models['E6A2']
ec=next(c for c in components(be) if len(c)==6)
ac=next(c for c in components(be) if len(c)==2)
cegram=S.Matrix.hstack(*be).T*S.Matrix.hstack(*be)
d4ids=next(list(js) for js in combinations(ec,4)
           if sorted(sum(cegram[i,j]==-1 for j in js) for i in js)==[1,1,1,3])
a5ids=None
for js in combinations(ec,5):
    kk=rootspan([be[j] for j in js])
    if root_type(kk)=='A5' and {key(fe*r) for r in kk}=={key(r) for r in kk}:
        a5ids=list(js);break
assert a5ids
settings=[
    ('A2_in_A8',sa,ca,[ba[chain[3]],ba[chain[4]]],fa,3),
    ('D4_in_E6',se,ce,[be[i] for i in d4ids],I,4),
    ('whole_E6',se,ce,[be[i] for i in ec],I,2),
    ('external_A2',se,ce,[be[i] for i in ac],fe,4),
    ('A5_in_E6',se,ce,[be[i] for i in a5ids],fe,2),
]
# The quartic source uses a different actual A2 embedding; check it literally.
sq=S.Matrix([1,1,-2,0,0,0,0,0])/3
fq=S.diag(S.eye(4),S.Matrix([[0,-1],[1,0]]),S.Matrix([[0,-1],[1,0]]))
cq=centralizer(sq)
kq=[r for r in cq if not any(r[4:])]
assert root_type(kq)=='A2'
settings.append(('external_A2_quartic',sq,cq,base(kq),fq,2))
for name,s,cs,kbasis,f,r in settings:
    kr=rootspan(kbasis); pp=perpendicular(kr); pb=base(pp)
    km=Bi*S.Matrix.hstack(*base(kr)); pm=Bi*S.Matrix.hstack(*pb)
    assert all(x.q==1 for x in km) and all(x.q==1 for x in pm)
    assert smith_nonzero(km)==[1]*km.cols
    assert smith_nonzero(pm)==[1]*pm.cols
    assert km.cols+pm.cols==8
    jr=[x for x in pp if (x.T*s)[0].q==1]; jb=base(jr)
    jm=Bi*S.Matrix.hstack(*jb)
    coordinates=pm.gauss_jordan_solve(jm)[0]
    assert all(x.q==1 for x in coordinates)
    pairing=(jm.T*gram*pm).applyfunc(int)
    z=saturated_kernel(pairing)
    joint=coordinates.row_join(z)
    index=abs(int(joint.det()))
    fpm=pm.gauss_jordan_solve(Bi*f*B*pm)[0]
    fz=z.gauss_jordan_solve(fpm*z)[0] if z.cols else S.zeros(0)
    out['centralizers'].append({
        'name':name,'s':s.tolist(),'K':root_type(kr),'K_basis':km.tolist(),
        'K_root_count':len(kr),'K_lattice_smith':smith_nonzero(km),
        'perp':root_type(pp),'perp_basis':pm.tolist(),
        'perp_lattice_smith':smith_nonzero(pm),'perp_root_count':len(pp),
        'J_derived':root_type(jr),'J_root_count':len(jr),
        'J_basis_in_perp':coordinates.tolist(),'center_basis_in_perp':z.tolist(),
        'derived_times_connected_center_index':index,
        'center_character_smith':smith_nonzero(pairing),
        'F_diagram_on_J':diag_action(f,jr),
        'F_on_connected_center':fz.tolist(),
        'F_center_charpoly':str(fz.charpoly().as_expr()),
        'F_opposite_center_charpoly':str((-fz).charpoly().as_expr()),
        'five_rank_in_relevant_source_row':r})

gap=['# Uniform E8 order-five Weyl actions, generated by geometry.py.']
root_index={key(r):i+1 for i,r in enumerate(roots)}
def gp(m):return 'PermList('+str([root_index[key(m*r)] for r in roots])+')'
gap.append('WE8:=Group(['+','.join(gp(reflection(r)) for r in simple)+']);;')
for family,(s,cs,bb,f) in models.items():
    ids=range(8) if family=='A8' else next(c for c in components(bb) if len(c)==6)
    a4=path(bb,ids,4);w=I
    for j in a4:w=w*reflection(bb[j])
    assert w**5==I and w!=I
    wi=Bi*w*B
    assert all(x.q==1 for x in wi)
    sf=[]
    for j in range(1,5):
        ss=smith_nonzero(wi**j-I);assert ss==[1]*4;sf.append(ss)
    perp4=perpendicular([bb[i] for i in a4])
    assert root_type(perp4)=='A4'
    out['five_cycles'].append({'family':family,'w':wi.tolist(),
        'charpoly':str(wi.charpoly().as_expr()),'smith_of_nonidentity_minus_one':sf,
        'fixed_module_all_n':'(Z/nZ)^4','perpendicular_A4_root_count':len(perp4),
        'P_type':'C5','fixed_A_order':'5^(4a)',
        'abelian_outside_A_bound':'5^(4a+1) < 5^(8a), for every a>=1'})
    vname='WA' if family=='A8' else 'WE'
    gap.append(vname+':=Group(['+','.join(gp(reflection(r)) for r in bb)+']);;')
    gap.append('p:='+gp(w)+';; P:=Group(p);;')
    gap.append('Print("'+family+'|Wsource|",Size('+vname+'),"|Nsource|",Size(Normalizer('+vname+',P)),"|NE8|",Size(Normalizer(WE8,P)),"|CE8|",Size(Centralizer(WE8,p)),"\\n");')
gap.append('WD4:=Group(['+','.join(gp(reflection(be[i])) for i in d4ids)+']);;')
gap.append('Print("D4_in_E6A2|relative|",Size(Normalizer(WE,WD4))/Size(WD4),"\\n");')
gap.append('WE6quartic:=Group(['+','.join(gp(reflection(r)) for r in base(perpendicular(kq)))+']);;')
gap.append('Print("E6_quartic|centralizer|",Size(Centralizer(WE6quartic,'+gp(fq)+')),"\\n");')
gap.append('QUIT;')
(HERE/'normalizers.g').write_text('\n'.join(gap)+'\n')
out['status']='passed'
(HERE/'geometry.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
print(json.dumps({'status':'passed','centralizers':[(r['name'],r['perp'],r['J_derived'],r['derived_times_connected_center_index'],r['center_character_smith']) for r in out['centralizers']], 'five_cycles':[(r['family'],r['smith_of_nonidentity_minus_one']) for r in out['five_cycles']]}))
