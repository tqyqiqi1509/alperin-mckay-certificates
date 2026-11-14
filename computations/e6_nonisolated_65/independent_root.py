"""Independent E6 calculation: enumerate norm-two vectors, not imported data.

The only mathematical input is the Bourbaki Cartan matrix.  Quotients are
represented faithfully on the rational centre and rationality is imposed
on those representations, avoiding assumptions about a chosen F lift.
"""
from collections import deque
from itertools import product
from pathlib import Path
from math import gcd
import json
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

C = 2 * np.eye(6, dtype=np.int64)
for i, j in [(0,2),(2,3),(3,4),(4,5),(1,3)]:
    C[i,j] = C[j,i] = -1
I = np.eye(6, dtype=np.int64)
CS = sp.Matrix(C)
# An a priori box is obtained from the highest-root coefficients.  We
# verify Weyl stability, so the box cannot have discarded another root.
positive = []
for v in product(range(2),range(3),range(3),range(4),range(3),range(2)):
    a = np.array(v, dtype=np.int64)
    if int(a @ C @ a) == 2:
        positive.append(a)
R = np.array(positive + [-a for a in positive], dtype=np.int64)
rootset = {tuple(a) for a in R}
assert len(R) == 72
simple = [I-np.outer(I[:,i], C[i,:]) for i in range(6)]
assert all(tuple(s@r) in rootset for s in simple for r in R)

def closure(gens):
    seen = {I.tobytes(): I}
    todo = deque([I])
    while todo:
        w = todo.popleft()
        for s in gens:
            v = w@s
            key = v.tobytes()
            if key not in seen:
                seen[key] = v
                todo.append(v)
    return list(seen.values())

W = np.stack(closure(simple))
assert len(W) == 51840
def rbasis(rs):
    pos = [r for r in rs if all(a>=0 for a in r)]
    pset = {tuple(r) for r in pos}
    return [r for r in pos if not any(tuple(r-s) in pset for s in pos)]
def reflection(r):
    return I-np.outer(r,r@C)
def sm(a):
    return sp.Matrix(a)
def matpoly(p,m):
    ans=sp.zeros(6)
    for (k,),a in sp.Poly(p).terms():
        ans += a*m**k
    return ans
def intbasis(cols):
    if not cols:
        return np.zeros((6,0),dtype=np.int64)
    arr=[]
    for v in cols:
        den=sp.ilcm(*[x.q for x in v])
        arr.append([int(den*x) for x in v])
    return np.array(arr,dtype=np.int64).T
def is_in_subspace(rs,B):
    return np.all(rs@C@B==0,axis=1)
def normalizers(roots):
    if not len(roots):
        return np.ones(len(W),dtype=bool)
    rset={tuple(x) for x in roots}
    mask=np.ones(len(W),dtype=bool)
    # A basis suffices; all images must be roots of the same subsystem.
    for r in rbasis(roots):
        imgs=W@r
        mask &= np.array([tuple(a) in rset for a in imgs])
    return mask
def quotient_keys(mats,B):
    return {tuple((w@B).flat) for w in mats}

theta=max(positive,key=lambda a:sum(a))
c=np.column_stack([I[:,5],I[:,2],I[:,4],I[:,3],I[:,1],-theta])
gamma=np.column_stack([I[:,5],I[:,1],I[:,4],I[:,3],I[:,2],I[:,0]])
assert np.array_equal(c@c@c,I)
assert np.array_equal(gamma@c@gamma,c@c)
v3=CS.inv()*(sm(I[:,0])+sm(I[:,5]))/3
v6=CS.inv()*(sm(I[:,1])+sm(I[:,2])+sm(I[:,4]))/6
omega1=CS.inv()*sm(I[:,0])
assert sm(c)*v3-v3 == -omega1
assert sm(c)*v6-v6 == -omega1
x=sp.Symbol('x')
phi={d:sp.cyclotomic_poly(d,x) for d in [1,2,3,4,6]}
result={'cartan':C.tolist(),'root_count':len(R),'weyl_order':len(W),
        'cycle':c.tolist(),'points':{'x3':list(v3),'x6':list(v6)},'cases':{}}
root_index={tuple(r):i+1 for i,r in enumerate(R)}
def perm_text(w):
    return 'PermList('+str([root_index[tuple(w@r)] for r in R])+')'
def shortgens(mats):
    gens=[]
    have={I.tobytes()}
    for w in mats:
        if w.tobytes() not in have:
            gens.append(w)
            have={a.tobytes() for a in closure(gens)}
        if len(have)==len(mats):
            break
    assert len(have)==len(mats)
    return gens or [I]
gap=[]
for label,v in [('x3',v3),('x6',v6)]:
    den=int(sp.ilcm(*[z.q for z in v]))
    vnum=np.array([int(den*z) for z in v],dtype=np.int64)
    cr=np.array([r for r in R if (sm(r).T*CS*v)[0].q==1])
    cb=rbasis(cr)
    wc=closure([reflection(r) for r in cb])
    wcset={a.tobytes() for a in wc}
    connected=np.array([a.tobytes() in wcset for a in W])
    residual=(W@vnum-vnum)@C
    natural=np.all(residual % den == 0,axis=1)
    assert natural.sum()==3*len(wc)
    d4v=v if label=='x3' else 2*v
    d4r=np.array([r for r in R if (sm(r).T*CS*d4v)[0].q==1])
    d4b=rbasis(d4r)
    d4B=sm(np.column_stack(d4b))
    wd4=closure([reflection(r) for r in d4b])
    m=sp.eye(6)-2*d4B*(d4B.T*CS*d4B).inv()*d4B.T*CS
    m=np.array(m,dtype=np.int64)
    assert m.tobytes() in {a.tobytes() for a in wd4}
    fs=[('split_d1',1,I,None),('split_d2',2,m,None),
        ('cubic_d1',1,c,None),('cubic_d2',2,m@c,None)]
    if label=='x3':
        targets={4:phi[1]**2*phi[4]**2,3:phi[3]**3,6:phi[3]*phi[6]**2}
        for d in [4,3,6]:
            f=next(w if d==4 else w@c for w in wd4
              if sp.Poly(sm(w if d==4 else w@c).charpoly(x).as_expr()-targets[d],x).is_zero)
            fs.append(('split_d4' if d==4 else 'cubic_d'+str(d),d,f,None))
        fs.append(('split_d1_cuspidal',1,I,intbasis((d4B.T*CS).nullspace())))
    else:
        fs.append(('cubic_d3',3,c,None))
    result[label]={'central_roots':len(cr),'central_weyl':len(wc),
                   'natural_weyl':int(natural.sum()),'central_basis':[r.tolist() for r in cb]}
    for name,d,f,active in fs:
        if active is None:
            active=intbasis(matpoly(phi[d],sm(f)).nullspace())
        lr=R[is_in_subspace(R,active)]
        lb=rbasis(lr)
        lB=np.column_stack(lb) if len(lb) else np.zeros((6,0),dtype=np.int64)
        centre=intbasis(sm(lB.T@C).nullspace())
        nw=normalizers(lr)
        fixed=np.all(((W@f-f@W)@centre)==0,axis=(1,2))
        nf=nw & fixed
        k0=nw & connected
        kb=nw & natural
        cn=quotient_keys(W[nf&connected],centre)
        bn=quotient_keys(W[nf&natural],centre)
        an=quotient_keys(W[nf],centre)
        wl=closure([reflection(r) for r in lb])
        wlset={a.tobytes() for a in wl}
        wlmask=np.array([a.tobytes() in wlset for a in W])
        Lc=int(np.count_nonzero(wlmask&connected))
        Lb=int(np.count_nonzero(wlmask&natural))
        source=[r for r in lr if (sm(r).T*CS*v)[0].q==1]
        # Restriction ideals for each primitive cyclotomic embedding.
        zz={1:1,2:-1,3:(-1+sp.sqrt(-3))/2,4:sp.I,6:(1+sp.sqrt(-3))/2}[d]
        es=(sm(f)-zz*sp.eye(6)).nullspace() if name!='split_d1_cuspidal' else [sm(a) for a in active.T]
        restrict=[]
        for r in R:
            if np.all(r@C@active==0):
                continue
            vals=[sp.simplify((sm(r).T*CS*z)[0]) for z in es]
            rat=[]
            for z in vals:
                if d in [1,2]:
                    rat.append((z,sp.Integer(0)))
                else:
                    rat.append((sp.simplify(sp.re(z)),sp.simplify(sp.im(z)/(1 if d==4 else sp.sqrt(3)))))
            denom=int(sp.ilcm(*[a.q for ab in rat for a in ab]))
            norms=[int(denom**2*(a*a+(1 if d==4 else 3)*b*b)) for a,b in rat]
            ng=gcd(*norms)
            restrict.append((denom,abs(ng)))
            assert set(sp.factorint(denom*ng)) <= {2,3}
        # The entire centre need not be active: keep them distinct.
        centre_action=[tuple((w@active).flat) for w in W[nf]]
        assert len(set(centre_action)) == len(an)
        smith=smith_normal_form(CS*sm(lB),domain=ZZ)
        datum={'d':d,'F':f.tolist(),'charpoly':str(sp.factor(sm(f).charpoly(x).as_expr())),
               'active_basis':active.tolist(),'rank_A':int(active.shape[1]/sp.totient(d)),
               'L_root_count':len(lr),'L_simple_roots':[r.tolist() for r in lb],
               'source_root_count':len(source),'W_L_order':len(wl),
               'ambient_relative_order':len(an),'connected_relative_order':len(cn),
               'natural_relative_order':len(bn),'connected_kernel_order':Lc,
               'natural_kernel_order':Lb,'source_branches':Lb//Lc,
               'centre_smith':[abs(int(smith[i,i])) for i in range(min(smith.rows,smith.cols)) if smith[i,i]],
               'root_restriction_denominator_norm_gcd':sorted(set(restrict)),
               'faithful_on_active_space':True}
        # Actual component permutation; a two-cycle means one paired
        # finite A2(q^2), not two arbitrary unitary factors.
        pending={tuple(r) for r in lr}
        components=[]
        while pending:
            comp={pending.pop()}
            while True:
                more={r for r in pending if any(np.array(r)@C@np.array(s) for s in comp)}
                if not more:
                    break
                pending-=more
                comp|=more
            components.append(comp)
        datum['F_component_permutation']=[components.index({tuple(f@np.array(r)) for r in comp}) for comp in components]
        if label=='x6' and name in ['cubic_d1','cubic_d2']:
            # Prefer this direct, same-F primary decomposition to a
            # comparison with a non-rational literal x3 representative.
            zprimary=3*v6
            sprimary=4*v6
            assert all(a.q==1 for a in sm(lB).T*CS*zprimary)
            assert all(a.q==1 for a in CS*(sm(f)*zprimary-zprimary))
            assert all(a.q==1 for a in CS*(sm(f)*sprimary-sprimary))
            pcroots=[r for r in R if (sm(r).T*CS*sprimary)[0].q==1]
            pcgroup=closure([reflection(r) for r in rbasis(pcroots)])
            pcset={a.tobytes() for a in pcgroup}
            pcmask=np.array([a.tobytes() in pcset for a in W])
            pnum=4*vnum
            pnmask=np.all(((W@pnum-pnum)@C)%den==0,axis=1)
            pcord=len(quotient_keys(W[nf&pcmask],centre))
            pnord=len(quotient_keys(W[nf&pnmask],centre))
            assert pcord==pnord==12
            datum['same_F_primary_decomposition']={
              's3_dual':list(sprimary),'central_involution_dual':list(zprimary),
              'central_Levi_root_values':list(sm(lB).T*CS*zprimary),
              'involution_F_return':list(CS*(sm(f)*zprimary-zprimary)),
              'order3_F_return':list(CS*(sm(f)*sprimary-sprimary)),
              'order3_centralizer_roots':len(pcroots),
              'order3_connected_source_order':pcord,'order3_natural_source_order':pnord,
              'uses_no_rational_transporter':True}
            ts=[]
            for w in W[nf]:
                y=sm(w)*v6
                if all(z.q==1 for z in CS*(4*y-v3)):
                    z=y-v3
                    if all(a.q==1 for a in sm(lB).T*CS*z) and all(a.q==1 for a in 2*CS*z):
                        if all(a.q==1 for a in CS*(sm(f)*z-z)):
                            ts.append((w,y,z))
            assert ts
            w,y,z=ts[0]
            datum['quadratic_transport']={'number':len(ts),'w':w.tolist(),
              'y':list(y),'zeta_dual':list(z),'twice_adjoint_coordinates':list(2*CS*z),
              'Levi_root_values':list(sm(lB).T*CS*z),
              'F_return_adjoint_coordinates':list(CS*(sm(f)*z-z))}
        for nm,mask,kernel in [('natural',nf&natural,wlmask&natural),('connected',nf&connected,wlmask&connected)]:
            gg=shortgens(list(W[mask]))
            kk=shortgens(list(W[kernel]))
            gap += ['G:=Group(['+','.join(perm_text(a) for a in gg)+']);;',
                    'K:=Subgroup(G,['+','.join(perm_text(a) for a in kk)+']);;',
                    'Q:=FactorGroup(G,K);;',
                    'Print("'+label+'_'+name+'|'+nm+'|order=",Size(Q),"|structure=",StructureDescription(Q),"|centre=",StructureDescription(Centre(Q)),"|derived=",StructureDescription(DerivedSubgroup(Q)),"|k=",NrConjugacyClasses(Q),"\\n");']
        result['cases'][label+'_'+name]=datum
        print(label,name,{a:datum[a] for a in ['rank_A','L_root_count','source_root_count','ambient_relative_order','connected_relative_order','natural_relative_order','source_branches']},flush=True)

# Independently recover the primitive three-A2 carrier from perpendicularity.
H=result['cases']['x3_cubic_d1']['L_simple_roots']
HB=np.array(H,dtype=np.int64).T
K=R[np.all(R@C@HB==0,axis=1)]
KB=np.array(rbasis(K),dtype=np.int64).T
carrier=sm(np.column_stack([HB,KB]))
assert abs(carrier.det())==3
assert abs((CS*carrier).det())==9
result['A2_three_carrier']={'H_simple_roots':H,'K_simple_roots':KB.T.tolist(),
                           'root_lattice_index':3,'pinning_character_determinant':9}
dest=Path(__file__).with_suffix('.json')
dest.write_text(json.dumps(result,indent=2,default=str)+'\n')
Path(__file__).with_suffix('.g').write_text('\n'.join(gap)+'\nQUIT;\n')
print('Saved',dest,flush=True)
