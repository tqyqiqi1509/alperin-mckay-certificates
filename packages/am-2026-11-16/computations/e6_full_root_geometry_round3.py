from pathlib import Path
import json, itertools, math
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

ROOT=Path(__file__).resolve().parents[1]
C=np.array([[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]],dtype=np.int16)
I=np.eye(6,dtype=np.int16)
S=[I-np.outer(I[:,i],C[i,:]) for i in range(6)]
gamma=I[[5,1,4,3,2,0]]
highest=np.array([1,2,2,3,2,1],dtype=np.int16)
component=np.column_stack([I[:,5],I[:,2],I[:,4],I[:,3],I[:,1],-highest])

def closure(gens):
    seen={I.tobytes():I}; queue=[I]
    for x in queue:
        for g in gens:
            y=g@x;k=y.tobytes()
            if k not in seen:seen[k]=y;queue.append(y)
    return np.stack(queue)

W=closure(S)
roots=sorted({tuple(w[:,0]) for w in W})
roots=np.array(roots,dtype=np.int16)

def element_point(num,den):
    num=np.array(num,dtype=np.int32)
    dif=np.einsum('nij,j->ni',W,num)-num
    reg=np.all(dif%den==0,axis=1)
    nat=np.all((dif@C.T)%den==0,axis=1)
    rt=roots[(roots@C@num)%den==0]
    return W[reg],W[nat],rt

def cp_key(m):
    return tuple(int(x) for x in sp.Matrix(m).charpoly().all_coeffs())

def kernel_basis(m,d):
    x=sp.Symbol('x');poly=sp.Poly(sp.cyclotomic_poly(d,x),x)
    a=sp.zeros(6)
    for coeff in poly.all_coeffs():a=a*sp.Matrix(m)+coeff*sp.eye(6)
    ns=a.nullspace()
    # This is a rational basis. Integral saturation is checked separately below.
    cols=[]
    for v in ns:
        q=sp.ilcm(*[z.q for z in v]); vv=v*q
        gg=math.gcd(*[abs(int(z)) for z in vv]);cols.append(vv/gg)
    return np.array(sp.Matrix.hstack(*cols)).astype(np.int64) if cols else np.zeros((6,0),dtype=np.int64)

def representative(gens,f0,d):
    # Classify by trace powers first; sympy charpolys only on distinct trace signatures.
    fs=gens@f0
    cur=np.broadcast_to(I,fs.shape).copy();signatures=[]
    tr=[]
    for k in range(1,7):
        cur=np.einsum('nij,njk->nik',cur,fs)
        tr.append(np.trace(cur,axis1=1,axis2=2))
    sig=np.stack(tr,axis=1)
    keys,inds=np.unique(sig,axis=0,return_index=True)
    candidates=[]
    for idx in inds:
        f=fs[idx];b=kernel_basis(f,d)
        candidates.append((b.shape[1],cp_key(f),int(idx),f,b))
    maxdim=max(v[0] for v in candidates)
    selected=sorted([v for v in candidates if v[0]==maxdim],key=lambda v:v[1])
    return selected[0],len(selected)

def relative_count(group,f,b):
    # Required even when the marked source torus is smaller than the full eigenspace.
    orth=sp.Matrix(b).T.nullspace()
    ncols=[]
    for v in orth:
        d=int(sp.ilcm(*[x.q for x in v]));ncols.append(np.array(v*d,dtype=np.int64).reshape(6))
    if ncols:
        n=np.stack(ncols)
        preserves=np.all(n@group@b==0,axis=(1,2))
    else:preserves=np.ones(len(group),dtype=bool)
    return int(np.sum(preserves & np.all((group@f-f@group)@b==0,axis=(1,2))))

def pointwise_count(group,b):
    return int(np.sum(np.all(group@b==b,axis=(1,2))))

def row(name,reg,nat,f0,d):
    (_,cp,idx,f,b),alternatives=representative(reg,f0,d)
    r0=roots[np.all(roots@C@b==0,axis=1)]
    den=pointwise_count(W,b);nr=relative_count(W,f,b)
    dr=pointwise_count(reg,b);dn=pointwise_count(nat,b)
    rr=relative_count(reg,f,b);rn=relative_count(nat,f,b)
    vals=roots@C@b
    gcds=sorted(set(math.gcd(*[abs(int(y)) for y in v]) for v in vals if np.any(v)))
    primitive=math.gcd(*[abs(int(sp.Matrix(b[list(ii),:]).det())) for ii in itertools.combinations(range(6),b.shape[1])]) if b.shape[1] else 1
    return dict(name=name,d=d,active_dimension=b.shape[1],active_lattice_basis=b.tolist(),basis_primitive=primitive,characteristic_polynomial=cp,source_root_count=len(r0),source_weyl_order=den,full_relative_order=nr//den,regular_relative_order=rr//dr,natural_relative_order=rn//dn,regular_pointwise=dr,natural_pointwise=dn,nonzero_root_gcds=gcds,number_maximal_cp_types=alternatives,f=f.tolist())

Cs=sp.Matrix(C);w1=Cs.inv()[:,0];w6=Cs.inv()[:,5];w4=Cs.inv()[:,3]
def rational_numden(v):
    den=int(sp.ilcm(*[x.q for x in v]));return [int(x*den) for x in v],den

cases=[];geometries={}
for name,v in [('isolated3',w4/3),('strict3',(w1+w6)/3),('strict6',(Cs.inv()[:,1]+Cs.inv()[:,2]+Cs.inv()[:,4])/6),('involution',sp.Matrix(highest)/2)]:
    num,den=rational_numden(v);reg,nat,rt=element_point(num,den)
    geometries[name]=dict(regular_W=len(reg),natural_W=len(nat),roots=len(rt),num=num,den=den)
    if name=='isolated3':
        for form,f0,ds in [('split',I,[1,2,3]),('cubic',component,[1,2,3]),('nonmatching',-gamma,[1,2,6])]:
            for d in ds:cases.append(row(name+'_'+form,reg,nat,f0,d))
    elif name in ('strict3','strict6'):
        for form,f0,ds in [('split',I,[1,2,3,4,6]),('cubic',component,[1,2,3,6])]:
            for d in ds:cases.append(row(name+'_'+form,reg,nat,f0,d))
    else:
        for d in [1,2,3,4,5,6]:cases.append(row(name,reg,nat,I,d))

z3=np.array([1,0,2,0,1,2],dtype=np.int16)
diff=(W-I)%3
quotient_kernel=np.all((diff-z3[None,:,None]*diff[:,0:1,:])%3==0,axis=(1,2))
nonmatching_basis=np.column_stack([I[:,0]+I[:,5],I[:,2]+I[:,4],I[:,1],I[:,3]])
nonmatching_zero_roots=int(np.sum(np.all((roots@C@nonmatching_basis)%3==0,axis=1)))
def prescribed_row(name,reg,nat,f,b):
    rz=roots[np.all(roots@C@b==0,axis=1)]
    a=pointwise_count(W,b); ar=pointwise_count(reg,b); an=pointwise_count(nat,b)
    return dict(name=name,source_roots=len(rz),active_basis=b.tolist(),f=f.tolist(),full_relative=relative_count(W,f,b)//a,regular_relative=relative_count(reg,f,b)//ar,natural_relative=relative_count(nat,f,b)//an,source_weyl=a)
reg3,nat3,_=element_point(*rational_numden(w4/3))
regx3,natx3,_=element_point(*rational_numden((w1+w6)/3))
bd4=np.array([[2,0],[1,1],[2,1],[2,2],[1,2],[0,2]],dtype=np.int16)
extra=[prescribed_row('isolated3_nonmatching_A1cube',reg3,nat3,S[1]@(-gamma),kernel_basis(S[1]@(-gamma),1)),prescribed_row('isolated3_nonmatching_markedD4',reg3,nat3,-gamma,kernel_basis(-gamma,1)),prescribed_row('isolated3_matching_markedA3',reg3,nat3,S[0]@S[4],np.column_stack([I[:,0],I[:,4]])),prescribed_row('strict3_markedD4',regx3,natx3,I,bd4)]
# Finite relative action of a proper D4 Levi: its full image is S3, not always cyclic.
wd4=closure([S[i] for i in [1,2,3,4]])
wd5=closure([S[i] for i in [0,1,2,3,4]])
def normalizer_image_on_roots(group,subroots):
    wanted={tuple(a) for a in subroots}
    inds=[i for i,w in enumerate(group) if {tuple(a) for a in subroots@w.T}==wanted]
    return group[inds]
rd4=roots[np.all(roots@C@bd4==0,axis=1)]
bd5=np.array(Cs.inv()[:,5]*3,dtype=np.int16).reshape(6,1)
rd5=roots[np.all(roots@C@bd5==0,axis=1)]
nd4=normalizer_image_on_roots(W,rd4);nd5=normalizer_image_on_roots(W,rd5)
autW=np.concatenate([W,W@gamma],axis=0)
ed4=normalizer_image_on_roots(autW,rd4);ed5=normalizer_image_on_roots(autW,rd5)
def restriction_image_order(group,idx):
    return len({tuple(w[:,idx].flatten()) for w in group})
levi_normalizers=dict(D4_W=len(wd4),D4_inner_normalizer=len(nd4),D4_inner_relative=len(nd4)//len(wd4),D4_extended_root_image=restriction_image_order(ed4,[1,2,3,4]),D4_extended_outer_image=restriction_image_order(ed4,[1,2,3,4])//len(wd4),D5_W=len(wd5),D5_inner_normalizer=len(nd5),D5_inner_relative=len(nd5)//len(wd5),D5_extended_root_image=restriction_image_order(ed5,[0,1,2,3,4]),D5_extended_outer_image=restriction_image_order(ed5,[0,1,2,3,4])//len(wd5))
result=dict(marked_rows=extra,levi_normalizers=levi_normalizers,W_order=len(W),root_count=len(roots),component_order3=bool(np.array_equal(component@component@component,I)),component_isometry=bool(np.array_equal(component.T@C@component,C)),geometries=geometries,rows=cases,matching_3_quotient_W_kernel=int(np.sum(quotient_kernel)),nonmatching_3_zero_root_restrictions=nonmatching_zero_roots)
Path('outputs/e6_full_root_geometry_round3.json').write_text(json.dumps(result,indent=2))
print(json.dumps(dict(W_order=len(W),roots=len(roots),geometries=geometries,rows=[{k:v for k,v in r.items() if k not in ('active_lattice_basis','f','characteristic_polynomial')} for r in cases]),indent=2))
