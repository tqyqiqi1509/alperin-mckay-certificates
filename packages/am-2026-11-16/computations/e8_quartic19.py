"""Independent E8 order-four source carriers. Exact root lattices and Tits lifts."""
from e8_tits_engine_root17 import *
import pickle, sys

OUT=Path(__file__).parent
I240=Permutation(list(range(240)))
def shadow(p):return Permutation([p(i)%240 for i in range(240)])
NS=[lift(a) for a in SIMPLE]; SS=[shadow(n) for n in NS]
GLOBALPOS={INDEX[r] for r in ROOTS if all(x>=0 for x in COORD[r])}

def vector2(*entries):return tuple(2*x for x in entries)
def root(i,j,si=-1,sj=1):
    v=[0]*8;v[i-1]=2*si;v[j-1]=2*sj;return tuple(v)

ROWS=[
 {'name':'D5A3_regular','kind':'D','base':[root(1,6),root(1,6,1,1)]},
 {'name':'A7A1_regular','kind':'A','base':[(-1,-1,1,1,-1,1,-1,1),(-1,1,-1,-1,1,1,-1,1),(1,1,1,1,1,-1,-1,1)]},
 {'name':'D5A3_A1q2','kind':'D','base':[root(5,7),root(5,7,1,1)],'signed':(3,4,-1,-2,5,-8,-7,6)},
 {'name':'D5A3_internal_A3','kind':'D','base':[root(3,4),root(3,4,1,1),root(4,5),root(5,7)],'signed':(2,-1,-3,4,5,-8,-7,6)},
 {'name':'A7A1_external_A1','kind':'A','base':[root(1,3),root(4,6,1,1),root(5,7,1,1),(1,-1,-1,-1,-1,-1,-1,1)],'signed':(-2,3,-8,5,6,7,4,1)}]

def rowsolve(A,b,mod):
    M=[[int(x)%mod for x in row]+[int(bb)%mod] for row,bb in zip(A,b)]
    nr=len(M);nc=len(A[0]);piv=[];r=0
    for c in range(nc):
        k=next((k for k in range(r,nr) if M[k][c]%2==1),None)
        if k is None:continue
        M[r],M[k]=M[k],M[r];z=pow(M[r][c],-1,mod);M[r]=[(x*z)%mod for x in M[r]]
        for k in range(nr):
            if k!=r:
                z=M[k][c];M[k]=[(x-z*y)%mod for x,y in zip(M[k],M[r])]
        piv.append(c);r+=1
        if r==nr:break
    assert r==nr,('nonprimitive input',A,mod)
    v=[0]*nc
    for i,c in enumerate(piv):v[c]=M[i][-1]
    kernel=[]
    for c in range(nc):
        if c not in piv:
            x=[0]*nc;x[c]=1
            for i,p in enumerate(piv):x[p]=-M[i][c]%mod
            kernel.append(x)
    return v,kernel

def coeff_y(v):return tuple(sum(v[j]*SIMPLE[j][i] for j in range(8)) for i in range(8))

def lift_weyl(p):
    target=p;out=ID;steps=0
    while p!=I240:
        i=next(i for i,a in enumerate(SIMPLE) if p(INDEX[a]) not in GLOBALPOS)
        p=SS[i]*p;out=out*NS[i];steps+=1
        assert steps<=120
    assert shadow(out)==target
    return out

def source_data(base):
    gram=Matrix([[dot(a,b)//4 for b in base] for a in base]);gi=gram.inv()
    data={}
    for r in ROOTS:
        co=gi*Matrix([dot(a,r)//4 for a in base])
        if all(sum(co[j]*base[j][i] for j in range(len(base)))==r[i] for i in range(8)):
            assert all(x>=0 for x in co) or all(x<=0 for x in co),(base,r,co)
            data[INDEX[r]]=tuple(co)
    pos={i for i,c in data.items() if all(x>=0 for x in c)}
    return set(data),pos

def based_projection(p,base,pos):
    refl=[shadow(lift(a)) for a in base];original=p;steps=0
    while True:
        inv=p**-1
        j=next((j for j,a in enumerate(base) if inv(INDEX[a]) not in pos),None)
        if j is None:break
        p=p*refl[j];steps+=1
        assert steps<=120
    assert {p(INDEX[a]) for a in base}=={INDEX[a] for a in base}
    return p

def pinned_lift(p,base):
    n=lift_weyl(p)
    image=[ROOTS[p(INDEX[a])] for a in base]
    A=[[dot(a,b)//4 for b in SIMPLE] for a in image]
    signs=[n(INDEX[a])//240 for a in base]
    v,k=rowsolve(A,signs,2)
    n=n*torus2(coeff_y(v))
    assert all(n(INDEX[a])==p(INDEX[a]) for a in base)
    return n

def signed_weyl(signed):
    def action(v):
        out=[0]*8
        for i,z in enumerate(signed):out[abs(z)-1]=v[i]*(1 if z>0 else -1)
        return tuple(out)
    return Permutation([INDEX[action(r)] for r in ROOTS])

def regular_weyl(base):
    return Permutation([INDEX[tuple(-r[i]+sum(dot(r,a)//4*a[i] for a in base) for i in range(8))] for r in ROOTS])

def get_C(kind):
    cr=[r for r in ROOTS if (sum(r[i] for i in [5,6,7])%4==0 if kind=='D' else dot(r,(-1,1,1,1,1,1,1,7))%16==0)]
    cb=positive_base(cr);W=PermutationGroup([shadow(lift(a)) for a in cb])
    return cr,cb,W

def build(row):
    start=time.time();base=row['base'];lroots,lpos=source_data(base)
    cr,cb,WC=get_C(row['kind']);print(row['name'],'C order',WC.order(),flush=True)
    based=[];normalizers=0
    for p in WC.generate_schreier_sims():
        if all(p(INDEX[a]) in lroots for a in base):
            normalizers+=1;based.append(based_projection(p,base,lpos))
    uniq={tuple(p.array_form):p for p in based};gens=[];B=PermutationGroup([I240])
    for p in uniq.values():
        if not B.contains(p):
            gens.append(p);B=PermutationGroup(gens)
    A=[[dot(a,b)//4 for b in SIMPLE] for a in base]
    _,ker2=rowsolve(A,[0]*len(base),2);_,ker4=rowsolve(A,[0]*len(base),4)
    egens=[pinned_lift(p,base) for p in gens]+[torus2(coeff_y(v)) for v in ker2]
    E0=PermutationGroup(egens)
    w=signed_weyl(row['signed']) if 'signed' in row else regular_weyl(base)
    assert all(w(INDEX[a]) in lroots for a in base)
    w=based_projection(w,base,lpos)
    assert B.contains(w)
    n=pinned_lift(w,base)
    out={'row':row['name'],'source_root_count':len(lroots),'C_root_count':len(cr),'C_weyl_order':int(WC.order()),'normalizer_in_C_order':normalizers,'based_image_order':len(uniq),'based_generator_count':len(gens),'E0_order':int(E0.order()),'center2_order':2**len(ker2),'w_order':int(w.order()),'n_order':int(n.order()),'n_central_E0':all(n*g==g*n for g in egens),'seconds':time.time()-start}
    path=OUT/('e8_quartic19_'+row['name'])
    path.with_suffix('.pickle').write_bytes(pickle.dumps({'row':row,'out':out,'base':base,'B':B,'E0':E0,'n':n,'w':w,'ker2':ker2,'ker4':ker4}))
    path.with_suffix('.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
    return out

if __name__=='__main__':
    for row in ROWS:
        if len(sys.argv)==1 or row['name'] in sys.argv[1:]:build(row)
