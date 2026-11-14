"""Uniform integral source geometry for E8 order 5/ell 3 and order 6/ell 5.

No previous certificate, root list, or output array is imported.  E8 is generated
from its Dynkin graph.  All outputs concern integral lattices, not finite-q counts.
"""
from pathlib import Path
from collections import deque
from fractions import Fraction
import json, math
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

OUT=Path(__file__).resolve().parent
n=8
C=2*S.eye(n)
for i,j in [(0,2),(2,3),(1,3),(3,4),(4,5),(5,6),(6,7)]:
    C[i,j]=C[j,i]=-1
assert C.det()==1
eye=S.eye(n)
simples=[eye[:,i] for i in range(n)]
def vt(v): return tuple(int(x) for x in v)
def refl(r): return eye-r*(r.T*C)
rr={vt(r) for r in simples}|{vt(-r) for r in simples}
todo=deque(rr)
while todo:
    r=S.Matrix(todo.popleft())
    for a in simples:
        v=vt(refl(a)*r)
        if v not in rr:
            rr.add(v); todo.append(v)
roots=[S.Matrix(t) for t in sorted(rr)]
assert len(roots)==240
assert all((r.T*C*r)[0]==2 for r in roots)
highest=max(roots,key=lambda r:sum(r))
assert list(highest)==[2,3,4,6,5,4,3,2]

def components(base):
    G=base.T*C*base
    unseen=set(range(base.cols)); out=[]
    while unseen:
        a=min(unseen); comp={a}; q=[a]; unseen.remove(a)
        while q:
            b=q.pop()
            for c in list(unseen):
                if G[b,c]: comp.add(c); q.append(c); unseen.remove(c)
        # All components in this certificate are chains.
        endpoints=[b for b in comp if sum(bool(G[b,c]) for c in comp if c!=b)<=1]
        a=min(endpoints); ordered=[a]
        while len(ordered)<len(comp):
            nxt=[b for b in comp if b not in ordered and G[ordered[-1],b]]
            assert len(nxt)==1
            ordered.append(nxt[0])
        out.append(ordered)
    return sorted(out,key=lambda x:(-len(x),x))

def subsystem(m,j):
    r=[x for x in roots if int(x[j])%m==0]
    pos=[x for x in r if sum(x)>0]
    ps={vt(x) for x in pos}
    b=[x for x in pos if not any(vt(x-y) in ps for y in pos)]
    B=S.Matrix.hstack(*b)
    cc=components(B)
    # Use the component ordering solely produced by the integral graph.
    B=B[:,sum(cc,[])]
    cc=components(B)
    return r,B,cc

def matrix(M): return [[str(x) for x in row] for row in M.tolist()]
def residue_order(v):
    z=1
    for x in v: z=S.ilcm(z,S.denom(x))
    return int(z)
def source(F,K):
    eq=F-eye
    if K: eq=eq.col_join(S.Matrix.hstack(*K).T*C)
    bs=eq.nullspace()
    U=S.Matrix.hstack(*bs) if bs else S.zeros(n,0)
    rs=[r for r in roots if r.T*C*U==S.zeros(1,U.cols)]
    return U,rs
def smith(A):
    d=smith_normal_form(A,domain=ZZ)
    return [abs(int(d[i,i])) for i in range(min(d.rows,d.cols)) if d[i,i]]

def based_action(F,rs):
    pos=[x for x in rs if sum(x)>0]
    ps={vt(x) for x in pos}
    b=[x for x in pos if not any(vt(x-y) in ps for y in pos)]
    B=S.Matrix.hstack(*b)
    gram=B.T*C*B
    rho=B*gram.inv()*S.ones(B.cols,1)
    x=F*rho
    w=eye; word=[]
    while any((r.T*C*x)[0]<0 for r in b):
        i=next(i for i,r in enumerate(b) if (r.T*C*x)[0]<0)
        u=refl(b[i]); x=u*x; w=u*w; word.append(i)
        assert len(word)<=120
    pp=[next(j for j,r in enumerate(b) if w*F*a==r) for a in b]
    return {"permutation":pp,"inner_correction_word":word,
            "simple_roots":matrix(B)}

report={"construction":"E8 Dynkin graph and reflection closure; no old arrays",
        "Cartan":matrix(C),"roots":len(roots),"highest_root":list(map(int,highest))}
for m,j,ell in [(5,4,3),(6,3,5)]:
    cr,B,cc=subsystem(m,j)
    assert abs(B.det())==m
    inv=B.inv(); vv=C.inv()[:,j]/m
    discr=[]
    for a in simples:
        z=inv*a
        if residue_order(z)==m:
            discr.append(z)
    gen=discr[0]
    proj=[residue_order([gen[k] for k in comp]) for comp in cc]
    w0=S.zeros(n)
    for comp in cc:
        for a,b in zip(comp,comp[::-1]): w0[b,a]=-1
    w0=B*w0*inv
    assert all(x.q==1 for x in w0)
    f=-w0
    assert f*f==eye and f.T*C*f==C
    assert all(x.q==1 for x in f*vv+vv)
    out={"order":m,"ell":ell,"centralizing_roots":len(cr),
         "component_ranks":[len(x) for x in cc],"root_basis":matrix(B),
         "root_index":m,"gluing_generator":matrix(gen),
         "gluing_projection_orders":proj,"v":matrix(vv),"f":matrix(f),
         "source_data":[]}
    if m==5:
        assert [len(x) for x in cc]==[4,4] and proj==[5,5]
        G=S.zeros(n)
        for a,b in zip(cc[0],cc[1][::-1]): G[b,a]=1
        for a,b in zip(cc[1],cc[0]): G[b,a]=1
        G=B*G*inv
        assert all(x.q==1 for x in G)
        assert G.T*C*G==C and G**4==eye and G**2==f
        out["Gamma"]=matrix(G)
        out["rational_residues"]={}
        for k in range(4):
            out["rational_residues"][str(k)]=[q for q in range(1,5)
                if all(x.q==1 for x in q*G**k*vv-vv)]
        # Action on Y(E8)/Q and on P(C)/Y(E8), separately.
        generator=next(a for a in simples if residue_order(inv*a)==5)
        out["kernel_action"] = next(k for k in range(5)
            if all(x.q==1 for x in inv*(G*generator-k*generator)))
        out["dual_center_action"] = next(k for k in range(5)
            if all(x.q==1 for x in G*vv-k*vv))
        mids=[[B[:,a] for a in comp[1:3]] for comp in cc]
        specs=[("13",G,[]),("14",G,sum(mids,[])),
               ("15",f,[]),("16a",f,mids[0]),
               ("16b",f,mids[1]),("17",f,sum(mids,[]))]
        for row,F,K in specs:
            U,rs=source(F,K)
            # Ennola uses the -1 eigenspace of -F, exactly the same U.
            assert (-F+eye)*U==S.zeros(n,U.cols)
            item={"row":row,"center_rank":U.cols,
                "source_roots":len(rs),"K_rank":len(K),
                "center_basis":matrix(U),"K_basis":matrix(S.Matrix.hstack(*K)) if K else [],
                "based_F":based_action(F,rs),"based_Ennola_F":based_action(-F,rs)}
            if row in ["13","14"]:
                Uf,rsf=source(F**3,K)
                assert set(map(vt,rsf))==set(map(vt,rs))
                item["based_inverse_F"]=based_action(F**3,rs)
                item["based_inverse_Ennola_F"]=based_action(-F**3,rs)
                assert item["based_inverse_F"]["permutation"]==item["based_F"]["permutation"]
                assert item["based_inverse_Ennola_F"]["permutation"]==item["based_Ennola_F"]["permutation"]
            out["source_data"].append(item)
        # Weyl 3-cycles, independently constructed as reflection products.
        p=[]
        for comp in cc:
            p.append(refl(B[:,comp[0]])*refl(B[:,comp[1]]))
        assert all(z**3==eye for z in p) and p[0]*p[1]==p[1]*p[0]
        out["toral_smith"]={}
        for label,mats in [("axis1",[p[0]]),("axis2",[p[1]]),
                           ("diagonal",[p[0]*p[1]]),
                           ("antidiagonal",[p[0]*p[1]**2]),
                           ("plane",p)]:
            stack=S.Matrix.vstack(*[z-eye for z in mats])
            sd=smith(stack)
            assert all(d%3 for d in sd)
            out["toral_smith"][label]={"smith_nonzero":sd,"fixed_rank":8-len(sd)}
    else:
        assert [len(x) for x in cc]==[5,2,1] and proj==[6,3,2]
        ks=[[],[B[:,a] for a in cc[1]],[B[:,a] for a in cc[0]],
            [B[:,a] for a in cc[0]+cc[1]]]
        for row,K in zip(["2","3","4","5"],ks):
            U,rs=source(f,K)
            assert (-f+eye)*U==S.zeros(n,U.cols)
            assert all((r.T*C*(3*vv))[0].q==1 for r in rs)
            bp=based_action(f,rs); bm=based_action(-f,rs)
            assert bp["permutation"]==list(range(len(bp["permutation"])))
            assert bm["permutation"]==list(range(len(bm["permutation"])))
            out["source_data"].append({"row":row,"center_rank":U.cols,
                "source_roots":len(rs),"K_rank":len(K),
                "s2_central":True,"center_basis":matrix(U),
                "based_F":bp,"based_Ennola_F":bm})
        p=eye
        for a in cc[0][:4]: p=p*refl(B[:,a])
        assert p**5==eye
        sd=smith(p-eye)
        assert all(d%5 for d in sd)
        out["toral_smith"]={"five_cycle":{"smith_nonzero":sd,"fixed_rank":8-len(sd)}}
        out["s2_centralizing_roots"]=sum((r.T*C*(3*vv))[0].q==1 for r in roots)
    report[str(m)]=out

def partitions(n,limit=None):
    if n==0:
        yield ()
        return
    for a in range(min(n,limit if limit is not None else n),0,-1):
        for rest in partitions(n-a,a): yield (a,)+rest
def hooks(p):
    return [p[i]-j+sum(r>j for r in p)-i-1
            for i in range(len(p)) for j in range(p[i])]
report["order6_e4_partition_marks"]={
    "partitions_of_6":[{"partition":p,"hooks":hooks(p),"defect_zero":4 not in hooks(p)}
                       for p in partitions(6)],
    "cyclic_source_marks":[[p,q] for p in partitions(3) for q in partitions(2)]}
assert [p for p in partitions(6) if 4 not in hooks(p)]==[(4,1,1),(3,2,1),(3,1,1,1)]
(OUT/"certificate.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps({k:{a:v for a,v in report[k].items() if a in
    ["centralizing_roots","component_ranks","root_index","gluing_projection_orders",
     "kernel_action","dual_center_action","rational_residues","toral_smith"]}
    for k in ["5","6"]},indent=2))
for k in ["5","6"]:
    print(k,[(x["row"],x["center_rank"],x["source_roots"],x["K_rank"])
             for x in report[k]["source_data"]])
