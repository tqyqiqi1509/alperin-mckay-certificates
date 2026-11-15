"""Independent E7 integral source and full-centralizer checks.

The computations are root-datum identities uniform in q, not sampled groups.
Coordinates are doubled Euclidean coordinates, plus exact coroot coordinates.
"""
from itertools import combinations, product
from collections import deque
from pathlib import Path
import json, time
import sympy as S
from sympy.combinatorics import Permutation, PermutationGroup
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ

e=S.eye(8)
simple=[S.Matrix([1,-1,-1,-1,-1,-1,-1,1]),2*(e[:,0]+e[:,1]),2*(e[:,1]-e[:,0])]
simple += [2*(e[:,i]-e[:,i-1]) for i in range(2,6)]
B=S.Matrix.hstack(*simple)
C=B.T*B/4
roots=set()
for i,j in combinations(range(6),2):
    for a,b in product([-2,2], repeat=2):
        r=[0]*8;r[i]=a;r[j]=b;roots.add(tuple(r))
roots.update({(0,0,0,0,0,0,-2,2),(0,0,0,0,0,0,2,-2)})
for v in product([-1,1],repeat=6):
    if sum(t==-1 for t in v)%2==1:
        for a in [-1,1]:roots.add(tuple(v)+(a,-a))
roots=sorted(roots);index={r:i for i,r in enumerate(roots)}
coords=[tuple(B.gauss_jordan_solve(S.Matrix(r))[0]) for r in roots]
coordidx={r:i for i,r in enumerate(coords)}
assert len(roots)==126 and C.det()==2 and all(all(x.q==1 for x in v) for v in coords)

def matrix_reflect(a):return S.eye(7)-S.Matrix(a)*(S.Matrix(a).T*C)
def perm_of(M):return Permutation([coordidx[tuple(M*S.Matrix(a))] for a in coords])
simplemat=[matrix_reflect(tuple(S.eye(7)[:,i])) for i in range(7)]
W=PermutationGroup([perm_of(m) for m in simplemat])
assert W.order()==2903040
def subsystem_base(R):
    R=set(R);P={r for r in R if next(x for x in r if x)!=-abs(next(x for x in r if x))}
    return sorted(r for r in P if not any(tuple(r[i]-s[i] for i in range(7)) in P for s in P))
def components(base):
    adj={i:{j for j in range(len(base)) if i!=j and (S.Matrix(base[i]).T*C*S.Matrix(base[j]))[0]!=0} for i in range(len(base))}
    unseen=set(adj);out=[]
    while unseen:
        q=[min(unseen)];seen=set(q)
        for i in q:
            for j in adj[i]:
                if j not in seen:seen.add(j);q.append(j)
        unseen-=seen;out.append([base[i] for i in sorted(seen)])
    return out
def chain_order(base):
    adj={i:[j for j in range(len(base)) if i!=j and (S.Matrix(base[i]).T*C*S.Matrix(base[j]))[0]!=0] for i in range(len(base))}
    end=next(i for i in adj if len(adj[i])==1);q=[end]
    while len(q)<len(base):q.append(next(j for j in adj[q[-1]] if j not in q))
    return [base[i] for i in q]
def centralizer_roots_of_space(V):
    return [r for r in coords if S.Matrix(r).T*C*V==S.zeros(1,V.cols)]
def perp_roots(R):return centralizer_roots_of_space(S.Matrix.hstack(*map(S.Matrix,R)))
def root_smith(R):
    M=S.Matrix.hstack(*map(S.Matrix,subsystem_base(R)))
    D=smith_normal_form(C*M,domain=ZZ)
    return [abs(int(D[i,i])) for i in range(min(D.shape))]
def pointwise_weyl_order(R):
    U=W
    for r in subsystem_base(R):U=U.stabilizer(coordidx[r])
    return int(U.order())
simpleperms=[[g(i) for i in range(126)] for g in W.generators]
def normalizer_order(R):
    """Orbit the embedded root subsystem, not all elements of W(E7)."""
    start=tuple(sorted(coordidx[r] for r in R));seen={start};todo=[start]
    for current in todo:
        for g in simpleperms:
            nxt=tuple(sorted(g[i] for i in current))
            if nxt not in seen:seen.add(nxt);todo.append(nxt)
    return int(W.order())//len(seen),len(seen)
def row(R):
    p=perp_roots(R)
    nw,orb=normalizer_order(R)
    wh=PermutationGroup([perm_of(matrix_reflect(a)) for a in subsystem_base(R)])
    pb=subsystem_base(p)
    contents={}
    for a in coords:
        vals=[int((S.Matrix(a).T*C*S.Matrix(b))[0]) for b in pb]
        nonzero=[abs(x) for x in vals if x]
        k=int(S.igcd(*nonzero)) if len(nonzero)>1 else (nonzero[0] if nonzero else 0)
        contents[str(k)]=contents.get(str(k),0)+1
    return {"root_count":len(R),"component_ranks":[len(c) for c in components(subsystem_base(R))],"root_smith":root_smith(R),"perpendicular_roots":len(p),"perpendicular_component_ranks":[len(c) for c in components(pb)],"pointwise_W_order":pointwise_weyl_order(R),"normalizer_order":nw,"subsystem_orbit_size":orb,"relative_normalizer_order":nw//int(wh.order()),"base":[list(r) for r in subsystem_base(R)],"root_restriction_contents_on_perp_coroots":contents}
out={"roots":126,"W_order":int(W.order())}
for name,I in [("D4",[1,2,3,4]),("special_A1cube",[1,4,6]),("A1square",[1,4])]:
    R=[r for r in coords if all(r[i]==0 for i in range(7) if i not in I)]
    out[name]=row(R)

# Isolated A5 A2 at omega3/3; build its actual two-packet Frobenius.
R=[r for r in coords if r[2]%3==0]
parts=components(subsystem_base(R));assert sorted(map(len,parts))==[2,5]
A5=chain_order(next(c for c in parts if len(c)==5));A2=chain_order(next(c for c in parts if len(c)==2))
v2=S.eye(7)
for a in [A5[0],A5[2],A5[4],A2[0]]:v2=v2*matrix_reflect(a)
V=S.Matrix.hstack(*(v2+S.eye(7)).nullspace())
L=centralizer_roots_of_space(V)
out["A5A2_two_packet_A1cube"]=row(L)
out["A5A2_two_packet_A1cube"].update({"active_rank":V.cols,"frob_charpoly":str(v2.charpoly().as_expr().factor()),"source_frob_on_roots":[list(v2*S.Matrix(a)) for a in subsystem_base(L)]})

# D4-source bad-three colour mu=(0,0,0,0,0,2,-1,1) in unscaled Euclidean coordinates.
mu=S.Matrix([0,0,0,0,0,4,-2,2])
hist={}
for r in roots:
    a=(S.Matrix(r).T*mu)[0]/4;hist[str(a)]=hist.get(str(a),0)+1
out["D4_equal_colour_root_pairing_histogram"]=hist

# Full central overlap of the D4 source and its actual A1^3 centralizer.
perp_euclidean=[2*(e[:,5]-e[:,4]),2*(e[:,5]+e[:,4]),2*(e[:,7]-e[:,6])]
perp_c=[B.gauss_jordan_solve(v)[0] for v in perp_euclidean]
central_points=[]
for bits in product([0,1],repeat=3):
    v=sum((bit*b for bit,b in zip(bits,perp_c)),S.zeros(7,1))/2
    in_H=all(v[i].q==1 for i in [0,5,6])
    in_ZG=all(x.q==1 for x in C*v)
    central_points.append({"bits":bits,"in_D4":in_H,"in_Z_E7":in_ZG})
out["D4_perpendicular_central_overlap"]={"coroots":[list(v) for v in perp_c],"points":central_points,"all_ones_is_zG":all(x.q==1 for x in sum(perp_c,S.zeros(7,1))/2-C.inv()[:,6])}

def primitive_integer(v):
    z=v*S.ilcm(*[x.q for x in v]);d=S.igcd(*[abs(int(x)) for x in z if x]);return z/d
def geometric_semisimple_row(v):
    R=[a for a in coords if (S.Matrix(a).T*C*v)[0].q==1]
    bas=subsystem_base(R);M=S.Matrix.hstack(*map(S.Matrix,bas))
    center=S.Matrix.hstack(*(M.T*C).nullspace()) if M.rank()<7 else S.zeros(7,0)
    Levi=centralizer_roots_of_space(center)
    lbase=subsystem_base(Levi);N=S.Matrix.hstack(*map(S.Matrix,lbase))
    proj=(N.T*C*N).inv()*N.T*C*v
    diag=smith_normal_form(N,domain=ZZ)
    ans={"point_coroot_coordinates":list(v),"centralizer_root_count":len(R),"centralizer_component_ranks":[len(x) for x in components(bas)],"connected_center_dimension":center.cols,"minimal_connected_Levi_root_count":len(Levi),"minimal_Levi_component_ranks":[len(x) for x in components(lbase)],"derived_projection_in_Levi_coroots":list(proj),"projection_order_in_natural_derived_torus":int(S.ilcm(*[z.q for z in proj])),"derived_coroot_embedding_Smith":[abs(int(diag[i,i])) for i in range(min(diag.shape))]}
    if center.cols==1:
        z=primitive_integer(center[:,0]);ans.update({"primitive_center_coroot":list(z),"zG_in_connected_center":all(a.q==1 for a in z/2-C.inv()[:,6])})
    return ans
fw=[C.inv()[:,i] for i in range(7)]
xD=B.gauss_jordan_solve(S.Matrix([1,1,1,1,0,0,S.Rational(1,2),S.Rational(-1,2)]))[0]
out["nonisolated_connected_QI_projections"]={
    "order2_E6T1":geometric_semisimple_row(fw[6]/2),
    "order6_A2cubeT1":geometric_semisimple_row((fw[2]+fw[4])/6),
    "order4_D4A1squareT1":geometric_semisimple_row(xD)}

# The actual D6 inherited central colour, including four rational returns.
# All torus coordinates here are unscaled Euclidean coordinates.
E=e
d6base=[E[:,i]-E[:,i-1] for i in range(5,0,-1)]+[E[:,0]+E[:,1]]
DB=S.Matrix.hstack(*d6base)
h0=E[:,0]
zs=sum((E[:,i] for i in range(6)),S.zeros(8,1))/2
zc=zs+h0
b=E[:,7]-E[:,6]
omega7=B*fw[6]/2
alpha1=B[:,0]/2
def in_QD6(z):
    return z[6]==0 and z[7]==0 and all(x.q==1 for x in z[:6]) and sum(z[:6])%2==0
def qE7(z):return B.gauss_jordan_solve(2*z)[0]
assert in_QD6(omega7-alpha1-zs)
assert in_QD6(b/2-alpha1-zc)
assert all(in_QD6(2*z) and not in_QD6(z) for z in [h0,zs,zc])
x=S.Matrix([S.Rational(1,2)]*4+[0,0,S.Rational(1,4),S.Rational(-1,4)])
x0=S.Matrix([S.Rational(1,2)]*4+[0,0,0,0])
u=S.eye(8);u[6,6]=u[7,7]=0;u[6,7]=u[7,6]=1
w=u.copy();w[0,0]=w[4,4]=-1
assert all(tuple(2*w*S.Matrix(r)/2) in index for r in roots)
returns=[]
for qmod4,delta in product([1,3],[0,1]):
    Fmat=(w if delta else S.eye(8))*(u if qmod4==3 else S.eye(8))
    rE=qmod4*Fmat*x-x
    rD=qmod4*Fmat*x0-x0
    assert all(a.q==1 for a in qE7(rE-delta*omega7))
    assert in_QD6(rD-delta*h0)
    # Increasing q by four changes both return classes by the root lattice.
    assert all(a.q==1 for a in qE7(4*Fmat*x))
    assert in_QD6(4*Fmat*x0)
    returns.append({"q_mod4":qmod4,"component":delta,
       "torus_sign":int((b.T*Fmat*b)[0]/2),
       "E7_return_coordinates_representative":list(rE),
       "D6_return_coordinates_representative":list(rD),
       "E7_return_mod_Q": "zG" if delta else "1",
       "D6_return_mod_Q": "h0" if delta else "1",
       "D6_central_values_h0_zs_zc":[1,-1,-1] if delta else [1,1,1]})
out["D6_inherited_colour_and_returns"]={
 "simple_roots":[list(z) for z in d6base],
 "h0":list(h0),"zs":list(zs),"zc_representative":list(zc),
 "omega7_minus_alpha1":list(omega7-alpha1),
 "Levi_central_glue": "(zc,-1) of order 2; zG=zs is the other half-spin element",
 "centre_pairings_mod_Z":[[str((a.T*c)[0]%1) for c in [h0,zs,zc]] for a in [h0,zs,zc]],
 "four_rational_returns":returns}

def exact_json(x):
    if isinstance(x,S.Integer):return int(x)
    if isinstance(x,S.Rational):return str(x)
    raise TypeError(type(x).__name__)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2,default=exact_json))
print(json.dumps(out,default=exact_json))
