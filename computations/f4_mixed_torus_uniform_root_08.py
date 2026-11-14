"""Integral F4 root checks for all 3-adic valuations, not field samples.

For a mixed torus use w of order 3, rank(w-I)=2.  Put t=epsilon*q,
t=1+3^a*u, 3 not dividing u.  A root is trivial on the full Sylow-3
torus iff every coordinate of (t^2*w^2+t*w+I)*root is divisible by
3^(a+1). Its reduction is N*root+3^a*u*K*root. The finitely many
small valuations depend only on u mod 3; above them use an explicit
bound. Pointwise Weyl stabilizers are checked by the same exact criterion.
"""
from itertools import combinations, product
from pathlib import Path
import json
import sympy as s

E=s.eye(4)
basis=s.Matrix.hstack(E[:,1]-E[:,2],E[:,2]-E[:,3],E[:,3],
                     (E[:,0]-E[:,1]-E[:,2]-E[:,3])/2)
inv=basis.inv()
roots=[]
for i in range(4):
    for sign in (-1,1):
        roots.append(inv*(sign*E[:,i]))
for i,j in combinations(range(4),2):
    for a,b in product((-1,1),repeat=2):
        roots.append(inv*(a*E[:,i]+b*E[:,j]))
for signs in product((-1,1),repeat=4):
    roots.append(inv*s.Matrix(signs)/2)
assert len(set(tuple(v) for v in roots))==48
assert all(x.q==1 for v in roots for x in v)
reflections=[]
for i in range(4):
    alpha=basis[:,i]
    reflections.append(inv*(E-2*alpha*alpha.T/(alpha.dot(alpha)))*basis)

def key(M):
    return tuple(map(int,M))
weyl={key(E):E}
todo=[E]
while todo:
    M=todo.pop()
    for ref in reflections:
        V=ref*M
        k=key(V)
        if k not in weyl:
            weyl[k]=V
            todo.append(V)
assert len(weyl)==1152

report={"basis_columns":[[str(x) for x in row] for row in basis.tolist()],
        "root_count":48,"weyl_order":1152,"cases":[]}
for name,w in [('long_A2',reflections[0]*reflections[1]),
               ('short_A2',reflections[2]*reflections[3])]:
    assert w**3==E and (w-E).rank()==2
    N=E+w+w**2
    K=2*w**2+w
    perpendicular=[r for r in roots if w*r==r]
    pair=next((b,c) for b in perpendicular for c in perpendicular
              if (basis*b).dot(basis*b)==(basis*c).dot(basis*c)
              and 2*(basis*b).dot(basis*c)==-(basis*b).dot(basis*b))
    def root_reflection(r):
        b=basis*r
        return inv*(E-2*b*b.T/b.dot(b))*basis
    v=root_reflection(pair[0])*root_reflection(pair[1])
    assert v**3==E and v*w==w*v and (v-E).rank()==2
    delta=(v-E)**2
    bound=max(abs(int(x)) for v in roots for x in N*v)
    wbound=max(abs(int(x)) for h in weyl.values() for x in N*(h-E))
    first_large=1
    while 3**first_large<=max(bound,wbound):
        first_large+=1
    small=[]
    for a,u in product(range(1,first_large),(1,2)):
        P=N+3**a*u*K
        modulus=3**(a+1)
        bad=[list(map(int,v)) for v in roots
             if all(int(x)%modulus==0 for x in P*v)]
        pointwise=[h.tolist() for h in weyl.values()
                   if all(int(x)%modulus==0 for x in P*(h-E))]
        small.append({"a":a,"u_mod_3":u,"trivial_roots":bad,
                      "pointwise_weyl":pointwise,
                      "D_commutator_noncentral":any(int(x)%modulus for x in P*delta)})
    zero=[v for v in roots if N*v==s.zeros(4,1)]
    large_bad=[list(map(int,v)) for v in zero
               if all(int(x)%3==0 for x in K*v)]
    large_pointwise=[h.tolist() for h in weyl.values()
                     if N*(h-E)==s.zeros(4) and
                     all(int(x)%3==0 for x in K*(h-E))]
    # For a>=first_large, N*v divisible by 3^a iff N*v=0.
    report['cases'].append({"name":name,"w":w.tolist(),
        "N":N.tolist(),"K":K.tolist(),"small_valuations":small,
        "max_abs_N_root":bound,"N_zero_roots":len(zero),
        "max_abs_N_weyl_difference":wbound,
        "perpendicular_A2_roots":[list(map(int,r)) for r in pair],
        "defect_weyl_action_v":v.tolist(),
        "N_delta_nonzero":N*delta!=s.zeros(4),
        "max_abs_N_delta":max(abs(int(x)) for x in N*delta),
        "first_large_a":first_large,"large_a_trivial_roots":large_bad,
        "large_a_pointwise_weyl":large_pointwise})
out=Path(__file__).with_suffix('.json')
out.write_text(json.dumps(report,default=int,indent=2)+'\n')
print(json.dumps(report,default=int,indent=2))
