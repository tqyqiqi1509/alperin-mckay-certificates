"""Exact root-lattice certificate for the entire E8 D4-source 3-block family.

No field size is substituted. Smith forms control all 3-adic valuations.
The matrices and mathematical identification are explained in the companion note.
"""
from collections import deque, Counter
from pathlib import Path
import hashlib
import json
from itertools import product
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

OUT = Path(__file__).resolve().parent

Q = sp.Matrix([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]])
u = sp.Matrix([[-1,0,1,0],[-1,0,1,-1],[-1,0,0,0],[-1,1,0,-1]])
v = sp.Matrix([[1,0,0,0],[1,0,1,-1],[0,0,1,0],[0,1,0,-1]])
t0 = sp.Matrix([0,1,0,2])

def key(M):
    return tuple(int(x) for x in M)

def mul(a,b,n=4):
    return tuple(sum(a[n*i+k]*b[n*k+j] for k in range(n))
                 for i in range(n) for j in range(n))

def mv(a,b,n=4):
    return tuple(sum(a[n*i+j]*b[j] for j in range(n)) for i in range(n))

def closure(gens,n=4):
    one=tuple(int(i==j) for i in range(n) for j in range(n))
    seen={one}; todo=deque([one])
    while todo:
        a=todo.popleft()
        for g in gens:
            z=mul(a,g,n)
            if z not in seen:
                seen.add(z);todo.append(z)
    return sorted(seen)

def snf(M):
    d=smith_normal_form(M,domain=ZZ)
    return [abs(int(d[i,i])) for i in range(min(d.rows,d.cols))]

def roots(C):
    n=C.rows
    start=[tuple(int(i==j) for i in range(n)) for j in range(n)]
    seen=set(start);todo=deque(start)
    while todo:
        a=todo.popleft()
        for i in range(n):
            b=list(a);b[i]-=sum(int(C[i,j])*a[j] for j in range(n));b=tuple(b)
            if b not in seen:seen.add(b);todo.append(b)
    return sorted(seen)

assert u**3 == v**3 == sp.eye(4)
assert u*v == v*u
assert u.T*Q*u == v.T*Q*v == Q
gens=[]
for i in range(4):
    s=sp.eye(4)
    for j in range(4):s[i,j]-=Q[i,j]
    gens.append(key(s))
for perm in ([2,1,0,3],[2,1,3,0]):
    a=sp.zeros(4)
    for j in range(4):a[perm[j],j]=1
    assert a.T*Q*a == Q
    gens.append(key(a))
W=closure(gens)
assert len(W)==1152 and key(u) in W and key(v) in W
S=[u**i*v**j for i,j in product(range(3),repeat=2)]
forms=[snf(w-sp.eye(4)) for w in S if w!=sp.eye(4)]
assert Counter(map(tuple,forms))==Counter({(1,1,3,3):4,(1,1,0,0):4})
stack=(u-sp.eye(4)).col_join(v-sp.eye(4))
assert snf(stack)==[1,1,1,3]
orbit={tuple(x%3 for x in mv(w,tuple(t0))) for w in W}
Wt=[w for w in W if tuple(x%3 for x in mv(w,tuple(t0)))==tuple(t0)]
assert len(orbit)==32 and len(Wt)==36
reflections=[w for w in Wt if (sp.Matrix(4,4,w)-sp.eye(4)).rank()==1
             and mul(w,w)==key(sp.eye(4))]
assert len(reflections)==6
todo=set(reflections);reflection_components=[]
while todo:
    start=todo.pop();component={start};queue=[start]
    while queue:
        a=queue.pop()
        neighbours=[b for b in todo if mul(a,b)!=mul(b,a)]
        for b in neighbours:todo.remove(b);component.add(b);queue.append(b)
    reflection_components.append(component)
assert sorted(map(len,reflection_components))==[3,3]
assert [len(closure(list(x))) for x in reflection_components]==[6,6]
assert len(closure(reflections))==36
Sset={key(x) for x in S}
assert len(Sset)==9
normalizer=[]
for w in W:
    inverse=key(sp.Matrix(4,4,w).inv())
    if {mul(mul(w,s),inverse) for s in Sset}==Sset:normalizer.append(w)
assert len(normalizer)==72

# Bourbaki E8: edges 1-3-4-5-6-7-8 and 2-4.
C=2*sp.eye(8)
for i,j in [(0,2),(2,3),(3,4),(4,5),(5,6),(6,7),(1,3)]:
    C[i,j]=C[j,i]=-1
assert C.det()==1
R=roots(C);assert len(R)==240
X=sp.Matrix([[0,0,0,2],[0,0,1,2],[0,0,1,3],[0,0,2,4],
             [0,0,2,3],[0,0,2,2],[1,0,1,1],[0,1,0,0]])
assert X.T*C*X==Q
restricted=[tuple(int(x) for x in (sp.Matrix(1,8,r)*C*X)) for r in R]
annihilator=[r for r,z in zip(R,restricted) if all(x%3==0 for x in z)]
assert len(annihilator)==24
central_roots=[r for r in R if int((sp.Matrix(1,8,r)*C*X*t0)[0])%3==0]
assert len(central_roots)==78
# Connected components using nonorthogonality, with rank and root count.
todo=set(central_roots);components=[]
while todo:
    start=todo.pop();component={start};queue=[start]
    while queue:
        a=queue.pop()
        neighbours=[b for b in todo if (sp.Matrix(1,8,a)*C*sp.Matrix(b))[0]!=0]
        for b in neighbours:todo.remove(b);component.add(b);queue.append(b)
    components.append({"roots":len(component),"rank":sp.Matrix(list(component)).rank()})
assert sorted((c["roots"],c["rank"]) for c in components)==[(6,2),(72,6)]

# Opposition on the whole central torus and identity on the source D4 space.
H=sp.eye(8)[:,[1,2,3,4]]
B=H.row_join(X)
m=B*sp.diag(*([1]*4+[-1]*4))*B.inv()
assert all(x.q==1 for x in m) and m.T*C*m==C and m*m==sp.eye(8)
assert m*H==H and m*X==-X
positive=[r for r in central_roots if all(x>=0 for x in r)]
pset=set(positive)
simple_c=[r for r in positive if not any(
    tuple(r[i]-a[i] for i in range(8)) in pset for a in positive)]
assert len(simple_c)==8
mc=m;steps=0
while True:
    inverse=mc.inv()
    bad=[r for r in simple_c if all(x<=0 for x in inverse*sp.Matrix(r))]
    if not bad:break
    a=sp.Matrix(bad[0])
    mc=(sp.eye(8)-a*(a.T*C))*mc
    steps+=1
    assert steps<=len(positive)
perm=[simple_c.index(tuple(int(x) for x in mc*sp.Matrix(r))) for r in simple_c]
assert sorted(perm)==list(range(8))
simple_graph=[]
todo=set(range(8))
while todo:
    i=todo.pop();cc={i};queue=[i]
    while queue:
        i=queue.pop()
        ns=[j for j in todo if (sp.Matrix(1,8,simple_c[i])*C*sp.Matrix(simple_c[j]))[0]]
        for j in ns:todo.remove(j);cc.add(j);queue.append(j)
    assert {perm[j] for j in cc}==cc and any(perm[j]!=j for j in cc)
    simple_graph.append({"rank":len(cc),"indices":sorted(cc),
                         "permutation":[perm[j] for j in sorted(cc)]})
assert sorted(c["rank"] for c in simple_graph)==[2,6]

R4=roots(Q);assert len(R4)==24
idx={r:i+1 for i,r in enumerate(R4)}
perms=[[idx[mv(g,r)] for r in R4] for g in gens]
def gaplist(x):return "["+",".join(map(str,x))+"]"
gap="W:=Group("+",".join("PermList("+gaplist(p)+")" for p in perms)+");;\n"
gap+="ct:=CharacterTable(W);; dd:=List(Irr(ct),chi->chi[1]);;\n"
gap+='Print("order=",Size(W)," degrees=",dd," three_prime_count=",Number(dd,d->d mod 3<>0),"\\n");\nQUIT;\n'
(OUT/"e8_d4_source_weyl.g").write_text(gap)
record={
 "scope":"All q prime to 3 and all a=v3(q-epsilon)>=1; exact integral lattices",
 "D4_gram":list(map(list,Q.tolist())),"u":u.tolist(),"v":v.tolist(),
 "relative_weyl_order":len(W),"sylow_three_order":len(Sset),
 "nonidentity_smith_forms":forms,"joint_smith_form":snf(stack),
 "fixed_line_generator_mod3":list(t0),"nonzero_colour_orbit":len(orbit),
 "colour_stabilizer_order":len(Wt),"sylow_normalizer_order":len(normalizer),
 "colour_stabilizer_structure":"S3 x S3, certified by commuting rank-two reflection factors",
 "E8_root_count":len(R),"source_centralizer_root_count":len(annihilator),
 "nonzero_colour_centralizer_components":components,
 "opposition_integral_matrix":m.tolist(),
 "opposition_on_centralizer_components":simple_graph,"chamber_correction_steps":steps,
 "proof_of_characteristic":"For abelian J outside A, E=JA/A has order3 or9; |J|<=|E||C_A(E)|<3^(4a).",
}
record=json.loads(json.dumps(record,default=int))
record["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(OUT/"e8_d4_source_structure.json").write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({k:record[k] for k in ["relative_weyl_order","sylow_three_order","joint_smith_form",
    "nonzero_colour_orbit","colour_stabilizer_order","sylow_normalizer_order",
    "E8_root_count","source_centralizer_root_count","nonzero_colour_centralizer_components"]}))
