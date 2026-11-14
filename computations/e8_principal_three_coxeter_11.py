"""Universal integral certificate for the E8 principal 3-block argument.
No finite-q character values or sampled prime powers are used.
"""
import json
from pathlib import Path
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
n=8
C=2*s.eye(n)
for a,b in [(1,3),(3,4),(4,5),(5,6),(6,7),(7,8),(2,4)]:
    C[a-1,b-1]=C[b-1,a-1]=-1
reflections=[]
for i in range(n):
    R=s.eye(n)
    for j in range(n):R[i,j]-=C[i,j]
    assert R.T*C*R==C and R*R==s.eye(n)
    reflections.append(R)
c=s.eye(n)
for R in reflections:c=c*R
x=s.symbols("x")
assert s.expand(c.charpoly(x).as_expr())==s.cyclotomic_poly(30,x)
assert c**30==s.eye(n)
z=c**10
assert z**3==s.eye(n) and z!=s.eye(n)
assert s.eye(n)+z+z**2==s.zeros(n)
assert s.expand(z.charpoly(x).as_expr())==s.expand(s.cyclotomic_poly(3,x)**4)
snf=smith_normal_form(s.eye(n)-z,domain=s.ZZ)
diagonal=[abs(int(snf[i,i])) for i in range(n)]
assert diagonal==[1]*4+[3]*4
z5=c**6
assert z5**5==s.eye(n) and z5!=s.eye(n)
assert sum((z5**j for j in range(5)),s.zeros(n))==s.zeros(n)
assert s.expand(z5.charpoly(x).as_expr())==s.expand(s.cyclotomic_poly(5,x)**2)
snf5=smith_normal_form(s.eye(n)-z5,domain=s.ZZ)
diagonal5=[abs(int(snf5[i,i])) for i in range(n)]
assert diagonal5==[1]*6+[5]*2
result={
 "purpose":"E8(q), 3 not dividing q, all a: elliptic order-three Weyl element",
 "basis":"Bourbaki simple-root basis, edges 1-3-4-5-6-7-8 and 2-4",
 "cartan":[list(map(int,C.row(i))) for i in range(n)],
 "coxeter_matrix":[list(map(int,c.row(i))) for i in range(n)],
 "coxeter_characteristic_polynomial":str(c.charpoly(x).as_expr()),
 "order_three_matrix":[list(map(int,z.row(i))) for i in range(n)],
 "order_three_charpoly":str(s.factor(z.charpoly(x).as_expr())),
 "norm_1_z_z2_zero":True,
 "smith_1_minus_z":diagonal,
 "torus_fixed_character_group":"(Z/3Z)^4",
 "order_five_supplement":{
  "matrix":[list(map(int,z5.row(i))) for i in range(n)],
  "characteristic_polynomial":str(s.factor(z5.charpoly(x).as_expr())),
  "norm_sum_zero":True,
  "smith_1_minus_z":diagonal5,
  "torus_fixed_character_group":"(Z/5Z)^2",
  "scope":"Same actual rational Tits-lift argument for principal 5-blocks at e=1,2; no claim about e=4 here."
 },
 "logical_use":"A Sylow 3-subgroup S of W(E8) containing z has T^S <= T^z, a 3-group. This is used with the actual rational Tits lift and C_G(T_3)=T, not as a stand-alone block proof.",
}
out=Path(__file__).with_suffix(".json")
out.write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({k:result[k] for k in ["order_three_charpoly","norm_1_z_z2_zero","smith_1_minus_z","order_five_supplement"]}))
