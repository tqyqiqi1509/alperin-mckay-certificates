from e8_quartic19_phase import *
from sympy import Symbol, factor, eye, zeros, diag, igcd
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
q=Symbol('q'); out={}
for row in ROWS:
 d=pickle.loads((OUT/('e8_quartic19_'+row['name']+'.pickle')).read_bytes()); E=PermutationGroup(list(d['E0'].args));B=PermutationGroup(list(d['B'].args));w=d['w']
 A=Matrix([[dot(a,b)//4 for b in SIMPLE] for a in row['base']]);snf=smith_normal_form(A,domain=ZZ)
 V=Matrix.hstack(*A.nullspace()); M=Matrix.hstack(*[Matrix(action_vector(w,b)) for b in SIMPLE]); M=BM.T.inv()*M
 W=V.gauss_jordan_solve(M*V)[0]
 bw= B.centralizer(w); baseidx=[INDEX[a] for a in row['base']]
 induced=[Permutation([baseidx.index(g(i)) for i in baseidx]) for g in B.generators]
 item={'source_equation_smith_diagonal':[int(snf[i,i]) for i in range(len(row['base']))],'center_characteristic_polynomial':str(factor(W.charpoly(q).as_expr())),'center_rational_order_polynomial':str(factor((q*W-eye(W.rows)).det())),'rational_based_weyl_centralizer_order':int(bw.order()),'based_source_graph_image_order':int(PermutationGroup(induced).order()),'based_source_graph_generators':[list(g.array_form) for g in induced]}
 print(row['name'],item,flush=True);out[row['name']]=item

# Standard hyperbolic SO8 model: no division by 3.
def eu(i,j):
 m=zeros(8);m[i,j]=1;return m
def xr(i,j,s):
 if s==-1:return eu(i,j)-eu(j+4,i+4)
 return eu(i,j+4)-eu(j,i+4)
x=xr(0,1,-1)+xr(1,2,-1)+xr(1,2,1)
eta=[4,2,0,0,-4,-2,0,0];etaM=diag(*eta)
units=[eye(8)[:,i] for i in range(8)]
chain=Matrix.hstack(units[4],-units[5],units[6]+units[2],2*units[1],2*units[0],units[2]-units[6],units[3],units[7])
rootmats=[];grades=[]
for i in range(4):
 for j in range(i+1,4):
  for s in [-1,1]:
   y=xr(i,j,s);rootmats.append(y);grades.append(eta[i]+s*eta[j]);rootmats.append(y.T);grades.append(-eta[i]-s*eta[j])
cartans=[eu(i,i)-eu(i+4,i+4) for i in range(4)]
basis=cartans+rootmats
flat=lambda m:Matrix(list(m))
g0=cartans+[m for m,a in zip(rootmats,grades) if a==0]
g2=[m for m,a in zip(rootmats,grades) if a==2]
T2=Matrix.hstack(*[flat(m) for m in g2]);bracket=Matrix.hstack(*[T2.gauss_jordan_solve(flat(y*x-x*y))[0] for y in g0])
from itertools import combinations
mins=[int(bracket[:,list(c)].det()) for c in combinations(range(len(g0)),len(g2))]
tracegram=Matrix([[(a*b).trace() for b in basis] for a in basis])
graph=eye(8);graph.col_swap(3,7)
cayley=(eye(8)+x)*(eye(8)-x).inv()
out['GGGR_integral']={'x_nilpotency_ranks':[int((x**k).rank()) for k in range(1,6)],'Jordan_chain_determinant':int(chain.det()),'Jordan_chain_action':[[int(t) for t in r] for r in (chain.inv()*x*chain).tolist()],'grade0_dimension':len(g0),'grade2_dimension':len(g2),'bracket_rank':int(bracket.rank()),'bracket_minor_gcd':int(igcd(*mins)),'bracket_has_minor_pm2':2 in mins or -2 in mins,'trace_gram_determinant':str(tracegram.det()),'all_grades_even':all(a%2==0 for a in grades),'graph_fixes_x':graph*x*graph==x,'graph_fixes_eta':graph*etaM*graph==etaM,'inverse_cayley_recovers_x':(cayley-eye(8))*(cayley+eye(8)).inv()==x,'bracket_matrix':[[int(t) for t in r] for r in bracket.tolist()]}
print(json.dumps(out['GGGR_integral']),flush=True)
(OUT/'e8_quartic19_structural.json').write_text(json.dumps(out,indent=2))
