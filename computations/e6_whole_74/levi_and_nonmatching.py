"""Exact E6 Levi and nonmatching strict-class calculations.

Mathematical input: only the Bourbaki Cartan matrix and the two
specified rational semisimple points.  No predecessor arrays are read.
The output checks lattice identities used in the accompanying proof.
"""
from itertools import combinations, product
from collections import Counter
import json
from pathlib import Path
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

C = 2*s.eye(6)
edges = [(0,2),(2,3),(3,4),(4,5),(1,3)]
for i,j in edges:
    C[i,j] = C[j,i] = -1
gamma = s.zeros(6)
for i,j in enumerate([5,1,4,3,2,0]):
    gamma[j,i] = 1
theta = [1,2,2,3,2,1]
component = s.zeros(6)
# Affine E6 arm rotation: (0,1,6)(2,3,5), fixing node 4.
for j, target in enumerate([5,2,4,3,1,-1]):
    component[:,j] = -s.Matrix(theta) if target == -1 else s.eye(6)[:,target]
assert component**3 == s.eye(6)
assert component.T*C*component == C
assert gamma*component*gamma == component**-1
pos = [s.Matrix(a) for a in product(*(range(t+1) for t in theta))
       if (s.Matrix(a).T*C*s.Matrix(a))[0] == 2]
roots = pos + [-a for a in pos]
assert len(pos) == 36 and len(roots) == 72

def snf(A):
    D=smith_normal_form(A,domain=s.ZZ)
    return [abs(int(D[i,i])) for i in range(min(D.shape)) if D[i,i]]

def types(I):
    left=set(I); ans=[]
    while left:
        todo=[min(left)]; comp=set()
        while todo:
            x=todo.pop()
            if x in comp: continue
            comp.add(x); left.discard(x)
            todo += [y for y in I if C[x,y] == -1 and y not in comp]
        n=len(comp)
        degrees=[sum(C[x,y]==-1 for y in comp) for x in comp]
        ans.append(('D' if max(degrees,default=0)==3 else 'A')+str(n))
    return '+'.join(sorted(ans)) if ans else 'torus'

levi=[]
for r in range(6):
  for I in combinations(range(6),r):
    if not I:
        levi.append(dict(nodes=[],type='torus',centre_components=[],multiplication_kernel=[]))
        continue
    A=C[:,I]
    cent=[x for x in snf(A) if x>1]
    # The Levi centre component group has order dividing det(E6)=3.
    assert cent in [[],[3]]
    Kmat=C.extract(I,I)
    if cent:
        v=None
        for m in product(range(3),repeat=r):
            if not any(m): continue
            x=A*s.Matrix(m)/3
            if all(z.q==1 for z in x):
                v=x; break
        assert v is not None
        # Saturation is exactly Q_I+Zv because the torsion index is 3.
        Kmat=Kmat.row_join(v.extract(I,[0]))
    ker=[x for x in snf(Kmat) if x>1]
    assert s.prod(ker)*s.prod(cent)==C.extract(I,I).det()
    levi.append(dict(nodes=[i+1 for i in I],type=types(I),
                     centre_components=cent,multiplication_kernel=ker))

def finite_group(gens):
    n=gens[0].rows; e=s.eye(n)
    key=lambda M:tuple(M)
    seen={key(e):e}; todo=[e]
    while todo:
        x=todo.pop()
        for a in gens:
            y=x*a; k=key(y)
            if k not in seen: seen[k]=y; todo.append(y)
    return seen

points={
    'x3':C.inv()*(s.eye(6)[:,0]+s.eye(6)[:,5])/3,
    'x6':C.inv()*(s.eye(6)[:,1]+s.eye(6)[:,2]+s.eye(6)[:,4])/6,
}
nonmatching={}
for name,x in points.items():
    cr=[a for a in roots if (a.T*C*x)[0].q==1]
    span=s.Matrix.hstack(*cr)
    B=s.Matrix.hstack(*span.columnspace())
    closure=[a for a in roots if B.row_join(a).rank()==B.cols]
    closure_pos=[a for a in pos if B.row_join(a).rank()==B.cols]
    keys={tuple(a) for a in closure_pos}
    simple=[a for a in closure_pos if not any(tuple(a-b) in keys for b in closure_pos)]
    assert len(simple)==4 and len(closure)==24
    L=s.Matrix.hstack(*simple)
    gram=L.T*C*L
    G=L.gauss_jordan_solve(gamma*L)[0]
    gens=[]
    for i in range(4):
        Ri=s.eye(4)
        for j in range(4): Ri[i,j]-=gram[i,j]
        gens.append(Ri)
    W=finite_group(gens)
    assert len(W)==192
    assert tuple(G) not in W and tuple(-G) not in W
    assert tuple(-s.eye(4)) in W and G**2==s.eye(4)
    CC=L.gauss_jordan_solve(component*L)[0]
    assert tuple(CC) not in W and CC**3 == s.eye(4)
    assert all(y.q == 1 for y in C*(component*x-x))
    # No root is perpendicular to the whole D4 space. Hence its Weyl
    # normalizer modulo W(D4) embeds in Out(D4)=S3; the displayed 3
    # and 2 give equality and the nonmatching Frobenius fixed C2.
    assert not any(a.T*C*L == s.zeros(1,4) for a in roots)
    assert G*CC*G == CC**-1
    Zbasis=s.Matrix.hstack(*(L.T*C).nullspace())
    GZ=Zbasis.gauss_jordan_solve(gamma*Zbasis)[0]
    assert GZ.charpoly().as_expr()==s.Symbol('lambda')**2-1
    # The geometric Levi closure is the centralizer of its full central torus.
    assert {tuple(a) for a in roots if (a.T*C*Zbasis)==s.zeros(1,2)} == {tuple(a) for a in closure}
    central_part=x if name=='x3' else 4*x
    assert all((a.T*C*central_part)[0].q==1 for a in closure)
    if name=='x6':
        # The primary involution has the same A1^4 centralizer in this Levi.
        cr2=[a for a in closure if (a.T*C*(3*x))[0].q==1]
        assert {tuple(a) for a in cr2}=={tuple(a) for a in cr}
    nonmatching[name]={
        'point':[str(y) for y in x], 'connected_centralizer_roots':len(cr),
        'levi_closure_roots':len(closure),
        'levi_simple_roots':[list(map(int,a)) for a in simple],
        'levi_cartan':[list(map(int,gram.row(i))) for i in range(4)],
        'gamma_on_levi':[list(map(int,G.row(i))) for i in range(4)],
        'gamma_is_outer_on_D4':tuple(G) not in W,
        'minus_gamma_is_outer_on_D4':tuple(-G) not in W,
        'gamma_central_characteristic_polynomial':str(GZ.charpoly().as_expr()),
        'central_three_primary_part_checked':True,
        'component_rotation_matrix':[list(map(int,CC.row(i))) for i in range(4)],
        'component_rotation_is_outer_on_D4':tuple(CC) not in W,
        'component_inverted_by_nonmatching_frobenius':True,
        'geometric_levi_normalizer_quotient':'S3',
        'rational_nonmatching_levi_normalizer_quotient':'C2',
    }

grouped=Counter((x['type'],tuple(x['centre_components']),tuple(x['multiplication_kernel'])) for x in levi)
out={
    'root_count':len(roots),'proper_standard_levi_count':len(levi),
    'summary':[dict(type=t,centre_components=list(c),multiplication_kernel=list(k),standard_subsets=n)
               for (t,c,k),n in sorted(grouped.items())],
    'proper_levis':levi,'nonmatching':nonmatching,
}
path=Path(__file__).with_suffix('.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'summary':out['summary'],'nonmatching':nonmatching},indent=2))
