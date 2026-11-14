"""Integral, all-field E7 Sylow actions for the good-prime principal blocks.

These matrices live on the E7 root lattice.  The simply connected and
adjoint lattices differ by index two, so their ell-primary fixed point
orders agree at ell=5,7.  No finite field-size sampling is used.
"""
import json
from pathlib import Path
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

C=2*s.eye(7)
for i,j in [(1,3),(3,4),(2,4),(4,5),(5,6),(6,7)]:
    C[i-1,j-1]=C[j-1,i-1]=-1
assert C.det()==2
ref=[]
for i in range(7):
    R=s.eye(7)
    R[i,:]-=C[i,:]
    assert R**2==s.eye(7) and R.T*C*R==C
    ref.append(R)
x=s.Symbol('x')
out={'scope':'all q with ell=5 or 7, q=+-1 mod ell, all a>=1',
     'E7_cartan':C.tolist(),'weyl_order':2**10*3**4*5*7,'cases':[]}
for ell,chain in [(5,[1,3,4,5]),(7,[1,3,4,5,6,7])]:
    w=s.eye(7)
    for i in chain:
        w=w*ref[i-1]
    assert w**ell==s.eye(7) and w!=s.eye(7)
    char=w.charpoly(x).as_expr().expand()
    assert char==((x-1)**(8-ell)*s.cyclotomic_poly(ell,x)).expand()
    smith=smith_normal_form(w-s.eye(7),domain=ZZ)
    inv=[abs(int(smith[i,i])) for i in range(7)]
    nonzero=[a for a in inv if a]
    assert nonzero==[1]*(ell-1)
    # |A^w|=ell^((8-ell)*a), while |A|=ell^(7a).
    assert 1+(8-ell)<7
    out['cases'].append({'ell':ell,'simple_reflection_chain':chain,
      'integral_matrix':w.tolist(),'characteristic_polynomial':str(char),
      'smith_w_minus_1':inv,'fixed_rank':8-ell,
      'fixed_order_formula':f'{ell}^({8-ell}*a)',
      'outside_abelian_subgroup_bound':f'{ell}^(1+{8-ell}*a)',
      'toral_abelian_subgroup_order':f'{ell}^(7*a)'})
p=Path(__file__).with_suffix('.json')
p.write_text(json.dumps(out,indent=2,default=int)+'\n')
print(json.dumps(out,default=int))
