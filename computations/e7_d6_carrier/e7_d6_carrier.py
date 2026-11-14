"""Exact E7 root realization of the nonisolated order-four D4 A1^2 T branch.

E7 is the E8 root subsystem perpendicular to e7+e8.  Roots are doubled;
the semisimple cocharacter x is quadrupled.  No CHEVIE table is imported.
"""
from itertools import combinations, product
from fractions import Fraction
from pathlib import Path
import json

roots2 = []
for i,j in combinations(range(8),2):
    for a,b in product((-2,2), repeat=2):
        r = [0]*8
        r[i],r[j] = a,b
        if r[6]+r[7] == 0:
            roots2.append(tuple(r))
for r in product((-1,1), repeat=8):
    if r.count(-1) % 2 == 0 and r[6]+r[7] == 0:
        roots2.append(r)

def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

def weight_lattice_quarters(v4):
    return all(dot(v4,r) % 8 == 0 for r in roots2)

def e7_root_lattice_quarters(v4):
    # E8 root lattice intersected with e7+e8=0.
    if v4[6]+v4[7] != 0:
        return False
    residue = {a % 4 for a in v4}
    if residue not in ({0},{2}):
        return False
    return sum(v4) % 8 == 0

x4 = (2,2,2,2,0,0,1,-1)
Croots2 = [r for r in roots2 if dot(x4,r) % 8 == 0]
Leviroots2 = [r for r in roots2 if r[6] == r[7] == 0]
assert len(roots2) == 126
assert len(Croots2) == 28
assert len(Leviroots2) == 60
assert all(r[6] == r[7] == 0 for r in Croots2)
assert sum(all(r[i]==0 for i in (4,5,6,7)) for r in Croots2) == 24
assert sum(all(r[i]==0 for i in (0,1,2,3,6,7)) for r in Croots2) == 4

# w: reverse one coordinate in each D support and invert the central torus.
# It is a product of two D6 root reflections and the A1 reflection.
def w(v):
    out = list(v)
    out[0] *= -1
    out[4] *= -1
    out[6],out[7] = out[7],out[6]
    return tuple(out)

assert {w(r) for r in roots2} == set(roots2)
assert w(w(x4)) == x4
return4 = tuple(a-b for a,b in zip(w(x4),x4))
assert weight_lattice_quarters(return4)
assert not e7_root_lattice_quarters(return4)
assert not weight_lattice_quarters(x4)
assert not weight_lattice_quarters(tuple(2*a for a in x4))
assert weight_lattice_quarters(tuple(4*a for a in x4))

# Natural Z(E7) is the half-spin D6 centre, not the vector centre h0.
spin4 = (2,2,2,2,2,2,0,0)
vector4 = (4,0,0,0,0,0,0,0)
assert weight_lattice_quarters(spin4)
assert not e7_root_lattice_quarters(spin4)
assert e7_root_lattice_quarters(tuple(2*a for a in spin4))
assert not weight_lattice_quarters(vector4)

def torus_reverse(v):
    out=list(v)
    out[6],out[7]=out[7],out[6]
    return tuple(out)

rational_forms=[]
for q in (1,3):
    for component in (0,1):
        fx=tuple(q*a for a in x4)
        if q==3:
            fx=torus_reverse(fx)
        if component:
            fx=w(fx)
        difference=tuple(a-b for a,b in zip(fx,x4))
        assert weight_lattice_quarters(difference)
        spin_class=not e7_root_lattice_quarters(difference)
        assert spin_class == bool(component)
        torus_sign=(-1)**((q==3)+component)
        rational_forms.append({
            "q_mod_4":q,
            "component_twist":component,
            "torus_factor":"Phi1" if torus_sign==1 else "Phi2",
            "derived_centralizer":"D4(q) A1(q)^2" if component==0 else "2D4(q) A1(q^2)",
            "minimal_Levi_derived_form":"split D6(q)",
            "SC_lift_Frobenius_return":"nontrivial E7 centre" if spin_class else "identity",
            "table4_geometric_structure_rows": [15,16] if q==1 and component==0 else
                [17,18] if q==3 and component==0 else [20] if q==1 else [19],
        })

out = {
    "root_model": "E8 roots perpendicular to e7+e8",
    "x_times_four": x4,
    "adjoint_E7_semisimple_order": 4,
    "E7_root_count": len(roots2),
    "connected_centralizer_root_count": len(Croots2),
    "connected_centralizer_type": "D4 D2 T1 = D4 A1^2 T1",
    "minimal_Levi_root_count": len(Leviroots2),
    "minimal_Levi_type": "D6 T1",
    "component_representative": "flip e1 and e5, swap e7/e8",
    "component_return_times_four": return4,
    "component_action_on_Levi_central_torus": "inversion",
    "D6_adjoint_projection_times_two": [1,1,1,1,0,0],
    "D6_adjoint_projection_order": 2,
    "E7_central_involution_times_two": [1,1,1,1,1,1,0,0],
    "E7_centre_in_D6": "half-spin, while vector h0 is not E7 central",
    "rational_forms":rational_forms,
    "rational_row_caution":"The q mod 4 assignment of rows 19 and 20 here is opposite the sentence preceding KM13 Table 4; independently confirmed by CHEVIE in e7_quasiisolated_probe.log. This is recorded as an internal correction requiring external author verification, not an externally accepted erratum.",
    "scope": "exact algebraic root-data identities; rational extension operators retain actual Tits/Frobenius returns",
}
Path(__file__).with_name("e7_d6_carrier.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
