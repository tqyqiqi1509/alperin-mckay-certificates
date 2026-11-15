from cubic_functional_exact import C,O
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ
import json
D=smith_normal_form(C,domain=ZZ)
v=[int(D[i,i]) for i in range(min(D.shape)) if D[i,i]]
r={'ad_x_smith_nonzero':v,'rank_every_odd_characteristic':len(v),'centralizer_dimension_every_odd_characteristic':28-len(v)}
assert len(v)==20 and all(abs(t) in (1,2,4,8) for t in v)
(O/'cubic_functional_smith.json').write_text(json.dumps(r,indent=2))
print(r)
