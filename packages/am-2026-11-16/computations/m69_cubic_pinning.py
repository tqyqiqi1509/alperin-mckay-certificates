"""Check the two actual pinning coefficients omitted from the older cubic JSON.

This imports the already separately certified integral Chevalley engine. It
does not recount its Weyl group or use a finite-field sample.
"""
import json
from pathlib import Path
from e8_tits_engine_root17 import INDEX, ROOTS, SIMPLE, longest

base = [(2,-2,0,0,0,0,0,0), (0,2,-2,0,0,0,0,0)]
n0, _ = longest(SIMPLE)
n_a2, _ = longest(base)
a0 = n0 * n_a2**-1
images = []
for root in base:
    target = a0(INDEX[root])
    images.append({"input": root, "image": ROOTS[target % 240],
                   "coefficient": (-1)**(target // 240)})
assert [x["image"] for x in images] == base[::-1]
assert [x["coefficient"] for x in images] == [1,1]
assert a0.order() == 2
out = {"model": "integral E8 root vectors", "a0_order": 2,
       "simple_pinning_images": images, "finite_field_samples": 0}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out))
