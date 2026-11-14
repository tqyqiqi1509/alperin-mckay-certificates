"""Exact valuations from the supplied official CHEVIE export.

The input contains rational coefficients and cyclotomic multiplicities, so the
output is a proof for every valuation a, not a finite set of q experiments.
"""
from collections import Counter
from fractions import Fraction
from pathlib import Path
import ast
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
path = HERE / "e8_unipotent_degrees_reexport_observed.txt"
text = path.read_text()
assert "Error," not in text
records = []
for line in text.splitlines():
    if re.match(r"^(TYPE|HC|DEG|ENNOLA)\|", line):
        records.append(line)
    elif records:
        records[-1] += line.strip().replace("\\", "")

data = {}
for line in records:
    bits = line.split("|")
    if bits[0] == "TYPE":
        typ = bits[1].strip()
        data[typ] = {"series": [], "rows": []}
    elif bits[0] == "HC":
        data[typ]["series"].append(
            {"name": bits[1] or "torus", "indices": ast.literal_eval(bits[2])}
        )
    elif bits[0] == "ENNOLA":
        data[typ]["ennola"] = ast.literal_eval(bits[1])
    else:
        data[typ]["rows"].append({
            "index": int(bits[1]), "name": bits[2],
            "coefficient": bits[3], "cyclotomic": ast.literal_eval(bits[4]),
        })

def integer_valuation(n, prime):
    ans = 0
    while n % prime == 0:
        ans += 1
        n //= prime
    return ans

def valuation_pair(row, e, prime=3):
    coefficient = Fraction(row["coefficient"])
    constant = integer_valuation(coefficient.numerator, prime)
    constant -= integer_valuation(coefficient.denominator, prime)
    powers = row["cyclotomic"]
    j = e * prime
    while j <= len(powers):
        constant += powers[j - 1]
        j *= prime
    return [powers[e - 1], constant]

assert len(data["E6"]["rows"]) == 30
assert len(data["E8"]["rows"]) == 166
for typ, entry in data.items():
    for row in entry["rows"]:
        row["v3_e1"] = valuation_pair(row, 1)
        row["v3_e2"] = valuation_pair(row, 2)
    lookup = {r["index"]: r for r in entry["rows"]}
    for series in entry["series"]:
        series["e1_valuations"] = {
            str(k): v for k, v in Counter(
                tuple(lookup[i]["v3_e1"]) for i in series["indices"]
            ).items()
        }
    if typ == "E8":
        assert sorted(abs(i) for i in entry["ennola"]) == list(range(1, 167))
        for i, signed_j in enumerate(entry["ennola"], 1):
            assert lookup[abs(signed_j)]["v3_e2"] == lookup[i]["v3_e1"]
        for series in entry["series"]:
            series["e2_actual_indices"] = [
                abs(entry["ennola"][i - 1]) for i in series["indices"]
            ]
            series["e2_valuations"] = {
                str(k): v for k, v in Counter(
                    tuple(lookup[i]["v3_e2"])
                    for i in series["e2_actual_indices"]
                ).items()
            }
            assert series["e2_valuations"] == series["e1_valuations"]

e8d4 = next(s for s in data["E8"]["series"] if s["name"] == "D_4")
assert e8d4["e1_valuations"] == {"(4, 3)": 18, "(4, 4)": 3, "(4, 5)": 4}
e6d4 = next(s for s in data["E6"]["series"] if s["name"] == "D_4")
assert e6d4["e1_valuations"] == {"(4, 3)": 3}

result = {
    "status": "exact valuation of the supplied official CHEVIE polynomial export",
    "official_commit": "64365f30757f5374c86511b377d7d98eb15c1e02",
    "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    "method": "v3(c product Phi_j(q)^m_j)=m_e*a+v3(c)+sum m_(e*3^j)",
    "E6": data["E6"], "E8": data["E8"],
}
(HERE / "e8_unipotent_degree_valuation.json").write_text(
    json.dumps(result, indent=2) + "\n"
)
print(json.dumps({
    "E6_rows": 30, "E8_rows": 166,
    "D4_E8_e1_and_actual_e2": e8d4["e1_valuations"],
    "D4_E6_e1": e6d4["e1_valuations"],
}))
