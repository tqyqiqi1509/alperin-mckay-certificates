"""Symbolic valuations from a fresh execution of the author's CHEVIE data.
No field value is substituted.  Pairs (r,c) mean r*a+c for every a>=1.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import ast, re, json, hashlib

here=Path(__file__).resolve().parent
p=here/'e8_nonprincipal22_generic.stdout.txt'
txt=p.read_text();assert 'Error,' not in txt and 'Syntax error' not in txt and 'DONE22' in txt
# Runtime diagnostics checked independently; GAP3 final quit warning follows DONE22.
rows={};sers=[];hc=[]
for match in re.finditer(r'^DEG\|(\d+)\|([^|]*)\|([^|]*)\|(\[[^\]]*\])',txt,re.M):
 i,n,c,ds=match.groups();rows[int(i)]={'name':n,'coefficient':c,'factors':ast.literal_eval(ds)}
assert len(rows)==166
for m in re.finditer(r'^SER\|(\d+)\|([^|]*)\|([^|]*)\|(\d+)\|(\[[^\]]*\])',txt,re.M):
 e,l,n,w,ids=m.groups();sers.append({'e':int(e),'levi':l,'mark':n,'W':int(w),'indices':ast.literal_eval(ids)})
for m in re.finditer(r'^HC\|([^|]*)\|(\[[^\]]*\])',txt,re.M):
 n,ids=m.groups();hc.append({'name':n,'indices':ast.literal_eval(ids)})
def vp(n,p):
 n=abs(n);assert n
 k=0
 while n%p==0:k+=1;n//=p
 return k
def val(i,p,e):
 row=rows[i];co=Fraction(row['coefficient']);ds=row['factors'];c=vp(co.numerator,p)-vp(co.denominator,p);j=e*p
 while j<=30:c+=ds[j-1];j*=p
 return (ds[e-1],c)
def vals(ids,p,e):return {str(k):n for k,n in sorted(Counter(val(i,p,e) for i in ids).items())}
for s in sers:
 e=s['e'];s['valuation_5']=vals(s['indices'],5,e) if e in [1,2,4] else None
 s['valuation_3']=vals(s['indices'],3,e) if e in [1,2] else None
for h in hc:h['valuation_3_e1']=vals(h['indices'],3,1);h['valuation_5_e1']=vals(h['indices'],5,1)
e3D4=[s for s in sers if s['e'] in [1,2] and 'D4<' in s['levi']]
assert all(s['valuation_3']=={'(4, 3)':18,'(4, 4)':3,'(4, 5)':4} for s in e3D4)
assert all(s['valuation_5']=={'(4, 2)':25} for s in e3D4)
e5E6=[s for s in sers if s['e'] in [1,2] and 'E6.' in s['levi']]
assert all(s['valuation_5']=={'(6, 2)':6} for s in e5E6)
quartic=[s for s in sers if s['e']==4];assert len(quartic)==4
assert all(s['valuation_5']=={'(2, 1)':16} for s in quartic)
cusp4=[i for i in rows if rows[i]['factors'][3]==4]
exceptional5=[i for i in cusp4 if val(i,5,4)==(4,-1)]
assert len(cusp4)==43 and len(exceptional5)==5
assert vals(cusp4,5,4)=={'(4, -1)':5,'(4, 1)':38}
out={'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'commit':'64365f30757f5374c86511b377d7d98eb15c1e02','method':'m_e*a + v_l(c) + sum_j>=1 m_(e*l^j)', 'ordinary_hc':hc,'non_toral_positive_series':sers,'quartic_cuspidal_values':vals(cusp4,5,4),'quartic_exceptional_five':[(i,rows[i]['name']) for i in exceptional5],'rows':rows}
(here/'e8_nonprincipal22_valuations.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ['rows','non_toral_positive_series']},indent=2))
