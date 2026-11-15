"""All-prime symbolic check of full-J/source handoffs; no field q is sampled.

Input order polynomials are the actual root-coset computations of source21.
The only possible disconnected J in these noncyclic records has component C2;
all source centre components have order 1 or 2. Thus odd valuations equal the
connected-order valuations. The finite-prime split is exhaustive: for ell>20,
only Phi_e contributes among the displayed cyclotomic factors.
"""
from fractions import Fraction
from pathlib import Path
import json

BASE=Path(__file__).resolve().parent
data=json.loads((BASE/'e8_strict21_actual_sources.json').read_text())
small_primes=(3,5,7,11,13,17,19)

def vp(n,p):
    a=0
    while n and n%p==0:
        n//=p;a+=1
    return a

def val(poly,e,p):
    # The pair (u,v) means u*a+v for a=v_ell(Phi_e(q)).
    c=Fraction(poly['coefficient'])
    b=0 if p is None else vp(c.numerator,p)-vp(c.denominator,p)
    a=0
    for key,m in poly['phi'].items():
        n=int(key)
        if n==e:a+=m
        elif p is not None and n%e==0:
            x=n//e
            while x>1 and x%p==0:x//=p
            if x==1:b+=m
    return (a,b)

records=[]; failures=[]; toral=[]; wl=[]
for f in data['forms']:
    for k,s in enumerate(f['sources']):
        assert s['M_centre_components'] in (1,2)
        e=s['e']; rid=f"{f['index']}/{f['q_residue']}/{k}"
        # Generic means every prime >20 satisfying e | (ell-1), ell∤|s|.
        for p in (*small_primes,None):
            if p is not None and ((p-1)%e or f['order']%p==0):continue
            z=val(s['M_connected_centre_order'],e,p)
            j=val(s['J_identity_order'],e,p)
            m=val(s['M_order'],e,p); d=val(s['kappa_degree'],e,p)
            central=(m[0]-d[0],m[1]-d[1])==z
            defect_excess=(m[0]-d[0]-z[0],m[1]-d[1]-z[1])
            wv=0 if p is None else vp(s['relative_order'],p)
            cv=val(f['order_polynomial'],e,p)
            item={'record':rid,'C':f['rational_C_type'],'M':s['M_rational_type'],
                  'mark':s['unipotent_mark'],'e':e,'ell':p or 'generic>20',
                  'central_defect':central,'source_defect_excess':defect_excess,
                  'J':j,'A':z,'W_valuation':wv,
                  'toral_dual':not s['K_simple_indices']}
            if central:
                if s['K_simple_indices'] and j!=z:failures.append(item)
                if not s['K_simple_indices']:
                    item['maximal_defect']=(z[0],z[1]+wv)==cv
                    toral.append(item)
                    assert item['maximal_defect'],item
            records.append(item)
        def ppoly(name):
            v=s[name]
            return [str(Fraction(v['coefficient'])),[[int(i),n] for i,n in v['phi'].items()]]
        wl.append([rid,f['order'],e,not bool(s['K_simple_indices']),s['relative_order'],
                   *[ppoly(k) for k in ('M_order','kappa_degree','J_identity_order','M_connected_centre_order')]])

unexpected=[x for x in failures if not (x['ell']==3 and x['e'] in (1,2) and
             x['M'].startswith('D4') and
             (x['C'].startswith('D8') or x['C'].startswith('D5') or x['C'].startswith('2D5')))]
if unexpected: print(json.dumps({'unexpected':unexpected},indent=2))
assert not unexpected
assert all(x['source_defect_excess']==(0,5) for x in records if not x['central_defect'])
summary={'records':len(wl),'prime_classes':len(records),
         'central_defect_classes':sum(x['central_defect'] for x in records),
         'non_toral_full_J_failures':len(failures),'full_J_failure_types':
             sorted(set((x['C'],x['M'],str(x['ell']),x['e']) for x in failures)),
         'toral_central_defect_classes':len(toral),'all_toral_maximal_defect':True,
         'all_primes_proof':'ell<=20 checked symbolically; for ell>20 only Phi_e contributes',
         'full_components':'J/J0 is trivial or C2; Z(M)/Z0(M) has order1 or2'}
(BASE/'e66_full_J_valuation.json').write_text(json.dumps({'summary':summary,'records':records},indent=2))
(BASE/'e66_full_J_wolfram_input.json').write_text(json.dumps(wl,separators=(',',':')))
print(json.dumps(summary,indent=2))
