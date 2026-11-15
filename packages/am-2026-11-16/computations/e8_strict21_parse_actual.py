"""Parse exact generic CHEVIE source records and verify E8 lattice data.

The parameters are root data, alcove vertices and rational Weyl cosets; no
finite field instances are used.  Every noncyclic candidate source is retained
with its actual ambient matrix and every ordinary unipotent mark.
"""
import ast
from collections import Counter, defaultdict
from fractions import Fraction
import json
from math import prod
from pathlib import Path
import re
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

HERE=Path(__file__).parent


def records(path):
    text=path.read_text()
    # CHEVIE emits a separate informational line while building some
    # rational relative-group character tables.  It is retained in raw output
    # and copied to the JSON audit field, but is not part of a source record.
    text=re.sub(r'(?m)^warning: powermap not bound\s*$', '', text)
    prefixes=("AM21_ACTUAL_ROOT_DATA","ISOLATED_COUNT","ISOLATED","FORM_BEGIN",
              "FORM_ORDER","ACTUAL_E","ACTUAL_SOURCE","FORM_END","DONE_ACTUAL_STRICT21")
    pattern=re.compile(r"(?:^|\n|gap> )("+"|".join(prefixes)+r")(?=\||\n)")
    matches=list(pattern.finditer(text))
    out=[]
    for i,m in enumerate(matches):
        end=matches[i+1].start() if i+1<len(matches) else len(text)
        value=text[m.start(1):end].strip()
        value=re.sub(r"\bgap>\s*$","",value).strip()
        out.append(value.split('|'))
    assert out[-1][0]=='DONE_ACTUAL_STRICT21'
    return out


def cyc(text):
    text=text.strip()
    phi={int(i):int(v) if v else 1 for i,v in re.findall(r'P(\d+)(?:\^(-?\d+))?',text)}
    pref=re.split(r'q|P',text)[0].strip()
    coeff=Fraction(pref or '1')
    return {'display':text,'coefficient':str(coeff),'phi':phi}


def phi_valuation(poly,ell,e):
    """Return slope*a+constant, a=v_ell(Phi_e(q)), for odd ell and e|ell-1."""
    coeff=Fraction(poly['coefficient'])
    def v(n):
        a=0
        while n and n%ell==0:n//=ell;a+=1
        return a
    const=v(abs(coeff.numerator))-v(coeff.denominator)
    slope=poly['phi'].get(e,0)
    n=e*ell
    while n<=max(poly['phi'],default=0):
        const+=poly['phi'].get(n,0);n*=ell
    return slope,const


def mat_eval(poly,m):
    ans=sp.zeros(m.rows)
    for a in sp.Poly(poly).all_coeffs():ans=ans*m+a*sp.eye(m.rows)
    return ans


def source_geometry(row,roots,cartan):
    fm=sp.Matrix(row['F_matrix'])
    assert fm*cartan*fm.T==cartan
    k=sp.Matrix([roots[i-1] for i in row['K_simple_indices']])
    if not row['K_simple_indices']:k=sp.zeros(0,8)
    snf=smith_normal_form(k,domain=sp.ZZ)
    factors=[abs(int(snf[i,i])) for i in range(min(snf.shape)) if snf[i,i]]
    assert k.rank()==len(row['K_simple_indices'])
    # Roots are row vectors.  The F-stable cocharacter space can be
    # represented with the invariant bilinear form because E8 is self-dual.
    p=mat_eval(sp.cyclotomic_poly(row['e']),fm)
    conditions=k*cartan
    vmat=p.T.col_join(conditions)
    space=vmat.nullspace()
    assert len(space)==row['active_phi']*sp.totient(row['e'])
    vv=sp.Matrix.hstack(*space)
    ambient_roots=[r for r in roots if sp.Matrix([r])*cartan*vv==sp.zeros(1,vv.cols)]
    def stype(rr):
        if not rr:return '1'
        rem=set(tuple(x) for x in rr);components=[]
        while rem:
            a=rem.pop(); todo=[a];co=[a]
            while todo:
                a=todo.pop()
                linked=[b for b in rem if (sp.Matrix([a])*cartan*sp.Matrix(b))[0]]
                for b in linked:rem.remove(b);todo.append(b);co.append(b)
            rank=sp.Matrix(co).rank();count=len(co)
            if count==rank*(rank+1):typ=f'A{rank}'
            elif count==2*rank*(rank-1):typ=f'D{rank}'
            elif (rank,count) in ((6,72),(7,126),(8,240)):typ=f'E{rank}'
            else:raise AssertionError((rank,count))
            components.append(typ)
        return '*'.join(sorted(components))
    # Whole source centre components are precisely the torsion in X/Q_K.
    return {'K_root_lattice_smith':factors,'M_centre_components':prod(factors),
            'ambient_Levi_root_count':len(ambient_roots),
            'ambient_Levi_type':stype(ambient_roots),
            'active_torus_dimension':len(space)}


def main():
    name='e8_strict21_actual_sources_with_defect.out'
    path=HERE/name
    if not path.exists() or 'DONE_ACTUAL_STRICT21' not in path.read_text():
        name='e8_strict21_actual_sources.out';path=HERE/name
    blocks=records(path)
    data={'input':name,'forms':[],
          'chevie_warnings':re.findall(r'(?m)^warning:.*$',path.read_text())}; current=None;iso={}
    for rec in blocks:
        tag=rec[0]
        if tag=='AM21_ACTUAL_ROOT_DATA':
            data['simple_roots']=ast.literal_eval(rec[1])
            data['simple_coroots']=ast.literal_eval(rec[2])
            data['roots']=ast.literal_eval(rec[3])
        elif tag=='ISOLATED_COUNT':assert int(rec[1])==9
        elif tag=='ISOLATED':
            iso[int(rec[1])]={'index':int(rec[1]),'order':int(rec[2]),
                 'semisimple_point':rec[3].strip(),'C_type':rec[4],
                 'C_simple_indices':ast.literal_eval(rec[5])}
        elif tag=='FORM_BEGIN':
            current={**iso[int(rec[1])],'q_residue':int(rec[2]),'rational_C_type':rec[3],
                     'F_matrix':ast.literal_eval(rec[4]),'sources':[],'cuspidal_pair_counts':{}}
            data['forms'].append(current)
        elif tag=='FORM_ORDER':current['order_polynomial']=cyc(rec[1])
        elif tag=='ACTUAL_E':current['cuspidal_pair_counts'][int(rec[1])]=int(rec[2])
        elif tag=='ACTUAL_SOURCE':
            row={'e':int(rec[1]),'active_phi':int(rec[2]),'M_rational_type':rec[3],
                 'unipotent_mark':rec[4],'relative_order':int(rec[5]),
                 'K_simple_indices':ast.literal_eval(rec[6]),'F_matrix':ast.literal_eval(rec[7]),
                 'M_connected_centre_order':cyc(rec[8]),'J_augmented_rational_type':rec[9],
                 'J_identity_order':cyc(rec[10]),'global_unipotent_indices':ast.literal_eval(rec[11])}
            if len(rec)>12:
                row['geometric_normalizer_quotient_order']=int(rec[12])
                row['kappa_degree']=cyc(rec[13]);row['M_order']=cyc(rec[14])
            current['sources'].append(row)
    roots=data['roots'];cartan=sp.Matrix(data['simple_coroots'])
    geometric_cache={};count=0;bad=[]
    for f in data['forms']:
        for row in f['sources']:
            key=(tuple(row['K_simple_indices']),str(row['F_matrix']),row['e'])
            if key not in geometric_cache:
                geometric_cache[key]=source_geometry(row,roots,cartan)
            row.update(geometric_cache[key]);count+=1
            row['prime_specific_checks']=[]
            for ell in (3,5,7,11,13,17,19,23,29,31):
                if f['order']%ell==0 or (ell-1)%row['e']:continue
                test={'ell':ell,'J_valuation':phi_valuation(row['J_identity_order'],ell,row['e']),
                      'M_centre_valuation':phi_valuation(row['M_connected_centre_order'],ell,row['e'])}
                if 'kappa_degree' in row:
                    mo=phi_valuation(row['M_order'],ell,row['e'])
                    kd=phi_valuation(row['kappa_degree'],ell,row['e'])
                    test['source_defect']=(mo[0]-kd[0],mo[1]-kd[1])
                    test['source_central_defect']=test['source_defect']==(row['active_phi'],0)
                row['prime_specific_checks'].append(test)
                if test['J_valuation']!=(row['active_phi'],0):
                    bad.append({'C':f['rational_C_type'],'q_residue':f['q_residue'],
                                'e':row['e'],'mark':row['unipotent_mark'],'r':row['active_phi'],**test})
    data['summary']={'form_count':len(data['forms']),'noncyclic_marked_source_count':count,
                     'distinct_geometry_count':len(geometric_cache),
                     'relevant_e':sorted({e for f in data['forms'] for e in f['cuspidal_pair_counts']}),
                     'M_component_orders':sorted({r['M_centre_components'] for f in data['forms'] for r in f['sources']}),
                     'non_sylow_J_cases':bad}
    target=HERE/'e8_strict21_actual_sources.json'
    target.write_text(json.dumps(data,indent=2,ensure_ascii=False,default=int))
    print(json.dumps(data['summary'],indent=2,ensure_ascii=False,default=int))


if __name__=='__main__':main()
