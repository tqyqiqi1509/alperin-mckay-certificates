"""Run the exact computations in an isolated, disposable output directory."""
from pathlib import Path
import argparse,sys,subprocess,json,shutil,hashlib,time,datetime,re,concurrent.futures
ROOT=Path(__file__).resolve().parent
if sys.flags.optimize:
    raise SystemExit('Run without -O/PYTHONOPTIMIZE: exact assertion checks must be enabled.')
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--suite',choices=['python','gap3','gap4','all'],default='python')
ap.add_argument('--output',type=Path)
ap.add_argument('--gap3-command',default='["gap3","-b"]',help='JSON array: GAP3 executable and runtime arguments')
ap.add_argument('--gap4-command',default='["gap","-q","-b","--quitonbreak"]',help='JSON array: GAP4 executable and runtime arguments')
a=ap.parse_args();target=(a.output or ROOT/'runs'/datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')).resolve()
if target.exists():raise SystemExit('Choose a new --output directory; existing results are preserved.')
target.mkdir(parents=True);shutil.copytree(ROOT/'computations',target/'computations');(target/'outputs').mkdir();(target/'logs').mkdir()
for p in (ROOT/'data/chevie').iterdir():shutil.copy2(p,target/'computations'/p.name)
g3=json.loads(a.gap3_command);g4=json.loads(a.gap4_command);receipts=[]
def run(rel,kind='python',alias=None,marker=None):
    p=target/'computations'/rel;name=rel.replace('/','__');start=time.monotonic()
    cmd=[sys.executable,str(p)] if kind=='python' else (g3 if kind=='gap3' else g4)
    if kind=='gap4':cmd=cmd+[str(p)]
    r=subprocess.run(cmd,input=p.read_bytes() if kind=='gap3' else None,cwd=target,capture_output=True,timeout=1800)
    (target/'logs'/(name+'.stdout')).write_bytes(r.stdout);(target/'logs'/(name+'.stderr')).write_bytes(r.stderr)
    text=r.stdout.decode('utf-8',errors='replace');err=r.stderr.decode('utf-8',errors='replace')
    ok=r.returncode==0 and not re.search(r'(^|\n)(?:Syntax error|Error[, :])',text+'\n'+err)
    if marker:ok=ok and marker in text
    row={'program':rel,'runtime':kind,'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'exit_code':r.returncode,'seconds':round(time.monotonic()-start,2),'completion_marker':marker,'passed':bool(ok),'stdout_sha256':hashlib.sha256(r.stdout).hexdigest()}
    receipts.append(row);print(json.dumps(row),flush=True)
    if not ok:raise RuntimeError('Computation failed: '+rel+'; inspect retained logs.')
    if alias:(target/'computations'/alias).write_bytes(r.stdout)
    return text
def python_suite():
    names=['odd_d_round4_check.py','e6_carrier_round3_checks.py','e6_D4_carrier_round3_checks.py','e6_full_root_geometry_round3.py','e6_cubic_index_round3_checks.py','e7_uniform_carrier_11.py','e7_global20_rootdata.py','e8_tits_engine_root17.py','m69_cubic_pinning.py','e67_source_regular_phases.py','cubic_functional_exact.py','cubic_functional_smith.py','e7_whole_76/levi_lattices.py','e8_d4_source_structure.py','e8_quartic19.py','e8_quartic19_phase.py','e8_quartic19_liftscan.py','e8_quartic19_structural.py','e8_quartic19_internal_source.py','e8_quartic19_characters.py','e8_quartic19_field.py','quartic_centre_overlap.py','r4_d4_identity.py','e8_strict21_parse_actual.py','e8_strict21_check_tails.py','e66_full_J_valuation.py','e8_nonprincipal22_valuations.py']
    for n in names:run(n)
    for n in ['e8_principal_three_coxeter_11.py','e8_three_five_75/geometry.py']:run(n)
def gap3_suite():
    jobs=[('e8_strict21_actual_sources.g','e8_strict21_actual_sources.out','DONE_ACTUAL_STRICT21'),('e8_strict21_cyclic_tails.g','e8_strict21_cyclic_tails.out','DONE_TAILS21'),('e8_nonprincipal22_generic.g','e8_nonprincipal22_generic.stdout.txt','DONE22'),('e8_strict21_exceptional_catalog.g','e8_strict21_exceptional_catalog.out','DONE_STRICT21'),('r4_d4_springer.g',None,'DONE_R4_D4')]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda x:run(x[0],'gap3',x[1],x[2]),jobs))
    for n in ['e8_strict21_parse_actual.py','e8_strict21_check_tails.py','e66_full_J_valuation.py','e8_nonprincipal22_valuations.py']:run(n)
def gap4_suite():
    completion={'c60_independent_block_and_kernel.g':'C60_INDEPENDENT_BLOCK_KERNEL_COMPLETE','round2_colours.g':'C60_PER_COLOUR_ANCHORS_COMPLETE','round2_smallgrp.g':'ROUND2_SMALLGRP_COMPLETE','round2_seven.g':'ROUND2_SEVEN_AFFINE_COMPLETE'}
    for n in ['c60_independent_block_and_kernel.g','round2_colours.g','round2_smallgrp.g','round2_seven.g']:
        run(n,'gap4',marker=completion[n])
    for p in sorted((target/'computations/exceptional_covers').glob('*.g')):
        if p.name!='actual_c2_pairing.g':
            markers=re.findall(r'FRESH_[A-Z0-9_]+_COMPLETE',p.read_text())
            run('exceptional_covers/'+p.name,'gap4',marker=markers[-1] if markers else None)
    run('e7_56/independent_coefficients.py')
    text=run('e7_56/independent_coefficients.g','gap4')
    assert re.search(r'Subgroup classes=22\s+pairs=47\s+ordinary entries=260',text)
    assert re.search(r'Failed=\[\s*\]',text)
    run('e8_phi4_sylow5.py')
    run('e8_phi4_sylow5.g','gap4',alias='e8_phi4_sylow5_gap_observed.txt')
    run('e8_phi4_sylow5.py')
    assert json.loads((target/'computations/e8_phi4_sylow5.json').read_text())['relative_group_order']==46080
try:
    if a.suite in ['gap3','all']:gap3_suite()
    if a.suite in ['python','all']:python_suite()
    if a.suite in ['gap4','all']:gap4_suite()
finally:
    (target/'receipt.json').write_text(json.dumps({'suite':a.suite,'runs':receipts},indent=2)+'\n')
print('Completed suite:',a.suite)
