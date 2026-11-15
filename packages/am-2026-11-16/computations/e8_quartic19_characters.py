"""Explicit Hom(G,C4), including prescribed central line and true returns."""
from e8_quartic19_phase import *
from collections import deque

def characters(G):
 D=G.derived_subgroup(); table=G.coset_table(D); N=len(table)
 qgens=[Permutation([table[i][2*j] for i in range(N)]) for j in range(len(G.generators))]
 assert int(PermutationGroup(qgens).order())==N
 actors={0:Permutation(list(range(N)))};todo=deque([0])
 while todo:
  a=todo.popleft()
  for gen in qgens:
   b=gen(a)
   if b not in actors:actors[b]=actors[a]*gen;todo.append(b)
 assert len(actors)==N
 mul=lambda a,b:actors[b](a)
 H={0}; chars=[{0:0}]
 for gen in qgens:
  g=gen(0)
  if g in H:continue
  powers=[0];x=g
  while x not in H:powers.append(x);x=mul(x,g)
  r=len(powers);newchars=[]
  for chi in chars:
   for v in range(4):
    if r*v%4!=chi[x]:continue
    c={mul(a,h):(k*v+chi[h])%4 for k,a in enumerate(powers) for h in H}
    newchars.append(c)
  H={mul(a,h) for a in powers for h in H};chars=newchars
 assert len(H)==N and len(chars)==N,('abelian exponent exceeds four',N,len(chars))
 genmap={tuple(g.array_form):q for g,q in zip(G.generators,qgens)}
 genmap.update({tuple((g**-1).array_form):q**-1 for g,q in zip(G.generators,qgens)})
 def imageid(x):
  out=0
  for g in G.generator_product(x,original=True):out=genmap[tuple(g.array_form)](out)
  return out
 return chars,imageid,qgens,table

out=[]
for row in ROWS:
 if 'signed' not in row:
  d=pickle.loads((OUT/('e8_quartic19_'+row['name']+'.pickle')).read_bytes());G=PermutationGroup(list(d['E0'].args));a,_=longest(SIMPLE);b,_=longest(row['base']);h=a*b**-1;n=h
  zlist=[(torus2(coeff_y(v)),2*zeta(v,row)) for v in combos(d['ker2'],2)]
  extras=[(n,0),(n*h,0)]
 else:
  fn='e8_quartic19_'+row['name']+('_corrected_phase.pickle' if row['name']=='A7A1_external_A1' else '_phase.pickle')
  d=pickle.loads((OUT/fn).read_bytes());G=PermutationGroup([Permutation(x) for x in d['Qgens']]);n=Permutation(d['n']);h=Permutation(d['h']);zlist=[(Permutation(x),v) for y,x,v in d['zlist']]
  extras=[(n,2 if row['name']=='A7A1_external_A1' else 0),(h,2)]
 chars,im,qgens,table=characters(G)
 constraints=[(im(g),v) for g,v in zlist+extras]
 good=[chi for chi in chars if all(chi[a]==v for a,v in constraints)]
 assert good
 chi=good[0]
 result={'row':row['name'],'group_order':int(G.order()),'abelian_order':len(table),'number_C4_characters':len(chars),'number_solutions_with_returns':len(good),'rho_generator_exponents_i':[chi[q(0)] for q in qgens],'rho_n_exponent_i':chi[im(n)],'rho_h_exponent_i':chi[im(h)],'rho_nh_exponent_i':chi[im(n*h)],'abelian_coset_character_exponents_i':[chi[i] for i in range(len(table))],'rho_iota_exponent_i':chi[im(IOTA)] if 'signed' in row else None,'all_quotient_edges_verified':all(chi[table[a][2*j]]==(chi[a]+chi[qgens[j](0)])%4 for a in range(len(table)) for j in range(len(qgens)))}
 print(json.dumps(result),flush=True);out.append(result)
(OUT/'e8_quartic19_explicit_characters.json').write_text(json.dumps(out,indent=2))
