import itertools,json
from fractions import Fraction as Q
# Gaussian rationals represented by Python complex are exact here: dyadic entries only.
def mm(a,b):return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2)) for i in range(2) for j in range(2))
I=(1,0,0,1); X=(0,1,1,0); Z=(1,0,0,-1); S=(1,0,0,1j); H=tuple((1+1j)/2*x for x in (1,1,1,-1))
def closure(gs):
 seen={I}; todo=[I]
 while todo:
  a=todo.pop()
  for g in gs:
   b=mm(a,g)
   if b not in seen:seen.add(b);todo.append(b)
 return seen
def powm(a,n):
 b=I
 for _ in range(n):b=mm(b,a)
 return b
W=closure([tuple(1j*x for x in I),X,Z]);R=closure([S,H])
ref=[a for a in W if a!=I and (a[0]-1)*(a[3]-1)-a[1]*a[2]==0]
perms=list(itertools.permutations(range(4)))
signs=[s for s in itertools.product([-1,1],repeat=4) if s[0]*s[1]*s[2]*s[3]==1]
def wt(v):return all(x.denominator==1 for x in v) or all((x-Q(1,2)).denominator==1 for x in v)
def stab(v):return sum(wt(tuple(s[i]*v[p[i]]-v[i] for i in range(4))) for p in perms for s in signs)
roots=[tuple(Q(a if k==i else b if k==j else 0) for k in range(4)) for i in range(4) for j in range(i+1,4) for a,b in itertools.product([-1,1],repeat=2)]
pts=[(Q(1,2),Q(1,2),Q(0),Q(0)),(Q(1,2),Q(1,4),Q(1,4),Q(0))]
out={'W_order':len(W),'R_order':len(R),'S2_eq_Z':powm(S,2)==Z,'H2_eq_iI':powm(H,2)==tuple(1j*x for x in I),'SH3_eq_minusI':powm(mm(S,H),3)==tuple(-x for x in I),'W_nonidentity_eigen1_count':len(ref),'all_eigen1_involutions':all(powm(a,2)==I for a in ref),'D4_points':[{'point':list(map(str,v)),'full_W_stabilizer':stab(v),'centralizing_roots':sum(sum(x*y for x,y in zip(v,r)).denominator==1 for r in roots)} for v in pts]}
print(json.dumps(out,indent=2))
open('outputs/odd_d_round4_exact_check.json','w').write(json.dumps(out,indent=2))
