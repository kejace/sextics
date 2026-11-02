# wave-4 amix-component: exact decomposition A = A_K3 + A_mix in physical coordinates, Jacobian non-membership, E1 kernels
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
from fractions import Fraction
LOG='/tmp/hodge/results/polarized-k3-remainder/wave4/log.jsonl'
def log(**kw):
    kw['ts']=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    open(LOG,'a').write(json.dumps(kw, default=str)+'\n')
F, A, q, rec = core.load_instance('polarized-k3-remainder')
X=[core.parse(v) for v in ['X0','X1','X2','Y0','Y1','Y2']]
U=[core.add(X[i],X[i+3]) for i in range(3)]      # U_i = X_i + Y_i
V=[core.sub(X[i],X[i+3]) for i in range(3)]      # V_i = X_i - Y_i
def sq(P): return core.mul(P,P)
U2=[sq(u) for u in U]; V2=[sq(v) for v in V]
c=Fraction(1,64)
# Ahat = 2^-6 prod (U_i^2 - V_i^2); expand by V-degree b
terms={}
import itertools
for choice in itertools.product([0,1],repeat=3):   # 1 = pick -V_i^2
    P=core.parse('1')
    for i,ch in enumerate(choice):
        P=core.mul(P, V2[i] if ch else U2[i])
    sign=(-1)**sum(choice); b=2*sum(choice)
    terms[b]=core.add(terms.get(b,{}), core.scale(P, c*sign))
AK3=core.add(terms[0],terms[6]); A2=terms[2]; A4=terms[4]; Amix=core.add(A2,A4)
assert core.sub(core.add(AK3,Amix),A)=={} , 'decomposition failed'
print('A = A_K3 + A_mix exact:', core.to_str(AK3)[:200])
print('A_2 =', core.to_str(A2)); print('A_4 =', core.to_str(A4))
print('A_mix =', core.to_str(Amix))
# which are Y-odd? parity check
for nm,P in [('A_K3',AK3),('A_2',A2),('A_4',A4),('A_mix',Amix)]:
    print(nm, 'nterms', len(P), 'parity', core.parity_of(P))
out={'A_K3':core.to_str(AK3),'A_2':core.to_str(A2),'A_4':core.to_str(A4),'A_mix':core.to_str(Amix)}
json.dump(out, open('/tmp/hodge/results/polarized-k3-remainder/wave4/amix_pieces.json','w'), indent=1)
log(what='amix-decomposition-exact', result='A = A_K3 + A_mix verified exactly over Q (sympy-free, engine arithmetic)', file='wave4/amix_pieces.json')
# Jacobian membership mod two primes
for p in [core.PRIMES[0], core.PRIMES[1]]:
    res={}
    for nm,P in [('A',A),('A_K3',AK3),('A_2',A2),('A_4',A4),('A_mix',Amix)]:
        res[nm]=bool(jacobian.in_jacobian(F,P,p))
    print('p',p,'in J(F)?',res)
    log(what='jacobian-membership-pieces', prime=p, result=res)
