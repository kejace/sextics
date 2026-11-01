import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import NV, monomials, shift, pderiv
from flint import nmod_mat
import sympy as sp
P = core.PRIMES[0]
F, A, q, rec = core.load_instance('reflection-coupled-residual')
K8 = core.parse(json.load(open('kernel_d8_X03Y03.json'))[0])
X = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
def tosym(Pl): return sum(sp.Rational(c.numerator, c.denominator)*sp.prod([X[i]**e[i] for i in range(6)]) for e, c in Pl.items())
k = tosym(K8)
print('K8 terms', len(K8), 'factor:', sp.factor(k))
# substitution checks
print('K8 sector: Y2-even', all(e[5]%2==0 for e in K8), ' (Y0+Y1)-parity even', all((e[3]+e[4])%2==0 for e in K8))
# normal forms mod J_8 for comparison
Fp = core.reduce_mod(F, P); Fi = [pderiv(Fp, i, P) for i in range(NV)]
mons = monomials(8, 0); idx = {m: i for i, m in enumerate(mons)}
Jd = []
for i in range(NV):
    pi = 0 if i < 3 else 1
    for m in monomials(3, (0 - pi) % 2): Jd.append(shift(Fi[i], m, 1, P))
R, C = len(Jd), len(mons); flat = [0]*(R*C)
for r, poly in enumerate(Jd):
    for e, v in poly.items(): flat[r*C+idx[e]] = v
Mr, rJ = nmod_mat(R, C, flat, P).rref()
piv = []
for i in range(rJ):
    j = 0
    while int(Mr[i, j]) == 0: j += 1
    piv.append(j)
def nf(poly):
    v = [0]*C
    for e, c in poly.items(): v[idx[e]] = (v[idx[e]] + int(c)) % P
    for i, j in enumerate(piv):
        if v[j]:
            c = v[j]
            for kk in range(C):
                a = int(Mr[i, kk])
                if a: v[kk] = (v[kk] - c*a) % P
    return v
vK = nf(core.reduce_mod(K8, P))
def prop(v, w):
    # is v proportional to w?
    ratio = None
    for a, b in zip(v, w):
        if a == 0 and b == 0: continue
        if a == 0 or b == 0: return None
        r = a * pow(b, -1, P) % P
        if ratio is None: ratio = r
        elif ratio != r: return None
    return ratio
Hu = core.parse('3*Y2**4 + 5*X0*X1*Y0**2/7 + 10*X0*X1*Y1**2/7')  # dF/du, u=Y2^2
cands = {'Hu^2': core.mul(Hu, Hu), 'Y2^4*Hu': core.mul(core.parse('Y2**4'), Hu), 'Y2^8': core.parse('Y2**8'),
         'X0^4*Y0^4': core.parse('X0**4*Y0**4'), 'Y0^4*Y1^4': core.parse('Y0**4*Y1**4'), 'X0^8': core.parse('X0**8')}
for name, G in cands.items():
    print(name, 'nf nonzero coords', sum(1 for x in nf(core.reduce_mod(G, P)) if x), 'proportional to K8:', prop(nf(core.reduce_mod(G, P)), vK))
print('K8 nf nonzero coords', sum(1 for x in vK if x), 'of', C - rJ)
# gradient / singular locus of K8 alone (affine cone): check whether K8 is singular along something big via random point test of rank of Hessian? skip
# evaluate K8 on a few structured loci
print('K8|Y2=0 factor:', sp.factor(k.subs(X[5], 0)))
print('K8|X0=X1=0:', sp.factor(k.subs({X[0]:0, X[1]:0})))
print('K8|Y0=Y1=0:', sp.factor(k.subs({X[3]:0, X[4]:0})))
