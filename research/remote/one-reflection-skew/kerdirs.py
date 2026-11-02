"""Test the E1-kernel-mod-J directions (5-var, rho-even) at degree d: each reduced basis vector and random combos."""
import sys, pickle, random, os
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
from five import *
from flint import nmod_mat
random.seed(2)
d = int(sys.argv[1]); ncombo = int(sys.argv[2]) if len(sys.argv) > 2 else 6
p = P1
K = pickle.load(open(f'e1ker5_d{d}_par0.pkl', 'rb'))
basis = K['basis']  # rref basis of the full kernel (includes J_d)
# J_d even basis
Jc = []
for i in range(5):
    pi = (0 - (0 if i < 3 else 1)) % 2
    for m in mons5(d - 5, pi):
        Jc.append(shift(fi(i, p), m, 1, p))
Dm = mons5(d, 0); idx = {m: k for k, m in enumerate(Dm)}; nD = len(Dm)
def vec(P):
    v = [0] * nD
    for e, c in P.items(): v[idx[e]] = int(c) % p
    return v
MJ = nmod_mat(len(Jc), nD, [x for P in Jc for x in vec(P)], p)
Jr, rJ = MJ.rref()
Jrows = [[int(Jr[i, j]) for j in range(nD)] for i in range(rJ)]
# reduce each kernel vector modulo J (row-reduce against J's rref rows) and keep independent remainders
extra = []
cur = list(Jrows)
for b in basis:
    M = nmod_mat(len(cur) + 1, nD, [x for r in cur for x in r] + vec(b), p)
    if M.rank() > len(cur):
        cur.append(vec(b)); extra.append(b)
print(f'd={d}: kernel {len(basis)}, J {rJ}, extra directions {len(extra)}', flush=True)
def lift(poly):
    return {e: Fraction(v if v <= p // 2 else v - p) for e, v in poly.items()}
def small(poly):
    # reduce coefficients mod J to a small-height representative is hard; use symmetric residues
    return lift(poly)
cands = []
for k, b in enumerate(extra):
    cands.append((f'<E1ker5 d{d} extra #{k} ({len(b)} terms)>', lift(b)))
for k in range(ncombo):
    comb = {}
    for b in extra:
        c = random.randint(-2, 2)
        if c: comb = add(comb, scale(b, c % p, p), p)
    # add a random J element too
    for P in random.sample(Jc, min(3, len(Jc))):
        comb = add(comb, scale(P, random.randint(1, 5), p), p)
    cands.append((f'<E1ker5 d{d} random combo #{k}>', lift(comb)))
os.makedirs('served', exist_ok=True)
for name, D in cands:
    try:
        r = test5(name, f'5var:E1-kernel-directions-d{d}', note='lifted symmetric residues of a mod-p kernel vector', Dpoly=D)
    except Exception as ex:
        print('ERR', name, ex); continue
    fn = name.replace('/', '_').replace('*', 'x').replace('|', 'I').replace('<','').replace('>','').replace(' ','_') + '.pkl'
    with open('served/' + fn, 'wb') as fh:
        pickle.dump({'name': name, 'family': f'E1-kernel-directions-d{d}', 'D': {e: str(c) for e, c in D.items()}, 'p': r['p'], 'd': r['d'],
                     'served_basis': r['served_basis'], 'served_dim': r['served_dim'], 'served_mod_J': r['served_mod_J'],
                     'serves_A': r['serves_A'], 'serves_A+': r['serves_A+'], 'serves_A-': r['serves_A-'], 'sigma_sector': r['sigma_sector'],
                     'rank_N': r['rank_N'], 'nullity_N': r['nullity_N'], 'rows': r['rows'], 'cols': r['cols'], 'seconds': r['seconds']}, fh)
print('DONE')
