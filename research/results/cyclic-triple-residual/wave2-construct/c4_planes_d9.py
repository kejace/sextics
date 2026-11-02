# C4: degree-9 mu3-invariant E1-admissible denominators vanishing to order >= k along coordinate planes (k = 3, 4),
# per G-sector c0x (x = 0..5), using the exact E1 bases of the search lens; j=0 diagnostics for planes inside {X2 = 0};
# then mod-p sector tests (E2) of random small-integer members of the largest families.
import os, sys, json, time, random, itertools
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
import flint
from w2core import *
OUT = '/tmp/hodge/results/cyclic-triple-residual/wave2-construct'
LOG = '/tmp/hodge/results/cyclic-triple-residual/wave2-construct.log.jsonl'
def log(e):
    e['ts'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()); e['engine'] = 2
    open(LOG, 'a').write(json.dumps(e) + '\n')
F, A, q = load()
REM = '/tmp/hodge/remote/cyclic-triple-residual'
names = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
NTEST = int(sys.argv[1]) if len(sys.argv) > 1 else 4
spaces = {c: load_space(f'{REM}/w2-e1space-d9-c0{c}.json') for c in range(6)}
def kernel_combos(basis, P, k):
    """exact integer kernel basis (columns) of the order->=k conditions along plane {Z_i = 0, i in P}."""
    mons = {}
    for j, (D, N) in enumerate(basis):
        for e, v in D.items():
            if sum(e[i] for i in P) < k:
                mons.setdefault(e, {})[j] = Fraction(v)
    n = len(basis)
    if not mons:
        return [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    rows = []
    for row in mons.values():
        den = 1
        for v in row.values(): den = den * v.denominator // __import__('math').gcd(den, v.denominator)
        rows.append({j: int(v * den) for j, v in row.items()})
    M = flint.fmpz_mat(len(rows), n)
    for r, row in enumerate(rows):
        for j, v in row.items(): M[r, j] = v
    X, nullity = M.nullspace()
    return [[int(X[i, c]) for i in range(n)] for c in range(nullity)]
results = []; families = []
t0 = time.time()
for P in itertools.combinations(range(6), 3):
    for c in range(6):
        basis = spaces[c]['basis']
        rec = {'plane': [names[i] for i in P], 'in_branch_locus': 2 in P, 'sector': '0%d' % c, 'dim_E1': len(basis)}
        for k in (3, 4):
            K = kernel_combos(basis, P, k)
            rec['dim_order%d' % k] = len(K)
            if k == 3 and K:
                families.append((P, c, K))
        results.append(rec)
print('families with order>=3: ', sum(1 for r in results if r['dim_order3'] > 0), ' time %.1fs' % (time.time() - t0), flush=True)
for r in results:
    if r['dim_order3'] > 0: print(r, flush=True)
json.dump(results, open(f'{OUT}/c4_planes_d9_dims.json', 'w'), indent=1)
log({'step': 'C4-dims', 'result': 'degree-9 mu3-invariant E1-admissible order>=3/4 families per coordinate plane and sector',
     'n_families_order3': sum(1 for r in results if r['dim_order3'] > 0), 'file': 'wave2-construct/c4_planes_d9_dims.json'})
# ---- j = 0 diagnostics and E2 tests on the largest families
def member(basis, K, rng):
    coeffs = [0] * len(basis)
    for col in K:
        r = rng.randint(-3, 3)
        if r: coeffs = [a + r * b for a, b in zip(coeffs, col)]
    D = {}
    for a, (Dj, N) in zip(coeffs, basis):
        if a: D = padd(D, scale(Dj, a))
    return D
def diag_j0(D, P):
    """for a plane inside {X2=0}: leading cubic along the plane = c3(n1,n2) + a*X2^3; return whether a and c3 are nonzero."""
    a = {}; c3 = {}
    for e, v in D.items():
        nd = sum(e[i] for i in P)
        if nd == 3:
            if e[2] == 3: a[e] = v
            elif e[2] == 0: c3[e] = v
    return bool(a), bool(c3)
rng = random.Random(7)
families.sort(key=lambda f: -len(f[2]))
pref = [f for f in families if 2 in f[0]] + [f for f in families if 2 not in f[0]]
tested = 0
for P, c, K in pref:
    if tested >= NTEST: break
    basis = spaces[c]['basis']
    D = member(basis, K, rng)
    if not D: continue
    ordP = min(sum(e[i] for i in P) for e in D)
    info = {'step': 'C4-test', 'plane': [names[i] for i in P], 'sector': '0%d' % c, 'family_dim': len(K), 'order_along_plane': ordP,
            'terms': len(D), 'd': core.degree(D)}
    if 2 in P:
        info['leading_cubic_has_X2^3'], info['leading_cubic_has_c3'] = diag_j0(D, P)
    t1 = time.time()
    try:
        res = sect.test(F, A, q, D, exact=False, verbose=False, lift=True)
        info.update({k: v for k, v in res.items() if k not in ('N', 'Q')})
    except Exception as ex:
        info['error'] = repr(ex)
    info['seconds'] = round(time.time() - t1, 1)
    info['D_head'] = core.to_str(D)[:300]
    print(json.dumps(info)[:700], flush=True)
    log(info); tested += 1
print('done', flush=True)
