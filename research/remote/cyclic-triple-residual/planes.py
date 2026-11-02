"""Intersect the E1-admissible degree-8 coset  <D8> + J_8^(sector)  with  I_Pi^k (vanishing to order k along a coordinate
plane Pi, 3 coordinates = 0) and I_L^k (coordinate line L, 4 coordinates = 0).  Exact linear algebra over Q; then test the
resulting denominators (sector solver, exact when feasible)."""
import sys, time, json, itertools, random
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from sect import *
WD = '/tmp/hodge/results/cyclic-triple-residual'
def log(e):
    with open(f'{WD}/log.jsonl', 'a') as f: f.write(json.dumps(e) + '\n')
F, A, q, rec = core.load_instance('cyclic-triple-residual')
V = core.VARS
D8 = core.parse(json.load(open(f'{WD}/kernel_d8_0_4_0.json'))['D'])
Fz = integral(F); G = semi_invariant_subgroup(elems(Fz), [integral(A), integral(D8)]); gens = generators(G); chi = make_chi(gens)
c8 = chi_of(integral(D8), chi)
Fi = [core.pderiv(F, i) for i in range(6)]
Dm = [m for m in mons(8) if chi(m) == c8]; idx = {m: i for i, m in enumerate(Dm)}
Jpolys = []
for i in range(6):
    for m in mons(3):
        if cadd(chi(m), tuple((-x) % 6 for x in chi(E6[i]))) == c8:
            Jpolys.append(core.mul({m: Fraction(1)}, Fi[i]))
gensP = [D8] + Jpolys          # D = mu*D8 + sum lambda_k J_k ; unknown vector (mu, lambda)
nU = len(gensP)
def coeff_vec(P):
    v = [Fraction(0)]*len(Dm)
    for m, c in P.items(): v[idx[m]] += Fraction(c)
    return v
Cols = [coeff_vec(P) for P in gensP]
def solve_vanishing(low_mons):
    """Solutions (mu, lambda) with coefficient of every monomial in low_mons = 0.  Returns list of basis vectors (Fractions)."""
    rows = [[Cols[j][idx[m]] for j in range(nU)] for m in low_mons]
    if not rows: return [[Fraction(int(i == j)) for j in range(nU)] for i in range(nU)]
    L = 1
    for r in rows:
        for x in r: L = L * x.denominator // __import__('math').gcd(L, x.denominator)
    M = fmpq_mat(len(rows), nU, [int(x*L) for r in rows for x in r], 1) if False else fmpq_mat([[fmpq(x.numerator, x.denominator) for x in r] for r in rows])
    # nullspace via rref
    R, rk = M.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < nU and R[i, j] == 0: j += 1
        piv.append(j)
    free = [j for j in range(nU) if j not in piv]
    basis = []
    for f in free:
        v = [Fraction(0)]*nU; v[f] = Fraction(1)
        for i, pj in enumerate(piv):
            x = R[i, f]; v[pj] = -Fraction(int(x.p), int(x.q))
        basis.append(v)
    return basis
def build(vecv):
    D = {}
    for c, P in zip(vecv, gensP):
        if c: D = core.add(D, core.scale(P, c))
    return D
results = []
tested = 0
def consider(name, low_mons, kind):
    global tested
    basis = solve_vanishing(low_mons)
    with_mu = [b for b in basis if b[0] != 0]
    dim_aff = len(basis) - 1 if with_mu else None
    print(f'{name}: solution space dim {len(basis)}, contains D8-direction: {bool(with_mu)}', flush=True)
    results.append({'name': name, 'kind': kind, 'dim': len(basis), 'has_D8': bool(with_mu)})
    if not with_mu: return
    # normalise mu = 1: particular solution; then a random affine point too
    b0 = with_mu[0]; b0 = [x / b0[0] for x in b0]
    pts = [('particular', b0)]
    if len(basis) > 1:
        random.seed(1)
        v = list(b0)
        for b in basis:
            if b is with_mu[0]: continue
            c = Fraction(random.randint(-2, 2))
            if b[0] != 0:  # keep mu = 1
                b = [x - b[0]*y for x, y in zip(b, b0)]
            v = [x + c*y for x, y in zip(v, b)]
        pts.append(('random-affine', v))
    for tag, vv in pts:
        D = build(vv)
        if not D: continue
        t0 = time.time()
        out = test(F, A, 1, D, exact=(tested < 6), verbose=False); tested += 1
        ent = {'engine': 2, 'D': f'{name} [{tag}]: ' + core.to_str(D)[:300], 'd': 8, 'family': f'E1-coset-cap-{kind}', 'e1': True,
               'consistent': out['consistent'], 'rank': out['rank'], 'nullity': out['nullity'], 'seconds': round(time.time()-t0, 1),
               'note': f"D in (<D8>+J_8) with vanishing order >= k along {name}; solution space dim {len(basis)}; rows={out['rows']} cols={out['cols']} rank_aug={out['rank_aug']} exact={out['exact_rank']} consistent2={out.get('consistent2')}"}
        if out.get('verified'):
            ent['witness'] = {'d': 8, 'D': core.to_str(D), 'N': out['N'], 'Q': out['Q']}; ent['note'] += ' | *** SOLVED ***'
            with open(f'{WD}/witness-candidate-{int(time.time())}.json', 'w') as f:
                json.dump({'instanceId': 'cyclic-triple-residual', 'family': ent['family'], 'verified_by_core': True, **ent['witness']}, f, indent=1)
            print('*** SOLVED ***', name, flush=True)
        log(ent)
        print(f'    test [{tag}]: consistent={out["consistent"]} rank={out["rank"]}/{out["rank_aug"]} cols={out["cols"]} ({ent["seconds"]}s)', flush=True)
# planes: 3 vanishing coordinates; find the maximal k with a solution containing D8, report k=3.. as well
for S in itertools.combinations(range(6), 3):
    nm = 'Pi{' + ','.join(V[i] for i in S) + '=0}'
    kmax = None
    for k in range(8, 2, -1):
        low = [m for m in Dm if sum(m[i] for i in S) < k]
        basis = solve_vanishing(low)
        if any(b[0] != 0 for b in basis): kmax = k; break
    print(f'{nm}: maximal vanishing order along the plane within the E1 coset: {kmax}', flush=True)
    if kmax is not None and kmax >= 3:
        consider(f'{nm} k={kmax}', [m for m in Dm if sum(m[i] for i in S) < kmax], 'plane')
        if kmax > 3:
            consider(f'{nm} k=3', [m for m in Dm if sum(m[i] for i in S) < 3], 'plane')
# lines: 4 vanishing coordinates, k >= 4
for S in itertools.combinations(range(6), 4):
    nm = 'L{' + ','.join(V[i] for i in S) + '=0}'
    kmax = None
    for k in range(8, 3, -1):
        low = [m for m in Dm if sum(m[i] for i in S) < k]
        basis = solve_vanishing(low)
        if any(b[0] != 0 for b in basis): kmax = k; break
    print(f'{nm}: maximal vanishing order along the line within the E1 coset: {kmax}', flush=True)
    if kmax is not None and kmax >= 4:
        consider(f'{nm} k={kmax}', [m for m in Dm if sum(m[i] for i in S) < kmax], 'line')
json.dump(results, open(f'{WD}/planes_results.json', 'w'), indent=1)
