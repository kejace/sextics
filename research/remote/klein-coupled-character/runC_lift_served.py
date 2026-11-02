import sys, time, json, pickle
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from kcc2 import *
from kcc import WSystem, served_targets
F, A, q = load()
Fi = [core.pderiv(F, i) for i in range(6)]
X = [core.parse(v) for v in core.VARS]
out = {}
# (1) rational lift of the d=9 E1-kernel representatives (canonical rref basis, 3 primes)
def kernel_reps_rational(d, k, parD, nprimes=4):
    acc = None
    for p in core.PRIMES[:nprimes]:
        ctx = Ctx(F, A, p)
        reps, _ = e1_kernel(ctx, d, k, parD)
        Dm = wmonomials(d, parD, k)
        vecs = [[D.get(m, 0) for m in Dm] for D in reps]
        if acc is None: acc = (vecs, p); continue
        xs, mmod = acc
        if len(xs) != len(vecs): return None, reps
        xs = [[core.crt_pair(a, mmod, b, p)[0] for a, b in zip(v1, v2)] for v1, v2 in zip(xs, vecs)]
        acc = (xs, mmod * p)
        rat = [[core.ratrecon(a, acc[1]) for a in v] for v in acc[0]]
        if all(a is not None for v in rat for a in v):
            return [{m: a for m, a in zip(Dm, v) if a} for v in rat], reps
    return None, reps
lifted = {}
for parD in [0, 1]:
    for k in range(7):
        t = time.time()
        rats, reps = kernel_reps_rational(9, k, parD)
        print(f'--- d=9 parD={parD} k={k}: {len(reps)} reps; rational lift {"ok" if rats else "FAILED"} [{time.time()-t:.1f}s]', flush=True)
        if rats:
            # verify over Q: A*D in J via e1_lift at a fresh prime
            ctx = Ctx(F, A, core.PRIMES[5])
            for i, D in enumerate(rats):
                ok, _, _ = e1_lift(ctx, D)
                den = max(abs(c.denominator) for c in D.values()); num = max(abs(c.numerator) for c in D.values())
                print(f'   rep {i}: {len(D)} terms, max |num| {num}, max den {den}, E1 ok at 6th prime: {ok}')
                print('      D =', core.to_str(D)[:600])
            lifted[(parD, k)] = rats
pickle.dump({str(kk): v for kk, v in lifted.items()}, open('kernel_d9_rational.pkl', 'wb'))
json.dump({str(kk): [core.to_str(D) for D in v] for kk, v in lifted.items()}, open('kernel_d9_rational.json', 'w'), indent=1)
# (2) served chi3 targets for various D (engine WSystem, weight-restricted, target weight 3)
Jhat = core.parse(open('/home/kejace/devel/sextics/data/klein-coupled-character/instance.rec').read().split('F: ')[1].split('\n')[0])
Fstr = open('/home/kejace/devel/sextics/data/klein-coupled-character/instance.rec').read().split('F: ')[1].split('\n')[0]
Jhat = core.parse(f'7*X0*X2**2*(9*(Y0**6 - Y1**6 + Y2**6) - 9*({Fstr})) + 15*X1*Y2*Y0*(2*Y0**6 - 4*Y1**6 + 8*Y2**6 - 2*({Fstr}))')
cands = [('Jhat', Jhat), ('F_0', Fi[0]), ('F_1', Fi[1]), ('F_2', Fi[2]), ('F_3', Fi[3]), ('X0*F_0', core.mul(X[0], Fi[0])), ('Y0*F_3', core.mul(X[3], Fi[3]))]
for (parD, k), rats in lifted.items():
    for i, D in enumerate(rats): cands.append((f'kernel d=9 parD={parD} k={k} #{i}', D))
for name, D in cands:
    t = time.time()
    try:
        dim, basis, S = served_targets(F, 1, D)
        note = f'served chi3-sector targets (Y-odd deg 6, char 3; sector dim 26 mod J): dim {dim}; served polys: ' + ' | '.join(core.to_str(b)[:80] for b in basis[:6])
        log({'D': f'served-targets {name}', 'd': core.degree(D), 'family': 'served-targets', 'engine': 2, 'k': sigma_weight(D), 'e1': None,
             'consistent': None, 'rank': dim, 'nullity': None, 'seconds': round(time.time()-t, 1), 'note': note})
        print(f'SERVED {name}: dim {dim} [{time.time()-t:.1f}s]  ' + ' | '.join(core.to_str(b)[:60] for b in basis[:4]), flush=True)
    except Exception as ex:
        print('served-targets failed for', name, ex, flush=True)
# (3) Hessian 3x3 block determinants (degree 12, char 0)
def det3(M):
    return core.add(core.add(core.sub(core.mul(M[0][0], core.sub(core.mul(M[1][1], M[2][2]), core.mul(M[1][2], M[2][1]))),
                    core.mul(M[0][1], core.sub(core.mul(M[1][0], M[2][2]), core.mul(M[1][2], M[2][0])))),
                    core.mul(M[0][2], core.sub(core.mul(M[1][0], M[2][1]), core.mul(M[1][1], M[2][0])))), {})
H = [[core.pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
HXX = det3([[H[i][j] for j in range(3)] for i in range(3)])
HYY = det3([[H[i][j] for j in range(3, 6)] for i in range(3, 6)])
HXY = det3([[H[i][j] for j in range(3, 6)] for i in range(3)])
for name, D in [('det Hess_XX', HXX), ('det Hess_YY', HYY), ('det Hess_XY', HXY)]:
    print(name, 'deg', core.degree(D), 'char', sigma_weight(D), 'parity', core.parity_of(D), 'terms', len(D), flush=True)
    e = full_test(F, A, D, name, 'hessian-block', note='3x3 block of the Hessian matrix of F', verbose=True)
    print(f'{name}: e1={e["e1"]} consistent={e["consistent"]} rank={e["rank"]} nullity={e["nullity"]} {e["seconds"]}s', flush=True)
# (4) served chi3 targets for pure-Y and pure-X cone denominators (all sigma-invariant; E1 for A itself fails at deg<=8, but other sector targets may be served)
cones = [('Y0**3+Y1**3+Y2**3', 'cubic Y-cone'), ('Y0*Y1*Y2', 'triangle Y-cone'), ('Y0**6+Y1**6+Y2**6', 'Fermat Y-block G'), ('Y0**2', 'double hyperplane Y0'),
         ('Y0**4+Y1**4+Y2**4', 'quartic Y-cone'), ('Y0**3*Y1 + Y1**3*Y2 + Y2**3*Y0', 'Klein-quartic Y-cone'),
         ('X0*X1*X2', 'triangle X-cone'), ('X0**5*X1 + X1**5*X2 + X2**5*X0', 'Klein X-block K'), ('X0**3*X2', 'invariant X quartic monomial'), ('X0**2*X1**2*X2**2', '(xyz)^2')]
for expr, desc in cones:
    D = core.parse(expr); t = time.time()
    try:
        dim, basis, S = served_targets(F, 1, D)
        log({'D': f'served-targets {expr}', 'd': core.degree(D), 'family': 'served-targets-cone', 'engine': 2, 'k': sigma_weight(D), 'e1': None,
             'consistent': None, 'rank': dim, 'nullity': None, 'seconds': round(time.time()-t, 1),
             'note': f'{desc}: served chi3-sector targets dim {dim}; ' + ' | '.join(core.to_str(b)[:80] for b in basis[:6])})
        print(f'SERVED-CONE {expr} ({desc}): dim {dim} [{time.time()-t:.1f}s] ' + ' | '.join(core.to_str(b)[:60] for b in basis[:4]), flush=True)
    except Exception as ex:
        print('served-targets failed for', expr, ex, flush=True)
