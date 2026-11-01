"""(B) geometrically motivated candidates, degrees 5..12, sector solver (exact when d<=9, else mod p with 2 primes)."""
import sys, time, json, random
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from sect import *
WD = '/tmp/hodge/results/cyclic-triple-residual'
def log(e):
    with open(f'{WD}/log.jsonl', 'a') as f: f.write(json.dumps(e) + '\n')
F, A, q, rec = core.load_instance('cyclic-triple-residual')
P = core.parse
t = Fraction(5, 7)
X0,X1,X2,Y0,Y1,Y2 = [P(v) for v in core.VARS]
C = P('X0**3 + X0*Y1*Y2 + X1*Y0*Y1')
H = P('X0**6 + X1**6 + Y0**6 + Y1**6 + Y2**6 + 5*X0**3*X1*Y0*Y1/7 + 5*X0**2*Y0*Y1*Y2**2/7 + 5*X0*X1*Y1**3*Y2/7')
Z = P('X2**3'); Zp = core.add(core.scale(Z, 2), core.scale(C, t)); Hp = core.sub(H, core.scale(core.mul(C, C), t*t/4))
Fi = [core.pderiv(F, i) for i in range(6)]
D8 = P(json.load(open(f'{WD}/kernel_d8_0_4_0.json'))['D'])
M = core.mul
def run(D, name, fam, note='', exact=None):
    d = core.degree(D)
    if exact is None: exact = (d <= 6)
    t0 = time.time()
    try:
        out = test(F, A, 1, D, exact=exact, verbose=False)
    except AssertionError as ex:
        print('skip', name, ex, flush=True); return
    ent = {'engine': 2, 'D': name[:400], 'd': d, 'family': fam, 'e1': None, 'consistent': out['consistent'], 'rank': out['rank'],
           'nullity': out['nullity'], 'seconds': round(time.time()-t0, 1),
           'note': (note + f" | sector solver |G|={out['G_order']} rows={out['rows']} cols={out['cols']} rank_aug={out['rank_aug']} exact={out['exact_rank']}"
                    + ('' if out['exact_rank'] else f" p={out.get('prime')} consistent2={out.get('consistent2')}")).strip(' |')}
    if out.get('verified'):
        ent['witness'] = {'d': d, 'D': core.to_str(D), 'N': out['N'], 'Q': out['Q']}; ent['note'] += ' | *** SOLVED ***'
        with open(f'{WD}/witness-candidate-{int(time.time())}.json', 'w') as f:
            json.dump({'instanceId': 'cyclic-triple-residual', 'family': fam, 'verified_by_core': True, **ent['witness']}, f, indent=1)
        print('*** SOLVED ***', name, flush=True)
    log(ent)
    print(f"d={d:2d} cons={str(out['consistent']):5s} rank={out['rank']:5d}/{out['rank_aug']:5d} cols={out['cols']:5d} rows={out['rows']:5d} |G|={out['G_order']:2d} {ent['seconds']:6.1f}s {fam:18s} {name[:70]}", flush=True)
random.seed(5)
def polar(a):
    D = {}
    for ai, Fi_ in zip(a, Fi):
        if ai: D = core.add(D, core.scale(Fi_, Fraction(ai)))
    return D
for k in range(1):
    a = [random.randint(-5, 5) for _ in range(6)]
    run(polar(a), f'polar F_a, a={a}', 'polar-pencil', 'ramification divisor of projection from a random point')
for a, nm in [((0,1,0,0,0,1), 'a on line {X0=Y1=Y0=0}'), ((0,0,0,1,0,1), 'a on line {X0=Y1=X1=0}'), ((0,0,0,0,1,0), 'a = e_Y1 (sing pt of C=0)'),
              ((1,0,0,1,0,0), 'a = e_X0+e_Y0'), ((1,0,0,-1,0,0), 'a = e_X0-e_Y0'), ((1,1,0,0,0,0), 'a=e_X0+e_X1'), ((0,0,1,1,0,0), 'a=e_X2+e_Y0')]:
    run(polar(a), f'polar F_a, {nm}', 'polar-pencil')
for i, j in [(0, 3), (2, 2), (0, 0), (3, 3), (1, 4), (2, 5)]:
    run(M(Fi[i], Fi[j]), f'F_{core.VARS[i]}*F_{core.VARS[j]}', 'polar-products', exact=False)
cands = [
 (M(M(X2, X2), H), 'X2^2*H', 'pieces-8'), (M(M(X2, X2), Hp), "X2^2*H'", 'pieces-8'), (M(M(X2, X2), M(C, Zp)), "X2^2*C*Z'", 'pieces-8'),
 (M(C, M(C, M(X2, X2))), 'X2^2*C^2', 'pieces-8'), (M(Zp, M(Zp, M(X2, X2))), "X2^2*Z'^2", 'pieces-8'),
 (M(X2, D8), 'X2*D8', 'D8-products'), (M(M(X2, X2), D8), 'X2^2*D8', 'D8-products'), (M(C, D8), 'C*D8', 'D8-products'), (M(Zp, D8), "Z'*D8", 'D8-products'),
 (M(Z, H), 'X2^3*H', 'pieces-9'), (M(Z, Hp), "X2^3*H'", 'pieces-9'), (M(C, Hp), "C*H'", 'pieces-9'), (M(Zp, Hp), "Z'*H'", 'pieces-9'), (M(C, M(C, C)), 'C^3', 'pieces-9'),
 (M(Zp, M(Zp, Zp)), "Z'^3", 'pieces-9'), (M(C, M(C, Zp)), "C^2*Z'", 'pieces-9'), (M(C, M(Zp, Zp)), "C*Z'^2", 'pieces-9'),
 (M(H, M(C, Zp)), "H*C*Z'", 'pieces-12'), (M(H, Hp), "H*H'", 'pieces-12'), (M(H, M(C, C)), 'H*C^2', 'pieces-12'), (M(Hp, M(C, Zp)), "H'*C*Z'", 'pieces-12'),
 (M(H, M(Zp, Zp)), "H*Z'^2", 'pieces-12'), (M(H, H), 'H^2', 'pieces-12'), (M(Hp, Hp), "H'^2", 'pieces-12'),
]
for D, name, fam in cands:
    run(D, name, fam)
