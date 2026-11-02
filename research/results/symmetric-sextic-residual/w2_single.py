"""Single-D character-reduced solves (mod p) of every E1-kernel class (and random combinations with J) at
d = 8..12 in the four 1-dim character sectors, plus structured exact candidates (products of kernel elements,
Hessian blocks, ...).  8 worker processes.  Output: w2/<job>.log.jsonl"""
import sys, os, time, json, argparse
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from w2lib import *
from charsolve import e1_kernel_char, J_char_basis
import multiprocessing as mp
ap = argparse.ArgumentParser(); ap.add_argument('--job', default='single'); ap.add_argument('--workers', type=int, default=8)
ap.add_argument('--part', default='kernels,struct'); ap.add_argument('--dmin', type=int, default=8); ap.add_argument('--dmax', type=int, default=12)
args = ap.parse_args()
JOBLOG = f'{RES}/w2/{args.job}.log.jsonl'; os.makedirs(f'{RES}/w2', exist_ok=True)
F, A, q = load(); p = P0
QX, QY, e1X, e1Y, e3X, e3Y = e_k(X, 2), e_k(Y, 2), e_k(X, 1), e_k(Y, 1), e_k(X, 3), e_k(Y, 3)
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
DX, DY = delta(X), delta(Y)
def e2sq(V):  # e2(V^2) = e2^2 - 2 e1 e3
    return sub(pw(e_k(V, 2), 2), sc(mul(e_k(V, 1), e_k(V, 3)), 2))
D7, D7Y = mul(DX, e2sq(X)), mul(DY, e2sq(Y))
def quint(V):  # e2 e3 - 7 e1^3 e2 + 4 e1^5
    e1, e2, e3 = e_k(V, 1), e_k(V, 2), e_k(V, 3)
    return add(sub(mul(e2, e3), sc(mul(pw(e1, 3), e2), 7)), sc(pw(e1, 5), 4))
T8, T8Y = mul(DX, quint(X)), mul(DY, quint(Y))
Fi = [core.pderiv(F, i) for i in range(6)]
def radial(idx):
    out = {}
    for i in idx: out = add(out, mul(var(i), Fi[i]))
    return out
EXF, EYF = radial((0, 1, 2)), radial((3, 4, 5))
def det3(M):
    t = {}
    for (a, b, c), s in (((0, 1, 2), 1), ((0, 2, 1), -1), ((1, 0, 2), -1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((2, 1, 0), -1)):
        t = add(t, sc(mul(mul(M[0][a], M[1][b]), M[2][c]), s))
    return t
H = [[core.pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
HXX = det3([[H[i][j] for j in range(3)] for i in range(3)])
HYY = det3([[H[i][j] for j in range(3, 6)] for i in range(3, 6)])
HXY = det3([[H[i][j] for j in range(3, 6)] for i in range(3)])
STRUCT = [
    ('D7', D7, 'sgnX'), ('D7Y', D7Y, 'sgnY'), ('T8=DX*(e2e3-7e1^3e2+4e1^5)', T8, 'sgnX'), ('T8Y', T8Y, 'sgnY'),
    ('D7*e1X', mul(D7, e1X), 'sgnX'), ('D7*e1Y', mul(D7, e1Y), 'sgnX'), ('D7*QX', mul(D7, QX), 'sgnX'), ('D7*QY', mul(D7, QY), 'sgnX'),
    ('D7*DX', mul(D7, DX), 'triv'), ('D7*DY', mul(D7, DY), 'sgnsgn'), ('D7Y*DX', mul(D7Y, DX), 'sgnsgn'),
    ('D7*D7Y', mul(D7, D7Y), 'sgnsgn'), ('T8*D7Y', mul(T8, D7Y), 'sgnsgn'), ('D7*T8Y', mul(D7, T8Y), 'sgnsgn'),
    ('D7*T8', mul(D7, T8), 'triv'), ('e2sqX*e2sqY', mul(e2sq(X), e2sq(Y)), 'triv'), ('quintX*quintY', mul(quint(X), quint(Y)), 'triv'),
    ('e2sqX*quintY', mul(e2sq(X), quint(Y)), 'triv'), ('quintX*e2sqY', mul(quint(X), e2sq(Y)), 'triv'),
    ('HXX', HXX, 'triv'), ('HYY', HYY, 'triv'), ('HXY', HXY, 'sgnsgn'), ('HXX+HYY', add(HXX, HYY), 'triv'), ('HXX-HYY', sub(HXX, HYY), 'triv'),
    ('EXF*EYF', mul(EXF, EYF), 'triv'), ('EXF*D7Y', mul(EXF, D7Y), 'sgnY'), ('EYF*D7', mul(EYF, D7), 'sgnX'),
    ('DX^2*DY^2', mul(pw(DX, 2), pw(DY, 2)), 'triv'), ('T8*T8Y', mul(T8, T8Y), 'sgnsgn'),
]

def run_task(task):
    kind, label, chi, D, meta = task
    G = GROUPS[chi]()
    t0 = time.time()
    try:
        if not G.is_semi_invariant(D):
            return dict(meta, job=args.job, family=kind, label=label, chi=chi, d=core.degree(D), note='not semi-invariant', seconds=0)
        if core.parity_of(D) is None:
            return dict(meta, job=args.job, family=kind, label=label, chi=chi, d=core.degree(D), note='mixed parity', seconds=0)
        r = solve_single(F, A, q, D, G, p, lift=True)
        rec = dict(meta, job=args.job, family=kind, label=label, chi=chi, d=core.degree(D), pD=core.parity_of(D),
                   e1=r['e1'], consistent=r['consistent'], rank=r['rank'], nullity=r['nullity'], rows=r['rows'], cols=r['cols'],
                   seconds=round(time.time() - t0, 1))
        if r['consistent']:
            rec['lift_status'] = r.get('lift_status')
            if 'N' in r:
                ok, info = core.verify_witness(F, A, q, D, r['N'], r['Q'])
                rec['verified'] = bool(ok)
                if ok:
                    w = witness_json(D, r['N'], r['Q']); w.update({'instance': LANE, 'family': f'{kind}:{label}', 'engine': ENGINE})
                    path = f'{RES}/witness-candidate.json'
                    json.dump(w, open(path, 'w'), indent=1)
                    rec['WITNESS'] = path
                    print('WITNESS FOUND', path, flush=True)
        return rec
    except Exception as ex:
        return dict(meta, job=args.job, family=kind, label=label, chi=chi, error=repr(ex)[:300], seconds=round(time.time() - t0, 1))

tasks = []
parts = args.part.split(',')
if 'struct' in parts:
    for label, D, chi in STRUCT:
        tasks.append(('structured-single', label, chi, D, {}))
if 'kernels' in parts:
    for chi in ('triv', 'sgnX', 'sgnY', 'sgnsgn'):
        G = GROUPS[chi]()
        for d in range(args.dmin, args.dmax + 1):
            for pD in (0, 1):
                t0 = time.time()
                info, polys = e1_kernel_char(F, A, q, d, pD, G, p)
                print('kernel', chi, d, pD, info, f'{time.time()-t0:.1f}s', flush=True)
                if info['ker'] == 0:
                    continue
                Jp = J_char_basis(F, d, pD, G)
                for i, P in enumerate(polys):
                    tasks.append(('kernel-class', f'{chi}/d{d}/p{pD}/class{i}', chi, P, {'d': d, 'pD': pD, 'class': i, 'n_classes': info['ker']}))
                rng = np.random.default_rng(100 * d + pD)
                for trial in range(3):
                    D = {}
                    for P in polys:
                        D = add(D, sc(P, Fraction(int(rng.integers(-3, 4)))))
                    if Jp:
                        for j in rng.choice(len(Jp), size=min(4, len(Jp)), replace=False):
                            D = add(D, sc(Jp[j], Fraction(int(rng.integers(-3, 4)))))
                    if D:
                        tasks.append(('kernel-random-combo', f'{chi}/d{d}/p{pD}/combo{trial}', chi, D, {'d': d, 'pD': pD, 'trial': trial}))
print(f'{len(tasks)} tasks', flush=True)
# longest first (highest degree) to balance
tasks.sort(key=lambda t: -core.degree(t[3]))
ctx = mp.get_context('fork')
with ctx.Pool(args.workers) as pool:
    for rec in pool.imap_unordered(run_task, tasks):
        wlog(rec, JOBLOG)
        print(json.dumps({k: v for k, v in rec.items() if k in ('family', 'label', 'chi', 'd', 'pD', 'e1', 'consistent', 'rank', 'nullity', 'rows', 'cols', 'seconds', 'note', 'error', 'WITNESS')}), flush=True)
print('ALL DONE', flush=True)
