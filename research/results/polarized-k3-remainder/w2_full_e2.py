"""Job: full reduced-E2 variable-projection search (Levenberg-Marquardt on the sphere) over the Jacobian-type
family D = sum_i a_i F_i of degree d and Y-parity pD, F-multiples projected out.

usage: python3 w2_full_e2.py D PD NRESTARTS SEED [--maxiter 60] [--tag T] [--init file.jsonl]
Writes w2-full-d{D}-p{PD}-{tag}.jsonl: one line per restart (f, rel, iters, dist to spurious g*h family, a-vector),
rational reconstruction + exact mod-p confirmation whenever rel < 1e-10."""
import sys, time, json, argparse, os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
if '--procs' in sys.argv:
    _np = int(sys.argv[sys.argv.index('--procs') + 1])
    if _np > 1:
        for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
            os.environ[_v] = '1'
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np
import multiprocessing as mp
import w2lib as L
from fractions import Fraction as Fr
from hodge import core, fixed_d
from hodge.core import NV, mul, add, sub, scale, shift, pderiv, monomials

ap = argparse.ArgumentParser()
ap.add_argument('d', type=int); ap.add_argument('pD', type=int); ap.add_argument('n', type=int); ap.add_argument('seed', type=int)
ap.add_argument('--maxiter', type=int, default=60); ap.add_argument('--tag', default='a'); ap.add_argument('--init', default=None)
ap.add_argument('--rel-confirm', type=float, default=1e-10); ap.add_argument('--deflate', action='store_true')
ap.add_argument('--procs', type=int, default=1); ap.add_argument('--rows', type=int, default=0, help='random subset of reduced rows (necessary-condition relaxation)'); ap.add_argument('--uvbideg', type=int, default=None, help='restrict to the swap-symmetric (U,V)-bidegree piece (a, d-a)+(d-a, a) of J(F)_d')
args = ap.parse_args()
d, pD = args.d, args.pD
out = f'w2-full-d{d}-p{pD}-{args.tag}.jsonl'
T0 = time.time()
F, A, q = L.load()
if args.uvbideg is None:
    S = L.ReducedE2.jacobian(d, pD)
    basis = S.basis
    pairs = [(Dk, [mul(A, {m: Fr(1)}) if j == i0 else {} for j in range(NV)]) for (i0, m, Dk) in basis]
    # spurious family g*h for classification
    f = L.f_ternary(); g = sub(L.uv_to_xy(f), L.uv_to_xy(L.swapUV(f)))
    Mb = S.Mbasis.toarray()
    spur = []
    for h in monomials(d - 6, (pD + 1) % 2):
        v = S.S_D.vec(mul(g, {h: Fr(1)}))
        a, *_ = np.linalg.lstsq(Mb, v, rcond=None)
        spur.append(S.P.T @ a)
    Qs = np.linalg.qr(np.array(spur).T)[0] if spur else np.zeros((S.dimP, 0))
else:
    fam = L.uv_bidegree_family(d, pD, args.uvbideg)
    pairs = [(Dk, N0) for (_, Dk, N0, _) in fam]
    S = L.ReducedE2(d, pD, pairs)
    S.P = np.eye(S.K); S.dimP = S.K
if args.rows:
    S.subsample_rows(args.rows, seed=args.seed)
    out = out.replace('.jsonl', f'-rows{args.rows}.jsonl')
    basis = [(lab, None, Dk) for (lab, Dk, _, _) in fam]
    Qs = np.zeros((S.dimP, 0))
    out = f'w2-full-d{d}-p{pD}-uv{args.uvbideg}-{args.tag}.jsonl'
def dist_spur(u):
    w = u - Qs @ (Qs.T @ u); return float(np.linalg.norm(w)), w
def deflate(u):
    rho, w = dist_spur(u); rho = max(rho, 1e-12); return rho, w / rho
# timing probe
u0 = np.random.default_rng(0).standard_normal(S.dimP); u0 /= np.linalg.norm(u0)
t1 = time.time(); asm = S.assemble(S.P @ u0); r0 = S.residual(S.P @ u0, asm, with_jac=True); t_eval = time.time() - t1
L.log_line(out, {'experiment': 'full-e2-setup', 'd': d, 'pD': pD, 'uvbideg': args.uvbideg, 'row_subsample': args.rows or None, 'K': S.K, 'dimP': S.dimP, 'rows': S.S_B.n, 'red_rows': S.nred,
                 'nc': S.nc, 'nc_eff': S.nc_eff, 'nQ': S.S_Q.n, 'build_seconds': round(S.build_time, 1),
                 'eval_seconds': round(t_eval, 2), 'probe_f': r0['f'], 'probe_rel': r0['rel'], 'probe_rank': r0['rank'],
                 'spurious_family_dim': int(Qs.shape[1]), 'maxiter': args.maxiter, 'seed': args.seed, 'deflate': args.deflate})
print(json.dumps({'K': S.K, 'dimP': S.dimP, 'red_rows': S.nred, 'nc_eff': S.nc_eff, 'eval_seconds': round(t_eval, 2)}), flush=True)

def resfun(u):
    a = S.P @ u
    asm = S.assemble(a)
    res = S.residual(a, asm, with_jac=True)
    return res['f'], res['r'], res['J'] @ S.P

def try_confirm(a, rel, restart):
    """Rational reconstruction of a (in the basis coordinates) and exact mod-p re-solve."""
    tried = 0
    for max_den in (100, 1000, 10000, 100000):
        for frac in L.rational_candidates(a, max_den=max_den, tol=1e-6 if max_den <= 1000 else 1e-8):
            tried += 1
            D = {}
            N0 = [{} for _ in range(NV)]
            for k, (Dk, N0k) in enumerate(pairs):
                if frac[k]:
                    D = add(D, scale(Dk, frac[k]))
                    for j in range(NV):
                        if N0k[j]: N0[j] = add(N0[j], scale(N0k[j], frac[k]))
            if not D:
                continue
            divF, _ = L.jacobian.divisible_by(F, D)
            rec = {'experiment': 'full-e2-exact-confirm', 'd': d, 'pD': pD, 'restart': restart, 'rel': rel, 'max_den': max_den,
                   'F_divides_D': bool(divF), 'D_terms': len(D)}
            if divF:
                rec['note'] = 'reconstructed D is an F-multiple (spurious)'; L.log_line(out, rec); continue
            o = S.exact_consistent(D, N0, p=L.P0)
            rec.update({'consistent': o['consistent'], 'rank': o['rank'], 'nullity': o['nullity'], 'p': o['p']})
            if o['consistent']:
                status, data = L.lift_solution(S, D, N0)
                rec['lift_status'] = status
                if status == 'solved':
                    w = {'d': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}
                    json.dump(w, open('witness-candidate.json', 'w'), indent=1)
                    rec['witness'] = w
                    L.log_line(out, rec); print('WITNESS FOUND', flush=True); return True
            L.log_line(out, rec)
            if tried >= 6: return False
    return False

inits = []
if args.init:
    for line in open(args.init):
        r = json.loads(line)
        if 'a' in r: inits.append(np.array(r['a'], float))
LOCK = mp.Lock()
def log_locked(rec):
    with LOCK:
        L.log_line(out, rec)

def one_restart(it):
    t0 = time.time()
    rng = np.random.default_rng(args.seed * 100000 + it)
    if it < len(inits):
        u0 = S.P.T @ inits[it]
    else:
        u0 = rng.standard_normal(S.dimP)
    u0 /= np.linalg.norm(u0)
    res = L.lm_sphere(resfun, u0, maxiter=args.maxiter, deflate=deflate if args.deflate else None, ftol=1e-30)
    u = res['u']; a = S.P @ u
    asm = S.assemble(a); fin = S.residual(a, asm)
    rel = fin['rel']; rho, _ = dist_spur(u)
    rec = {'experiment': 'full-e2-restart', 'd': d, 'pD': pD, 'uvbideg': args.uvbideg, 'row_subsample': args.rows or None, 'seed': args.seed, 'restart': it,
           'init': 'file' if it < len(inits) else 'random', 'f': fin['f'], 'rel': rel, 'bnorm2': fin['bnorm2'], 'rank_Mc': fin['rank'],
           'iters': res['iters'], 'nfev': res['nfev'], 'dist_to_spurious': rho, 'seconds': round(time.time() - t0, 1),
           'hist': [float(x) for x in res['hist'][-5:]], 'a': [float(x) for x in a]}
    try:
        rec['uv_bidegree_profile'] = L.uv_bidegree_profile({e: float(c) for e, c in L.a_to_D(basis, a).items()}, d) if args.uvbideg is None else {str(args.uvbideg): 1.0}
    except Exception as ex:
        rec['uv_bidegree_profile'] = str(ex)
    log_locked(rec)
    print(json.dumps({k: rec[k] for k in ('restart', 'f', 'rel', 'iters', 'nfev', 'dist_to_spurious', 'seconds')}), flush=True)
    found = False
    if rel < args.rel_confirm:
        with LOCK:
            found = try_confirm(a, rel, it)
    return rel, found

if args.procs > 1:
    with mp.get_context('fork').Pool(args.procs) as pool:
        best = None
        for rel, found in pool.imap_unordered(one_restart, range(args.n)):
            best = rel if best is None else min(best, rel)
            if found:
                pool.terminate(); break
else:
    best = None
    for it in range(args.n):
        rel, found = one_restart(it)
        best = rel if best is None else min(best, rel)
        if found: break
print('done', round(time.time() - T0, 1), 's; best rel', best, flush=True)
