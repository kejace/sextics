"""Hessian (mod-J) relaxation of E2 per sigma-sector, numerical landscape + exact Singular dimension check.

  python3 -u w2_hess.py --tasks "6:0:0-6,6:1:0-6,7:0:0-6,7:1:0-6,8:0:0-6,8:1:0-6" --restarts 2000 --procs 8 [--singular] [--family J|KJ]

For D = sum a_k D_k with E1 lifts N_k, the necessary condition for E2 is  B_{D(a)}(N0(a)) in J(F)_{2d}
(for Jacobian-type D this is  A * (w^T Hess(F) w) in J).  Numerically: f(a) = |L^T b(a)|^2 with L an orthonormal
basis of the left null space of the J-block (degree 2d, Y-odd, weight 3+2k) and b(a) = -B_{D(a)}(N0(a)),
minimized over the unit sphere in the complement of the F-multiples.  Exact: the same quadrics mod 32003 in a
rational chart complementary to the F-multiples; Singular dim(std(I)) == 0 means no sigma-eigen denominator of
this degree/parity/weight admits any certificate (mod-p emptiness of a homogeneous system implies emptiness over Q).
"""
import os, sys, time, json, math, argparse, pickle, subprocess, traceback
os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize
import w2lib as L
from hodge import core
from hodge.core import NV, shift, mul, add, sub, scale, pderiv
import kcc2
from kcc import W, wt, wmonomials
from kcc2 import Ctx, PI
from multiprocessing import Pool
from flint import nmod_mat

LANE = L.LANE
F, A, q = L.load()
PSING = 32003
import os as _os
PSING = int(_os.environ.get('W2_SING_PRIME', PSING))


def build_family(d, parD, k, family, pkl):
    fam = []
    if family in ('KJ', 'K'):
        store = pickle.load(open(pkl, 'rb'))
        info = store[str((d, parD, k))]
        reps, lifts = info['reps'], info['lifts']
        if (len(reps) < info.get('nontrivial', 0) or any(N is None for N in lifts)) and info.get('nontrivial', 0) > 0:
            # incomplete rational representatives (large height): use the orthonormal floating-point kernel basis
            if info.get('reps_float') and len(info['reps_float']) == info['nontrivial']:
                reps, lifts = info['reps_float'], info['lifts_float']
            else:
                reps, lifts, finfo = L.float_kernel_reps(F, A, d, parD, k)
            if len(reps) != info['nontrivial']:
                raise RuntimeError(f'float kernel reps {len(reps)} != nontrivial {info["nontrivial"]}')
        for idx, (D, N) in enumerate(zip(reps, lifts)):
            if N is None:
                N, relf = L.float_e1_lift(F, A, D, restrict=True)
            fam.append((D, N, ('K', idx)))
    if family in ('KJ', 'J'):
        fam += L.jacobian_family(F, A, d, parD, k, restrict=True)
    return fam


def ring_dim_modp(d, k):
    a = 3; pA = 1; tk = (a + 2 * k) % 7
    ctx = Ctx(F, A, L.P0)
    tgt = wmonomials(2 * d, pA, tk)
    jk, jp = kcc2.jac_block_cols(ctx, 2 * d, pA, tk)
    MJ, _ = kcc2.dense_from_cols(jp, None, L.P0, {m: i for i, m in enumerate(tgt)})
    return len(tgt) - MJ.rank()


def singular_check(d, parD, k, fam, tag, timeout=1800):
    """Exact mod-32003 quadric system of the relaxation in a chart complementary to F-multiples; Singular dim."""
    p = PSING
    ctx = Ctx(F, A, p)
    a_, pA = 3, 1
    tk = (a_ + 2 * k) % 7
    tgt = wmonomials(2 * d, pA, tk)
    rows = {m: i for i, m in enumerate(tgt)}
    jk, jp = kcc2.jac_block_cols(ctx, 2 * d, pA, tk)
    MJ, _ = kcc2.dense_from_cols(jp, None, p, rows)
    Lm, nl = MJ.transpose().nullspace()
    Lt = nmod_mat(nl, len(rows), [int(Lm[r, c]) for c in range(nl) for r in range(len(rows))], p)
    K = len(fam)
    Dp = [core.reduce_mod(Dk, p) for Dk, _, _ in fam]
    for _, Nk, _ in fam:
        for Ni in Nk:
            for v in Ni.values():
                if isinstance(v, float):
                    raise RuntimeError('Singular check needs exact (rational) lifts')
    Np = [[core.reduce_mod(Ni, p) for Ni in Nk] for _, Nk, _ in fam]
    divN = []
    for l in range(K):
        acc = {}
        for i in range(NV):
            if Np[l][i]:
                acc = add(acc, pderiv(Np[l][i], i, p), p)
        divN.append(acc)
    Dgrad = [[pderiv(Dp[kk], i, p) for i in range(NV)] for kk in range(K)]
    # coefficient tensor c[m][k][l]
    C = {}
    for kk in range(K):
        for l in range(K):
            v = mul(Dp[kk], divN[l], p)
            for i in range(NV):
                if Np[l][i]:
                    v = sub(v, mul(Np[l][i], Dgrad[kk][i], p), p)
            vec = [0] * len(rows)
            for e, c in v.items():
                vec[rows[e]] = c % p
            colv = nmod_mat(len(rows), 1, vec, p)
            q = Lt * colv
            for m in range(nl):
                c = int(q[m, 0])
                if c:
                    C.setdefault(m, {})[(kk, l)] = c
    # chart: drop a_{(0, x*g)} for every g of degree d-6 in the sector (F-multiples), only for J-type members
    drop = set()
    if d >= 6:
        tagidx = {fam[i][2]: i for i in range(K)}
        for g in wmonomials(d - 6, parD, k):
            xg = tuple(g[j] + (1 if j == 0 else 0) for j in range(6))
            key = ('J', 0, xg)
            if key in tagidx:
                drop.add(tagidx[key])
    keep = [i for i in range(K) if i not in drop]
    var = {i: f'u{j}' for j, i in enumerate(keep)}
    polys = []
    for m, terms in C.items():
        parts = []
        for (kk, l), c in terms.items():
            if kk in drop or l in drop:
                continue
            parts.append(f'{c}*{var[kk]}*{var[l]}')
        if parts:
            polys.append('+'.join(parts))
    path = f'{LANE}/w2-hess-{tag}-d{d}-p{parD}-k{k}.sing'
    with open(path, 'w') as f:
        f.write(f'ring r = {p}, ({",".join(var[i] for i in keep)}), dp;\n')
        f.write('option(redSB);\n')
        f.write('ideal I = ' + ',\n'.join(polys) + ';\n')
        f.write('ideal G = std(I);\n')
        f.write('int dm = dim(G);\n')
        f.write('print("DIM " + string(dm) + " GBSIZE " + string(size(G)) + " VARS " + string(nvars(r)) + " NEQ " + string(size(I)));\n')
        f.write('if (dm == 0) { print("VDIM " + string(vdim(G))); }\n')
        f.write('quit;\n')
    t = time.time()
    try:
        r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=timeout)
        outp = r.stdout.strip().splitlines()
        res = {'singular': outp[-2:] if outp else r.stderr[-300:], 'seconds': round(time.time() - t, 1), 'nvars': len(keep), 'neq': len(polys), 'nl': nl}
        for line in outp:
            if line.startswith('DIM'):
                res['dim'] = int(line.split()[1]); res['gbsize'] = int(line.split()[3])
            if line.startswith('VDIM'):
                res['vdim'] = int(line.split()[1])
    except subprocess.TimeoutExpired:
        res = {'singular': 'timeout', 'seconds': timeout, 'nvars': len(keep), 'neq': len(polys), 'nl': nl}
    return res


def run_task(args):
    d, parD, k, family, restarts, seed, do_sing, tag, pkl, sing_timeout = args
    t0 = time.time()
    try:
        fam = build_family(d, parD, k, family, pkl)
        if not fam:
            L.wlog({'family': f'hess-relax-{family}', 'd': d, 'parD': parD, 'k': k, 'K': 0, 'tag': tag, 'note': 'empty family'})
            return {'task': (d, parD, k), 'empty': True}
        kv = L.KoszulVarPro(F, A, fam, restrict=True, k=k, reduce_F=True, exclude_F_multiples=True)
        rx = L.Relax(kv)
        nl = rx.nl; rank = rx.jrank; s = rx.sv_gap
        ring_dim_p = ring_dim_modp(d, k)
        out = rx.search(restarts, seed)
        rels = np.array([c['rel'] for c in out])
        rec = {'family': f'hess-relax-{family}', 'd': d, 'parD': parD, 'k': k, 'K': kv.K, 'dimP': kv.dimP, 'rows': kv.rows,
               'nl_float': nl, 'ring_dim_modp': ring_dim_p, 'jrank_float': rank, 'sv_gap': list(s),
               'n_equations': nl, 'n_unknowns': kv.dimP, 'restarts': restarts, 'best_rel': float(rels.min()),
               'median_rel': float(np.median(rels)), 'n_lt_1e-6': int(np.sum(rels < 1e-6)), 'n_lt_1e-10': int(np.sum(rels < 1e-10)),
               'seed': seed, 'tag': tag, 'seconds': round(time.time() - t0, 1),
               'note': 'mod-J (Hessian) relaxation of E2: |proj_(R/J) B_D(N0)|^2 / |B_D(N0)|^2 minimized over the sphere in the F-complement'}
        if do_sing:
            try:
                rec['singular'] = singular_check(d, parD, k, fam, tag, timeout=sing_timeout)
            except Exception as ex:
                rec['singular'] = {'error': repr(ex)}
        L.wlog(rec)
        st_path = f'{LANE}/w2-hess-{tag}-d{d}-p{parD}-k{k}.json'
        json.dump({'rec': rec, 'top': [{'rel': c['rel'], 'a': [float(x) for x in c['a']]} for c in out[:5]]}, open(st_path, 'w'))
        return {'task': (d, parD, k), 'best_rel': rec['best_rel'], 'n_eq': nl, 'n_unk': kv.dimP, 'singular': rec.get('singular', {}).get('dim', None) if do_sing else None, 'seconds': rec['seconds']}
    except Exception as ex:
        L.wlog({'family': f'hess-relax-{family}', 'd': d, 'parD': parD, 'k': k, 'tag': tag, 'error': repr(ex), 'traceback': traceback.format_exc()[-1500:]})
        return {'task': (d, parD, k), 'error': repr(ex)}


def main():
    from w2_varpro import parse_tasks
    ap = argparse.ArgumentParser()
    ap.add_argument('--tasks', required=True)
    ap.add_argument('--family', default='J')
    ap.add_argument('--restarts', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--singular', action='store_true')
    ap.add_argument('--sing-timeout', type=int, default=1800)
    ap.add_argument('--tag', default='hess')
    ap.add_argument('--pkl', default=f'{LANE}/w2_sectors.pkl')
    a = ap.parse_args()
    tasks = parse_tasks(a.tasks)
    jobs = [(d, parD, k, a.family, a.restarts, a.seed + 17 * i, a.singular, a.tag, a.pkl, a.sing_timeout) for i, (d, parD, k) in enumerate(tasks)]
    t0 = time.time()
    with Pool(a.procs) as pool:
        for r in pool.imap_unordered(run_task, jobs):
            print(json.dumps(r, default=str), f'[{time.time()-t0:.0f}s]', flush=True)
    print('HESS_DONE', flush=True)


if __name__ == '__main__':
    main()
