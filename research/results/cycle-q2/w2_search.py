"""Wave-2 joint search driver (lane cycle-q2): sector-restricted Koszul-reduced complex varpro over linear
families of E1-admissible denominators, with incremental logging, rational reconstruction and exact confirmation.

Usage:
  python3 w2_search.py --d 8 --family kernel --keys all --restarts 40 --maxiter 150 --nproc 8 --seed 0 --tag k8
Families:
  kernel         K_kappa = {D in (S_d)_kappa : A D in J}           (E1 kernel + J, the whole admissible sector space)
  jac            (J_d)_kappa = span{ n F_i }                        (Jacobian-type w-space)
  cone4:i,j,k,l  K_kappa ∩ forms in the 4 variables i,j,k,l        (cone with vertex the complementary line)
  cone3:i,j,k    K_kappa ∩ forms in 3 variables                    (cone with vertex the complementary plane)
  mult4:i,j      K_kappa ∩ {mult >= 4 at the 6 points of X ∩ line(i,j)}   (also mult5:i,j)
  sing2:i,j,k    K_kappa ∩ {vanishing to order 2 along the Fermat curve X ∩ plane(i,j,k)}  (I_C^2)
  full           all of (S_d)_kappa (no E1 restriction; relaxation check only)
Options: --mode complex|real, --key-budget SECONDS, --min-K, --max-K, --rel-thr (candidate threshold), --nontrivial
(only sectors with dimK > dimJ), --exclude-keys.
Every batch of restarts appends one line to wave2-search.log.jsonl; candidates with rel < rel-thr are reconstructed
and re-solved exactly at once; witnesses go to witness-candidate.json.
"""
import os as _os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "FLINT_NUM_THREADS"):
    _os.environ[_v] = "1"
import sys, os, time, json, argparse, traceback
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import numpy as np
import multiprocessing as mp
import w2lib as W
from hodge import core
from fractions import Fraction

LANE = '/tmp/hodge/results/cycle-q2'
VARS = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']


def jet_matrix(Dm, point, order):
    """complex matrix (n_alpha x nD): d^alpha m (point) for |alpha| <= order"""
    import itertools
    alphas = [a for a in itertools.product(range(order + 1), repeat=6) if sum(a) <= order]
    rows = []
    for al in alphas:
        row = np.zeros(len(Dm), dtype=complex)
        for j, m in enumerate(Dm):
            coef = 1.0
            ok = True
            ex = []
            for k in range(6):
                if al[k] > m[k]:
                    ok = False; break
                # falling factorial m_k (m_k-1) ... (m_k - al_k + 1)
                for t in range(al[k]):
                    coef *= (m[k] - t)
                ex.append(m[k] - al[k])
            if not ok:
                continue
            val = coef
            for k in range(6):
                if ex[k]:
                    val *= point[k] ** ex[k]
            row[j] = val
        rows.append(row)
    return np.array(rows)


def family_basis(S, d, kappa, kind, q=2):
    """returns (Dm, Bmat, info) with Bmat an orthonormal real basis (nD x K) of the family inside (S_d)_kappa"""
    Dm, Kb, info = W.e1_kernel_numeric(S, d, kappa, q)
    info = dict(info)
    nD = len(Dm)
    if nD == 0:
        return Dm, np.zeros((0, 0)), info
    if kind == 'full':
        return Dm, np.eye(nD), info
    if kind == 'kernel':
        return Dm, Kb, info
    if kind.startswith('prod'):
        # products of E1-kernel elements of degree d1 (sector kappa1) and d2 = d - d1 (sector kappa2 = kappa - kappa1):
        # relaxation of {D1*D2} and in particular of the squares {D1^2} (d1 = d2, kappa1 = kappa2 = kappa/2)
        d1 = int(kind[4:]) if len(kind) > 4 else d // 2
        d2 = d - d1
        cols = []
        for kappa1 in S.all_keys(d1):
            kappa2 = S.add(kappa, S.neg(kappa1))
            if d1 == d2 and kappa2 < kappa1:
                continue
            Dm1, K1, _ = W.e1_kernel_numeric(S, d1, kappa1, q)
            if K1.size == 0 or K1.shape[1] == 0:
                continue
            Dm2, K2, _ = W.e1_kernel_numeric(S, d2, kappa2, q)
            if K2.size == 0 or K2.shape[1] == 0:
                continue
            idx = {tuple(int(x) for x in m): j for j, m in enumerate(Dm)}
            P1 = [W.dict_from_vec(Dm1, K1[:, a]) for a in range(K1.shape[1])]
            P2 = [W.dict_from_vec(Dm2, K2[:, b]) for b in range(K2.shape[1])]
            for Pa in P1:
                for Pb in P2:
                    prod = core.mul(Pa, Pb)
                    v = np.zeros(nD)
                    for e, c in prod.items():
                        v[idx[e]] += float(c)
                    cols.append(v)
        if not cols:
            return Dm, np.zeros((nD, 0)), info
        B = W.orth_basis(np.array(cols).T)
        info['n_products'] = len(cols)
        return Dm, B, info
    if kind == 'jac':
        jb = W.jac_basis(S, d, kappa)
        if not jb:
            return Dm, np.zeros((nD, 0)), info
        V = W.orth_basis(W.vecs_from_dicts(Dm, jb))
        info['dimJ'] = V.shape[1]
        return Dm, V, info
    if kind.startswith('cone4:') or kind.startswith('cone3:'):
        vs = [int(x) for x in kind.split(':')[1].split(',')]
        mask = np.ones(6, dtype=bool); mask[vs] = False
        sel = np.all(Dm[:, mask] == 0, axis=1)
        if not np.any(sel):
            return Dm, np.zeros((nD, 0)), info
        E = np.eye(nD)[:, sel]
        B = W.intersect_subspaces(Kb, E)
        info['support_dim'] = int(sel.sum())
        return Dm, B, info
    if kind.startswith('mult4:') or kind.startswith('mult5:'):
        order = 3 if kind.startswith('mult4') else 4
        i, j = [int(x) for x in kind.split(':')[1].split(',')]
        eta = np.exp(1j * np.pi / 6)   # eta^6 = -1
        pt = np.zeros(6, dtype=complex); pt[i] = 1.0; pt[j] = eta
        Jm = jet_matrix(Dm, pt, order)
        M = Jm @ Kb                      # conditions on kernel coordinates
        U, s, Vt = np.linalg.svd(M, full_matrices=True)
        rk = int(np.sum(s > 1e-9 * max(s[0], 1.0))) if len(s) else 0
        Nc = Vt[rk:, :].conj().T          # complex null space in kernel coords (K x k)
        if Nc.shape[1] == 0:
            return Dm, np.zeros((nD, 0)), info
        Rb = W.orth_basis(np.hstack([Nc.real, Nc.imag]))
        info['jet_rank'] = rk; info['complex_null'] = int(Nc.shape[1]); info['real_null'] = int(Rb.shape[1])
        return Dm, Kb @ Rb, info
    if kind.startswith('sing2:'):
        # vanishing to order 2 along C = X ∩ plane(i,j,k): I_C = (other vars, G) with G = x_i^6 + x_j^6 + x_k^6.
        # I_C^2 in degree d: (others)^2 S + (others) G S + G^2 S.  Equivalent linear conditions: write D = D0 + D1 + D2..
        # where D_t has total degree t in the other variables; need D0 in (G^2), D1 in (G) (coefficient-wise).
        vs = [int(x) for x in kind.split(':')[1].split(',')]
        others = [v for v in range(6) if v not in vs]
        odeg = Dm[:, others].sum(axis=1)
        G = {tuple(6 if t == v else 0 for t in range(6)): Fraction(1) for v in vs}
        conds = []
        # D0: monomials with odeg == 0 must form a multiple of G^2 (degree d, in 3 vars): condition = orthogonal to complement
        for t, Gt in ((0, core.mul(G, G)), (1, G)):
            sel = np.where(odeg == t)[0]
            if len(sel) == 0:
                continue
            # group by the 'other' monomial part
            keyo = [tuple(Dm[s, o] for o in others) for s in sel]
            groups = {}
            for s, ko in zip(sel, keyo):
                groups.setdefault(ko, []).append(s)
            for ko, idxs in groups.items():
                # the 3-variable part must be divisible by Gt: multiples of Gt of the right degree span a subspace;
                # condition = projection onto the complement of that subspace vanishes
                deg3 = int(Dm[idxs[0]][vs].sum())
                dg = deg3 - 6 * (2 if t == 0 else 1)
                sub = np.zeros((len(idxs), 0))
                if dg >= 0:
                    import itertools
                    mons3 = [m for m in itertools.product(range(dg + 1), repeat=3) if sum(m) == dg]
                    cols = []
                    idxmap = {tuple(int(x) for x in Dm[s][vs]): r for r, s in enumerate(idxs)}
                    for m3 in mons3:
                        shiftm = [0] * 6
                        for a, v in enumerate(vs):
                            shiftm[v] = m3[a]
                        P = core.shift(Gt, tuple(shiftm), 1)
                        col = np.zeros(len(idxs))
                        for e, v in P.items():
                            key3 = tuple(e[v_] for v_ in vs)
                            if key3 in idxmap:
                                col[idxmap[key3]] = float(v)
                        cols.append(col)
                    sub = W.orth_basis(np.array(cols).T) if cols else sub
                comp = np.eye(len(idxs)) - (sub @ sub.T if sub.shape[1] else 0)
                compb = W.orth_basis(comp) if np.linalg.norm(comp) > 1e-12 else np.zeros((len(idxs), 0))
                for c in range(compb.shape[1]):
                    row = np.zeros(nD); row[idxs] = compb[:, c]; conds.append(row)
        if not conds:
            return Dm, Kb, info
        C = np.array(conds)
        M = C @ Kb
        U, s, Vt = np.linalg.svd(M, full_matrices=True)
        rk = int(np.sum(s > 1e-9 * max(s[0], 1.0))) if len(s) else 0
        N = Vt[rk:, :].T
        info['sing_conditions'] = int(rk)
        return Dm, Kb @ N, info
    raise ValueError(kind)


def reconstruct_and_confirm(S, kvp, cand, args, keyinfo):
    """rational reconstruction of a candidate and exact confirmation; returns dict"""
    F, A, q = S.F, S.A, kvp.q
    dv = kvp.D_coeffs(cand['a'])
    v, imrel = W.normalize_phase(dv)
    out = {'imag_rel': imrel, 'tried': []}
    if imrel > 1e-6:
        out['status'] = 'complex-candidate'
        out['D_float'] = {core.to_str({tuple(int(x) for x in m): 1}): [float(c.real), float(c.imag)] for m, c in zip(kvp.Dm, v) if abs(c) > 1e-12}
        return out
    vr = v.real
    cands = W.rational_candidates(vr)
    if not cands:
        rel = W.lll_integer_relation(vr)
        if rel:
            cands = [rel]
    out['n_rational'] = len(cands)
    for fr in cands[:6]:
        D = W.exact_D_from_fracs(kvp.Dm, fr)
        if not D:
            continue
        Ds = core.to_str(D)
        try:
            res = W.exact_confirm(F, A, q, D, do_lift=True)
        except Exception as ex:
            res = {'status': f'error: {ex}'}
        out['tried'].append({'D': Ds[:400], 'status': res.get('status'), 'e1': res.get('e1'), 'consistent_p1': res.get('consistent_p1'),
                             'consistent_p2': res.get('consistent_p2'), 'F_divides_D': res.get('F_divides_D'), 'rank': res.get('rank'), 'nullity': res.get('nullity')})
        if res.get('status') == 'WITNESS':
            wit = {'instance': 'cycle-q2', 'd': core.degree(D), 'D': Ds, 'N': res['N'], 'Q': res['Q'], 'engine': W.ENGINE,
                   'family': args.family, 'key': keyinfo['key'], 'tag': args.tag, 'rel': cand['rel']}
            wfile = 'witness-candidate.json' if not getattr(args, 'target', None) else f'witness-candidate-sibling-q{args.q}.json'
            wit['target'] = core.to_str(A); wit['q'] = q
            with open(os.path.join(LANE, wfile), 'w') as f:
                json.dump(wit, f, indent=1)
            W.log(step='WITNESS', **{k: v for k, v in wit.items() if k not in ('N', 'Q')}, note='WITNESS FOUND: exact certificate verified, F does not divide D')
            out['status'] = 'WITNESS'
            return out
        if res.get('status') in ('consistent-but-lift-failed', 'spurious-F-multiple'):
            out['status'] = res['status']
            return out
    out['status'] = 'no-rational-match' if not cands else 'rational-candidates-inconsistent'
    return out


def load_target(args):
    F, A, q, rec = core.load_instance('cycle-q2')
    if getattr(args, 'target', None):
        A = core.parse(args.target)
        q = args.q
    return F, A, q


def subfamily_specs(args):
    specs = []
    import itertools
    for kind in args.subfamilies.split(','):
        kind = kind.strip()
        if kind == 'cone4':
            specs += ['cone4:' + ','.join(map(str, c)) for c in itertools.combinations(range(6), 4)]
        elif kind == 'cone3':
            specs += ['cone3:' + ','.join(map(str, c)) for c in itertools.combinations(range(6), 3)]
        elif kind == 'mult4':
            specs += ['mult4:' + ','.join(map(str, c)) for c in itertools.combinations(range(6), 2)]
        elif kind == 'mult5':
            specs += ['mult5:' + ','.join(map(str, c)) for c in itertools.combinations(range(6), 2)]
        elif kind == 'sing2':
            specs += ['sing2:' + ','.join(map(str, c)) for c in itertools.combinations(range(6), 3)]
        elif kind:
            specs.append(kind)
    return specs


def run_key(task):
    args, kappa, part = task if len(task) == 3 else (task[0], task[1], 0)
    if args.family == 'multi':
        return run_key_multi(task)
    F, A, q = load_target(args)
    S = W.Sectors(F, A)
    d = args.d
    t0 = time.time()
    try:
        Dm, Bmat, info = family_basis(S, d, kappa, args.family, q)
    except Exception as ex:
        W.log(step='search', family=args.family, d=d, key=kappa, tag=args.tag, error=f'family: {ex}', note=traceback.format_exc()[-800:])
        return {'key': kappa, 'error': str(ex)}
    K = Bmat.shape[1] if Bmat.size else 0
    par = S.parity_of_key(d, kappa)
    base = {'step': 'search', 'family': args.family, 'd': d, 'key': kappa, 'part': part, 'parity': par, 'tag': args.tag, 'K': int(K), 'q': q,
            'target': core.to_str(A) if getattr(args, 'target', None) else 'A',
            'dimS': int(info.get('dimS', 0)), 'dimK': int(info.get('dimK', 0)), 'mode': args.mode, 'method': args.method}
    if K == 0 or K < args.min_K or (args.max_K and K > args.max_K):
        W.log(**base, restarts=0, note='family empty or outside K range; skipped', seconds=round(time.time() - t0, 1))
        return {'key': kappa, 'K': K, 'skipped': True}
    try:
        kvp = W.KoszulVarPro(S, d, kappa, Dm, Bmat, q=q, verbose=False)
    except Exception as ex:
        W.log(**base, error=f'build: {ex}', note=traceback.format_exc()[-800:])
        return {'key': kappa, 'error': str(ex)}
    base.update({'rows': kvp.R, 'cols': kvp.ncols, 'nC': kvp.nC, 'nQ': kvp.nQ, 'e1_residual': kvp.e1_residual, 'build_s': round(kvp.build_time, 1)})
    rng = np.random.default_rng(args.seed * 100003 + kappa * 101 + part)
    done = 0
    best = None
    rels = []
    cand_file = os.path.join(LANE, f'w2cand_{args.tag}_d{d}.jsonl')
    n_cand = 0
    tk0 = time.time()
    batch = max(1, min(args.batch, args.restarts))
    while done < args.restarts:
        if args.key_budget and time.time() - tk0 > args.key_budget:
            break
        n = min(batch, args.restarts - done)
        seed = int(rng.integers(0, 2 ** 31 - 1))
        tb = time.time()
        if args.method == 'lm':
            cands = kvp.search_lm(n_restarts=n, seed=seed, max_nfev=args.maxiter, complex_mode=(args.mode == 'complex'), stop_rel=args.rel_thr * 1e-3)
        else:
            cands = kvp.search(n_restarts=n, seed=seed, maxiter=args.maxiter, complex_mode=(args.mode == 'complex'), stop_rel=args.rel_thr * 1e-3)
        done += len(cands)
        for c in cands:
            rels.append(c['rel'])
        cands.sort(key=lambda c: (c['spurious'], c['rel']))
        cb = cands[0]
        if best is None or (not cb['spurious'] and cb['rel'] < best['rel']) or (best['spurious'] and not cb['spurious']):
            best = cb
        arr = np.array(rels)
        hist = {'lt1e-4': int(np.sum(arr < 1e-4)), 'lt1e-6': int(np.sum(arr < 1e-6)), 'lt1e-8': int(np.sum(arr < 1e-8)), 'lt1e-10': int(np.sum(arr < 1e-10))}
        rec = dict(base)
        rec.update({'restarts': done, 'batch_seed': seed, 'batch_n': n, 'batch_s': round(time.time() - tb, 1), 'best_rel': float(best['rel']),
                    'best_f': float(best['f']), 'best_spurious': bool(best['spurious']), 'median_rel': float(np.median(arr)), 'min_rel_batch': float(cb['rel']),
                    'hist': hist, 'seconds': round(time.time() - t0, 1), 'nfev': kvp._nfev})
        # candidates
        for c in cands:
            if (c['rel'] < args.rel_thr or c['f'] < 1e-18) and not c['spurious']:
                n_cand += 1
                try:
                    rc = reconstruct_and_confirm(S, kvp, c, args, {'key': kappa})
                except Exception as ex:
                    rc = {'status': f'error {ex}', 'tb': traceback.format_exc()[-600:]}
                dv = kvp.D_coeffs(c['a'])
                with open(cand_file, 'a') as f:
                    f.write(json.dumps({'family': args.family, 'd': d, 'key': kappa, 'rel': c['rel'], 'f': c['f'], 'a_re': [float(x) for x in c['a'].real],
                                        'a_im': [float(x) for x in c['a'].imag], 'Dm': [tuple(int(x) for x in m) for m in kvp.Dm],
                                        'D_re': [float(x) for x in dv.real], 'D_im': [float(x) for x in dv.imag], 'recon': rc}, default=str) + '\n')
                W.log(step='candidate', family=args.family, d=d, key=kappa, tag=args.tag, rel=c['rel'], f=c['f'], recon_status=rc.get('status'),
                      imag_rel=rc.get('imag_rel'), tried=rc.get('tried'), note='numerical candidate below threshold; see w2cand file')
                if rc.get('status') == 'WITNESS':
                    rec['note'] = 'WITNESS FOUND'
                    W.log(**rec)
                    return {'key': kappa, 'K': K, 'best_rel': best['rel'], 'witness': True}
        rec['n_candidates'] = n_cand
        rec['note'] = f'{args.family} d={d} key={kappa} K={K}: {done} restarts, best rel {best["rel"]:.2e}'
        W.log(**rec)
    return {'key': kappa, 'K': K, 'best_rel': float(best['rel']) if best else None, 'restarts': done, 'rows': kvp.R, 'cols': kvp.ncols,
            'seconds': round(time.time() - t0, 1)}


def run_key_multi(task):
    """build the sector kernel system once, then search every structured sub-family inside it"""
    args, kappa, part = task if len(task) == 3 else (task[0], task[1], 0)
    F, A, q = load_target(args)
    S = W.Sectors(F, A)
    d = args.d
    t0 = time.time()
    Dm, Kb, info = W.e1_kernel_numeric(S, d, kappa, q)
    K = Kb.shape[1] if Kb.size else 0
    par = S.parity_of_key(d, kappa)
    if K == 0:
        W.log(step='search', family='multi', d=d, key=kappa, tag=args.tag, K=0, note='empty kernel; skipped')
        return {'key': kappa, 'K': 0, 'skipped': True}
    specs = subfamily_specs(args)
    fams = []
    for spec in specs:
        try:
            Dm2, Bsub, inf2 = family_basis(S, d, kappa, spec, q)
        except Exception as ex:
            W.log(step='search', family=spec, d=d, key=kappa, tag=args.tag, error=f'family: {ex}')
            continue
        k = Bsub.shape[1] if Bsub.size else 0
        if k == 0:
            continue
        C = Kb.T @ Bsub
        # check Bsub lies in the kernel (it should, by construction)
        resid = float(np.linalg.norm(Kb @ C - Bsub) / max(np.linalg.norm(Bsub), 1e-300))
        fams.append((spec, C, k, inf2, resid))
    if not fams:
        W.log(step='search', family='multi', d=d, key=kappa, tag=args.tag, K=int(K), n_subfamilies=0, note='no nonempty structured sub-family in this sector',
              seconds=round(time.time() - t0, 1))
        return {'key': kappa, 'K': K, 'skipped': True}
    try:
        kvp = W.KoszulVarPro(S, d, kappa, Dm, Kb, q=q, verbose=False)
    except Exception as ex:
        W.log(step='search', family='multi', d=d, key=kappa, tag=args.tag, error=f'build: {ex}', note=traceback.format_exc()[-800:])
        return {'key': kappa, 'error': str(ex)}
    rng = np.random.default_rng(args.seed * 100003 + kappa * 101 + part)
    cand_file = os.path.join(LANE, f'w2cand_{args.tag}_d{d}.jsonl')
    results = []
    best_overall = None
    for spec, C, k, inf2, resid in fams:
        sub = kvp.subfamily(C)
        tb = time.time()
        seed = int(rng.integers(0, 2 ** 31 - 1))
        n = args.restarts
        if args.method == 'lm':
            cands = sub.search_lm(n_restarts=n, seed=seed, max_nfev=args.maxiter, complex_mode=(args.mode == 'complex'), stop_rel=args.rel_thr * 1e-3)
        else:
            cands = sub.search(n_restarts=n, seed=seed, maxiter=args.maxiter, complex_mode=(args.mode == 'complex'), stop_rel=args.rel_thr * 1e-3)
        cands.sort(key=lambda c: (c['spurious'], c['rel']))
        cb = cands[0]
        arr = np.array([c['rel'] for c in cands])
        hist = {'lt1e-4': int(np.sum(arr < 1e-4)), 'lt1e-6': int(np.sum(arr < 1e-6)), 'lt1e-8': int(np.sum(arr < 1e-8)), 'lt1e-10': int(np.sum(arr < 1e-10))}
        rec = {'step': 'search', 'family': spec, 'd': d, 'key': kappa, 'part': part, 'parity': par, 'tag': args.tag, 'K': int(k), 'K_sector': int(K),
               'dimS': int(info.get('dimS', 0)), 'dimK': int(info.get('dimK', 0)), 'mode': args.mode, 'method': args.method, 'q': q,
               'rows': kvp.R, 'cols': kvp.ncols, 'nC': kvp.nC, 'nQ': kvp.nQ, 'e1_residual': kvp.e1_residual, 'subfamily_resid': resid,
               'family_info': {kk: (vv if isinstance(vv, (int, float, str)) else str(vv)) for kk, vv in inf2.items() if kk not in ('dimS', 'dimK', 'rankJ', 'gap')},
               'restarts': len(cands), 'batch_seed': seed, 'batch_s': round(time.time() - tb, 1), 'best_rel': float(cb['rel']), 'best_f': float(cb['f']),
               'best_spurious': bool(cb['spurious']), 'median_rel': float(np.median(arr)), 'hist': hist, 'seconds': round(time.time() - t0, 1)}
        n_cand = 0
        for c in cands:
            if (c['rel'] < args.rel_thr or c['f'] < 1e-18) and not c['spurious']:
                n_cand += 1
                try:
                    rc = reconstruct_and_confirm(S, sub, c, args, {'key': kappa})
                except Exception as ex:
                    rc = {'status': f'error {ex}', 'tb': traceback.format_exc()[-600:]}
                dv = sub.D_coeffs(c['a'])
                with open(cand_file, 'a') as f:
                    f.write(json.dumps({'family': spec, 'd': d, 'key': kappa, 'rel': c['rel'], 'f': c['f'], 'a_re': [float(x) for x in c['a'].real],
                                        'a_im': [float(x) for x in c['a'].imag], 'Dm': [tuple(int(x) for x in m) for m in sub.Dm],
                                        'D_re': [float(x) for x in dv.real], 'D_im': [float(x) for x in dv.imag], 'recon': rc}, default=str) + '\n')
                W.log(step='candidate', family=spec, d=d, key=kappa, tag=args.tag, rel=c['rel'], f=c['f'], recon_status=rc.get('status'),
                      imag_rel=rc.get('imag_rel'), tried=rc.get('tried'), note='numerical candidate below threshold; see w2cand file')
                if rc.get('status') == 'WITNESS':
                    rec['note'] = 'WITNESS FOUND'
                    W.log(**rec)
                    return {'key': kappa, 'K': K, 'best_rel': cb['rel'], 'witness': True}
        rec['n_candidates'] = n_cand
        rec['note'] = f'{spec} d={d} key={kappa} K={k} (sector kernel {K}): {len(cands)} restarts, best rel {cb["rel"]:.2e}'
        W.log(**rec)
        results.append((spec, k, float(cb['rel'])))
        if best_overall is None or cb['rel'] < best_overall:
            best_overall = float(cb['rel'])
    return {'key': kappa, 'K': K, 'best_rel': best_overall, 'restarts': args.restarts, 'rows': kvp.R, 'cols': kvp.ncols,
            'seconds': round(time.time() - t0, 1), 'n_subfamilies': len(results)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--family', default='kernel')
    ap.add_argument('--keys', default='all')
    ap.add_argument('--exclude-keys', default='')
    ap.add_argument('--nontrivial', action='store_true', help='only sectors with dimK > dimJ')
    ap.add_argument('--restarts', type=int, default=20)
    ap.add_argument('--batch', type=int, default=10)
    ap.add_argument('--maxiter', type=int, default=150)
    ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--mode', default='complex')
    ap.add_argument('--tag', default='run')
    ap.add_argument('--key-budget', type=float, default=0)
    ap.add_argument('--min-K', type=int, default=1)
    ap.add_argument('--max-K', type=int, default=0)
    ap.add_argument('--rel-thr', type=float, default=1e-10)
    ap.add_argument('--method', default='lm', help='lm (Levenberg-Marquardt varpro, default) or lbfgs')
    ap.add_argument('--split', type=int, default=1, help='split the restarts of every key into this many parallel parts')
    ap.add_argument('--target', default=None, help='alternative target polynomial (e.g. the q=1 sibling B = X0*X1**2*Y0*Y1**2)')
    ap.add_argument('--subfamilies', default='cone4,mult4,cone3', help='for --family multi: comma list of structured kinds')
    ap.add_argument('--q', type=int, default=2, help='pole order for --target')
    args = ap.parse_args()
    F, A, q = load_target(args)
    S = W.Sectors(F, A)
    if args.keys == 'all':
        keys = S.all_keys(args.d)
    else:
        keys = [int(k) for k in args.keys.split(',')]
    if args.exclude_keys:
        ex = {int(k) for k in args.exclude_keys.split(',')}
        keys = [k for k in keys if k not in ex]
    if args.nontrivial:
        dims_file = os.path.join(LANE, f'w2_e1dims_d{args.d}.json')
        if os.path.exists(dims_file):
            dd = json.load(open(dims_file))['summary']['per_key']
            keys = [k for k in keys if str(k) in dd and dd[str(k)][2] > dd[str(k)][3]]
    t0 = time.time()
    W.log(step='search-start', family=args.family, d=args.d, tag=args.tag, n_keys=len(keys), restarts=args.restarts, maxiter=args.maxiter,
          mode=args.mode, nproc=args.nproc, seed=args.seed, q=q, target=core.to_str(A) if args.target else 'A',
          note=f'starting {args.family} d={args.d} over {len(keys)} sectors (q={q}, target {core.to_str(A) if args.target else "A"})')
    print(f'[w2_search] {args.family} d={args.d}: {len(keys)} sectors, {args.restarts} restarts each, {args.nproc} procs', flush=True)
    results = []
    with mp.get_context('fork').Pool(args.nproc) as pool:
        tasks = [(args, k, part) for k in keys for part in range(max(1, args.split))]
        for r in pool.imap_unordered(run_key, tasks):
            results.append(r)
            print(f'  key {r.get("key")}: K={r.get("K")} best_rel={r.get("best_rel")} restarts={r.get("restarts")} '
                  f'{r.get("rows")}x{r.get("cols")} {r.get("seconds")}s {"WITNESS" if r.get("witness") else ""} {r.get("error", "")}', flush=True)
            if r.get('witness'):
                print('WITNESS FOUND', flush=True)
    bests = [r['best_rel'] for r in results if r.get('best_rel') is not None]
    W.log(step='search-end', family=args.family, d=args.d, tag=args.tag, n_keys=len(keys), n_searched=len(bests),
          best_rel=min(bests) if bests else None, median_best_rel=float(np.median(bests)) if bests else None,
          seconds=round(time.time() - t0, 1), witness=any(r.get('witness') for r in results),
          note=f'finished {args.family} d={args.d}: best rel over sectors {min(bests) if bests else None}')
    with open(os.path.join(LANE, f'w2_summary_{args.tag}_d{args.d}.json'), 'w') as f:
        json.dump({'args': vars(args), 'results': results}, f, default=str)
    print(f'[w2_search] done in {time.time()-t0:.0f}s; best rel {min(bests) if bests else None}', flush=True)


if __name__ == '__main__':
    main()
