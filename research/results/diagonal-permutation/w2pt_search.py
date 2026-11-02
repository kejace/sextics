"""Driver: multi-process pointwise variable-projection search (see w2pt.py).

usage: python3 w2pt_search.py --d 8 --tag j1 --workers 8 --restarts 100 --seed 1000 [--extra D8] [--sub families.json:key]
                              [--maxiter 400] [--margin 60] [--hit 1e-9]
Each worker uses its own random point set and RNG; one JSON line per restart is appended to w2-pt-d<d>-<tag>.jsonl;
hits (rel residual < --hit) are refined, re-tested on an independent point set, rationally reconstructed and
confirmed exactly with the engine (fixed_d).  A verified witness is written to witness-candidate.json.
"""
import os, sys, time, json, argparse
os.environ['OPENBLAS_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'; os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np
import multiprocessing as mp
from scipy.linalg import null_space as sla_null_space
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2pt import *
from w2sym import w_sector_matrix, c_sector_matrix, symmetrize_pair, sig_poly
from hodge import core, fixed_d
from fractions import Fraction
LANE = os.path.dirname(os.path.abspath(__file__))
IID = 'diagonal-permutation'


def e1_lift_float(F, A, D):
    """Least-squares N (6 polys, float coefficients) with N(F) = A*D (D must satisfy E1).  Returns (N, relres)."""
    import scipy.linalg as sla
    d = core.degree(D)
    pars = fixed_d.numerator_parities(D, A)
    degN = d + 1
    Fi = [pderiv(F, i) for i in range(NV)]
    cols = []
    for i in range(NV):
        for m in monomials(degN, pars[i]):
            cols.append((i, m))
    rhs = mul(A, D)
    rows = {}
    for (i, m) in cols:
        for e in shift(Fi[i], m):
            rows.setdefault(e, len(rows))
    for e in rhs:
        rows.setdefault(e, len(rows))
    M = np.zeros((len(rows), len(cols)))
    for j, (i, m) in enumerate(cols):
        for e, v in shift(Fi[i], m).items():
            M[rows[e], j] = float(v)
    b = np.zeros(len(rows))
    for e, v in rhs.items():
        b[rows[e]] = float(v)
    x, *_ = sla.lstsq(M, b, cond=1e-12, lapack_driver='gelsd', check_finite=False)
    rel = float(np.linalg.norm(M @ x - b) / np.linalg.norm(b))
    N = [{} for _ in range(NV)]
    for j, (i, m) in enumerate(cols):
        if abs(x[j]) > 1e-14:
            N[i][m] = float(x[j])
    return N, rel


args_d = [None]
args_pD = [0]

def load_extra(F, A, names):
    extra = []
    for nm in names:
        if nm == 'D8':
            D = core.parse(json.load(open(os.path.join(LANE, 'D8_e1kernel.json')))['D'])
            N, rel = e1_lift_float(F, A, D)
            print(f'[extra] D8 float E1 lift residual {rel:.2e}', flush=True)
            extra.append((D, N))
        elif nm.startswith('kerfam'):
            fn = nm.split(':')[1] if ':' in nm else f'w2-kerfam-d{args_d[0]}-p{args_pD[0]}.npz'
            z = np.load(os.path.join(LANE, fn))
            Dmons = [tuple(int(x) for x in e) for e in z['Dmons']]
            Ni = z['Ncols_i']; Nm = [tuple(int(x) for x in e) for e in z['Ncols_m']]
            for k in range(z['Dvecs'].shape[0]):
                D = {Dmons[j]: float(v) for j, v in enumerate(z['Dvecs'][k]) if abs(v) > 1e-15}
                N = [{} for _ in range(NV)]
                for t, v in enumerate(z['Nvecs'][k]):
                    if abs(v) > 1e-15: N[int(Ni[t])][Nm[t]] = float(v)
                extra.append((D, N))
            print(f'[extra] loaded {len(extra)} kernel members from {fn}', flush=True)
        else:
            raise ValueError(nm)
    return extra


def sub_w_matrix(F, d, wbasis_index, Dstrs):
    """Express exact Jacobian-type polys D (strings) as w-vectors (exact fmpq solve), return K0 x s float matrix."""
    from flint import fmpq_mat, fmpq
    Fi = [pderiv(F, i) for i in range(NV)]
    cols = []
    for (i, m) in wbasis_index:
        cols.append(shift(Fi[i], m))
    rows = {}
    Ds = [core.parse(s) for s in Dstrs]
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for D in Ds:
        for e in D: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    M = fmpq_mat(R, C + len(Ds))
    for j, c in enumerate(cols):
        for e, v in c.items(): M[rows[e], j] = fmpq(v.numerator, v.denominator)
    for t, D in enumerate(Ds):
        for e, v in D.items(): M[rows[e], C + t] = fmpq(v.numerator, v.denominator)
    Mr, rk = M.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < C + len(Ds) and Mr[i, j] == 0: j += 1
        piv.append(j)
    W = np.zeros((C, len(Ds)))
    for t in range(len(Ds)):
        assert (C + t) not in piv, 'D not in the Jacobian family'
        for i, j in enumerate(piv):
            if j < C:
                v = Mr[i, C + t]
                W[j, t] = -float(Fraction(int(v.p), int(v.q)))
    return W


def worker(args, wid, conn):
    try:
        F, A, q, rec = core.load_instance(IID)
        d = args.d
        args_d[0] = d; args_pD[0] = args.pD
        extra = load_extra(F, A, args.extra.split(',')) if args.extra else None
        if extra and args.sector:
            extra = [symmetrize_pair(D, N, args.sector) for (D, N) in extra]
            extra = [(D, N) for (D, N) in extra if D]
        t0 = time.time()
        # dimensions of the search variety (sphere in a) x (Koszul unknowns), sector-aware, to choose the point count:
        # with complex points each point gives 2 real equations; need 2r > dim(a) - 1 + dim(c) (no false zeros), plus margin
        wb0 = [(i, m) for i in range(NV) for m in monomials(d - 5, (args.pD + PAR[i]) % 2)]
        K0 = len(wb0)
        Kx = len(extra) if extra else 0
        pA = core.parity_of(A)
        cols0 = [(i, j, m) for i in range(NV) for j in range(i + 1, NV) for m in monomials(d - 4, ((pA - args.pD - PAR[i]) % 2 - PAR[j]) % 2)]
        nK = len(cols0)
        if args.sector:
            dim_a = w_sector_matrix(wb0, args.sector).shape[1] + Kx
            dim_c = c_sector_matrix(cols0, -args.sector).shape[1]
        else:
            dim_a, dim_c = K0 + Kx, nK
        r = (dim_a + dim_c + 1) // 2 + args.margin
        pts = Points(F, r, seed=args.seed * 7919 + wid)
        fam = jacobian_family(F, A, d, pts, extra=extra, pD=args.pD)
        kb = KoszulBlock(F, A, d, pts, pD=args.pD)
        subspace = None
        if args.sub:
            fn, key = args.sub.split(':')
            Dstrs = json.load(open(os.path.join(LANE, fn)))[key]
            if isinstance(Dstrs, dict): Dstrs = Dstrs.get('all') or sum(Dstrs.values(), [])
            wb = [m['w'] for m in fam.members if m['w'] is not None]
            W = sub_w_matrix(F, d, wb, Dstrs)
            subspace = np.zeros((fam.K, W.shape[1])); subspace[:W.shape[0]] = W
        # directions with D(a) = 0 (linear dependence among members) are spurious: add them
        allm = {}
        for m in fam.members:
            for e in m['Dpoly']: allm.setdefault(e, len(allm))
        Dmat = np.zeros((fam.K, len(allm)))
        for k, m in enumerate(fam.members):
            for e, v in m['Dpoly'].items(): Dmat[k, allm[e]] = float(v)
        ns = sla_null_space(Dmat.T)
        if ns.shape[1]:
            fam.spurious = np.hstack([fam.spurious, ns]) if fam.spurious.size else ns
        if args.vanish:
            # linear conditions on a: coefficients of D(a) on monomials outside the ideal I^k vanish
            spec = args.vanish
            if spec.startswith('X'): idxs, k = [0, 1, 2], int(spec[1:])
            elif spec.startswith('Y'): idxs, k = [3, 4, 5], int(spec[1:])
            elif spec.startswith('L01x'): idxs, k = [0, 1, 3, 4], int(spec[4:])
            elif spec.startswith('L02x'): idxs, k = [0, 2, 3, 5], int(spec[4:])
            elif spec.startswith('P0x'): idxs, k = [0, 3], int(spec[3:])
            else: raise ValueError(spec)
            bad = [j for e, j in allm.items() if sum(e[i] for i in idxs) < k]
            Cmat = Dmat[:, bad]                       # K x nbad ; need a with a @ Cmat = 0
            Wv = sla_null_space(Cmat.T)
            subspace = Wv if subspace is None else Wv @ np.linalg.lstsq(Wv, subspace, rcond=None)[0]
            conn.send({'type': 'info', 'wid': wid, 'vanish': spec, 'subspace_dim': int(Wv.shape[1])})
        csector = None
        if args.sector:
            wb = [m['w'] for m in fam.members if m['w'] is not None]
            Ws = w_sector_matrix(wb, args.sector)
            Wfull = np.zeros((fam.K, Ws.shape[1] + (fam.K - len(wb))))
            Wfull[:len(wb), :Ws.shape[1]] = Ws
            for t in range(fam.K - len(wb)):
                Wfull[len(wb) + t, Ws.shape[1] + t] = 1.0
            if subspace is None:
                subspace = Wfull
            else:
                # intersection of span(Wfull) and span(subspace): null space of [Wfull, -subspace]
                Nn = sla_null_space(np.hstack([Wfull, -subspace]))
                subspace = Wfull @ Nn[:Wfull.shape[1], :]
            csector = c_sector_matrix(kb.cols, -args.sector)
        vp = VarPro(fam, kb, pts, csector=csector)
        if subspace is not None:
            # restrict to span(subspace) orthogonal to spurious directions
            S = fam.spurious
            Wn = subspace
            if S.size and S.shape[1]:
                Qs, _ = np.linalg.qr(S)
                Wn = Wn - Qs @ (Qs.T @ Wn)
            U, s, _ = np.linalg.svd(Wn, full_matrices=False)
            rk = int(np.sum(s > 1e-10 * max(1.0, s[0])))
            vp.P = U[:, :rk]; vp.dimP = rk
        # reference scale of |B_D(N^0)|^2 on the sphere (median over random directions), for an absolute hit criterion
        rng0 = np.random.default_rng(12345 + wid)
        refs = []
        for _ in range(12):
            a0 = vp.P @ rng0.standard_normal(vp.dimP); a0 /= np.linalg.norm(a0)
            refs.append(vp.bnorm2(a0))
        ref = float(np.median(refs))
        setup = time.time() - t0
        out_path = os.path.join(LANE, f'w2-pt-d{d}-{args.tag}.jsonl')
        label = fam.label + (f'|pD{args.pD}' if args.pD else '') + (f'|vanish:{args.vanish}' if args.vanish else '') + (f'|sub:{args.sub}' if args.sub else '') + (f'|sector{args.sector:+d}' if args.sector else '')
        conn.send({'type': 'setup', 'wid': wid, 'r': r, 'K': fam.K, 'nK': kb.nK, 'dimP': vp.dimP, 'nspur': vp.nspur,
                   'setup_s': round(setup, 1), 'Fres': pts.Fres, 'label': label, 'nKsector': int(vp.nKc), 'ref_B0norm2': ref})
        rels = []
        rng_seed = args.seed * 1000003 + wid * 101
        def cb(c):
            rels.append(c['rel'])
            c['rel_abs'] = c['f'] / ref
            line = {'d': d, 'family': label, 'tag': args.tag, 'wid': wid, 'restart': c['restart'], 'rel': c['rel'], 'f': c['f'], 'rel_abs': c['rel_abs'],
                    'iters': c['iters'], 'seconds': c['seconds'], 'amax': float(np.max(np.abs(c['a']))),
                    'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
            with open(out_path, 'a') as fh: fh.write(json.dumps(line) + '\n')
            if min(c['rel'], c['rel_abs']) < args.hit:
                try:
                    handle_hit(args, wid, F, A, q, fam, kb, pts, vp, c, label)
                except Exception as ex:
                    import traceback
                    with open(os.path.join(LANE, f'w2-hits-d{d}.jsonl'), 'a') as fh:
                        fh.write(json.dumps({'d': d, 'tag': args.tag, 'wid': wid, 'restart': c['restart'], 'rel': c['rel'],
                                             'a': [float(x) for x in c['a']], 'error': repr(ex), 'tb': traceback.format_exc()[-2000:]}) + '\n')
            if len(rels) % 25 == 0:
                conn.send({'type': 'progress', 'wid': wid, 'n': len(rels), 'best': min(rels), 'median': float(np.median(rels)),
                           'elapsed': round(time.time() - t0, 1)})
        vp.search(args.restarts, seed=rng_seed, maxiter=args.maxiter, callback=cb)
        hist = np.histogram(np.log10(np.maximum(rels, 1e-300)), bins=[-300, -14, -12, -10, -8, -6, -4, -3, -2, -1, 0, 10])[0].tolist()
        conn.send({'type': 'done', 'wid': wid, 'n': len(rels), 'best': min(rels) if rels else None,
                   'median': float(np.median(rels)) if rels else None, 'hist_log10': hist, 'elapsed': round(time.time() - t0, 1)})
    except Exception as ex:
        import traceback
        conn.send({'type': 'error', 'wid': wid, 'error': repr(ex), 'tb': traceback.format_exc()[-3000:]})
    conn.close()


def handle_hit(args, wid, F, A, q, fam, kb, pts, vp, c, label):
    d = args.d
    tag = f'{args.tag}-w{wid}-r{c["restart"]}'
    rep = {'d': d, 'family': label, 'tag': tag, 'rel_initial': c['rel'], 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
    a = vp.refine(c['a'])
    f, rel, cc, rho = vp.value(a)
    rep['rel_refined'] = rel
    # independent point set
    pts2 = Points(F, pts.r, seed=args.seed * 31 + 17 * wid + c['restart'] + 999)
    extra = [(m['D'], None) for m in fam.members if m['w'] is None]
    fam2 = jacobian_family(F, A, d, pts2, extra=None, pD=args.pD) if not extra else None
    if fam2 is not None:
        kb2 = KoszulBlock(F, A, d, pts2, pD=args.pD)
        vp2 = VarPro(fam2, kb2, pts2)
        f2, rel2, _, _ = vp2.value(a)
        rep['rel_independent_points'] = rel2
    rep['a'] = [float(x) for x in a]
    json.dump(rep, open(os.path.join(LANE, f'w2-cand-d{d}-{tag}.json'), 'w'))
    # rational reconstruction + exact confirmation (only for exact families)
    status = None
    if fam.exact:
        for fr in rational_candidates(a):
            D = fam.D_exact(fr)
            try:
                status, data = confirm_exact(F, A, q, D)
            except Exception as ex:
                status, data = 'error', {'error': repr(ex)}
            rep.setdefault('exact_attempts', []).append({'den_max': max(x.denominator for x in fr), 'status': status,
                                                           'info': {k: v for k, v in data.items() if k in ('rank', 'nullity', 'note', 'e1', 'error')}})
            if status == 'solved':
                N, Q = data['N'], data['Q']
                ok, info = core.verify_witness(F, A, q, D, N, Q)
                w = {'d': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q), 'verified': ok}
                json.dump(w, open(os.path.join(LANE, 'witness-candidate.json'), 'w'), indent=1)
                json.dump(w, open(os.path.join(LANE, f'witness-candidate-{tag}.json'), 'w'), indent=1)
                rep['WITNESS'] = True
                break
    if not fam.exact:
        # reconstruct the polynomial D(a) directly (normalize largest coefficient to 1)
        Df = {}
        for k, m in enumerate(fam.members):
            for e, v in m['Dpoly'].items(): Df[e] = Df.get(e, 0.0) + a[k] * float(v)
        emax = max(Df, key=lambda e: abs(Df[e])); scale_ = Df[emax]
        for md in (100, 10 ** 4, 10 ** 6):
            Dr = {e: Fraction(v / scale_).limit_denominator(md) for e, v in Df.items()}
            Dr = {e: v for e, v in Dr.items() if v}
            err = max(abs(float(Dr.get(e, 0)) - v / scale_) for e, v in Df.items())
            if err > 1e-6: continue
            try:
                status, data = confirm_exact(F, A, q, Dr)
            except Exception as ex:
                status, data = 'error', {'error': repr(ex)}
            rep.setdefault('exact_attempts', []).append({'den_max': md, 'status': status, 'recon_err': err,
                                                           'info': {k: v for k, v in data.items() if k in ('rank', 'nullity', 'note', 'e1', 'error')}})
            if status == 'solved':
                N, Q = data['N'], data['Q']
                ok, info = core.verify_witness(F, A, q, Dr, N, Q)
                w = {'d': d, 'D': core.to_str(Dr), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q), 'verified': ok}
                json.dump(w, open(os.path.join(LANE, 'witness-candidate.json'), 'w'), indent=1)
                json.dump(w, open(os.path.join(LANE, f'witness-candidate-{tag}.json'), 'w'), indent=1)
                rep['WITNESS'] = True
                break
    rep['exact_status'] = status
    with open(os.path.join(LANE, f'w2-hits-d{d}.jsonl'), 'a') as fh: fh.write(json.dumps(rep) + '\n')
    print('[HIT]', json.dumps({k: v for k, v in rep.items() if k != 'a'}), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--tag', required=True)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--restarts', type=int, default=100, help='restarts per worker')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--extra', default='', help='comma list of extra members, e.g. D8')
    ap.add_argument('--sub', default='', help='file.json:key with exact Jacobian-type D strings spanning a sub-family')
    ap.add_argument('--maxiter', type=int, default=400)
    ap.add_argument('--margin', type=int, default=100, help='extra complex points beyond (dim a + dim c)/2')
    ap.add_argument('--hit', type=float, default=1e-9)
    ap.add_argument('--vanish', default='', help='restrict D(a) to an ideal power: X3, X4, Y3, Y4, L01x4, P0x2 ...')
    ap.add_argument('--pD', type=int, default=0, help='Y-parity of D (0 even, 1 odd)')
    ap.add_argument('--sector', type=int, default=0, help='(01)-sector: +1 D invariant (N anti-equivariant), -1 D anti-invariant, 0 none')
    args = ap.parse_args()
    ctx = mp.get_context('fork')
    procs = []
    for wid in range(args.workers):
        parent, child = ctx.Pipe(duplex=False)
        pr = ctx.Process(target=worker, args=(args, wid, child)); pr.start(); child.close()
        procs.append((pr, parent))
    done = 0; summary = []
    t0 = time.time()
    while done < len(procs):
        for pr, conn in procs:
            if conn.poll(0.5):
                try:
                    msg = conn.recv()
                except EOFError:
                    continue
                print(f'[{time.strftime("%H:%M:%S")}] ' + json.dumps(msg)[:700], flush=True)
                if msg['type'] in ('done', 'error'):
                    done += 1; summary.append(msg)
        procs = [(pr, conn) for pr, conn in procs if pr.is_alive() or conn.poll(0)]
        if not procs: break
    tot = sum(m.get('n', 0) for m in summary)
    best = min([m['best'] for m in summary if m.get('best') is not None], default=None)
    hist = [sum(m.get('hist_log10', [0]*11)[i] for m in summary) for i in range(11)]
    out = {'exp': 'w2pt-search', 'd': args.d, 'tag': args.tag, 'extra': args.extra, 'sub': args.sub, 'sector': args.sector, 'pD': args.pD, 'vanish': args.vanish, 'workers': args.workers,
           'restarts_total': tot, 'best_rel': best, 'hist_log10_bins': '[<-14,-12,-10,-8,-6,-4,-3,-2,-1,0,>0]', 'hist': hist,
           'maxiter': args.maxiter, 'seconds': round(time.time() - t0, 1), 'seed': args.seed, 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
    with open(os.path.join(LANE, 'w2-pt-summary.jsonl'), 'a') as fh: fh.write(json.dumps(out) + '\n')
    print('[SUMMARY]', json.dumps(out), flush=True)


if __name__ == '__main__':
    main()
