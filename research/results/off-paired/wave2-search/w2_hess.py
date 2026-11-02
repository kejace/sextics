#!/usr/bin/env python3
"""Hessian relaxation on the Jacobian-type family: D = sum_i w_i F_i with N_i := A w_i (no Koszul corrections), so
E1 holds exactly and the certificate reduces to E(w) := D * sum_i d_i(A w_i) - A * sum_i w_i d_i D  in  F * R_{2d-6}.
Quadratic in w; residual = distance of E(w) to span(F*R_{2d-6}) (sector-restricted), minimized by LM over the unit
sphere of w-space (w_i of degree d-5, character c + chi(z_i)).  Sufficient, not necessary.

  python3 w2_hess.py --d 10 --c 0,0 --restarts 2000 --workers 8 --job hess10c00
"""
import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); os.chdir(os.path.dirname(os.path.abspath(__file__)))
import argparse, json, time, math, glob
import multiprocessing as mp
import numpy as np, scipy.sparse as sp, scipy.linalg as sla
from fractions import Fraction
import w2lib
from hodge import core
from hodge.core import NV, mul, add, sub, scale, pderiv, shift

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True); ap.add_argument('--c', default='0,0'); ap.add_argument('--restarts', type=int, default=500)
ap.add_argument('--workers', type=int, default=8); ap.add_argument('--seed0', type=int, default=0); ap.add_argument('--maxiter', type=int, default=100)
ap.add_argument('--job', required=True); ap.add_argument('--time-limit', type=float, default=3000); ap.add_argument('--group', default='G')
args = ap.parse_args()
c = tuple(int(x) for x in args.c.split(','))
w2lib.set_group(args.group)
F, A, q = w2lib.load(); d = args.d
cA = w2lib.oplib.char_of(A)
T0 = time.time()
Fi = [pderiv(F, i) for i in range(NV)]
# w-space basis: (i, m) with deg m = d-5, char(m) = c + VARCHAR[i]
W = []
for i in range(NV):
    for m in w2lib.monos_char(d - 5, w2lib.cadd(c, w2lib.VARCHAR[i])):
        W.append((i, m))
K = len(W)
crow = w2lib.cadd(cA, w2lib.cadd(c, c))
rowmons = w2lib.monos_char(2 * d, crow)
rk = w2lib.keys_of(rowmons); order = np.argsort(rk); rks = rk[order]; R = len(rowmons)
def row_index(ks):
    pos = np.minimum(np.searchsorted(rks, ks), len(rks) - 1)
    assert np.all(rks[pos] == ks), 'row outside sector'
    return order[pos]
# quadratic form: E(w) = sum_{k,l} w_k w_l E_kl, E_kl = (m_k F_ik) * d_j(A m_l)[j = i_l]  - A m_l * d_{i_l}(m_k F_{i_k})
# (from D = sum_k w_k m_k F_{i_k}, N_j = A w_j = A sum_{l: i_l = j} w_l m_l)
Dk = [shift(Fi[i], m, 1) for (i, m) in W]
dAm = [pderiv(mul(A, {m: Fraction(1)}), i) for (i, m) in W]           # d_{i_l}(A m_l)
Am = [mul(A, {m: Fraction(1)}) for (i, m) in W]
ri_l, ci_l, va_l = [], [], []
for k in range(K):
    Dk_di = [pderiv(Dk[k], j) for j in range(NV)]
    for l in range(K):
        il = W[l][0]
        E = sub(mul(Dk[k], dAm[l]), mul(Am[l], Dk_di[il]))
        if not E: continue
        ks, vs = w2lib.poly_arrays(E)
        rows = row_index(ks)
        ri_l.append(rows); ci_l.append(np.full(len(ks), k * K + l, dtype=np.int64)); va_l.append(vs)
TE = sp.csr_matrix((np.concatenate(va_l), (np.concatenate(ri_l), np.concatenate(ci_l))), shape=(R, K * K))
# F * R_{2d-6} block and its orthonormal basis
Qm = w2lib.monos_char(2 * d - 6, crow)
fk, fv = w2lib.poly_arrays(F)
MQ = np.zeros((R, len(Qm)))
if Qm:
    qk = w2lib.keys_of(Qm)
    kk = (fk[:, None] + qk[None, :]).ravel(); rows = row_index(kk)
    cols = np.broadcast_to(np.arange(len(Qm))[None, :], (len(fk), len(Qm))).ravel()
    np.add.at(MQ, (rows, cols), np.broadcast_to(fv[:, None], (len(fk), len(Qm))).ravel())
    QQ, _ = np.linalg.qr(MQ, mode='reduced')
else:
    QQ = np.zeros((R, 0))
# remove the trivial directions: w with D(w) = 0 (Koszul: w = (F_j e_i - F_i e_j) m) and gauge D in F*R: handled by
# restricting to the complement of {w : D(w) in F*R_{d-6}} -- compute the linear map w -> D coefficients
Dm = w2lib.monos_char(d, c); Didx = {m: i for i, m in enumerate(Dm)}
Bw = np.zeros((K, len(Dm)))
for k in range(K):
    for e, v in Dk[k].items(): Bw[k, Didx[e]] = float(v)
G = np.zeros((max(1, len(w2lib.monos_char(d - 6, c))), len(Dm)))
gp = w2lib.gauge_polys(F, d, c)
for k, P in enumerate(gp):
    for e, v in P.items(): G[k, Didx[e]] = float(v)
# w-directions whose D lies in span(gauge): solve Bw^T w = G^T g  -> nullspace of [Bw^T | -G^T]
if gp:
    Mbig = np.hstack([Bw.T, -G.T]); ns = sla.null_space(Mbig); S = ns[:K, :]
else:
    S = sla.null_space(Bw.T)   # D(w) = 0 directions only
if S.shape[1]:
    U, s, Vt = np.linalg.svd(S, full_matrices=True); r = int(np.sum(s > 1e-10 * s[0])); P = U[:, r:]
else:
    P = np.eye(K)
dimP = P.shape[1]
print(f'[hess] d={d} c={c} K={K} dimP={dimP} rows={R} nQ={len(Qm)} nnz(TE)={TE.nnz} build {time.time()-T0:.1f}s', flush=True)
w2lib.append_json(w2lib.REMOTE_LOG, {'ts': w2lib.now(), 'step': 'hess-build', 'job': args.job, 'd': d, 'sector': list(c),
                  'K_w': K, 'dimP': dimP, 'rows': R, 'nQ': len(Qm)})

def resid(u):
    a = P @ u; a = a / np.linalg.norm(a)
    E = TE @ np.outer(a, a).ravel()
    r = E - QQ @ (QQ.T @ E) if QQ.shape[1] else E
    return r, a, float(E @ E)

def resid_jac(u):
    a_full = P @ u; na = np.linalg.norm(a_full); a = a_full / na
    E = TE @ np.outer(a, a).ravel()
    r = E - QQ @ (QQ.T @ E) if QQ.shape[1] else E
    # dE/da_k = sum_l a_l (E_kl + E_lk): columns of TE indexed (k*K + l)
    Sw = sp.csr_matrix((np.concatenate([np.tile(a, K), np.repeat(a, K)]),
                        (np.concatenate([np.arange(K * K), np.arange(K * K)]),
                         np.concatenate([np.repeat(np.arange(K), K), np.tile(np.arange(K), K)]))), shape=(K * K, K))
    J = np.asarray((TE @ Sw).todense())
    if QQ.shape[1]: J = J - QQ @ (QQ.T @ J)
    J = (J - np.outer(J @ a, a)) / na
    return r, J @ P, float(E @ E)

def lm(seed, maxiter):
    rng = np.random.default_rng(seed); u = rng.standard_normal(dimP); u /= np.linalg.norm(u)
    r, J, en = resid_jac(u); f = float(r @ r); lam = 1e-3; st = 'maxiter'; it = 0
    while it < maxiter:
        it += 1; g = J.T @ r; H = J.T @ J; dH = np.diag(H).copy(); dH[dH < 1e-300] = 1e-300; acc = False
        for _ in range(30):
            try: delta = np.linalg.solve(H + lam * np.diag(dH), -g)
            except np.linalg.LinAlgError: lam *= 10; continue
            un = u + delta; un /= np.linalg.norm(un); rn, an, enn = resid(un); fn = float(rn @ rn)
            if fn < f: acc = True; break
            lam *= 4
            if lam > 1e14: break
        if not acc: st = 'stalled'; break
        impr = (f - fn) / max(f, 1e-300); u, f = un, fn; lam = max(lam / 3, 1e-15)
        r, J, en = resid_jac(u)
        if math.sqrt(f / max(en, 1e-300)) < 1e-13: st = 'converged'; break
        if impr < 1e-12 and lam > 1e6: st = 'flat'; break
    a = P @ u; a /= np.linalg.norm(a)
    return {'seed': int(seed), 'relres': math.sqrt(f / max(en, 1e-300)), 'f': f, 'Enorm2': en, 'iters': it, 'status': st, 'a': a}

def worker(k, seeds):
    out = f'wave2-{args.job}.w{k}.jsonl'
    for s in seeds:
        if time.time() - T0 > args.time_limit: break
        t = time.time(); r = lm(s, args.maxiter)
        rec = {kk: v for kk, v in r.items() if kk != 'a'}; rec['time'] = round(time.time() - t, 2)
        if r['relres'] < 1e-4:
            Dc = Bw.T @ r['a']; rec['Dc'] = [float(x) for x in Dc]
        w2lib.append_json(out, rec)

seeds = list(range(args.seed0, args.seed0 + args.restarts)); chunks = [seeds[k::args.workers] for k in range(args.workers)]
procs = [mp.get_context('fork').Process(target=worker, args=(k, chunks[k])) for k in range(args.workers)]
for p_ in procs: p_.start()
for p_ in procs: p_.join()
recs = []
for fn in sorted(glob.glob(f'wave2-{args.job}.w*.jsonl')):
    for line in open(fn):
        try: recs.append(json.loads(line))
        except Exception: pass
rels = sorted(r['relres'] for r in recs)
summary = {'ts': w2lib.now(), 'step': 'hess', 'job': args.job, 'd': d, 'sector': list(c), 'K_w': K, 'dimP': dimP, 'rows': R, 'nQ': len(Qm),
           'restarts': len(rels), 'best_relres': rels[0] if rels else None, 'median_relres': rels[len(rels)//2] if rels else None,
           'minima': sorted({f'{x:.3e}' for x in rels})[:10], 'wall': round(time.time() - T0, 1),
           'mean_time': round(float(np.mean([r['time'] for r in recs])), 2) if recs else None}
# exact confirmation of small residuals via the general machinery (D = Bw^T a; N free in the exact solve)
cands = sorted([r for r in recs if 'Dc' in r and r['relres'] < 1e-9], key=lambda r: r['relres'])
if cands:
    basis = [{m: Fraction(1)} for m in Dm]
    fam = w2lib.SectorVarPro(F, A, q, basis, c, gauge=w2lib.gauge_polys(F, d, c))
    atts = []
    for r in cands[:10]:
        att = w2lib.try_reconstruct_and_confirm(fam, np.array(r['Dc']))
        atts.append({'seed': r['seed'], 'relres': r['relres'], 'attempts': att})
        if any(a.get('status') == 'solved' for a in att):
            sol = [a for a in att if a.get('status') == 'solved'][0]
            json.dump({'d': d, 'D': sol['D'], 'N': sol['N'], 'Q': sol['Q'], 'job': args.job}, open(f'witness-candidate-{args.job}.json', 'w'))
            json.dump({'d': d, 'D': sol['D'], 'N': sol['N'], 'Q': sol['Q'], 'job': args.job}, open('witness-candidate.json', 'w'))
            print('!!!!!!!! WITNESS FOUND (exact, verified) !!!!!!!!', flush=True); summary['witness'] = True; break
    summary['confirm_attempts'] = atts
json.dump(summary, open(f'wave2-{args.job}.summary.json', 'w'), indent=1)
w2lib.append_json(w2lib.REMOTE_LOG, {k: v for k, v in summary.items() if k != 'confirm_attempts'})
print(json.dumps({k: v for k, v in summary.items() if k != 'confirm_attempts'}), flush=True)
