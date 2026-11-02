"""E1 kernels K_d^{par} = {D in S_d^par : A D in J(F)} at new primes, both parities, with Hilbert certification.

usage: python3 wave2_excl_kernel.py <jobname> <workers> <spec>
  spec = comma-separated items  d:par:primeindex[:basis]   e.g.  10:0:1:b,10:1:1:b,9:0:1,...
Output: wave2-excl.<jobname>.jsonl (one line per task) and kernel files wave2-excl-K_d<d>_par<par>_p<p>.json
"""
import sys, os, json, time, traceback
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import w2x
from w2x import core, monomials, shift, pderiv, nmod_mat, F, A, hilb_R, hilb_S, dim_J, logline
from hodge import jacobian
import multiprocessing as mp

JOB = sys.argv[1]
WORKERS = int(sys.argv[2])
SPEC = sys.argv[3]
LOG = os.path.join(w2x.LANE, f'wave2-excl.{JOB}.jsonl')


def kernel_task(d, parity, p, want_basis):
    t0 = time.time()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(6)]
    tpar = (parity + 1) % 2                      # A is Y-odd
    cols_list = monomials(d + 6, tpar)
    col = {e: j for j, e in enumerate(cols_list)}
    C = len(cols_list)
    Jrows = []
    for i in range(6):
        par_i = 0 if i < 3 else 1
        for m in monomials(d + 1, (tpar - par_i) % 2):
            Jrows.append(shift(Fi[i], m, 1, p))
    Amons = monomials(d, parity)
    Arows = [shift(Ap, m, 1, p) for m in Amons]
    RJ, RA = len(Jrows), len(Arows)
    flat = [0] * ((RJ + RA) * C)
    for r, poly in enumerate(Jrows + Arows):
        base = r * C
        for e, v in poly.items():
            flat[base + col[e]] = v
    MJ = nmod_mat(RJ, C, flat[:RJ * C], p)
    rJ = MJ.rank()
    del MJ
    M = nmod_mat(RJ + RA, C, flat, p)
    del flat
    rM = M.rank()
    t1 = time.time()
    dimK = RA - (rM - rJ)
    dimJd = dim_J(d, parity) if d >= 5 else 0
    hilbJ = dim_J(d + 6, tpar)
    rec = dict(exp='e1kernel', d=d, parity=parity, prime=p, cols=C, rowsJ=RJ, rowsA=RA, rankJ=rJ, rank_stack=rM,
               dimK=dimK, dimJd=dimJd, dimK_mod_J=dimK - dimJd, dimR_target_modp=C - rJ, hilbR_target=hilb_R(d + 6, tpar),
               hilbJ_target=hilbJ, J_rank_certified=(rJ == hilbJ), hilbR_d=hilb_R(d, parity),
               gorenstein_lower=max(0, hilb_R(d, parity) - hilb_R(d + 6, tpar)), seconds_rank=round(t1 - t0, 1), job=JOB)
    rec['rigorous'] = rec['J_rank_certified']
    rec['note'] = ('K_d = J_d over Q (rigorous: modular kernel dim is an upper bound and equals dim J_d)' if (dimK == dimJd and rec['rigorous'])
                   else ('kernel exceeds J_d by %d (upper bound over Q; equals Gorenstein lower bound %d? %s)' % (dimK - dimJd, rec['gorenstein_lower'], dimK - dimJd == rec['gorenstein_lower'])))
    if want_basis and dimK > dimJd:
        Mt = M.transpose(); del M
        X, nul = Mt.nullspace(); del Mt
        cvecs = []
        for kk in range(nul):
            v = [int(X[RJ + r, kk]) for r in range(RA)]
            if any(v):
                cvecs.append(v)
        del X
        Cm = nmod_mat(len(cvecs), RA, [x for v in cvecs for x in v], p)
        Cr, rk = Cm.rref(); del Cm
        Kb = [[int(Cr[i, r]) for r in range(RA)] for i in range(rk)]
        del Cr
        # reduce modulo J_d^{parity}: rref of J_d block over the same monomial basis Amons
        Jd, _ = jacobian.jacobian_matrix(F, d, p, parity, Amons)
        Jr, rJd = Jd.rref(); del Jd
        pivots = []
        for i in range(rJd):
            j = 0
            while j < RA and int(Jr[i, j]) == 0:
                j += 1
            pivots.append(j)
        Kp = nmod_mat(rk, rJd, [Kb[i][j] for i in range(rk) for j in pivots], p)
        Kmat = nmod_mat(rk, RA, [x for row in Kb for x in row], p)
        Kred = Kmat - Kp * Jr
        Kr, rk2 = Kred.rref()
        basis_modJ = [{' '.join(map(str, Amons[j])): int(Kr[i, j]) for j in range(RA) if int(Kr[i, j])} for i in range(rk2)]
        Kfull = [{' '.join(map(str, Amons[j])): Kb[i][j] for j in range(RA) if Kb[i][j]} for i in range(rk)]
        out = {'d': d, 'parity': parity, 'p': p, 'dimK': rk, 'dimJd_modp': rJd, 'dimK_mod_J': rk2,
               'monomial_key': 'exponents e0..e5 of X0 X1 X2 Y0 Y1 Y2 separated by spaces; values are residues mod p',
               'standard_monomials': [' '.join(map(str, Amons[j])) for j in range(RA) if j not in set(pivots)],
               'K_mod_J_basis': basis_modJ, 'K_basis': Kfull}
        path = os.path.join(w2x.LANE, f'wave2-excl-K_d{d}_par{parity}_p{p}.json')
        with open(path, 'w') as f:
            json.dump(out, f)
        rec.update(basis_file=path, dimK_basis=rk, dimK_mod_J_basis=rk2, dimJd_modp=rJd)
    rec['seconds'] = round(time.time() - t0, 1)
    return rec


def run(item):
    parts = item.split(':')
    d, par, pi = int(parts[0]), int(parts[1]), int(parts[2])
    basis = len(parts) > 3 and parts[3].startswith('b')
    p = core.PRIMES[pi]
    try:
        rec = kernel_task(d, par, p, basis)
    except Exception as ex:
        rec = dict(exp='e1kernel', d=d, parity=par, prime=p, result='ERROR', error=traceback.format_exc()[-800:], job=JOB)
    logline(LOG, rec)
    print(json.dumps({k: v for k, v in rec.items() if k not in ('note',)}), flush=True)
    return rec


if __name__ == '__main__':
    items = [s for s in SPEC.split(',') if s]
    # heavy first
    items.sort(key=lambda s: -int(s.split(':')[0]))
    with mp.get_context('fork').Pool(WORKERS) as pool:
        for _ in pool.imap_unordered(run, items):
            pass
    print('ALL DONE', flush=True)
