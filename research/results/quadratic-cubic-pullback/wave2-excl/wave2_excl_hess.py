"""Hessian relaxation of the polar family (rigorous, Hilbert-certified).

For D = w.grad F the E2 condition Phi(w) + L(c,w) = F Q with L(c,w), F Q in J(F) forces Phi(w) in J(F)_{2d}, and
Phi(w) = -A w.grad D (mod J) = -A sum_{k,l} w_k w_l F_kl (mod J).  Hence every admissible w satisfies
      A . (w^T Hess(F) w)  in  J(F)_{2d}.
Relaxing W_kl = w_k w_l to independent unknowns gives a LINEAR condition on W whose kernel K_H is computed mod p;
the modular dimension is an upper bound for the rational one as soon as the Jacobian block of degree 2d has its
Hilbert rank (same certification as the E1 sweeps).  Admissible w satisfy w w^T in K_H.

usage: python3 wave2_excl_hess.py <jobname> <d> <parity 0|1|all> <primeindex or raw prime> [msolve=1]
"""
import sys, os, json, time, subprocess, random
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import w2x
from w2x import core, F, A, monomials, mul, add, sub, scale, shift, pderiv, logline, hilb_R, hilb_S, dim_J
from flint import nmod_mat
JOB, D_DEG, PAR = sys.argv[1], int(sys.argv[2]), sys.argv[3]
pp = int(sys.argv[4]); P = core.PRIMES[pp] if pp < 100 else pp
MSOLVE = any(a == 'msolve=1' for a in sys.argv[5:])
LOG = os.path.join(w2x.LANE, f'wave2-excl.{JOB}.jsonl')
parz = [0, 0, 0, 1, 1, 1]
pA = core.parity_of(A)


def main():
    t0 = time.time()
    d = D_DEG
    pD = None if PAR in ('all', 'cross') else int(PAR)
    Fp = core.reduce_mod(F, P); Ap = core.reduce_mod(A, P)
    Fi = [pderiv(Fp, i, P) for i in range(6)]
    Fij = [[pderiv(Fi[i], j, P) for j in range(6)] for i in range(6)]
    wb = []
    if PAR == 'cross':
        for i in range(6):
            for m in monomials(d - 5, (0 - parz[i]) % 2):
                wb.append((i, m))
        n_even = len(wb)
        for i in range(6):
            for m in monomials(d - 5, (1 - parz[i]) % 2):
                wb.append((i, m))
        nw = len(wb)
        Wkeys = [(a, b) for a in range(n_even) for b in range(n_even, nw)]
        tpar = (pA + 1) % 2
    else:
        for i in range(6):
            pw = None if pD is None else (pD - parz[i]) % 2
            for m in monomials(d - 5, pw):
                wb.append((i, m))
        nw = len(wb)
        Wkeys = [(a, b) for a in range(nw) for b in range(a, nw)]
        # target degree 2d, parity: A (odd) * even(D^2-like) -> parity of A*w_a*w_b*F_ab = pA + par(w_a)+par(w_b)+par(z_a)+par(z_b) = pA + 2pD = pA
        tpar = None if pD is None else pA
    Wpolys = []
    for (a, b) in Wkeys:
        ia, ma = wb[a]; ib, mb = wb[b]
        t = mul(Ap, shift(Fij[ia][ib], tuple(x + y for x, y in zip(ma, mb)), 1, P), P)
        if a != b:
            t = scale(t, 2, P)
        Wpolys.append(t)
    cols_list = monomials(2 * d, tpar)
    col = {e: j for j, e in enumerate(cols_list)}
    C = len(cols_list)
    Jrows = []
    for i in range(6):
        mpar = None if tpar is None else (tpar - parz[i]) % 2
        for m in monomials(2 * d - 5, mpar):
            Jrows.append(shift(Fi[i], m, 1, P))
    RJ, RW = len(Jrows), len(Wpolys)
    flat = [0] * ((RJ + RW) * C)
    for r, poly in enumerate(Jrows + Wpolys):
        base = r * C
        for e, v in poly.items():
            flat[base + col[e]] = v
    MJ = nmod_mat(RJ, C, flat[:RJ * C], P); rJ = MJ.rank(); del MJ
    M = nmod_mat(RJ + RW, C, flat, P); del flat
    rM = M.rank()
    dimK = RW - (rM - rJ)
    hilbJ = dim_J(2 * d, tpar) if tpar is not None else dim_J(2 * d)
    rec = dict(exp='hessrelax', d=d, parity=PAR, prime=P, nw=nw, nW=RW, cols=C, rowsJ=RJ, rankJ=rJ, hilbJ=hilbJ,
               J_rank_certified=(rJ == hilbJ), rank_stack=rM, dim_KH_modp=dimK, job=JOB, seconds_rank=round(time.time() - t0, 1))
    rec['rigorous_upper_bound'] = rec['J_rank_certified']
    print(json.dumps(rec), flush=True)
    names = [f"{i}:{''.join((v + ('^%d' % e if e > 1 else '')) for v, e in zip('X0 X1 X2 Y0 Y1 Y2'.split(), m) if e) or '1'}" for (i, m) in wb]
    if dimK > 0:
        Mt = M.transpose(); del M
        X, nul = Mt.nullspace(); del Mt
        cvecs = []
        for k in range(nul):
            v = [int(X[RJ + r, k]) for r in range(RW)]
            if any(v):
                cvecs.append(v)
        Cm = nmod_mat(len(cvecs), RW, [x for v in cvecs for x in v], P)
        Cr, rk = Cm.rref()
        basis = []
        for i in range(rk):
            row = {}
            for a in range(RW):
                v = int(Cr[i, a])
                if v:
                    row[a] = v if v <= P // 2 else v - P
            basis.append(row)
        rec['dim_KH_basis'] = rk
        json.dump({'d': d, 'parity': PAR, 'p': P, 'wb': [[i, list(m)] for (i, m) in wb], 'Wkeys': Wkeys,
                   'basis': [{str(a): v for a, v in row.items()} for row in basis]},
                  open(os.path.join(w2x.LANE, f'wave2-excl-hess-{JOB}-d{d}-par{PAR}-p{P}-KH.json'), 'w'))
        rec['KH_basis'] = [{f'{names[Wkeys[a][0]]}*{names[Wkeys[a][1]]}': v for a, v in row.items()} for row in basis]
        rec['KH_max_abs_coeff'] = max(abs(v) for row in basis for v in row.values())
        for row in basis[:min(rk, 40)]:
            print('  KH basis:', {f'{names[Wkeys[a][0]]}*{names[Wkeys[a][1]]}': v for a, v in row.items()} if len(row) <= 10 else f'{len(row)} terms, maxabs {max(abs(v) for v in row.values())}')
        # annihilator quadrics and msolve (characteristic P must be msolve-compatible, < 2^30)
        Bm = nmod_mat(rk, RW, [int(Cr[i, a]) for i in range(rk) for a in range(RW)], P)
        Lam, nlam = Bm.nullspace()
        quadrics = []
        for k in range(nlam):
            terms = []
            for idx, (a, b) in enumerate(Wkeys):
                c = int(Lam[idx, k])
                if c:
                    terms.append(f'{c}*w{a}*w{b}')
            quadrics.append('+'.join(terms))
        rec['n_quadrics'] = nlam
        if MSOLVE and P < 2 ** 30:
            base = os.path.join(w2x.LANE, f'wave2-excl-hess-{JOB}-d{d}-par{PAR}-p{P}')
            out = {}
            def run(tag, extra):
                path = base + f'-{tag}.ms'
                with open(path, 'w') as f:
                    f.write(','.join(f'w{a}' for a in range(nw)) + f'\n{P}\n' + ',\n'.join(quadrics + extra) + '\n')
                opath = path.replace('.ms', '.out')
                try:
                    r = subprocess.run(['msolve', '-v', '0', '-t', '2', '-P', '1', '-f', path, '-o', opath], capture_output=True, text=True, timeout=900)
                    o = open(opath).read() if os.path.exists(opath) else ''
                    return {'rc': r.returncode, 'out_head': o[:200], 'out_len': len(o), 'stderr': r.stderr[-200:]}
                except subprocess.TimeoutExpired:
                    return {'timeout': 900}
            out['homogeneous'] = run('hom', [])
            random.seed(3); r1 = [random.randint(1, 97) for _ in range(nw)]
            out['chart1'] = run('chart1', ['+'.join(f'{c}*w{a}' for a, c in enumerate(r1)) + f'+{P-1}'])
            rec['msolve'] = out
    else:
        rec['result'] = f'K_H = 0: NO polar denominator of degree {d} (sector {PAR}) -- rigorous (Hilbert-certified)'
    rec['seconds'] = round(time.time() - t0, 1)
    logline(LOG, rec)
    print(json.dumps({k: v for k, v in rec.items() if k != 'KH_basis'}), flush=True)


if __name__ == '__main__':
    main()
