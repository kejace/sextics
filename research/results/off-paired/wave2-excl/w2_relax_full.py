"""Full degree-2 relaxation of a bilinear sector system: lambda in the left kernel of [M_0 | ... | M_{nt-1}] (s-parts), quadrics
lambda^T q(t).  Homogeneous in t; m-primary mod p (pure power of each variable among GB leading monomials) => exclusion of the
sector over Q, provided the rational lambda, D_a, N_a, S_k are verified (see w2_relax.py docstring)."""
import sys, os, json, time, subprocess, re
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import w2lib, w2_bilinear as WB, w2_relax as WR
from flint import nmod_mat

def full_relaxation(sysd):
    p = sysd['p']; nt, ns = sysd['nt'], sysd['ns']; eqs = sysd['eqs']; neq = len(eqs)
    # C_y^T : (nt*ns) x neq
    flat = [0] * (nt * ns * neq)
    for r, eq in enumerate(eqs):
        for (k1, i1, k2, i2), cf in eq.items():
            if k2 == 's':
                flat[(i1 * ns + i2) * neq + r] = cf
    M = nmod_mat(nt * ns, neq, flat, p) if nt * ns else None
    if M is None:
        lams = [[1 if j == r else 0 for j in range(neq)] for r in range(neq)]
    else:
        X, nul = M.nullspace()
        lams = [[int(X[r, k]) for r in range(neq)] for k in range(nul)]
    quads = []
    for lam in lams:
        Qd = {}
        for r, eq in enumerate(eqs):
            if not lam[r]: continue
            for (k1, i1, k2, i2), cf in eq.items():
                if k2 == 't':
                    key = (i1, i2)
                    Qd[key] = (Qd.get(key, 0) + lam[r] * cf) % p
        Qd = {k: v for k, v in Qd.items() if v}
        if Qd: quads.append((lam, Qd))
    # independent quadrics only (rref over the monomial basis)
    keys = sorted({k for _, Qd in quads for k in Qd})
    kidx = {k: i for i, k in enumerate(keys)}
    if quads:
        A = nmod_mat(len(quads), len(keys), [Qd.get(k, 0) for _, Qd in quads for k in keys], p)
        rank = A.rank()
    else:
        rank = 0
    return lams, quads, rank

if __name__ == '__main__':
    d = int(sys.argv[1]); chars = [(int(s[0]), int(s[1])) for s in sys.argv[2].split(',')]; job = sys.argv[3]
    timeout = int(sys.argv[4]) if len(sys.argv) > 4 else 120
    F, A, q = w2lib.load(); p = WB.PM
    outdir = f'{w2lib.PK}/relaxfull_{job}'; os.makedirs(outdir, exist_ok=True)
    for c in chars:
        t0 = time.time()
        sysd, info = WB.build_fast(F, A, q, d, c, p, 'sector', 'all', addJ=False, verbose=False)
        if sysd is None:
            print(json.dumps({'d': d, 'char': c, 'result': 'empty-t-space'})); continue
        lams, quads, rank = full_relaxation(sysd)
        nt = sysd['nt']
        names = [f't{i}' for i in range(nt)]
        polys = []
        for _, Qd in quads:
            polys.append('+'.join(f'{cf}*t{i1}*t{i2}' if i1 != i2 else f'{cf}*t{i1}^2' for (i1, i2), cf in Qd.items()))
        path = f'{outdir}/d{d}_c{c[0]}{c[1]}.ms'
        open(path, 'w').write(','.join(names) + '\n' + str(p) + '\n' + ',\n'.join(polys) + '\n')
        out = path.replace('.ms', '.gb')
        st = 'ok'
        try:
            subprocess.run(['msolve', '-g', '2', '-f', path, '-o', out, '-t', '1', '-v', '0'], capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            st = 'timeout'
        rec = {'exp': 'relax-full-gb', 'job': job, 'd': d, 'char': c, 'prime': p, 'nt': nt, 'ns': sysd['ns'], 'neq': len(sysd['eqs']),
               'relax_dim': len(lams), 'n_quadrics': len(quads), 'quadric_rank': rank}
        if st == 'ok' and os.path.exists(out):
            gpolys, lead = WR.parse_gb(out, nt)
            pure = {}
            for l in lead:
                nz = [i for i in range(nt) if l[i]]
                if len(nz) == 1: pure[nz[0]] = min(pure.get(nz[0], 10**9), l[nz[0]])
            mprimary = all(i in pure for i in range(nt))
            rec.update(gb_size=len(gpolys), pure_powers=[pure.get(i) for i in range(nt)], mprimary=mprimary,
                       result='mprimary-modp' if mprimary else 'not-mprimary')
        else:
            rec.update(result='msolve-' + st)
        rec['seconds'] = round(time.time() - t0, 1)
        w2lib.log(rec); print(json.dumps(rec), flush=True)
