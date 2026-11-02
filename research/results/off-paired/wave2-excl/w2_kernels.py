"""Non-Jacobian E1-kernel dimensions per character sector and degree: colon_p(d,c) - dimJ_exact(d,c) (modular colon dim is
an upper bound for the rational one; dim J exact by Koszul).  Also the sector sizes of the bilinear systems.
usage: python3 w2_kernels.py dmin dmax [prime]"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import w2lib
from w2lib import char, char_of, cadd, VARCHAR, monos_char, dim_J_sector, equivariant_hilbert
from hodge import core
from hodge.core import NV, monomials, mul, add, sub, scale, pderiv, shift, reduce_mod
from flint import nmod_mat

def colon_dim_sector(F, A, q, d, c, p):
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p); cA = char_of(Ap)
    degN = d + 6 * q - 5
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    Dmons = monos_char(d, c)
    if not Dmons: return 0, 0
    cols = [shift(Fi[i], m, q, p) for i in range(NV) for m in monos_char(degN, cadd(cadd(c, VARCHAR[i]), cA))]
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    for m in Dmons:
        for e in shift(Ap, m, 1, p): rows.setdefault(e, len(rows))
    R = len(rows); C = len(cols)
    flat = [0] * (R * (C + len(Dmons)))
    W = C + len(Dmons)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * W + j] = v
    for k, m in enumerate(Dmons):
        for e, v in shift(Ap, m, 1, p).items(): flat[rows[e] * W + C + k] = v
    M = nmod_mat(R, W, flat, p)
    rJ = nmod_mat(R, C, [flat[r * W + j] for r in range(R) for j in range(C)], p).rank() if C else 0
    rA = M.rank()
    # colon dim = #D + rank(J part) - rank(all)  (kernel projected to D-coordinates)
    return len(Dmons) + rJ - rA, rJ

if __name__ == '__main__':
    dmin, dmax = int(sys.argv[1]), int(sys.argv[2]); p = int(sys.argv[3]) if len(sys.argv) > 3 else w2lib.P1
    F, A, q = w2lib.load()
    cA = char_of(A)
    out = {}
    for d in range(dmin, dmax + 1):
        t0 = time.time()
        row = {}
        tot = {'colon': 0, 'J': 0, 'nonjac': 0, 'gauge': 0}
        for a in range(6):
            for b in range(6):
                c = (a, b)
                cd, rJ = colon_dim_sector(F, A, q, d, c, p)
                dJ = dim_J_sector(d, c)
                dJ6 = dim_J_sector(d + 6, cadd(c, cA))
                gauge = len(monos_char(d - 6, c)) if d >= 6 else 0
                row[f'{a}{b}'] = {'colon': cd, 'dimJ_d': dJ, 'nonjac': cd - dJ, 'gauge': gauge, 'rankJ_d+6_p': rJ, 'dimJ_d+6_exact': dJ6,
                                  'ring_d': len(monos_char(d, c)) - dJ, 'ring_d+6': equivariant_hilbert(d + 6, cadd(c, cA))}
                tot['colon'] += cd; tot['J'] += dJ; tot['nonjac'] += cd - dJ; tot['gauge'] += gauge
        out[d] = {'sectors': row, 'total': tot, 'seconds': round(time.time() - t0, 1), 'prime': p}
        print(d, tot, f'{time.time()-t0:.1f}s', flush=True)
        print('  nonjac per sector:', {k: v['nonjac'] for k, v in row.items() if v['nonjac']}, flush=True)
        rec = {'exp': 'kernel-dims', 'd': d, 'prime': p, 'total': tot, 'nonjac_per_sector': {k: v['nonjac'] for k, v in row.items()},
               'result': 'table', 'seconds': round(time.time() - t0, 1), 'note': 'colon dim mod p (upper bound), dim J exact (Koszul)'}
        w2lib.log(rec)
    json.dump(out, open(f'{w2lib.PK}/kernel_dims_p{p}.json', 'w'), indent=1)
