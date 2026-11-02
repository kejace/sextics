"""Hessian relaxation for MIXED-parity Jacobian denominators D = w(F), w = w_+ + w_- (both parity sectors at once).
The parity WLOG (D -> D*tau D) is a WLOG for existence only, not for degree-by-degree exclusion (pointed out by the verify agent),
so the degree bound must also cover mixed w.  Condition: A*Hess_F(w,w) in J_{2d}; with Hess(w,w) = Hess(w+,w+) + Hess(w-,w-) (even)
+ 2 Hess(w+,w-) (odd).  For d <= 7 (mult. by A injective on R_{2d-6}^{+-}): even part in J_{2d-6}^+ and odd part in J_{2d-6}^-.
Variables: a complement of T_d = {h Euler} in W_d^+ (+) W_d^-.  Usage: python3 w2x_hess_mixed.py d [p] [tag] [msolve_threads]
"""
import sys, time, json, subprocess, re
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl')
from w2x_lib import *
from w2x_hess import field_basis, euler_mult_basis, complement_basis, write_msolve

def build_mixed(d, p):
    t0 = time.time()
    Fp, Fi = partials(p)
    Fij = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
    Bp, Bm = field_basis(d, 0), field_basis(d, 1)
    B = Bp + Bm; nplus = len(Bp)
    # T_d in the joint coordinates: h*Euler with h even -> lives in the + block, h odd -> in the - block
    T = []
    for par, off in ((0, 0), (1, nplus)):
        for v in euler_mult_basis(d, par):
            T.append({k + off: c for k, c in v.items()})
    comp = complement_basis(B, T, p)
    n = len(comp)
    deg = 2 * d - 6
    red_e = Reducer(Fi, deg, 0, p); red_o = Reducer(Fi, deg, 1, p)
    pairs = [(a, b) for a in range(n) for b in range(a, n)]
    def hess(a, b):
        i, m = B[comp[a]]; j, m2 = B[comp[b]]
        return shift(Fij[i][j], tuple(x + y for x, y in zip(m, m2)), 1, p)
    blocks = []
    for red in (red_e, red_o):
        C = len(red.mons)
        flat = [0] * (len(pairs) * C)
        for r, (a, b) in enumerate(pairs):
            h = hess(a, b)
            if a != b: h = scale(h, 2, p)
            for e, v in h.items():
                if e in red.idx: flat[r * C + red.idx[e]] = v     # only the matching-parity part lands in this reducer
        V = nmod_mat(len(pairs), C, flat, p); del flat
        sub = lambda M, rows, cols: nmod_mat(len(rows), len(cols), [int(M[i, j]) for i in rows for j in cols], p)
        NF = sub(V, range(len(pairs)), red.standard) - sub(V, range(len(pairs)), red.pivots) * nmod_mat(red.rank, len(red.standard), [c for row in red.rows_std for c in row], p)
        blocks.append(NF)
    # concatenate columns: even coordinates then odd coordinates
    NFe, NFo = blocks
    tot = nmod_mat(len(pairs), NFe.ncols() + NFo.ncols(), [int(NFe[r, c]) for r in range(len(pairs)) for c in range(NFe.ncols())] if False else [0] * (len(pairs) * (NFe.ncols() + NFo.ncols())), p)
    ce, co = NFe.ncols(), NFo.ncols()
    flat = []
    for r in range(len(pairs)):
        flat.extend(int(NFe[r, c]) for c in range(ce)); flat.extend(int(NFo[r, c]) for c in range(co))
    tot = nmod_mat(len(pairs), ce + co, flat, p)
    rk = tot.rank()
    info = {'d': d, 'parD': 'mixed', 'p': p, 'dimW': len(B), 'dimT': len(T), 'dimWprime': n, 'dimSym2': len(pairs), 'deg_hess': deg,
            'rank_p_J_even': red_e.rank, 'exact_J_even': J_dim(deg, 0), 'rank_p_J_odd': red_o.rank, 'exact_J_odd': J_dim(deg, 1),
            'dim_R_even': ce, 'dim_R_odd': co, 'rank_linearization': rk, 'linearization_proves_empty': rk == len(pairs), 't_build': round(time.time() - t0, 1)}
    return info, tot, pairs, n

if __name__ == '__main__':
    d = int(sys.argv[1]); p = int(sys.argv[2]) if len(sys.argv) > 2 else P_NEW
    tag = sys.argv[3] if len(sys.argv) > 3 else 'local'; thr = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    info, NF, pairs, n = build_mixed(d, p)
    qdim = None
    if not info['linearization_proves_empty']:
        path = f'{PK}/hessmixed_d{d}_p{p}.ms'
        info['n_quadrics_written'] = write_msolve(d, 'mixed', p, NF, pairs, n, path); info['msolve_input'] = path
        if d <= 6:
            lg = path.replace('.ms', '.msolve.log')
            subprocess.run(f'msolve -v 1 -g 1 -t {thr} -f {path} -o {path.replace(".ms", ".lm")} > {lg} 2>&1', shell=True, timeout=900)
            mm = re.search(r'Dimension of quotient:\s*(\d+)', open(lg).read())
            qdim = int(mm.group(1)) if mm else None
            info['msolve_quotient_dim'] = qdim
    print(json.dumps(info), flush=True)
    with open(f'{PK}/hess_mixed' + ('' if tag == 'local' else f'.{tag}') + '.jsonl', 'a') as f: f.write(json.dumps(info) + '\n')
    wlog({'exp': 'b-hess-mixed-parity', 'family': f'J_{d} mixed-parity fields w = w+ + w-', 'd': d, 'prime': p,
          'result': ("V'={0} by linearization (rigorous)" if info['linearization_proves_empty'] else f"msolve quotient dim {qdim}" + (" => V'={0} (rigorous)" if qdim is not None else '; msolve input written')),
          'rank': info['rank_linearization'], 'seconds': info['t_build'], 'job': tag,
          'note': f"dim W={info['dimW']}, dim T={info['dimT']}, dim W'={n}, Sym2={len(pairs)}, quadrics: {info['dim_R_even']} even + {info['dim_R_odd']} odd coordinates"},
         fn=None if tag == 'local' else f'{PK}/wave2-excl-remote-{tag}.jsonl')
