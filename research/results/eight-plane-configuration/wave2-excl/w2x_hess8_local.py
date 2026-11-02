"""d = 8 Hessian relaxation on structured subfamilies, using the wave-1 E1-kernel basis at degree 10 (p = 2147483629):
condition  A*Hess_F(w,w) in J_16  <=>  Hess_F(w,w) in K := {D in S_10^+ : A D in J} (the projected kernel, which contains J_10^+).
K is given by 1044 rref vectors (wave-1 pickle); the quotient S_10^+/K has dimension 468 = number of quadrics in w.
Usage: python3 w2x_hess8_local.py <par> <family> [msolve_threads]
"""
import sys, time, json, pickle, re, subprocess, os
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl')
from w2x_lib import *
from w2x_hess import field_basis, euler_mult_basis
from w2x_hess_sub import subfamily, to_vec, write_ms
P1 = 2147483629
MIRROR = '/tmp/hodge/remote/eight-plane-configuration'

def kernel_quotient(par, p):
    """returns (mons, idx, pivots, rows_std, standard) for the subspace K subset S_10^{par} from the wave-1 rref basis."""
    K = pickle.load(open(f'{MIRROR}/e1kernel_basis_d10_par{par}.pkl', 'rb'))
    mons = monomials(10, par); idx = {m: i for i, m in enumerate(mons)}
    R, C = len(K), len(mons)
    flat = [0] * (R * C)
    for r, v in enumerate(K):
        for e, c in v.items(): flat[r * C + idx[e]] = c % p
    M = nmod_mat(R, C, flat, p); Mr, rank = M.rref()
    assert rank == R, (rank, R)
    piv = []
    for i in range(rank):
        j = 0
        while j < C and int(Mr[i, j]) == 0: j += 1
        piv.append(j)
    pivset = set(piv); std = [j for j in range(C) if j not in pivset]
    rows_std = [[int(Mr[i, j]) for j in std] for i in range(rank)]
    return mons, idx, piv, rows_std, std

def build(par, fam, p=P1):
    t0 = time.time()
    d = 8
    Fp, Fi = partials(p)
    Fij = [[pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
    fields = subfamily(fam, d, par, p)
    B = field_basis(d, par); bidx = {b: k for k, b in enumerate(B)}; n = len(B)
    T = euler_mult_basis(d, par)
    U = [to_vec(w, B, bidx, p) for w in fields]
    def mat(vs):
        flat = [0] * (len(vs) * n)
        for r, v in enumerate(vs):
            for k, c in v.items(): flat[r * n + k] = c
        return nmod_mat(len(vs), n, flat, p)
    rT = mat(T).rank() if T else 0
    keep = []; cur = rT; rows = list(T)
    for v in U:
        r = mat(rows + [v]).rank()
        if r > cur: keep.append(v); rows.append(v); cur = r
    m = len(keep)
    comp = []
    for v in keep:
        w = [{} for _ in range(NV)]
        for k, c in v.items():
            i, mo = B[k]; w[i][mo] = c
        comp.append(w)
    mons, idx, piv, rows_std, std = kernel_quotient(0, p)   # Hess(w,w) is always Y-even, whatever the parity of D
    pairs = [(a, b) for a in range(m) for b in range(a, m)]
    C = len(mons)
    flat = [0] * (len(pairs) * C)
    for r, (a, b) in enumerate(pairs):
        h = {}
        for i in range(NV):
            if not comp[a][i]: continue
            for j in range(NV):
                if not comp[b][j]: continue
                h = add(h, mul(mul(comp[a][i], comp[b][j], p), Fij[i][j], p), p)
        if a != b: h = scale(h, 2, p)
        for e, v in h.items(): flat[r * C + idx[e]] = v
    V = nmod_mat(len(pairs), C, flat, p); del flat
    def sub(M, rws, cols):
        return nmod_mat(len(rws), len(cols), [int(M[i, j]) for i in rws for j in cols], p)
    Vstd = sub(V, range(len(pairs)), std); Vpiv = sub(V, range(len(pairs)), piv); del V
    RS = nmod_mat(len(piv), len(std), [c for row in rows_std for c in row], p)
    NF = Vstd - Vpiv * RS      # coordinates in S_10^+ / K  (468)
    rk = NF.rank()
    info = {'d': 8, 'parD': par, 'p': p, 'family': fam, 'dimW': n, 'dimT': rT, 'dimU': mat(U).rank(), 'dimUprime': m, 'dimSym2': len(pairs),
            'n_quadrics_max': len(std), 'rank_linearization': rk, 'linearization_proves_empty': rk == len(pairs), 't_build': round(time.time() - t0, 1)}
    return info, NF, pairs, m

if __name__ == '__main__':
    par = int(sys.argv[1]); fam = sys.argv[2]; thr = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    info, NF, pairs, m = build(par, fam)
    safe = fam.replace(':', '-').replace('^', 'p')
    qdim = None
    if m and not info['linearization_proves_empty']:
        path = f'{PK}/hsub_d8_par{par}_{safe}_p{P1}.ms'
        info['n_quadrics_written'] = write_ms(NF, pairs, m, P1, path); info['msolve_input'] = path
        t1 = time.time()
        lg = path.replace('.ms', '.msolve.log')
        subprocess.run(f'msolve -v 1 -g 1 -t {thr} -f {path} -o {path.replace(".ms", ".lm")} > {lg} 2>&1', shell=True, timeout=int(sys.argv[4]) if len(sys.argv) > 4 else 600)
        mm = re.search(r'Dimension of quotient:\s*(\d+)', open(lg).read())
        qdim = int(mm.group(1)) if mm else None
        info['msolve_quotient_dim'] = qdim; info['t_msolve'] = round(time.time() - t1, 1)
    print(json.dumps(info), flush=True)
    with open(f'{PK}/hess_sub8_local.jsonl', 'a') as f: f.write(json.dumps(info) + '\n')
    wlog({'exp': 'b-hess-subfamily-d8', 'family': f'{fam} subfamily of J_8, parity {par}', 'd': 8, 'prime': P1,
          'result': ('empty: U inside T' if m == 0 else "V'={0} by linearization (rigorous)" if info['linearization_proves_empty'] else f"msolve quotient dim {qdim}: " + ("V'={0} (rigorous)" if qdim is not None else 'msolve did not finish')),
          'rank': info['rank_linearization'], 'seconds': info['t_build'] + info.get('t_msolve', 0), 'job': 'local',
          'note': f"dimU={info['dimU']}, dimU'={m}, Sym2={info['dimSym2']}, quadrics<=468 via wave-1 kernel K_10 (p=2147483629)"})
