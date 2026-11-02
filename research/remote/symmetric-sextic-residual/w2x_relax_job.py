"""Hessian-relaxation exclusion job.
usage: python3 w2x_relax_job.py --kind full|sector --d D --par 0|1|mixed [--chi triv|sgnX|sgnY|sgnsgn]
                                --prime P --threads T --timeout S [--singular]
Writes wave2-excl/relax_<tag>.{ms,lead,json,sing} and appends a receipt to wave2-excl.log.jsonl."""
import sys, time, json, argparse, subprocess, os
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import numpy as np
import w2x_lib as L
from hodge import core

ap = argparse.ArgumentParser()
ap.add_argument('--kind', default='full'); ap.add_argument('--d', type=int, required=True)
ap.add_argument('--par', default='0'); ap.add_argument('--chi', default='triv')
ap.add_argument('--prime', type=int, default=L.P1); ap.add_argument('--threads', type=int, default=16)
ap.add_argument('--timeout', type=int, default=7200); ap.add_argument('--singular', action='store_true')
ap.add_argument('--nomsolve', action='store_true')
a = ap.parse_args()
F, A, q = L.load(); p = a.prime; d = a.d
tag = f'{a.kind}_{a.chi if a.kind=="sector" else "all"}_d{d}_par{a.par}_p{p}'
out = f'{L.PK}/relax_{tag}'
os.makedirs(L.PK, exist_ok=True)
t0 = time.time()
rec = {'exp': 'hessian-relaxation', 'kind': a.kind, 'd': d, 'par': a.par, 'chi': a.chi if a.kind == 'sector' else None, 'prime': p, 'job': tag}
parities = [0, 1] if a.par == 'mixed' else [int(a.par)]
pairs_all = []; fam_info = []
for pD in parities:
    if a.kind == 'full':
        fam = L.E1Family(F, A, d, pD, p, verbose=True)
        r, ci = fam.rank_check()
        pairs, nr = fam.complement_basis()
        fam_info.append({'pD': pD, 'rows': fam.rows, 'nA': fam.nA, 'nJ': fam.nJ, 'dimK': fam.dimK, 'dimJd': fam.dimJd, 'ci_dimJd': fam.ci_dimJd,
                         'ker_R': fam.ker_R, 'rankJtop_modp': r, 'rankJtop_CI': ci, 'radial': nr, 'complement': len(pairs)})
    else:
        fam = L.SectorFamily(F, A, d, pD, a.chi, p, verbose=True)
        pairs, nr = fam.complement_basis()
        fam_info.append({'pD': pD, 'rows': fam.rows, 'nN': fam.nN, 'nD': fam.nD, 'dimK': fam.dimK, 'dimJd': fam.dimJd, 'ker_R': fam.ker_R,
                         'radial': nr, 'complement': len(pairs), 'rankM': fam.rankM})
    pairs_all += [(pD, D, u) for D, u in pairs]
n = len(pairs_all)
rec['family'] = fam_info; rec['n_unknowns'] = n
print(json.dumps(fam_info), flush=True)
# products, grouped by parity of R_kl: parity(D_k) + parity(D_l) + 1
pairs = [(D, u) for _, D, u in pairs_all]
pars = [pD for pD, _, _ in pairs_all]
indices = {par: L.MonIndex(core.monomials(2 * d, par)) for par in set((pars[k] + pars[l] + 1) % 2 for k in range(n) for l in range(n))}
keys = []; blocks = {par: [] for par in indices}
t1 = time.time()
Darr = [L.arr(D, p) for D, _ in pairs]
gradD = [[L.arr(core.pderiv(D, i, p), p) for i in range(6)] for D, _ in pairs]
divu = [L.arr(L.divergence(u, p), p) for _, u in pairs]
uarr = [[L.arr(u[i], p) for i in range(6)] for _, u in pairs]
def Rkl(k, l, index):
    v = L.mulvec(Darr[k][0], Darr[k][1], divu[l][0], divu[l][1], index, p)
    for i in range(6):
        v = v - L.mulvec(uarr[l][i][0], uarr[l][i][1], gradD[k][i][0], gradD[k][i][1], index, p)
    return v % p
rowpar = []
for k in range(n):
    for l in range(k, n):
        par = (pars[k] + pars[l] + 1) % 2
        idx = indices[par]
        v = Rkl(k, k, idx) if k == l else (Rkl(k, l, idx) + Rkl(l, k, idx)) % p
        blocks[par].append(v); keys.append((k, l)); rowpar.append(par)
print(f'products done {time.time()-t1:.1f}s', flush=True)
quads = []
ranks = {}
for par, rows in blocks.items():
    if not rows: continue
    Pm = np.array(rows, dtype=np.int64)
    kk = [keys[i] for i in range(len(keys)) if rowpar[i] == par]
    if a.kind == 'full':
        red = L.Reducer(F, 2 * d, par, p, verbose=True)
        C = L.coords_from_matrix(red, Pm)
        ranks[par] = {'rank_modp': red.rank, 'rank_CI': red.ci_rank, 'ring_dim': red.ring_dim}
    else:
        fred = L.FoldedReducer(F, 2 * d, par, p, verbose=True)
        C = L.folded_coords_from_matrix(fred, Pm, indices[par])
        ranks[par] = {'rank_modp': fred.rank, 'rank_CI': fred.ci_rank, 'ring_dim': fred.ring_dim, 'ring_CI': fred.ci_ring_dim}
    quads += L.quadrics_from_coords(kk, C)
rec['reducers'] = ranks; rec['n_quadrics'] = len(quads)
print(f'quadrics: {len(quads)} in {n} unknowns ({time.time()-t0:.1f}s)', flush=True)
names = L.write_msolve(out + '.ms', n, quads, p)
L.write_singular(out + '.sing', n, quads, p)
rec['files'] = {'msolve': out + '.ms', 'singular': out + '.sing'}
res = {}
if n == 0:
    res['verdict'] = 'empty family'
elif not a.nomsolve:
    t2 = time.time()
    cmd = ['timeout', str(a.timeout), 'msolve', '-f', out + '.ms', '-o', out + '.lead', '-g', '1', '-t', str(a.threads), '-v', '1']
    r = subprocess.run(cmd, capture_output=True, text=True)
    res['msolve_rc'] = r.returncode; res['msolve_seconds'] = round(time.time() - t2, 1)
    open(out + '.msolve.log', 'w').write(r.stdout + '\n' + r.stderr)
    if r.returncode == 0 and os.path.exists(out + '.lead'):
        pp = L.gb_pure_powers(out + '.lead', names)
        missing = [v for v, k in pp.items() if k is None]
        res['pure_powers'] = pp; res['missing_pure_power'] = missing
        res['verdict'] = 'V={0} (only trivial solution)' if not missing else f'nontrivial solutions possible: {len(missing)} variables without pure power'
        res['max_pure_power'] = max(k for k in pp.values() if k is not None) if any(k is not None for k in pp.values()) else None
    else:
        res['verdict'] = 'msolve failed/timeout'
if a.singular and n > 0:
    t3 = time.time()
    r = subprocess.run(['timeout', str(a.timeout), 'Singular', '-q', out + '.sing'], capture_output=True, text=True)
    res['singular_out'] = r.stdout.strip().replace('\n', ' ')[:500]; res['singular_seconds'] = round(time.time() - t3, 1)
rec['result'] = res; rec['seconds'] = round(time.time() - t0, 1)
json.dump(rec, open(out + '.json', 'w'), indent=1)
L.log(rec)
print(json.dumps(res), flush=True)
