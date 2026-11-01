"""Experiment 0 + 5: sanity (A in J? D=1), sector dimensions, and the E1 kernel sweep
dim{D in (S_d)_kappa : A D in J} vs dim (J_d)_kappa for d = 0..12 and every character kappa."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, numpy as np
from hodge import core
from flint import nmod_mat

F, A, q, rec = core.load_instance('cycle-q2')
p = eq.PRIME
Gw = eq.symmetry_group(F)
gens = eq.generators(Gw); ng = gens.shape[0]
kA = eq.key_of_poly(A, gens)
kx = [int(eq.keys_of(np.eye(6, dtype=np.int64)[i:i+1], gens)[0]) for i in range(6)]
Fp = core.reduce_mod(F, p)
Fi = [core.pderiv(Fp, i, p) for i in range(6)]
Ap = core.reduce_mod(A, p)
out = {'instance': 'cycle-q2', 'group_order': int(len(Gw)), 'generators': gens.tolist(), 'keyA': kA}

# --- sanity: D = 1 ---
D1 = core.parse('1')
r = eq.test_D(F, A, q, D1, 'sanity', note='D=1: A in J(F)? full system', verbose=True)
out['D1'] = {k: v for k, v in r.items() if not k.startswith('_')}
print('D=1:', out['D1'], flush=True)

# --- sector dimensions of S_k and J_k in A's class and in the invariant class ---
def jac_rank_in_sector(k, kappa):
    """rank of {m F_i : deg m = k-5, key(m)+key(x_i)^{-1}... = kappa} inside (S_k)_kappa"""
    if k < 5: return 0, 0
    allm = eq.monomial_exps(k - 5); km = eq.keys_of(allm, gens)
    cols = []
    for i in range(6):
        # key(m F_i) = key(m) - key(x_i)  (F_i has character chi(x_i)^{-1})
        need = eq.key_add(kappa, kx[i], ng)   # key(m) = kappa + key(x_i)
        sel = allm[km == need]
        for m in sel:
            cols.append(core.shift(Fi[i], tuple(int(x) for x in m), 1, p))
    allk = eq.monomial_exps(k); kk = eq.keys_of(allk, gens)
    rows_exps = allk[kk == kappa]
    nS = len(rows_exps)
    if not cols: return 0, nS
    ridx = {tuple(int(x) for x in e): j for j, e in enumerate(rows_exps)}
    M = np.zeros((nS, len(cols)), dtype=np.int64)
    for j, c in enumerate(cols):
        for e, v in c.items():
            M[ridx[e], j] = v
    Mn = nmod_mat(nS, len(cols), M.ravel().tolist(), p)
    return Mn.rank(), nS

for name, kappa in [('A-class', kA), ('invariant', 0)]:
    for k in [6, 12]:
        rJ, nS = jac_rank_in_sector(k, kappa)
        out[f'sector_{name}_deg{k}'] = {'dim_S': nS, 'dim_J': rJ, 'dim_R': nS - rJ}
        print(f'sector {name} deg {k}: dim S={nS} dim J={rJ} dim R={nS-rJ}', flush=True)

# --- E1 kernel sweep ---
def e1_kernel(d, kappa):
    """K = {D in (S_d)_kappa : A D in J}; returns (dim K, dim (J_d)_kappa, dim (S_d)_kappa, basis of K mod p as exps list)"""
    allD = eq.monomial_exps(d); kD = eq.keys_of(allD, gens)
    Dsel = allD[kD == kappa]
    nD = len(Dsel)
    if nD == 0: return 0, 0, 0, []
    degN = d + 7
    allN = eq.monomial_exps(degN); kN = eq.keys_of(allN, gens)
    cols = []  # polys mod p
    meta = []
    for m in Dsel:
        cols.append(core.shift(Ap, tuple(int(x) for x in m), 1, p)); meta.append(('D', tuple(int(x) for x in m)))
    for i in range(6):
        need = eq.key_add(eq.key_add(kA, kappa, ng), kx[i], ng)
        for m in allN[kN == need]:
            cols.append(core.shift(Fi[i], tuple(int(x) for x in m), p - q, p)); meta.append(('N', i))
    rows = {}
    for c in cols:
        for e in c:
            rows.setdefault(e, len(rows))
    M = np.zeros((len(rows), len(cols)), dtype=np.int64)
    for j, c in enumerate(cols):
        for e, v in c.items():
            M[rows[e], j] = v
    Mn = nmod_mat(len(rows), len(cols), M.ravel().tolist(), p)
    X, nullity = Mn.nullspace()
    if nullity == 0:
        dimK = 0; basis = []
    else:
        proj = nmod_mat(nD, nullity, [int(X[j, kk]) for j in range(nD) for kk in range(nullity)], p)
        T = proj.transpose(); Tr, rk = T.rref()
        dimK = rk
        basis = []
        for i in range(rk):
            poly = {}
            for jj in range(nD):
                v = int(Tr[i, jj])
                if v: poly[tuple(int(x) for x in Dsel[jj])] = v
            basis.append(poly)
    rJ, _ = jac_rank_in_sector(d, kappa)
    return dimK, rJ, nD, basis

sweep = {}
t0 = time.time()
for d in range(0, 13):
    allD = eq.monomial_exps(d); kD = eq.keys_of(allD, gens)
    keys = sorted(set(kD.tolist()))
    tot_nontriv = 0
    per = {}
    for kappa in keys:
        dimK, rJ, nD, basis = e1_kernel(d, kappa)
        nontriv = dimK - rJ
        per[kappa] = {'dimK': dimK, 'dimJ': rJ, 'dimS': nD, 'nontrivial': nontriv}
        if nontriv > 0:
            per[kappa]['basis_mod_p'] = [core.to_str({e: int(v) for e, v in b.items()}) for b in basis[:6]]
        tot_nontriv += nontriv
    sweep[d] = {'num_keys': len(keys), 'total_nontrivial_kernel': tot_nontriv,
                'total_dimS': int(len(allD)), 'total_dimJ': sum(v['dimJ'] for v in per.values()),
                'total_dimK': sum(v['dimK'] for v in per.values()), 'per_key': per}
    print(f'd={d}: keys={len(keys)} dim S_d={len(allD)} dim J_d={sweep[d]["total_dimJ"]} dim K_d={sweep[d]["total_dimK"]} nontrivial kernel (K mod J) = {tot_nontriv}  [{time.time()-t0:.0f}s]', flush=True)
    for kappa, v in per.items():
        if v['nontrivial'] > 0:
            print(f'    key {kappa}: dimS={v["dimS"]} dimJ={v["dimJ"]} dimK={v["dimK"]} -> nontrivial {v["nontrivial"]}; e.g. {v.get("basis_mod_p", [])[:2]}', flush=True)
out['e1_kernel_sweep'] = sweep
json.dump(out, open('/tmp/hodge/results/cycle-q2/exp_a_kernel.json', 'w'), indent=1, default=str)
print('DONE', flush=True)
