"""Multi-target L-relaxation for the K3 block.  If D serves one K3-sector (3,1)-class with nonzero Λ²T-component then,
by the dichotomy (End_Hdg(T) = Q; and ρ = 1 so that the K3 block IS (Λ²T)^{3,1}), D serves EVERY K3-block class
B(U) - B(V), B in R_{f,6}.  Each such target T_B gives (for the same D = sum P_i a_i) the necessary condition
B_D(T_B a) in L_{2d}.  We impose all of them simultaneously: constraint rows stacked over a basis of R_{f,6}
(19 standard monomials mod p; full rank mod p => a basis over Q as well).  Targets with B in J(f)_6 are zero classes
and must give ZERO constraints -- used as a correctness check of the implementation.
Columns (quadratic monomials c_a c_b) are mu6-graded by char(a)+char(b); constraints from quotient block r only touch
columns of character r, so everything is done block by block.
usage: relax_multi.py d p [--out=PREFIX] [--squares] [--char=b]"""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder/wave2-excl'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import numpy as np, scipy.sparse as sp
from flint import nmod_mat
from hodge.core import monomials
from relax import L_generators, parameters, quadratic_columns, vdeg, square_elimination
from uvring import TernaryJacobianRing, load_f

def build_multi(d, p, char=None, out=None, verbose=True):
    t0 = time.time()
    deg2 = 2 * d
    pars = parameters(d, True, char)
    npar = len(pars)
    gens = L_generators(deg2)
    allmons = monomials(deg2)
    # quotient maps per block
    Ys = {}
    blockinfo = []
    for r in range(6):
        mons = [m for m in allmons if vdeg(m) % 6 == r]
        idx = {m: i for i, m in enumerate(mons)}
        gr = [g for (c, g) in gens if c == r]
        mr, nr = len(mons), len(gr)
        flat = [0] * (mr * nr)
        for j, g in enumerate(gr):
            for e, v in g.items(): flat[idx[e] * nr + j] = v % p
        Lam = nmod_mat(mr, nr, flat, p)
        Yn, nul = Lam.transpose().nullspace()
        Y = np.array([[int(Yn[i, k]) for k in range(nul)] for i in range(mr)], dtype=np.int64)
        Ys[r] = (idx, Y)
        blockinfo.append({'r': r, 'rows': mr, 'L_cols': nr, 'L_rank_p': mr - nul, 'quotient_dim': nul})
    if verbose: print('[multi] quotient maps built', json.dumps(blockinfo), round(time.time() - t0, 1), 's', flush=True)
    # K3 targets: B(U) - B(V) for B monomial of degree 6; standard monomials of R_{f,6} mod p and J-monomials
    f = load_f(p); Rf = TernaryJacobianRing(f, p=p, maxdeg=6)
    std = Rf.data[6]['std']; piv = [Rf.data[6]['mons'][j] for j in Rf.data[6]['pivots']]
    def TB(m):
        if isinstance(m, dict):   # ternary poly with integer coefficients
            T = {}
            for e, v in m.items():
                T[(e[0], e[1], e[2], 0, 0, 0)] = T.get((e[0], e[1], e[2], 0, 0, 0), 0) + v
                T[(0, 0, 0, e[0], e[1], e[2])] = T.get((0, 0, 0, e[0], e[1], e[2]), 0) - v
            return {e: v for e, v in T.items() if v}
        return {(m[0], m[1], m[2], 0, 0, 0): 1, (0, 0, 0, m[0], m[1], m[2]): -1}
    pair_char = None
    results = {'d': d, 'p': p, 'char': char, 'params': npar, 'blocks': blockinfo, 'targets': []}
    stacked = {r: None for r in range(6)}     # per column-character r: reduced row basis (numpy int64)
    cols_by_char = None
    def process(m, label):
        nonlocal cols_by_char
        T = TB(m)
        cols = quadratic_columns(d, T, pars, verbose=False)
        npairs = len(cols)
        if cols_by_char is None:
            pc = [((0 if pars[a][0] < 3 else 5) + vdeg(pars[a][1]) + (0 if pars[b][0] < 3 else 5) + vdeg(pars[b][1])) % 6 for (a, b, q) in cols]
            cols_by_char = {r: [j for j in range(npairs) if pc[j] == r] for r in range(6)}
        ranks = {}
        for r in range(6):
            idx, Y = Ys[r]
            mr = len(idx)
            ri, ci, vals = [], [], []
            colsr = cols_by_char[r]; pos = {j: t for t, j in enumerate(colsr)}
            for j in colsr:
                a, b, q = cols[j]
                for e, v in q.items():
                    if vdeg(e) % 6 == r:
                        ri.append(idx[e]); ci.append(pos[j]); vals.append(int(v) % p)
            if not colsr: continue
            Qm = sp.csr_matrix((np.array(vals, dtype=np.int64), (np.array(ri), np.array(ci))), shape=(mr, len(colsr)))
            QT = Qm.T.tocsr()
            lo = (QT @ (Y & 0xFFFF)) % p; hi = (QT @ (Y >> 16)) % p
            Cr = np.ascontiguousarray(((lo + (hi * (1 << 16)) % p) % p).T)   # nul x ncols_r
            M = nmod_mat(Cr.shape[0], Cr.shape[1], Cr.ravel().tolist(), p)
            rk = M.rank(); ranks[r] = rk
            if stacked[r] is None:
                base = Cr
            else:
                base = np.vstack([stacked[r], Cr])
            Mb = nmod_mat(base.shape[0], base.shape[1], base.ravel().tolist(), p)
            Rb, rb = Mb.rref()
            stacked[r] = np.array([[int(Rb[i, j]) for j in range(base.shape[1])] for i in range(rb)], dtype=np.int64) if rb else np.zeros((0, base.shape[1]), dtype=np.int64)
        cum = {r: (int(stacked[r].shape[0]) if stacked[r] is not None else 0) for r in range(6)}
        rec = {'target': label, 'B': (list(m) if not isinstance(m, dict) else 'poly'), 'rank_per_char': ranks, 'cumulative_rank_per_char': cum, 'cumulative_total': sum(cum.values()),
               'columns_per_char': {r: len(cols_by_char[r]) for r in range(6)}, 'seconds': round(time.time() - t0, 1)}
        results['targets'].append(rec)
        if verbose: print('[multi]', json.dumps(rec), flush=True)
        if out: json.dump(results, open(out + '.json', 'w'))
        return cols
    # correctness check: genuine elements of J(f)_6 (7 z_j d_i f, integer coefficients) are zero classes => zero constraints
    from relax import poly7f
    from hodge.core import pderiv as _pd, mul as _mul
    f7 = {(e[0], e[1], e[2]): v for e, v in poly7f(0).items()}
    def tmul(P, Q):
        out = {}
        for e1, c1 in P.items():
            for e2, c2 in Q.items():
                e = (e1[0]+e2[0], e1[1]+e2[1], e1[2]+e2[2]); out[e] = out.get(e, 0) + c1 * c2
        return {e: v for e, v in out.items() if v}
    def tder(P, i):
        return {tuple(e[k] - (k == i) for k in range(3)): v * e[i] for e, v in P.items() if e[i]}
    for (j, i) in [(0, 0), (1, 2)]:
        Bj = tmul({tuple(int(k == j) for k in range(3)): 1}, tder(f7, i))
        rec0 = process(Bj, f'J-check z{j}*d{i}(7f)')
        last = results['targets'][-1]
        assert sum(last['rank_per_char'].values()) == 0, ('J-check FAILED: zero class gave constraints', last)
    results['J_check'] = 'passed (zero constraints for genuine J(f)_6 targets)'
    for m in std:
        cols = process(m, 'std')
    pairs = [(a, b) for (a, b, q) in cols]
    results['pairs'] = len(pairs)
    # assemble full stacked matrix in original column order for square elimination
    npairs = len(pairs)
    rows = []
    for r in range(6):
        if stacked[r] is None or stacked[r].shape[0] == 0: continue
        full = np.zeros((stacked[r].shape[0], npairs), dtype=np.int64)
        full[:, cols_by_char[r]] = stacked[r]
        rows.append(full)
    Call = np.vstack(rows) if rows else np.zeros((0, npairs), dtype=np.int64)
    results['stacked_rank_total'] = int(Call.shape[0])
    results['seconds'] = round(time.time() - t0, 1)
    if out:
        np.savez_compressed(out + '.npz', C=Call, pairs=np.array(pairs, dtype=np.int32))
        json.dump(results, open(out + '.json', 'w'))
    return results, Call, pairs, pars

if __name__ == '__main__':
    d = int(sys.argv[1]); p = int(sys.argv[2]); out = None; char = None
    for a in sys.argv[3:]:
        if a.startswith('--out='): out = a.split('=', 1)[1]
        if a.startswith('--char='): char = int(a.split('=', 1)[1])
    res, C, pairs, pars = build_multi(d, p, char, out)
    print('[multi-result]', json.dumps({k: v for k, v in res.items() if k not in ('targets', 'blocks')}), flush=True)
    if '--squares' in sys.argv:
        killed, steps = square_elimination(C, pairs, len(pars), p=p)
        res['square_elimination_mod_p'] = {'killed': sorted(killed), 'steps': steps}
        res['parameter_basis'] = [(i, list(m)) for (i, m) in pars]
        if out: json.dump(res, open(out + '.json', 'w'))
        print('[multi-squares]', json.dumps({'killed': len(killed), 'params': len(pars), 'steps': steps}), flush=True)
