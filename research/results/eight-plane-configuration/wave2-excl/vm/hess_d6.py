"""d=6: Hessian relaxation on (i) pure-parity sectors (independent re-derivation of the other agent's msolve result,
plus flint Macaulay-rank cross-check), (ii) the MIXED-parity 35-dimensional space W' = (linear fields)/Euler."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl/vm')
from hess import *
p = int(sys.argv[1]) if len(sys.argv) > 1 else 2147483587
mode = sys.argv[2] if len(sys.argv) > 2 else 'all'

def linear_fields(par=None):
    """basis of linear fields w (w_i in S_1) with w(F) of parity par (None: all), excluding one Euler direction."""
    out = []
    for i in range(6):
        pi = 0 if i < 3 else 1
        for m in monomials(1):
            if par is not None and (yparity(m) + pi) % 2 != par: continue
            w = [{} for _ in range(6)]; w[i] = {m: 1}
            out.append(w)
    return out

def complement_mod_euler(basis):
    """Euler = sum_i Z_i d_i is in span(basis) iff all six (i, Z_i) are present; drop (5, Z_5) then."""
    keep = []
    dropped = None
    e5 = tuple(0 if k < 5 else 1 for k in range(6))
    for w in basis:
        if w[5] == {e5: 1} and dropped is None:
            dropped = w; continue
        keep.append(w)
    return keep, dropped is not None

for par, label in [(0, 'even'), (1, 'odd'), (None, 'mixed')]:
    if mode != 'all' and mode != label: continue
    basis = linear_fields(par)
    basis, dropped = complement_mod_euler(basis)
    t0 = time.time()
    H = HessianSystem(basis, 6, p, verbose=True)
    rec = {'exp': 'hess-d6', 'd': 6, 'sector': label, 'p': p, 'dimW': len(basis) + (1 if dropped else 0), 'dimT': 1 if dropped else 0,
           'n_vars': H.n, 'deg_hess': H.deg, 'rank_p_J': H.rank_J_p, 'J_exact': H.J_exact, 'R_deg': len(H.red.standard),
           'rank_linearization': H.rank_lin, 'dimSym2': H.dimSym2, 'agent': 'excl-vm'}
    ms = f'{PK}/hess_d6_{label}_p{p}.ms'; out = ms.replace('.ms', '.msout')
    nq = H.write_msolve(ms)
    rec['n_quadrics_written'] = nq
    rec['t_msolve'] = run_msolve(ms, out, threads=1)
    hf, lms = hilbert_from_leading(out, H.n, maxdeg=8)
    rec['hilbert_S_mod_I'] = hf; rec['n_leading_monomials'] = len(lms)
    rec['V_prime_is_zero_mod_p'] = bool(hf and hf[-1] == 0)
    rec['rigorous_over_Q'] = bool(rec['V_prime_is_zero_mod_p'] and H.rank_J_p == H.J_exact)
    if label != 'mixed':
        # flint cross-check: Macaulay ranks at degree 3 and 4
        for k in (3, 4):
            r, C, R = H.macaulay_rank(k)
            rec[f'macaulay_deg{k}'] = {'rank': r, 'dim_S_k': C, 'rows': R, 'full': r == C}
    rec['seconds'] = round(time.time() - t0, 1)
    print(json.dumps(rec), flush=True); log(rec)
