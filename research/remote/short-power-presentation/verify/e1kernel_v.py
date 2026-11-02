"""Independent E1-kernel recomputation: dim {D in S_d^{par}: A D in J_{d+6}} mod J_d, with a prime
different from wave 1.  Reports rank J_{d+6} vs Hilbert expectation so the result is rigorous over Q."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/verify')
from vlib import *
from flint import nmod_mat, ctx
p = int(sys.argv[1]); jobs = sys.argv[2:]
ctx.threads = 2
Fs, As, q = sympy_F()
F = poly_from_sympy(Fs); A = poly_from_sympy(As)
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
H = hilb()
for j in jobs:
    d = int(j[:-1]); par = 0 if j[-1] == 'e' else 1
    t0 = time.time()
    monsD = monomials(d, par)
    tpar = (par + 1) % 2
    J = jac_gens(Fp, d + 6, tpar, p)
    AD = [shift(Ap, m, 1, p) for m in monsD]
    rows = {}
    for P in J + AD:
        for e in P: rows.setdefault(e, len(rows))
    MJ, _ = nmod_from_polys(J, p, rows); rJ = MJ.rank(); del MJ
    MA, _ = nmod_from_polys(J + AD, p, rows); rA = MA.rank(); del MA
    kdim = len(monsD) - (rA - rJ)
    JD = jac_gens(Fp, d, par, p)
    rJD = nmod_from_polys(JD, p)[0].rank() if JD else 0
    expJ = H[d + 6]['J_minus' if tpar else 'J_plus']; expJD = H[d]['J_minus' if par else 'J_plus']
    out = dict(check='e1kernel', d=d, parD=par, p=p, dimS=len(monsD), rankJ_target=rJ, expected_rankJ_target=expJ,
               rank_match=(rJ == expJ), dimJ_d=rJD, expected_dimJ_d=expJD, kernel_incl_J=kdim, kernel_mod_J=kdim - rJD,
               rows=len(rows), cols=len(J) + len(monsD), seconds=round(time.time() - t0, 1))
    print(json.dumps(out), flush=True)
    logjson(claim=f'E1 kernel K_{d}^{"-" if par else "+"} mod J (independent code, prime {p})', verdict='DATA', **out)
