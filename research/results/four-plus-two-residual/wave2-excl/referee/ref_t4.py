"""Direct fixed-D tests with the independent library at the unused prime 10^9+9, with the exact Koszul nullity ([RANK]):
B, H (Thom-Sebastiani summands, d = 6), random members of the 10-param and 8-param eigenvector Jacobian families (d = 6),
and the degree-6 core-kernel representatives #0 of each parity (regenerated with the agent's deterministic small_kernel_reps,
then re-verified exactly here)."""
import sys, time, random, json; sys.path.insert(0, '/tmp/hodge/work/referee42'); sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import reflib as R
from fractions import Fraction
LOG = '/tmp/hodge/work/referee42/out/ref_t4.jsonl'
PNEW = 1000000009
F, A = R.load_instance()
H = {e: c for e, c in F.items() if e[2] == 0 and e[5] == 0}; B = {e: c for e, c in F.items() if e[2] + e[5] == 6}
Fi = [R.pderiv(F, i) for i in range(6)]
def run(tag, D, pf=False):
    r = R.fixed_D_test(F, A, D, PNEW, parity_free=pf); r['D'] = tag; r['exp'] = 'fixedD-independent'; R.log(LOG, r)
run('B(X2,Y2)', B); run('H(X0,X1,Y0,Y1)', H)
rng = random.Random(2026)
# 10-param family: w_a = sum_{c in class(a)} w_ac z_c, classes {X0,X1},{Y0,Y1},{X2},{Y2}; avoid the Euler line
pairs10 = [(a, c) for cl in ((0, 1), (3, 4), (2,), (5,)) for a in cl for c in cl]
for t in range(2):
    D = {}
    co = {P: rng.randint(-4, 4) for P in pairs10}
    for (a, c), x in co.items():
        if x: D = R.padd(D, R.pscale(R.pmul(R.z(c), Fi[a]), Fraction(x)))
    run(f'10-param member w={co}', D)
pairs8 = [(a, c) for a in (0, 1) for c in (3, 4)] + [(a, c) for a in (3, 4) for c in (0, 1)]
for t in range(2):
    D = {}
    co = {P: rng.randint(-4, 4) for P in pairs8}
    for (a, c), x in co.items():
        if x: D = R.padd(D, R.pscale(R.pmul(R.z(c), Fi[a]), Fraction(x)))
    run(f'8-param member w={co}', D)
# degree-6 core-kernel reps via the agent's deterministic generator, re-verified exactly here
import w2lib as W
for par in (0, 1):
    t0 = time.time()
    r = W.small_kernel_reps(W.H, W.CORE, W.A, 6, par)
    reps = r['reps_modJ']
    D = {e: Fraction(c) for e, c in reps[0].items()}
    J12 = R.jac_cols(H, R.CORE, 12, (par + 1) % 2); rJ = R.rank_exact(J12)
    inker = R.rank_exact(J12 + [R.pmul(A, D)]) == rJ
    J6 = R.jac_cols(H, R.CORE, 6, par); notinJ = R.rank_exact(J6 + [D]) == R.rank_exact(J6) + 1
    R.log(LOG, {'exp': 'ker6-rep-check', 'parity': par, 'n_reps': len(reps), 'rep0_terms': len(D), 'max_abs_coeff': max(abs(int(c)) for c in D.values()),
                'A*rep0_in_J(H)_12_exact': inker, 'rep0_not_in_J(H)_6_exact': notinJ, 'seconds': round(time.time() - t0, 1)})
    run(f'ker6 rep #0 par{par} (core-only kernel element)', D)
    if len(reps) > 1:
        D2 = {}
        for c_ in reps[:5]:
            D2 = R.padd(D2, R.pscale({e: Fraction(v) for e, v in c_.items()}, Fraction(rng.randint(-3, 3))))
        if D2: run(f'random combo of first 5 ker6 reps par{par}', D2)
