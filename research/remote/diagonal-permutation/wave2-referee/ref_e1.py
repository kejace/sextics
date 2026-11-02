"""E1 kernels ker([A]. : R_d^par -> R_{d+6}^{1-par}) at a fresh prime, with good-reduction checks in degrees d and d+6.
Also (d = 8, par = 0): is the wave-1 rational D8 (D8_e1kernel.json) E1-admissible mod this prime?
Appends one JSON line per (d, par) to ref_e1_p<p>.jsonl (incremental log)."""
import sys, json, time; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation/wave2-referee')
import ref_common as R
from flint import nmod_mat
import sympy
p = int(sys.argv[1]) if len(sys.argv) > 1 else 999999937
dmin = int(sys.argv[2]) if len(sys.argv) > 2 else 1
dmax = int(sys.argv[3]) if len(sys.argv) > 3 else 8
F, A, q = R.load_instance()
plus, minus = R.koszul_series(30)
kos = lambda n, par: (plus if par == 0 else minus)[n]
parts = [R.pderiv(F, i) for i in range(6)]
Fp = {e: R.modp(c, p) for e, c in F.items()}
Ap = {e: R.modp(c, p) for e, c in A.items()}
partsp = [{e: R.modp(c, p) for e, c in P.items()} for P in parts]
LOG = f'ref_e1_p{p}.jsonl'

def jac_columns(n, par):
    """generators m*F_i of J_n with Y-parity par (F_i Y-even for i<3, Y-odd for i>=3); n >= 5."""
    cols = []
    if n < 5: return cols
    for i in range(6):
        mpar = par if i < 3 else 1 - par
        for m in R.monomials(n - 5, mpar):
            cols.append(R.mul({m: 1}, partsp[i]))
    return cols

def build(rows_mons, cols):
    idx = {m: k for k, m in enumerate(rows_mons)}
    M = nmod_mat(len(rows_mons), len(cols), p)
    for j, P in enumerate(cols):
        for e, c in P.items():
            M[idx[e], j] = c % p
    return M

D8 = None
try:
    d8 = json.load(open('/tmp/hodge/results/diagonal-permutation/D8_e1kernel.json'))
    s = d8['D'] if isinstance(d8.get('D'), str) else None
    if s:
        P = sympy.Poly(sympy.sympify(s), *R.VARS)
        D8 = {tuple(int(e) for e in mon): R.modp(sympy.Rational(c), p) for mon, c in zip(P.monoms(), P.coeffs())}
except Exception as ex:
    print('D8 load failed:', ex)

for d in range(dmin, dmax + 1):
    for par in (0, 1):
        t0 = time.time()
        tpar = 1 - par
        rows = R.monomials(d + 6, tpar)
        colsJ = jac_columns(d + 6, tpar)
        Smons = R.monomials(d, par)
        colsA = [R.mul(Ap, {m: 1}) for m in Smons]
        MJ = build(rows, colsJ); rJ = MJ.rank(); del MJ
        MA = build(rows, colsJ + colsA); rAug = MA.rank(); del MA
        good6 = (len(rows) - rJ == kos(d + 6, tpar))
        colsJd = jac_columns(d, par)
        rJd = build(R.monomials(d, par), colsJd).rank() if colsJd else 0
        goodd = (len(Smons) - rJd == kos(d, par))
        adm = len(Smons) - (rAug - rJ)          # dim {D in S_d^par : A D in J}
        ker = adm - rJd                          # kernel of [A] on R_d^par
        rec = {'exp': 'e1_kernel', 'p': p, 'd': d, 'par': par, 'dim_S_d': len(Smons), 'dim_J_d': rJd, 'ring_dim': len(Smons) - rJd,
               'koszul_d': kos(d, par), 'rows_d6': len(rows), 'rank_J_d6': rJ, 'ring_dim_d6': len(rows) - rJ, 'koszul_d6': kos(d + 6, tpar),
               'good_reduction_d': goodd, 'good_reduction_d6': good6, 'rank_aug': rAug, 'kernel_dim': ker, 'seconds': round(time.time() - t0, 1)}
        if d == 8 and par == 0 and D8 is not None:
            MD = build(rows, colsJ + [R.mul(Ap, D8)]); rD = MD.rank(); del MD
            rec['wave1_D8_terms'] = len(D8); rec['A_D8_in_J_mod_p'] = (rD == rJ)
        print(json.dumps(rec), flush=True)
        open(LOG, 'a').write(json.dumps(rec) + '\n')
print('DONE')
