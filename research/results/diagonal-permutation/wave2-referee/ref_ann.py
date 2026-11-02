"""Sum-principle catalogue replay: Ann_6(D) = ker([D]. : R_6^- -> R_{6+e}) at a fresh prime for a sample of factors,
and for Hessian forms H(w,w) of Jacobian-type denominators.  Zero kernel mod p (with good reduction in degree 6+e)
is a zero kernel over Q.  Appends to ref_ann_p<p>.jsonl."""
import sys, json, time; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation/wave2-referee')
import ref_common as R
from flint import nmod_mat
import sympy
p = int(sys.argv[1]) if len(sys.argv) > 1 else 999999937
F, A, q = R.load_instance()
plus, minus = R.koszul_series(30)
parts = [R.pderiv(F, i) for i in range(6)]
partsp = [{e: R.modp(c, p) for e, c in P.items()} for P in parts]
X0, X1, X2, Y0, Y1, Y2 = R.VARS
cat = {
 'p2X': X0**2 + X1**2 + X2**2, 'pXY': X0*Y0 + X1*Y1 + X2*Y2, 'VdmX': (X0-X1)*(X1-X2)*(X2-X0), 'Hesse': X0**3+X1**3+X2**3-3*X0*X1*X2,
 'sX2Y2': X0**2*Y0**2 + X1**2*Y1**2 + X2**2*Y2**2, 'A': X0**3*Y0**3 - X1**3*Y1**3, 'sumT': X0**3*Y0**3 + X1**3*Y1**3 + X2**3*Y2**3,
 'B2': X0**3*Y0**3 + X1**3*Y1**3 - 2*X2**3*Y2**3, 'p6X': X0**6 + X1**6 + X2**6,
}
def todict(expr):
    P = sympy.Poly(sympy.expand(expr), *R.VARS)
    return {tuple(int(e) for e in mon): R.modp(sympy.Rational(c), p) for mon, c in zip(P.monoms(), P.coeffs())}
items = [(k, todict(v)) for k, v in cat.items()]
# Hessian forms of Jacobian-type fields: w = d_0 (D = F_X0, deg 5) -> H = F_00 (deg 4); w = X-Euler -> H = sum X_i X_j F_ij (deg 6)
H00 = {e: R.modp(c, p) for e, c in R.pderiv(parts[0], 0).items()}
HE = {}
for i in range(3):
    for j in range(3):
        HE = R.add(HE, R.mul(R.mul({tuple(1 if k == i else 0 for k in range(6)): 1}, {tuple(1 if k == j else 0 for k in range(6)): 1}), R.pderiv(parts[i], j)))
HE = {e: R.modp(c, p) for e, c in HE.items()}
items += [('H(d0,d0)=F_00 [D=F_X0]', H00), ('H(EX,EX) [D=X-Euler(F)]', HE)]
def jac_columns(n, par):
    cols = []
    for i in range(6):
        mpar = par if i < 3 else 1 - par
        for m in R.monomials(n - 5, mpar):
            cols.append(R.mul({m: 1}, partsp[i]))
    return cols
def build(rows_mons, cols):
    idx = {m: k for k, m in enumerate(rows_mons)}
    M = nmod_mat(len(rows_mons), len(cols), p)
    for j, P in enumerate(cols):
        for e, c in P.items(): M[idx[e], j] = c % p
    return M
S6m = R.monomials(6, 1)
rJ6 = build(S6m, jac_columns(6, 1)).rank()
for label, D in items:
    t0 = time.time()
    e = R.degree(D); dpar = (list(D.keys())[0][3] + list(D.keys())[0][4] + list(D.keys())[0][5]) % 2
    assert all((k[3]+k[4]+k[5]) % 2 == dpar for k in D)
    tpar = (1 + dpar) % 2
    rows = R.monomials(6 + e, tpar); colsJ = jac_columns(6 + e, tpar)
    rJ = build(rows, colsJ).rank()
    rAug = build(rows, colsJ + [R.mul(D, {m: 1}) for m in S6m]).rank()
    good = (len(rows) - rJ == (plus if tpar == 0 else minus)[6 + e])
    adm = len(S6m) - (rAug - rJ); ker = adm - rJ6
    rec = {'exp': 'ann6', 'p': p, 'D': label, 'deg_D': e, 'D_parity': dpar, 'dim_R6_minus': len(S6m) - rJ6, 'good_reduction_6e': good,
           'dim_Ann6': ker, 'seconds': round(time.time() - t0, 1)}
    print(json.dumps(rec), flush=True); open(f'ref_ann_p{p}.jsonl', 'a').write(json.dumps(rec) + '\n')
print('DONE')
