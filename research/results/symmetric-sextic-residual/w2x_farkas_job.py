"""Exact rational inconsistency (Farkas) certificates for named denominators via the exact twisted system,
unfolded to a certificate for the FULL unreduced fixed-D system, plus an independent exact verification."""
import sys, time, json, pickle, math
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import numpy as np
import w2x_lib as L
import charsolve as CS
from hodge import core
from fractions import Fraction
F, A, q = L.load()
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
X, Y = L.X, L.Y
DX, DY = L.dlt(X), L.dlt(Y)
QX, QY, e1X, e1Y, e3X, e3Y = L.e_k(X,2), L.e_k(Y,2), L.e_k(X,1), L.e_k(Y,1), L.e_k(X,3), L.e_k(Y,3)
e2sqX = sub(pw(QX, 2), sc(mul(e1X, e3X), 2)); e2sqY = sub(pw(QY, 2), sc(mul(e1Y, e3Y), 2))
D7, D7Y = mul(DX, e2sqX), mul(DY, e2sqY)
Fi = [core.pderiv(F, i) for i in range(6)]
def radial(idx):
    out = {}
    for i in idx: out = add(out, mul(L.var(i), Fi[i]))
    return out
EXF, EYF = radial((0,1,2)), radial((3,4,5))
cands = [('D7=DX*e2(X^2)', D7, 'sgnX'), ('D7Y=DY*e2(Y^2)', D7Y, 'sgnY'), ('EX(F)', EXF, 'triv'), ('EY(F)', EYF, 'triv'),
         ('EX(F)-EY(F)', sub(EXF, EYF), 'triv'), ('DX*DY', mul(DX, DY), 'sgnsgn'), ('DX*DY*QX', mul(mul(DX,DY),QX), 'sgnsgn'), ('DX*DY*e1X^2', mul(mul(DX,DY),pw(e1X,2)), 'sgnsgn'),
         ('DX*e1X*e2(X^2)=D7*e1X', mul(D7, e1X), 'sgnX'), ('DX*(e2e3-7e1^3e2+4e1^5)', mul(DX, add(sub(mul(QX,e3X), sc(mul(pw(e1X,3),QX),7)), sc(pw(e1X,5),4))), 'sgnX')]
only = sys.argv[1:]
def full_check(y_orbit, S, Dz, label):
    """unfold the folded Farkas vector to the full monomial space and verify y.M_full = 0, y.b_full != 0 exactly, column by column,
    without using the folded matrix: columns of the full system are the coefficient polynomials of each unknown."""
    G = S.G; Fz, Az = S.F, S.A
    d = S.d; degN = d + 1; degQ = 2*d - 6; degE = 2*d + 6; pA = 1; pD = core.parity_of(Dz)
    # y_full[monomial] = y_orbit[orbit(monomial)] * sign (twisted fold uses no sign for the identity E because E is H-invariant: plain fold)
    rowmons = core.monomials(degE, pA)
    E = np.array(rowmons, dtype=np.int64)
    idx = np.searchsorted(S.orbit_codes, G.canon(E))
    yfull = {m: y_orbit[i] for m, i in zip(rowmons, idx.tolist()) if y_orbit[i]}
    def dot(P): return sum(yfull.get(e, 0) * int(c) for e, c in P.items())
    Fi = [core.pderiv(Fz, i) for i in range(6)]; Di = [core.pderiv(Dz, i) for i in range(6)]; FD = core.mul(Fz, Dz)
    bad = 0; ncols = 0
    for s in range(6):
        G_s = core.sub(core.scale(core.mul(Dz, Fi[s]), -1), core.mul(Fz, Di[s]))
        for m in core.monomials(degN, CS.pN_of(pA, pD, s)):
            P = core.shift(G_s, m, 1)
            if m[s]:
                m2 = list(m); m2[s] -= 1
                P = core.add(P, core.shift(FD, tuple(m2), m[s]))
            if dot(P) != 0: bad += 1
            ncols += 1
    if degQ >= 0:
        F2 = core.scale(core.mul(Fz, Fz), -1)
        for m in core.monomials(degQ, pA):
            if dot(core.shift(F2, m, 1)) != 0: bad += 1
            ncols += 1
    yb = dot(core.scale(core.mul(Az, core.mul(Dz, Dz)), -1))
    return bad, ncols, yb, len(yfull)
for label, D, chi in cands:
    if only and not any(o in label for o in only): continue
    t0 = time.time()
    G = L.group(chi)
    S = L.ExactCharSystem(F, A, q, D, G, verbose=True)
    r, r_aug = S.ranks()
    rec = {'exp': 'farkas-exact', 'D': label, 'd': S.d, 'chi': chi, 'rows': S.nrows, 'cols': S.ncols, 'rank_Q': r, 'rank_aug_Q': r_aug, 'prime': None, 'job': 'farkas'}
    if r_aug > r:
        y = S.farkas()
        ok, yb = S.check_farkas(y)
        bad, ncols, yb_full, nnz = full_check(y, S, S.D, label)
        rec.update({'result': 'inconsistent over Q (rigorous)', 'farkas_folded_ok': ok, 'y_dot_b': str(yb), 'full_system_cols_checked': ncols, 'full_cols_violating': bad, 'full_y_dot_b': str(yb_full), 'y_support': nnz,
                    'max_abs_y': max(abs(v) for v in y)})
        pickle.dump({'label': label, 'chi': chi, 'D': {str(k): str(v) for k, v in D.items()}, 'scaleF': S.lf, 'scaleD': S.ld, 'y_orbit': y, 'orbit_codes': S.orbit_codes.tolist(), 'degE': 2*S.d+6}, open(f'{L.PK}/farkas_{label.split("=")[0].replace("*","x").replace("(","").replace(")","").replace("+","p").replace("-","m").replace("^","")}.pkl', 'wb'))
    else:
        rec['result'] = 'consistent over Q (!!)'
    rec['seconds'] = round(time.time() - t0, 1)
    L.log(rec); print(json.dumps(rec), flush=True)
