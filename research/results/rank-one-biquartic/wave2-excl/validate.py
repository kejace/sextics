"""Positive and negative controls for w2.py (independent E1/E2/relaxation code)."""
import sys, time, json
sys.path.insert(0, '.')
import w2
from hodge import core, fixed_d
from fractions import Fraction
p = 1000000007
out = {}
# --- control 1: split sextic F0 = f + g, A = X0^3 Y0^3, witness D = f, N = A Z_X / 6 (strict) ---
F, A, q, rec = core.load_instance(w2.IID)
f = core.parse('X0**6 + 5*X0**5*X1/7 + 5*X0**2*X1**3*X2/7 + 15*X0*X2**5/7 + X1**6 + 10*X1**5*X2/7 + X2**6')
g = core.parse('Y0**6 + 10*Y0**5*Y1/7 + 10*Y0**3*Y1*Y2**2/7 + 5*Y0*Y2**5/7 + Y1**6 + 15*Y1**5*Y2/7 + Y2**6')
L = core.parse('X0**2 + X0*X1 + 2*X1*X2 + 3*X2**2'); M = core.parse('Y0**2 + Y0*Y1 + 2*Y1*Y2 + 3*Y2**2')
P = core.mul(core.mul(L, L), M)
assert core.sub(F, core.add(core.add(f, g), core.scale(P, Fraction(5, 7)))) == {}
out['decomposition_exact'] = True
F0 = core.add(f, g)
R, F0p, Ap, F0i = w2.load_any(p, F0, A)
fp = R.from_dict(f)
n = w2.e1_lift(R, F0p, Ap, F0i, fp)
assert n is not None and R.apply(n, F0p) == Ap * fp
B = w2.E2Builder(R, F0p, F0i, 6, 0)
r = B.solve(fp, n)
print('control1 split F0, D=f:', r, flush=True)
out['control1_split_D_f'] = r['consistent']
assert r['consistent']
# relaxation: B_f(n) in J_12 ?
red = w2.Reducer(R, F0i, 12, 1)
NF = red.nf_matrix([R.B(fp, n)])
z = all(int(NF[0, j]) == 0 for j in range(red.ring_dim))
print('control1 relaxation vanishes at f:', z, flush=True)
out['control1_relaxation_zero'] = z
assert z
# the relaxation on the whole d=6 even family of F0 must have f as a solution: evaluate quadrics at f's coordinates
fam = w2.e1_family(R, F0p, Ap, F0i, 6, 0, dpred=w2.not_div_X06)
# f reduced mod F0 (f = F0 - g; F0-reduced representative of f is -g + (f - F0) ... just check -g is in the family span):
quad, red2 = w2.relaxation_quadrics(R, F0i, fam, 6, 0, red=red)
# coordinates of D = -g (F0-reduced form of f, since f = F0 - g) in the family basis:
Dmons = fam['Dmons']; basis = fam['basis']
import flint
gp = R.from_dict(g)
target = [( -c) % p for c in w2.Index(6, 0, w2.not_div_X06).vec(gp)]
# solve basis combination
idx = w2.Index(6, 0, w2.not_div_X06)
Mb = flint.nmod_mat(idx.n, len(basis) + 1, [0] * (idx.n * (len(basis) + 1)), p)
flat = [0] * (idx.n * (len(basis) + 1))
for k, (D, _) in enumerate(basis):
    for e, c in D.to_dict().items():
        flat[idx.idx[e] * (len(basis) + 1) + k] = c
for i, c in enumerate(target):
    flat[i * (len(basis) + 1) + len(basis)] = c
Mb = flint.nmod_mat(idx.n, len(basis) + 1, flat, p)
Mr, rk = Mb.rref(); piv = w2.pivots_of_rref(Mr, rk)
assert piv[-1] != len(basis), 'f not in family?!'
a = [0] * len(basis)
for r_, j in enumerate(piv): a[j] = int(Mr[r_, len(basis)])
vals = []
for qd in quad:
    s = 0
    for (k, l), v in qd.items():
        s = (s + v * a[k] * a[l]) % p
    vals.append(s)
print('control1 quadrics at a(f): nonzero count', sum(1 for v in vals if v), 'of', len(quad), flush=True)
out['control1_quadrics_vanish_at_f'] = (sum(1 for v in vals if v) == 0)
assert out['control1_quadrics_vanish_at_f']
# --- control 2: accepted q=1 witnesses of degree 4 and 6 ---
for iid in ['coupled-pair-rational-d4-strict', 'f3-block-x0-quartic-strict', 'completed-square-nambu-holdout', 'block-q2']:
    try:
        Fw, Aw, qw, recw = core.load_instance(iid)
        dw, Dw, Nw, Qw = core.load_witness(iid)
    except Exception as ex:
        print(iid, 'skip', ex); continue
    if qw != 1:
        print(iid, 'q =', qw, 'skip'); continue
    Rw, Fwp, Awp, Fwi = w2.load_any(p, Fw, Aw)
    Dp = Rw.from_dict(Dw)
    pD = Rw.parity(Dp)
    nn = w2.e1_lift(Rw, Fwp, Awp, Fwi, Dp)
    assert nn is not None
    Bw = w2.E2Builder(Rw, Fwp, Fwi, dw, pD)
    rw = Bw.solve(Dp, nn)
    redw = w2.Reducer(Rw, Fwi, 2 * dw, 1 if pD is not None else None, verbose=False)
    NFw = redw.nf_matrix([Rw.B(Dp, nn)])
    zw = all(int(NFw[0, j]) == 0 for j in range(redw.ring_dim))
    print('control2', iid, 'd=', dw, 'parity', pD, 'E2 consistent:', rw['consistent'], 'relaxation zero:', zw, flush=True)
    out[f'control2_{iid}'] = dict(consistent=rw['consistent'], relaxation_zero=zw)
    assert rw['consistent'] and zw
# --- control 3: negative, this instance: D = E_X(F) (in J), compare with engine v2 fast solver ---
R, Fp, Ap, Fi, F, A = w2.load(p)
RX = R.zero()
for i in range(3):
    RX += R.gens[i] * Fi[i]
n = w2.e1_lift(R, Fp, Ap, Fi, RX)
B6 = w2.E2Builder(R, Fp, Fi, 6, 0)
r = B6.solve(RX, n)
print('control3 this instance D = E_X(F): w2 normalized E2 ->', r, flush=True)
out['control3_EX_F_w2'] = r
json.dump(out, open('validate.json', 'w'), indent=1, default=str)
w2.log({'exp': 'validate', 'result': out, 'p': p, 'note': 'positive controls (split F0 with D=f; accepted d=4/6 witnesses) and negative control D=E_X(F)'})
print('VALIDATION OK')
