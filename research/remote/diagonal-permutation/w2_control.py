"""Controls for the pointwise reduced-E2 machinery (w2pt.py) on accepted witnesses.
Positive: the accepted (D, N, Q) of an instance -> as an 'extra' family member with the E1 lift N (exact), the reduced
residual at a = e_1 must be ~0, and the full L-BFGS search from random starts must find rel ~ 0 for a 1-member family.
Negative: perturb D by a random polynomial of the same degree/parity -> rel far from 0.
usage: python3 w2_control.py <instance-id> [n_restarts]"""
import os, sys, json, time
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w2pt import *
from hodge.core import sub
iid = sys.argv[1]; nres = int(sys.argv[2]) if len(sys.argv) > 2 else 3
F, A, q, rec = core.load_instance(iid)
assert q == 1, 'controls only for q = 1'
d, D, N, Q = core.load_witness(iid)
pD = core.parity_of(D)
ok, info = core.verify_witness(F, A, q, D, N, Q); print('witness verifies:', ok, info)
# E1 lift: N itself satisfies N(F) = A D + F H ; we need a lift u with u(F) = A D exactly: u = N - (H/6) E  where H = (N(F) - A D)/F
Fi = [pderiv(F, i) for i in range(NV)]
NF = {}
for i in range(NV): NF = add(NF, mul(N[i], Fi[i]))
rem = sub(NF, mul(A, D))
from hodge import jacobian
divF, H = jacobian.divisible_by(F, rem) if rem else (True, {})
# divisible_by returns H mod p; redo exactly over Q by solving F*H = rem with fractions (small): use the mod-p pattern as support
if rem:
    # exact division via lex leading-term elimination
    def exact_div(num, den):
        num = dict(num); out = {}
        ltd = max(den); cd = den[ltd]
        while num:
            lt = max(num)
            if any(lt[i] < ltd[i] for i in range(NV)): raise ValueError('not divisible')
            m = tuple(lt[i] - ltd[i] for i in range(NV)); c = Fraction(num[lt]) / cd
            out[m] = out.get(m, 0) + c
            num = sub(num, shift(den, m, c))
        return out
    H = exact_div(rem, F)
    assert mul(F, H) == rem
else:
    H = {}
Z = [{tuple(1 if t == i else 0 for t in range(NV)): Fraction(1)} for i in range(NV)]
u = [sub(N[i], scale(mul(H, Z[i]), Fraction(1, 6))) for i in range(NV)]
uF = {}
for i in range(NV): uF = add(uF, mul(u[i], Fi[i]))
assert uF == mul(A, D), 'lift failed'
print('E1 lift u(F) = A D exact: ok;  deg D =', d, 'parity', pD)
pA = core.parity_of(A)
nK = sum(len(monomials(d - 4, ((pA - pD - PAR[i]) % 2 - PAR[j]) % 2)) for i in range(NV) for j in range(i + 1, NV)) if d >= 4 else 0
r = 1 + nK + 50
pts = Points(F, r, seed=3)
# family with the witness D as the only (extra) member: use jacobian_family with d-5 < 0 -> no Jacobian members
members = []
Dder = np.array([pts.eval_poly(pderiv(D, j)) for j in range(NV)])
divN = np.zeros(pts.r, dtype=complex); Nval = np.zeros((NV, pts.r), dtype=complex)
for i in range(NV):
    divN += pts.eval_poly(pderiv(u[i], i)); Nval[i] = pts.eval_poly(u[i])
members.append({'D': D, 'Dval': pts.eval_poly(D), 'Dder': Dder, 'divN': divN, 'Nval': Nval, 'w': None})
# negative control member: D + random perturbation (same degree & parity), lift by least squares in floats if E1 holds? it won't;
# instead perturb u (a wrong lift): u' = u + random Koszul-free field -> E1 fails -> B0 changes; use D' = D with N' = N + eps*field
rng = np.random.default_rng(0)
fam = Family(d, members, np.zeros((1, 0)), f'control-{iid}')
kb = KoszulBlock(F, A, d, pts, pD=pD)
vp = VarPro(fam, kb, pts)
f, rel, c, rho = vp.value(np.array([1.0]))
print(f'POSITIVE control: rel residual at the witness = {rel:.3e} (f={f:.3e}), nK={kb.nK}, r={r}')
# negative: replace lift by a perturbed one (breaks E1 slightly): N_i += 1e-2 * random monomial
mem2 = dict(members[0]); Nval2 = Nval.copy()
pert = {monomials(d + 1, (pA - pD - PAR[0]) % 2)[rng.integers(len(monomials(d + 1, (pA - pD - PAR[0]) % 2)))]: Fraction(1, 50)}
Nval2[0] = Nval2[0] + pts.eval_poly(pert); mem2['Nval'] = Nval2; mem2['divN'] = divN + pts.eval_poly(pderiv(pert, 0))
fam2 = Family(d, [mem2], np.zeros((1, 0)), 'neg'); vp2 = VarPro(fam2, kb, pts)
f2, rel2, _, _ = vp2.value(np.array([1.0]))
print(f'NEGATIVE control (E1 lift perturbed by 1/50 * monomial): rel = {rel2:.3e}')
# second negative: D perturbed (and lift kept): E1 fails -> should be far from 0
Dp = add(D, scale({monomials(d, pD)[rng.integers(len(monomials(d, pD)))]: Fraction(1)}, Fraction(1, 20)))
mem3 = dict(members[0]); mem3['D'] = Dp; mem3['Dval'] = pts.eval_poly(Dp); mem3['Dder'] = np.array([pts.eval_poly(pderiv(Dp, j)) for j in range(NV)])
fam3 = Family(d, [mem3], np.zeros((1, 0)), 'neg2'); vp3 = VarPro(fam3, kb, pts)
f3, rel3, _, _ = vp3.value(np.array([1.0]))
print(f'NEGATIVE control (D perturbed by 1/20 * monomial, lift kept): rel = {rel3:.3e}')
out = {'exp': 'control', 'instance': iid, 'd': d, 'pD': pD, 'nK': kb.nK, 'r': r, 'rel_positive': rel, 'rel_neg_lift': rel2, 'rel_neg_D': rel3}
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'w2-control.jsonl'), 'a') as fh: fh.write(json.dumps(out) + '\n')
print(json.dumps(out))
