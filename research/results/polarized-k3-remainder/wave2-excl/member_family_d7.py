"""Exact (Q) certification that the whole family D = f(U)*ell (ell any linear form) lies in the degree-7 L-relaxation of
the K3-sector target: with a(ell) = (U ell/6, 0), B_D(T a) is a quadratic form in the 6 coefficients of ell; we test the
membership in L_14 of its 21 coefficient polynomials (6 squares + 15 cross terms) exactly over Q (fmpz ranks), per
character block.  Membership of all 21 <=> membership for every ell in C^6."""
import sys, json, time, math
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder/wave2-excl')
from fractions import Fraction as Fr
from flint import fmpz_mat, nmod_mat
from hodge.core import monomials, mul, add, sub, scale, pderiv
from relax import G, P, target, vdeg, z, L_generators
from relax_member import B_of
kind = sys.argv[1] if len(sys.argv) > 1 else 'k3'
T = target(kind)
U = [z(i) for i in range(3)]
ells = [z(i) for i in range(6)]
def a_of(ell): return [scale(mul(U[k], ell), Fr(1, 42)) for k in range(3)] + [{}, {}, {}]
# coefficient polynomials of the quadratic form: B(ell_i + ell_j) - B(ell_i) - B(ell_j) for cross terms
Bs = {}
for i in range(6):
    Bs[(i, i)] = B_of(a_of(ells[i]), T)[1]
for i in range(6):
    for j in range(i + 1, 6):
        Bij = B_of(a_of(add(ells[i], ells[j])), T)[1]
        Bs[(i, j)] = sub(sub(Bij, Bs[(i, i)]), Bs[(j, j)])
gens = L_generators(14)
t0 = time.time()
results = {}
mons_by_r = {r: [m for m in monomials(14) if vdeg(m) % 6 == r] for r in range(6)}
Lmat = {}
for r in range(6):
    mons = mons_by_r[r]; idx = {m: i for i, m in enumerate(mons)}
    gr = [g for (c, g) in gens if c == r]
    rows = [[0] * len(gr) for _ in mons]
    for j, g in enumerate(gr):
        for e, v in g.items(): rows[idx[e]][j] = int(v)
    Lmat[r] = (rows, idx, len(gr))
rankL = {r: fmpz_mat(Lmat[r][0]).rank() for r in range(6) if Lmat[r][2]}
print('L ranks over Q:', rankL, round(time.time() - t0, 1), 's', flush=True)
allin = True
for key, B in Bs.items():
    chs = sorted({vdeg(e) % 6 for e in B})
    ok = True
    for r in chs:
        rows, idx, ng = Lmat[r]
        den = 1
        for e, v in B.items():
            if vdeg(e) % 6 == r: fv = Fr(v); den = den * fv.denominator // math.gcd(den, fv.denominator)
        aug = [row + [0] for row in rows]
        for e, v in B.items():
            if vdeg(e) % 6 == r: aug[idx[e]][ng] = int(Fr(v) * den)
        ra = fmpz_mat(aug).rank()
        ok &= (ra == rankL[r])
    results[str(key)] = {'chars': chs, 'in_L_Q': ok}
    allin &= ok
    print(key, chs, 'in_L_Q =', ok, round(time.time() - t0, 1), 's', flush=True)
out = {'target': kind, 'd': 7, 'family': 'D = f(U) * ell, a = (U ell/6, 0)', 'all_21_coefficients_in_L_Q': allin, 'per_coefficient': results,
       'L_ranks_Q': rankL, 'seconds': round(time.time() - t0, 1)}
json.dump(out, open(f'member-family-d7-{kind}.json', 'w'), indent=1)
print('ALL IN L over Q:', allin)
