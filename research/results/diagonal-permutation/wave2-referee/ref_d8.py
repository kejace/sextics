"""Fresh-prime replay of the D8 relaxation test (lane diagonal-permutation), independent code.
D8 := the wave-1 rational representative (../D8_e1kernel.json, 414 terms) reduced mod p.
Steps: (1) A*D8 in J_14 mod p?  (2) solve u(F) = A*D8 (u_i of degree 9, parity-restricted) mod p;
(3) Psi = D8*div(u) - u(D8) (degree 16, Y-odd); (4) Psi in J_16 mod p?  (with good reduction of J_16^- checked).
By the regular-sequence/Cramer transfer (packet README §1.3), Psi notin J_16 mod p at a good prime excludes every
certificate whose denominator is a p-integral E1-admissible lift of D8 mod p, up to multiples of F.
Receipt: ref_d8_p<p>.json (written incrementally)."""
import sys, json, time, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ref_common as R
from flint import nmod_mat
import sympy
p = int(sys.argv[1]) if len(sys.argv) > 1 else 999999937
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'ref_d8_p{p}.json')
rec = {'p': p, 'steps': {}}
def save(): json.dump(rec, open(OUT, 'w'), indent=1)
t0 = time.time()
F, A, q = R.load_instance()
plus, minus = R.koszul_series(30)
parts = [R.pderiv(F, i) for i in range(6)]
partsp = [{e: R.modp(c, p) for e, c in P.items()} for P in parts]
Ap = {e: R.modp(c, p) for e, c in A.items()}
d8 = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'D8_e1kernel.json')))
P8 = sympy.Poly(sympy.sympify(d8['D']), *R.VARS)
D8 = {tuple(int(e) for e in mon): R.modp(sympy.Rational(c), p) for mon, c in zip(P8.monoms(), P8.coeffs())}
D8 = {e: c for e, c in D8.items() if c}
rec['D8_terms_mod_p'] = len(D8); rec['D8_degree'] = R.degree(D8); save()
ypar = lambda e: (e[3] + e[4] + e[5]) % 2
assert all(ypar(e) == 0 for e in D8)

def jac_columns_with_index(n, par):
    """columns m*F_i of J_n^par together with the unknown labels (i, m)."""
    cols, labels = [], []
    for i in range(6):
        mpar = par if i < 3 else 1 - par
        for m in R.monomials(n - 5, mpar):
            cols.append(R.mul({m: 1}, partsp[i])); labels.append((i, m))
    return cols, labels

def fill(M, rows_idx, cols, offset=0):
    for j, P in enumerate(cols):
        for e, c in P.items():
            M[rows_idx[e], offset + j] = c % p

# ---- step 1+2: E1 system  u(F) = A*D8  in degree 14, Y-odd
rows14 = R.monomials(14, 1); idx14 = {m: k for k, m in enumerate(rows14)}
cols14, labels14 = jac_columns_with_index(14, 1)
b = R.mul(Ap, D8)
M = nmod_mat(len(rows14), len(cols14) + 1, p)
fill(M, idx14, cols14); fill(M, idx14, [b], offset=len(cols14))
t1 = time.time()
Rr, rank_aug = M.rref(); del M
ncol = len(cols14)
pivots = []
for r in range(rank_aug):
    for j in range(ncol + 1):
        if int(Rr[r, j]):
            pivots.append(j); break
consistent = (ncol not in pivots)
rank_J14 = rank_aug - (0 if consistent else 1)
good14 = (len(rows14) - rank_J14 == minus[14])
rec['steps']['e1'] = {'rows': len(rows14), 'unknowns': ncol, 'rank_J14_odd': rank_J14, 'koszul_J14_odd': len(rows14) - minus[14],
                      'good_reduction_14': good14, 'A_D8_in_J_mod_p': consistent, 'expected_nullity_koszul': ncol - rank_J14, 'seconds': round(time.time() - t1, 1)}
save()
print(json.dumps(rec['steps']['e1']), flush=True)
if not consistent:
    rec['verdict'] = 'A*D8 not in J mod p: the wave-1 rational D8 is NOT E1-admissible at this prime (or bad reduction)'; save(); sys.exit(0)
# particular solution: pivot variables take the last-column value
u_coeffs = {}
for r, j in enumerate(pivots):
    v = int(Rr[r, ncol])
    if v: u_coeffs[labels14[j]] = v
del Rr
u = [dict() for _ in range(6)]
for (i, m), v in u_coeffs.items():
    u[i][m] = v
# verify u(F) == A*D8 mod p
lhs = {}
for i in range(6): lhs = R.add(lhs, R.mul(u[i], partsp[i]))
lhs = {e: c % p for e, c in lhs.items() if c % p}
diff = R.add(lhs, {e: (-c) % p for e, c in b.items()}); diff = {e: c % p for e, c in diff.items() if c % p}
rec['steps']['u_check'] = {'u_F_equals_A_D8': (len(diff) == 0), 'nonzero_u_components': sum(1 for ui in u if ui)}
save(); print(json.dumps(rec['steps']['u_check']), flush=True)
assert not diff
# ---- step 3: Psi = D8 * div(u) - u(D8)
divu = {}
for i in range(6): divu = R.add(divu, R.pderiv(u[i], i))
Psi = R.mul(D8, divu)
for i in range(6): Psi = R.add(Psi, R.mul(u[i], R.pderiv(D8, i)), s=-1)
Psi = {e: c % p for e, c in Psi.items() if c % p}
rec['steps']['psi'] = {'degree': R.degree(Psi), 'terms': len(Psi), 'Y_odd': all(ypar(e) == 1 for e in Psi)}
save(); print(json.dumps(rec['steps']['psi']), flush=True)
# ---- step 4: Psi in J_16^- ?
rows16 = R.monomials(16, 1); idx16 = {m: k for k, m in enumerate(rows16)}
cols16, _ = jac_columns_with_index(16, 1)
M = nmod_mat(len(rows16), len(cols16) + 1, p)
fill(M, idx16, cols16); fill(M, idx16, [Psi], offset=len(cols16))
t2 = time.time()
Rr, rank_aug16 = M.rref(); del M
last_is_pivot = False
for r in range(rank_aug16):
    for j in range(len(cols16) + 1):
        if int(Rr[r, j]):
            if j == len(cols16): last_is_pivot = True
            break
del Rr
rank_J16 = rank_aug16 - (1 if last_is_pivot else 0)
good16 = (len(rows16) - rank_J16 == minus[16])
rec['steps']['relax'] = {'rows': len(rows16), 'generators': len(cols16), 'rank_J16_odd': rank_J16, 'koszul_J16_odd': len(rows16) - minus[16],
                         'good_reduction_16': good16, 'Psi_in_J_mod_p': (not last_is_pivot), 'seconds': round(time.time() - t2, 1)}
rec['verdict'] = ('Psi(D8;u) NOT in J mod p (good reduction): no certificate for any p-integral E1-admissible lift of D8 mod p, up to F-multiples'
                  if (last_is_pivot and good16) else ('Psi in J mod p: relaxation passes (no exclusion)' if not last_is_pivot else 'bad reduction'))
rec['seconds'] = round(time.time() - t0, 1); save()
print(json.dumps(rec['steps']['relax'])); print(rec['verdict']); print('DONE')
