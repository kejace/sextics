"""(a) Exact re-verification over Q of the structural claims of attack1 / dossier.
Writes receipt verify_exact.json and one log line.  Everything here is exact (python Fractions / fmpq_mat)."""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl')
from qx import *
from hodge import symmetry
t0 = time.time()
F, A, q = load()
Q = quadrics()
out = {'claims': {}}
# 1. F = sum Q_i^3 exactly
S3 = {}
for q_ in Q: S3 = add(S3, power(q_, 3))
out['claims']['F_eq_sum_Qi3'] = (sub(F, S3) == {})
out['claims']['F_terms'] = len(F); out['claims']['F_Yeven'] = (core.parity_of(F) == 0)
out['claims']['A'] = core.to_str(A); out['claims']['q'] = q
# 2. counts
S, R = hilbert()
out['claims']['S6_odd'] = S[6][1]; out['claims']['S6_even'] = S[6][0]
out['claims']['R6_odd_h31_plus'] = R[6][1]; out['claims']['R6_even_h31_minus'] = R[6][0]
# 3. J_6^odd generators exact rank (18), A not in J_6, delta not in J_6, A,delta independent mod J_6 (exact over Q)
Fi = [pderiv(F, i) for i in range(NV)]
gens = []
for i in range(NV):
    par_i = 0 if i < 3 else 1
    for m in monomials(1, (1 - par_i) % 2):
        gens.append(shift(Fi[i], m, 1))
dl = delta()
mons6o = monomials(6, 1)
M = fmpq_from_polys(gens, mons6o); rJ = M.rank()
MA = fmpq_from_polys(gens + [A], mons6o); rA = MA.rank()
Md = fmpq_from_polys(gens + [dl], mons6o); rd = Md.rank()
MAd = fmpq_from_polys(gens + [A, dl], mons6o); rAd = MAd.rank()
out['claims']['rank_J6_odd_exact'] = rJ
out['claims']['A_notin_J6'] = (rA == rJ + 1)
out['claims']['delta_notin_J6'] = (rd == rJ + 1)
out['claims']['A_delta_independent_mod_J6'] = (rAd == rJ + 2)
out['claims']['delta_terms'] = len(dl); out['claims']['delta_Yodd'] = (core.parity_of(dl) == 1); out['claims']['delta_deg'] = core.degree(dl)
# delta mod 5 = 4 X0X1X2Y0Y1Y2 ?  (coefficients are rationals with denominators prime to 5)
d5 = core.reduce_mod(dl, 5)
out['claims']['delta_mod5'] = {str(k): v for k, v in d5.items()}
# 4. F mod 5 = Fermat (smoothness by good reduction)
F5 = core.reduce_mod(F, 5)
ferm = core.reduce_mod(core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6'), 5)
out['claims']['F_mod5_is_Fermat'] = (F5 == ferm)
# 5. second partials: rank 21 (F is not a join); first partials rank 6
sec = []
for i in range(NV):
    for j in range(i, NV):
        sec.append(pderiv(Fi[i], j))
mons4 = monomials(4)
out['claims']['rank_second_partials_exact'] = fmpq_from_polys(sec, mons4).rank()
out['claims']['rank_first_partials_exact'] = fmpq_from_polys(Fi, monomials(5)).rank()
# 6. diagonal symmetry group
G = symmetry.diagonal_group(F)
out['claims']['diagonal_group_order'] = G.order
out['claims']['diagonal_group_elements'] = [list(map(int, w)) for w in G.elements]
# 7. bigraded pieces of F
big = {}
for e, c in F.items():
    k = e[0] + e[1] + e[2]
    big[k] = big.get(k, 0) + 1
out['claims']['bigraded_term_counts_Xdeg'] = big
# 8. A is Y-odd of degree 6; residue tau-invariant
out['claims']['A_Yodd'] = (core.parity_of(A) == 1)
# 9. Fermat control witness re-verify (exact)
try:
    w = json.load(open(LANE + '/witness-candidate-fermat-92047.json'))
    F0 = core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
    D0 = core.parse(w['D']); N0 = [core.parse(s) for s in w['N']]; Q0 = core.parse(w['Q'])
    ok, info = core.verify_witness(F0, A, 1, D0, N0, Q0)
    out['claims']['fermat_control_witness_verifies'] = bool(ok)
except Exception as ex:
    out['claims']['fermat_control_witness_verifies'] = f'error: {ex}'
out['seconds'] = round(time.time() - t0, 1)
print(json.dumps(out, indent=1))
receipt('verify_exact.json', out)
log({'exp': 'verify_exact', 'family': 'audit', 'result': {k: v for k, v in out['claims'].items() if k not in ('diagonal_group_elements', 'delta_mod5')}, 'seconds': out['seconds']})
