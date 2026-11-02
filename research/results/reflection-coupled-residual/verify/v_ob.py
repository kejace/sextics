"""E2-mod-J obstruction ob(D) = [D div N - N(D)] in (S/J)_{2d} for D = K8 (or any Y-even D with A*D in J), computed mod p
with independent code.  Steps: particular E1 solution N0 (rref), rho0, membership test in J_{2d}.
Rigor: ob is independent of the E1 solution (Koszul); ker of the combined (N,c) system is structural (exact Koszul
complex), so 'rho0 not in J mod p' (with the Hilbert check) implies ob_Q(D) != 0, i.e. NO certificate with this D.
'rho0 in J mod p' is only evidence for ob_Q = 0.   Usage: python3 v_ob.py <p> [D-json-file index | D-expr]"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/verify')
from vlib import *
from flint import nmod_mat
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000007
Fs, As, qv = load_rec(); F = modp(parse(Fs), p); A = modp(parse(As), p); cA = char_of(A)
if len(sys.argv) > 2 and not sys.argv[2].endswith('.json'):
    DQ = parse(sys.argv[2]); tag = sys.argv[2][:30]
else:
    DQ = parse(json.load(open('/tmp/hodge/results/reflection-coupled-residual/kernel_d8_X03Y03.json'))[0]); tag = 'K8'
D = modp(DQ, p); d = max(sum(e) for e in D); cD = char_of(DQ); cT = chi_mul(cA, cD)
Fi = [deriv(F, i, p) for i in range(6)]; cFi = [char_of(f) for f in Fi]
t0 = time.time()
# ---- E1 particular solution: sum_i N_i F_i = A*D, N_i in S_{d+1} with char cT*cFi
monsT = monos_chi(d + 6, cT); idxT = {m: j for j, m in enumerate(monsT)}
cols = [(i, m) for i in range(6) for m in monos_chi(d + 1, chi_mul(cT, cFi[i]))]
R, C = len(monsT), len(cols)
M = nmod_mat(R, C + 1, p)
for j, (i, m) in enumerate(cols):
    for e, v in shift(Fi[i], m).items(): M[idxT[e], j] = v
AD = mul(A, D, p)
for e, v in AD.items(): M[idxT[e], C] = v
Mr, rk = M.rref()
N0 = [dict() for _ in range(6)]; consistent = True
for r in range(rk):
    j = 0
    while Mr[r, j] == 0: j += 1
    if j == C: consistent = False; break
    val = int(Mr[r, C])
    if val: i, m = cols[j]; N0[i][m] = val
print(f'[{tag}] E1: {R} eqs x {C} unknowns, rank {rk}, consistent mod p = {consistent} ({time.time()-t0:.0f}s)', flush=True)
if not consistent:
    logjson(VLOG, {'check': 'ob', 'D': tag, 'p': p, 'e1_consistent': False}); sys.exit()
acc = {}
for i in range(6): acc = add(acc, mul(N0[i], Fi[i], p), p)
assert add(acc, scale(AD, p - 1, p), p) == {}, 'E1 verification failed'
# ---- rho0 = D div N0 - N0(D)
divN = {}
for i in range(6): divN = add(divN, deriv(N0[i], i, p), p)
ND = {}
for i in range(6): ND = add(ND, mul(N0[i], deriv(D, i, p), p), p)
rho0 = add(mul(D, divN, p), scale(ND, p - 1, p), p)
assert all(chi(e) == cT and sum(e) == 2 * d for e in rho0)
# ---- membership in J_{2d}^{cT}
mons2 = monos_chi(2 * d, cT); idx2 = {m: j for j, m in enumerate(mons2)}
Jrows = [shift(Fi[i], m) for i in range(6) for m in monos_chi(2 * d - 5, chi_mul(cT, cFi[i]))]
t1 = time.time()
rJ = nmod_from_rows(Jrows, idx2, len(mons2), p).rank()
rJr = nmod_from_rows(Jrows + [rho0], idx2, len(mons2), p).rank()
quot = len(mons2) - rJ; fc = fermat_count(2 * d, cT)
inJ = (rJr == rJ)
print(f'[{tag}] J_{2*d}^{cT}: {len(Jrows)} rows x {len(mons2)} cols, rank {rJ}, dim (S/J) = {quot} (Fermat {fc}, hilbert_ok={quot==fc}); '
      f'rho0 in J mod p: {inJ}  -> ob(D) {"= 0 mod p (evidence only)" if inJ else "!= 0 mod p  => RIGOROUS: no certificate with this D over Q"} ({time.time()-t1:.0f}s)', flush=True)
logjson(VLOG, {'check': 'ob', 'D': tag, 'd': d, 'p': p, 'e1_rank': rk, 'rankJ2d': rJ, 'hilbert_ok': quot == fc, 'ob_zero_modp': inJ})
# save N0 and rho0 for reuse
json.dump({'p': p, 'N0': [{','.join(map(str, e)): v for e, v in Ni.items()} for Ni in N0],
           'rho0': {','.join(map(str, e)): v for e, v in rho0.items()}}, open(f'ob_{tag}_p{p}.json', 'w'))
