"""K8 properties at a new prime: (i) A*K8 in J_14 mod p (third prime); (ii) non-proportionality mod J_8 to
H_u^2, Y2^4 H_u, Y2^8, X0^4 Y0^4, X0^8 and the (0,0)-sector J_8 dimension; (iii) character."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/verify')
from vlib import *
from flint import nmod_mat
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000007
Fs, As, qv = load_rec(); FQ = parse(Fs); F = modp(FQ, p); A = modp(parse(As), p)
K8Q = parse(json.load(open('/tmp/hodge/results/reflection-coupled-residual/kernel_d8_X03Y03.json'))[0]); K8 = modp(K8Q, p)
Fi = [deriv(F, i, p) for i in range(6)]; cFi = [char_of(f) for f in Fi]
print('char K8', char_of(K8Q), 'terms', len(K8Q))
# (i) membership of A*K8 in J_14^{(0,1)}
t0 = time.time()
cT = (0, 1); monsT = monos_chi(14, cT); idxT = {m: j for j, m in enumerate(monsT)}
Jrows = [shift(Fi[i], m) for i in range(6) for m in monos_chi(9, chi_mul(cT, cFi[i]))]
AK = mul(A, K8, p)
rJ = nmod_from_rows(Jrows, idxT, len(monsT), p).rank()
rJA = nmod_from_rows(Jrows + [AK], idxT, len(monsT), p).rank()
print(f'(i) p={p}: rank J_14^(0,1) = {rJ}, dim S_14^(0,1) = {len(monsT)}, dim (S/J) = {len(monsT)-rJ} (Fermat {fermat_count(14,cT)}); rank with A*K8 appended = {rJA} -> A*K8 in J mod p: {rJA == rJ} ({time.time()-t0:.0f}s)')
# (ii) proportionality tests in (S/J)_8^{(0,0)}
mons8 = monos_chi(8, (0, 0)); idx8 = {m: j for j, m in enumerate(mons8)}
J8 = [shift(Fi[i], m) for i in range(6) for m in monos_chi(3, cFi[i])]
rJ8 = nmod_from_rows(J8, idx8, len(mons8), p).rank()
print(f'(ii) dim S_8^(0,0) = {len(mons8)}, rank J_8^(0,0) = {rJ8} (expected 112), dim (S/J)_8^(0,0) = {len(mons8)-rJ8} (Fermat {fermat_count(8,(0,0))})')
Hu = parse('3*Y2**4 + 5*X0*X1*(Y0**2 + 2*Y1**2)/7')  # dF/du with u = Y2^2
cands = {'H_u^2': mul(Hu, Hu), 'Y2^4*H_u': mul(parse('Y2**4'), Hu), 'Y2^8': parse('Y2**8'), 'X0^4*Y0^4': parse('X0**4*Y0**4'), 'X0^8': parse('X0**8'),
         'K8 itself (sanity: must be dependent)': K8Q}
rK = nmod_from_rows(J8 + [K8], idx8, len(mons8), p).rank()
print('rank J8+K8 =', rK, '(K8 not in J_8 mod p:', rK == rJ8 + 1, ')')
for name, C in cands.items():
    Cp = modp(C, p)
    rC = nmod_from_rows(J8 + [Cp], idx8, len(mons8), p).rank()
    rKC = nmod_from_rows(J8 + [K8, Cp], idx8, len(mons8), p).rank()
    print(f'   {name:40s}: in J_8? {rC == rJ8}; proportional to K8 mod J_8? {rKC == rJ8 + 1}')
# check that H_u is what the attack note says: F = H(x, u), dF/dY2 = 2 Y2 H_u
dF5 = deriv(FQ, 5)
Y2Hu2 = scale(mul(parse('Y2'), Hu), Fraction(2))
print('dF/dY2 == 2*Y2*H_u exactly:', add(dF5, scale(Y2Hu2, Fraction(-1))) == {})
logjson(VLOG, {'check': 'k8props', 'p': p, 'AK8_in_J_modp': rJA == rJ, 'rankJ14': rJ, 'rankJ8_00': rJ8})
