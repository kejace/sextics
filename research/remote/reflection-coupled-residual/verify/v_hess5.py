"""Independent check of wave2-excl R1: for D = sum c_i F_i (c constants, the whole E1 space at d=5),
ob(c) = -[A sum c_i c_j F_ij] in (S/J)_10.  Verify F_{X2Y2} = 0 and that the other 20 classes [A F_ij] are linearly
independent mod J_10 at a new prime (independence mod p implies independence over Q: rank_p <= rank_Q)."""
import sys, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/verify')
from vlib import *
from flint import nmod_mat
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000007
Fs, As, qv = load_rec(); FQ = parse(Fs); F = modp(FQ, p); A = modp(parse(As), p); cA = char_of(A)
Fi = [deriv(F, i, p) for i in range(6)]; cFi = [char_of(f) for f in Fi]
Fij = {(i, j): deriv(Fi[i], j, p) for i in range(6) for j in range(i, 6)}
zero = [k for k, v in Fij.items() if not v]
print('vanishing Hessian entries:', zero)
t0 = time.time()
# all 21 classes A*F_ij in (S/J)_10, all characters at once
mons = monos(10); idx = {m: j for j, m in enumerate(mons)}
Jrows = [shift(Fi[i], m) for i in range(6) for m in monos(5)]
rJ = nmod_from_rows(Jrows, idx, len(mons), p).rank()
print(f'rank J_10 = {rJ}, dim (S/J)_10 = {len(mons)-rJ} (Fermat {fermat_count(10)})')
cls = {k: mul(A, v, p) for k, v in Fij.items() if v}
keys = sorted(cls)
rAll = nmod_from_rows(Jrows + [cls[k] for k in keys], idx, len(mons), p).rank()
print(f'{len(keys)} nonzero classes [A F_ij]: rank mod J = {rAll - rJ} -> linearly independent: {rAll - rJ == len(keys)} ({time.time()-t0:.0f}s)')
# by character blocks
blocks = {}
for k in keys: blocks.setdefault(chi_mul(cA, chi_mul(cFi[k[0]], cFi[k[1]])), []).append(k)
print('block sizes by character of A F_ij:', {str(c): len(v) for c, v in sorted(blocks.items())})
# explicit conclusion: c_i c_j = 0 for all 20 pairs => c_i^2 = 0 => c = 0
ok = (rAll - rJ == len(keys)) and all((i, i) in keys for i in range(6))
print('R1 conclusion (no degree-5 denominator) follows:', ok)
logjson(VLOG, {'check': 'hess5', 'p': p, 'zero_entries': [list(k) for k in zero], 'n_classes': len(keys), 'rank_modJ': rAll - rJ, 'independent': rAll - rJ == len(keys), 'hilbert_ok': len(mons) - rJ == fermat_count(10)})
