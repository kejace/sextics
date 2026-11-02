"""E2: smoothness re-verification. (i) Singular at p=1000003: jacob(2F) zero-dimensional, vdim.
(ii) own sector-wise check at p=2147483629: dim (S/J)_k per G-sector for k=0..25; R_25 = 0 and the
Hilbert function (1+t+..+t^4)^6 (=> F_i regular sequence => smooth; smooth mod p => smooth over Q)."""
import sys, json, subprocess, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import xlib
from xlib import *
F2, A, q = load_cycle_q2()
t0 = time.time()
sing = f"""ring r = 1000003,(X0,X1,X2,Y0,Y1,Y2),dp;
poly F2 = {to_str(F2).replace('**','^')};
ideal J = jacob(F2);
ideal G = std(J);
print("dim=" + string(dim(G)));
print("vdim=" + string(vdim(G)));
quit;
"""
open('/tmp/hodge/results/cycle-q2/wave2-excl/receipts/smooth_1000003.sing', 'w').write(sing)
out = subprocess.run(['Singular', '-q', '/tmp/hodge/results/cycle-q2/wave2-excl/receipts/smooth_1000003.sing'], capture_output=True, text=True).stdout
print(out)
rec = dict(exp='E2-smoothness-singular', p=1000003, singular_output=out.strip().replace('\n', ' | '), seconds=round(time.time() - t0, 1))
xlib.log(rec)
# (ii) own check
t0 = time.time()
JS = JacobianSectors(F2, P2)
hilb = [0] * 26
sector_dims = {}
for k in range(0, 26):
    tot = 0
    for key in all_sector_keys(JS.gens, k):
        dmk = JS.ring_dim(k, key)
        tot += dmk
        sector_dims[(k, key)] = dmk
    hilb[k] = tot
    print('k', k, 'dim R_k =', tot, flush=True)
# expected (1+t+t^2+t^3+t^4)^6
import numpy as np
poly = np.array([1])
for _ in range(6):
    poly = np.convolve(poly, np.array([1, 1, 1, 1, 1]))
exp = [int(poly[k]) if k < len(poly) else 0 for k in range(26)]
# target sector of A (residue character = chi of A*Omega: key(A) + key(volume) ) : compute dims of R_6, R_12, R_18 in the sector of A shifted
kA = key_of_poly(A, JS.gens)
kvol = key_of_poly({(1, 1, 1, 1, 1, 1): 1}, JS.gens)
kF = JS.kF
# Griffiths: numerators of degree 6k-6 at pole order k: class of B Omega/F^k has character key(B)+kvol-k*key(F) (mod scalars this is what matters)
# the sector of alpha: numerator A at k=3. Hodge numbers of that sector: numerators B of degree 6k-6 with key(B)+kvol-k*kF == kA+kvol-3*kF
res = {}
target = key_add(kA, kvol, key_mul(key_neg(kF), 3))
for k in range(1, 6):
    deg = 6 * k - 6
    need = key_add(target, key_mul(kF, k), key_neg(kvol))
    res[f'h^{5-k},{k-1}'] = JS.ring_dim(deg, need) if deg <= 24 else 0
rec = dict(exp='E2-smoothness-own', engine=ENGINE_TAG, p=P2, hilbert=hilb, expected=exp, hilbert_matches=(hilb == exp), R25_zero=(hilb[25] == 0),
           target_sector_hodge=res, seconds=round(time.time() - t0, 1),
           note='R_25 = 0 mod p => J_25 = S_25 over Q => F smooth over C; Hilbert function matches the regular-sequence series')
print(json.dumps(rec))
xlib.log(rec)
