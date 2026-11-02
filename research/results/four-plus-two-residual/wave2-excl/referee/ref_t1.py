import sys, time; sys.path.insert(0, '/tmp/hodge/work/referee42')
import reflib as R
from fractions import Fraction
LOG = '/tmp/hodge/work/referee42/out/ref_t1.jsonl'
PNEW = 1000000009   # prime not used by any lens so far (they used 2^31-1, 10^9+7, 998244353, 10^9+289, 61-bit primes)
F, A = R.load_instance()
H = {e: c for e, c in F.items() if e[2] == 0 and e[5] == 0}; B = {e: c for e, c in F.items() if e[2] + e[5] == 6}
# (1) Hilbert functions of R(H), R(B) mod PNEW vs complete-intersection series; dim_p >= dim_Q >= CI  => exact
def ci(n):
    poly = [1]
    for _ in range(n):
        new = [0] * (len(poly) + 4)
        for i, c in enumerate(poly):
            for j in range(5): new[i + j] += c
        poly = new
    return poly
for name, P, idx in (('B', B, R.BIN), ('H', H, R.CORE)):
    t0 = time.time(); hf = []
    for k in range(0, 4 * len(idx) + 2):
        Jc = R.jac_cols(P, idx, k)
        hf.append(len(R.monomials(k, None, idx)) - (R.rank_modp(Jc, PNEW) if Jc else 0))
    exp = ci(len(idx)) + [0]
    R.log(LOG, {'exp': 'hilbert-modp', 'poly': name, 'p': PNEW, 'hilbert': hf, 'ci': exp, 'equal_CI': hf == exp, 'seconds': round(time.time() - t0, 1),
                'note': 'dim_p >= dim_Q >= CI coefficientwise; equality mod p proves the CI Hilbert function over Q => partials are a regular sequence over Q and mod p'})
# (2) core kernels ker(a. : R(H)_k -> R(H)_{k+6}) by parity: mod PNEW for k=0..12, exact over Q for k=5,6,7,8
for k in range(0, 13):
    for par in (0, 1):
        t0 = time.time()
        r = R.mult_kernel_dims(H, R.CORE, A, k, par, p=PNEW)
        r.update({'exp': 'core-kernel-modp', 'p': PNEW, 'seconds': round(time.time() - t0, 1)}); R.log(LOG, r)
for k in (5, 6, 7, 8):
    for par in (0, 1):
        t0 = time.time()
        r = R.mult_kernel_dims(H, R.CORE, A, k, par, p=None)
        r.update({'exp': 'core-kernel-exact-Q', 'seconds': round(time.time() - t0, 1)}); R.log(LOG, r)
# (3) direct 6-variable E1 kernels {D : A D in J(F)} for d = 1..6, both parities, mod PNEW
for d in range(1, 7):
    for par in (0, 1):
        t0 = time.time()
        r = R.mult_kernel_dims(F, tuple(range(6)), A, d, par, p=PNEW)
        r.update({'exp': '6var-kernel-modp', 'p': PNEW, 'seconds': round(time.time() - t0, 1)}); R.log(LOG, r)
