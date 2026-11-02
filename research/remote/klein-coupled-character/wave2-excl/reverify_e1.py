"""(a) adversarial re-verification of the E1-kernel table (wave 1: p = 2^31-1 / dossier: 32003) at a new prime,
with an independent implementation (kx.e1_kernel_block, functional-based), for d = 1..10, both parities, all 7
characters.  Phase 2: EXACT rational kernels for d = 9 (and 10) per block via fmpz_mat nullspace (rigorous dims).
Output: e1_table_p<p>.json, e1_exact_d<d>.json; one log line per (d, parity)."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
import kx
from hodge import core
from hodge.core import shift, pderiv, scale, mul
from flint import fmpz_mat, fmpq_mat, nmod_mat
p = int(sys.argv[1]); DMAX = int(sys.argv[2]); EXACT = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 and sys.argv[3] else []
F, A, q = kx.load()
table = {}
t0 = time.time()
for d in range(1, DMAX + 1):
    for par in (0, 1):
        dims = []; tot = 0
        for ch in range(7):
            r = kx.e1_kernel_block(F, A, d, p, parD=par, k=ch, want_basis=(d >= 9))
            dims.append(r['kernel_modJ']); tot += r['kernel_modJ']
            table[f'{d},{par},{ch}'] = {k: v for k, v in r.items() if k != 'reps'}
            if r.get('reps'):
                table[f'{d},{par},{ch}']['reps'] = [core.to_str(core.symmetric(x, p) if hasattr(core, 'symmetric') else x) for x in r['reps']]
        print(f'd={d} par={par} ker[A] mod J by char: {dims} total {tot}  ({round(time.time()-t0,1)}s)', flush=True)
        kx.log(exp='e1-kernel-table', d=d, parD=par, p=p, kernel_by_char=dims, total=tot,
               dimR_d=[table[f'{d},{par},{c}'].get('dimR_d') for c in range(7)],
               dimR_target=[table[f'{d},{par},{c}'].get('dimR_t') for c in range(7)], seconds=round(time.time() - t0, 1))
json.dump(table, open(f'e1_table_p{p}.json', 'w'))

# ---- phase 2: exact rational kernels -------------------------------------------------------------
def exact_block(d, par, ch):
    """exact {D in S_d^(par,ch): A D in J_{d+6}} / J_d^(par,ch) over Q.  Returns (kernel_total, dimJ_d, kernel_modJ, reps as integer polys)."""
    F7 = scale(F, 7)
    Fi = [pderiv(F7, i) for i in range(6)]
    Dm = kx.wmonomials(d, par, ch)
    if not Dm: return 0, 0, 0, []
    tpar = (par + core.parity_of(A)) % 2; tk = (ch + kx.char_of(A)) % 7
    tgt = kx.wmonomials(d + 6, tpar, tk); rows = {m: i for i, m in enumerate(tgt)}
    cols = [shift(A, m, 1) for m in Dm]
    for i in range(6):
        for m in kx.wmonomials(d + 1, (tpar + kx.PI[i]) % 2, (tk + kx.W[i]) % 7):
            cols.append(shift(Fi[i], m, 1))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * C + j] = int(v)
    M = fmpz_mat(R, C, flat)
    X, nul = M.nullspace()            # integer basis of the right kernel (first nul columns)
    nD = len(Dm)
    # projection onto D-coordinates: rank over Q
    P = fmpq_mat(nD, nul, [X[r, c] for r in range(nD) for c in range(nul)]) if nul else None
    if nul == 0: ker_total = 0; kervecs = []
    else:
        Pr, rk = P.rref()
        ker_total = rk
        # take a basis: rows of rref of P^T
        PT = fmpq_mat(nul, nD, [X[r, c] for c in range(nul) for r in range(nD)])
        PTr, rk2 = PT.rref()
        assert rk2 == ker_total
        kervecs = [[PTr[i, r] for r in range(nD)] for i in range(rk2)]
    # J_d block
    Jcols = []
    for i in range(6):
        for m in kx.wmonomials(d - 5, (par + kx.PI[i]) % 2, (ch + kx.W[i]) % 7):
            Jcols.append(shift(Fi[i], m, 1))
    Drows = {m: i for i, m in enumerate(Dm)}
    if Jcols:
        flatJ = [0] * (nD * len(Jcols))
        for j, col in enumerate(Jcols):
            for e, v in col.items(): flatJ[Drows[e] * len(Jcols) + j] = int(v)
        JM = fmpq_mat(nD, len(Jcols), flatJ); dimJ = JM.rank()
    else:
        JM = None; dimJ = 0
    # reduce kernel vectors modulo J_d: stack [J | ker] and rref
    reps = []
    if ker_total > dimJ:
        Cc = len(Jcols) + ker_total
        flatA = [0] * (nD * Cc)
        if JM is not None:
            for r in range(nD):
                for j in range(len(Jcols)):
                    flatA[r * Cc + j] = JM[r, j]
        for c, vec in enumerate(kervecs):
            for r in range(nD): flatA[r * Cc + len(Jcols) + c] = vec[r]
        MA = fmpq_mat(nD, Cc, flatA)
        MAr, rk = MA.rref()
        piv = []
        for i in range(rk):
            j = 0
            while j < Cc and MAr[i, j] == 0: j += 1
            piv.append(j)
        for j in piv:
            if j >= len(Jcols):
                vec = kervecs[j - len(Jcols)]
                # make integral
                from fractions import Fraction
                import math
                fr = [Fraction(int(v.p), int(v.q)) for v in vec]
                den = 1
                for f in fr: den = den * f.denominator // math.gcd(den, f.denominator)
                reps.append({Dm[r]: int(f * den) for r, f in enumerate(fr) if f != 0})
    return ker_total, dimJ, ker_total - dimJ, reps

for d in EXACT:
    res = {}
    for par in (0, 1):
        dims = []
        for ch in range(7):
            t1 = time.time()
            kt, dj, km, reps = exact_block(d, par, ch)
            dims.append(km)
            res[f'{d},{par},{ch}'] = {'kernel_total': kt, 'dimJ_d': dj, 'kernel_modJ': km, 'reps': [core.to_str(r) for r in reps]}
            print(f'EXACT d={d} par={par} ch={ch}: kernel_total={kt} dimJ_d={dj} modJ={km} ({round(time.time()-t1,1)}s)', flush=True)
        kx.log(exp='e1-kernel-exact-Q', d=d, parD=par, kernel_by_char=dims, total=sum(dims), seconds=round(time.time() - t0, 1),
               note='exact rational nullspace (fmpz_mat) of [A*S_d | J_{d+6} generators] per (parity,character) block, reduced mod J_d')
    json.dump(res, open(f'e1_exact_d{d}.json', 'w'))
print('DONE', round(time.time() - t0, 1))
