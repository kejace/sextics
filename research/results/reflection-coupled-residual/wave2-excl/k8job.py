"""K8 completion job: (1) rigorous fixed-D test of K8 (restricted + full), (2) K8 not in J (two primes),
(3) exact rational E1 lift (character-restricted E1 system, CRT over primes, exact verification),
(4) wave-1 e2lin.py (E2-mod-J linearisation over the E1 family) as a subprocess."""
import sys, json, time, subprocess
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
import rigorfix
from hodge import core, jacobian
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv
from flint import nmod_mat
from fractions import Fraction

K8s = json.load(open('/tmp/hodge/results/reflection-coupled-residual/kernel_d8_X03Y03.json'))[0]
K8 = core.parse(K8s)
d = core.degree(K8)
assert d == 8 and all(L.chi(e) == (0, 0) for e in K8)
t0 = time.time()
# (1) rigorous fixed-D test
rec = rigorfix.run('K8-e1kernel-d8', K8s, full=True)
# (2) K8 not in J (reduction mod p of K8 lies outside J_p => K8 outside J_Q)
for p in L.PRIMES[2:4]:
    ok, info = jacobian.in_jacobian(L.F, K8, p)
    L.receipt('k8job', exp='K8-in-J', prime=p, in_J_modp=ok, info=info,
              note='K8 not in J mod p implies K8 not in J over Q (denominators of K8 prime to p)')
# (3) exact E1 lift: N_i in S_9 of character chi(A)+chi(z_i), sum N_i F_i = A K8
Fq = L.F; Aq = L.A
AD = mul(Aq, K8)
cols = []
for i in range(NV):
    ci = L.cadd(L.CHI_A, L.CHI_Z[i])
    for m in L.mons_chi(9, ci):
        cols.append((i, m))
rows = L.mons_chi(14, L.CHI_A); ridx = {m: k for k, m in enumerate(rows)}
assert all(e in ridx for e in AD)
acc = None; mod = 1; pivots0 = None; lifted = None
for p in L.PRIMES[2:12]:
    t1 = time.time()
    Fp = core.reduce_mod(Fq, p); Fi = [pderiv(Fp, i, p) for i in range(NV)]
    ADp = core.reduce_mod(AD, p)
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, (i, m) in enumerate(cols):
        for e, v in shift(Fi[i], m, 1, p).items():
            flat[ridx[e] * (C + 1) + j] = v
    for e, v in ADp.items():
        flat[ridx[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    Mr, rk = M.rref()
    piv = []
    for r in range(rk):
        j = 0
        while int(Mr[r, j]) == 0: j += 1
        piv.append(j)
    consistent = not (piv and piv[-1] == C)
    x = [0] * C
    if consistent:
        for r, j in enumerate(piv):
            x[j] = int(Mr[r, C])
    print(f'[E1 lift p={p}] rows={R} cols={C} rank={rk} consistent={consistent} ({time.time()-t1:.0f}s)', flush=True)
    L.receipt('k8job', exp='K8-E1-modp', prime=p, rows=R, cols=C, rank=rk, consistent=consistent)
    if not consistent:
        break
    if pivots0 is None:
        pivots0 = piv; acc = x; mod = p
    elif piv != pivots0:
        L.receipt('k8job', exp='K8-E1-lift', note='pivot structure changed; restarting CRT', prime=p)
        pivots0 = piv; acc = x; mod = p
    else:
        acc = [core.crt_pair(a, mod, b, p)[0] for a, b in zip(acc, x)]; mod *= p
    rat = [core.ratrecon(a, mod) for a in acc]
    if all(v is not None for v in rat):
        N = [{} for _ in range(NV)]
        for (i, m), v in zip(cols, rat):
            if v: N[i][m] = v
        s = {}
        Fiq = [pderiv(Fq, i) for i in range(NV)]
        for i in range(NV):
            s = add(s, mul(N[i], Fiq[i]))
        exact = (sub(s, AD) == {})
        L.receipt('k8job', exp='K8-E1-lift', primes_used=len([q for q in L.PRIMES[2:12] if q <= p]) , modulus_bits=mod.bit_length(),
                  exact_identity=exact, n_terms=[len(n) for n in N],
                  max_den_bits=max(Fraction(v).denominator.bit_length() for v in rat))
        if exact:
            lifted = N
            json.dump({'D': K8s, 'N': [core.to_str(n) for n in N]}, open(f'{L.OUT}/K8_E1_exact.json', 'w'))
            break
print('E1 lift done, exact =', lifted is not None, f'({time.time()-t0:.0f}s)', flush=True)
# (4) wave-1 e2lin.py (never completed in wave 1)
t2 = time.time()
r = subprocess.run(['python3', '/tmp/hodge/results/reflection-coupled-residual/e2lin.py',
                    '/tmp/hodge/results/reflection-coupled-residual/kernel_d8_X03Y03.json', '0'],
                   cwd='/tmp/hodge/results/reflection-coupled-residual', capture_output=True, text=True)
open(f'{L.OUT}/e2lin_d8.out', 'w').write(r.stdout + '\n--- stderr ---\n' + r.stderr[-3000:])
try:
    e2 = json.load(open('/tmp/hodge/results/reflection-coupled-residual/e2lin_d8.json'))
except Exception as ex:
    e2 = {'error': repr(ex)}
L.receipt('k8job', exp='e2lin-d8-family', seconds=round(time.time() - t2, 1), result=e2)
print('K8JOB DONE', flush=True)
