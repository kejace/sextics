"""Structure of the degree-7 E1 kernel: (a) E1 membership of structured candidates, (b) divisibility of kernel directions by linear/quadric forms mod J, (c) fixed-D consistency (mod p) of kernel basis elements and random combinations."""
import sys, json, random, time
import lib
from lib import *
from flint import nmod_mat
from kernelQ import canon_kernel
p = core.PRIMES[0]
names = ['X0','X1','X2','Y0','Y1','Y2']
lin = [P(n) for n in names]
Fi = [core.pderiv(F, i) for i in range(6)]
# (a) structured candidates
cands = []
for nm, K in [('f', f), ('g', g), ('P', Pc), ('RX', RX), ('RY', RY), ('f+g', core.add(f, g)), ('LM^2', core.mul(core.mul(L, M), M))]:
    for j, l in enumerate(lin):
        cands.append((core.mul(l, K), f'{names[j]}*{nm}'))
for i in range(6):
    for qn, Qf in [('L', L), ('M', M), ('X0^2', P('X0**2')), ('Y0^2', P('Y0**2')), ('X0*Y0', P('X0*Y0'))]:
        cands.append((core.mul(Qf, Fi[i]), f'{qn}*F_{i}'))
print('(a) E1 membership of structured degree-7 candidates:')
for D, nm in cands:
    if core.parity_of(D) is None: continue
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    print(f'   {nm:14s} parity={core.parity_of(D)} E1={ok}')
    log({'D': nm, 'd': 7, 'family': 'e1-structured-d7', 'e1': bool(ok), 'consistent': None if ok else False,
         'rank': info['rank'], 'nullity': info['cols'] - info['rank'], 'seconds': 0, 'note': 'E1 membership only', 'engine': 2})
# (b) kernel mod p and divisibility
Dmons, ker, rkJ, rkE = canon_kernel(7, 0, p)
print('(b) kernel dim', len(ker), 'J', rkJ, 'E1', rkE)
idxD = {m: k for k, m in enumerate(Dmons)}; nD = len(Dmons)
Fp = core.reduce_mod(F, p)
Jcols = []
for i in range(6):
    pi = 0 if i < 3 else 1
    for m in core.monomials(2, (0 - pi) % 2):
        Jcols.append(core.shift(core.pderiv(Fp, i, p), m, 1, p))
def dim_span(polys):
    flat = [0] * (len(polys) * nD)
    for j, c in enumerate(polys):
        for e, v in c.items(): flat[j * nD + idxD[e]] = v
    return nmod_mat(len(polys), nD, flat, p).rank()
kerpolys = [{Dmons[j]: v for j, v in enumerate(vec) if v} for vec in ker]
base = dim_span(Jcols)
for nm, l in [(n, P(n)) for n in names] + [('L', L), ('M', M), ('f', f), ('g', g), ('P', Pc), ('RX', RX), ('RY', RY)]:
    dl = core.degree(l); lp = core.reduce_mod(l, p)
    par_need = (0 - core.parity_of(l)) % 2
    mult = [core.shift(lp, m, 1, p) for m in core.monomials(7 - dl, par_need)]
    V = Jcols + mult
    dV = dim_span(V)
    dVK = dim_span(V + kerpolys)
    print(f'   multiplier {nm:3s}: dim(J+l*S)={dV}, dim(J+l*S+K)={dVK}, kernel directions divisible by {nm} mod J: {3 - (dVK - dV)}')
# (c) fixed-D consistency mod p for kernel reps (integer lifts) and random combos
print('(c) fixed-D consistency (mod p) of kernel representatives:')
random.seed(3)
tests = [(kerpolys[i], f'ker7 basis {i}') for i in range(len(kerpolys))]
for k in range(3):
    cs = [random.randint(1, p - 1) for _ in kerpolys]
    D = {}
    for c, kp in zip(cs, kerpolys): D = core.add(D, core.scale(kp, c, p), p)
    tests.append((D, f'ker7 random combo {k}'))
for D, nm in tests:
    Dint = {e: int(v) for e, v in D.items()}
    t0 = time.time()
    S = fixed_d.FixedDSystem(F, A, q, Dint, p)
    r = fixed_d.solve_compressed(S)
    print(f'   {nm}: consistent={r["consistent"]} rank={r["rank"]} nullity={r["nullity"]} rows={S.nrows} cols={S.ncols} ({time.time()-t0:.1f}s)', flush=True)
    log({'D': nm + ' (mod-p representative, canonical rref form)', 'd': 7, 'family': 'e1-kernel-d7', 'e1': True,
         'consistent': bool(r['consistent']), 'rank': int(r['rank']), 'nullity': int(r['nullity']), 'rows': S.nrows, 'cols': S.ncols,
         'seconds': round(time.time() - t0, 1), 'note': 'mod-p test only; representative modulo J is the canonical one', 'engine': 2})
