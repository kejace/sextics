"""Generic degree-d kernel structure + mod-p fixed-D tests. Usage: python3 e1struct.py d [ntests]"""
import sys, json, random, time
import lib
from lib import *
from flint import nmod_mat
from kernelQ import canon_kernel
d = int(sys.argv[1]); ntests = int(sys.argv[2]) if len(sys.argv) > 2 else 4
p = core.PRIMES[0]
names = ['X0','X1','X2','Y0','Y1','Y2']
Dmons, ker, rkJ, rkE = canon_kernel(d, 0, p)
print(f'd={d}: kernel dim', len(ker), 'J', rkJ, 'E1', rkE, flush=True)
idxD = {m: k for k, m in enumerate(Dmons)}; nD = len(Dmons)
Fp = core.reduce_mod(F, p)
Jcols = []
for i in range(6):
    pi = 0 if i < 3 else 1
    for m in core.monomials(d - 5, (0 - pi) % 2):
        Jcols.append(core.shift(core.pderiv(Fp, i, p), m, 1, p))
def dim_span(polys):
    flat = [0] * (len(polys) * nD)
    for j, c in enumerate(polys):
        for e, v in c.items(): flat[j * nD + idxD[e]] = v
    return nmod_mat(len(polys), nD, flat, p).rank()
kerpolys = [{Dmons[j]: v for j, v in enumerate(vec) if v} for vec in ker]
nk = len(kerpolys)
mults = [(n, P(n)) for n in names] + [('L', L), ('M', M), ('f', f), ('g', g), ('P', Pc), ('RX', RX), ('RY', RY)]
mults += [(f'F_{i}', core.pderiv(F, i)) for i in range(6)]
for nm, l in mults:
    dl = core.degree(l)
    if dl > d: continue
    lp = core.reduce_mod(l, p)
    par_need = (0 - core.parity_of(l)) % 2
    mult = [core.shift(lp, m, 1, p) for m in core.monomials(d - dl, par_need)]
    V = Jcols + mult
    dV = dim_span(V); dVK = dim_span(V + kerpolys)
    print(f'   multiplier {nm:4s}: kernel directions divisible by {nm} mod J: {nk - (dVK - dV)}', flush=True)
# combined: divisible by any X_i (ideal (X0,X1,X2)) or any Y_j
for nm, gens in [('(X0,X1,X2)', names[:3]), ('(Y0,Y1,Y2)', names[3:])]:
    mult = []
    for n in gens:
        lp = core.reduce_mod(P(n), p)
        mult += [core.shift(lp, m, 1, p) for m in core.monomials(d - 1, (0 - core.parity_of(P(n))) % 2)]
    V = Jcols + mult
    dV = dim_span(V); dVK = dim_span(V + kerpolys)
    print(f'   ideal {nm}: kernel directions inside ideal mod J: {nk - (dVK - dV)}', flush=True)
with open(f'{WORK}/kernel_modp_d{d}.json', 'w') as fh:
    json.dump({'p': p, 'basis': [core.to_str(k) for k in kerpolys]}, fh)
print(f'(c) fixed-D consistency (mod p) of degree-{d} kernel representatives:', flush=True)
random.seed(5)
tests = [(kerpolys[i], f'ker{d} basis {i}') for i in range(nk)]
for k in range(ntests):
    cs = [random.randint(1, p - 1) for _ in kerpolys]
    D = {}
    for c, kp in zip(cs, kerpolys): D = core.add(D, core.scale(kp, c, p), p)
    tests.append((D, f'ker{d} random combo {k}'))
for D, nm in tests:
    Dint = {e: int(v) for e, v in D.items()}
    t0 = time.time()
    S = fixed_d.FixedDSystem(F, A, q, Dint, p)
    r = fixed_d.solve_compressed(S)
    print(f'   {nm}: consistent={r["consistent"]} rank={r["rank"]} nullity={r["nullity"]} rows={S.nrows} cols={S.ncols} ({time.time()-t0:.1f}s)', flush=True)
    log({'D': nm + ' (mod-p representative, canonical rref form)', 'd': d, 'family': f'e1-kernel-d{d}', 'e1': True,
         'consistent': bool(r['consistent']), 'rank': int(r['rank']), 'nullity': int(r['nullity']), 'rows': S.nrows, 'cols': S.ncols,
         'seconds': round(time.time() - t0, 1), 'note': 'mod-p test only; representative modulo J is the canonical one', 'engine': 2})
