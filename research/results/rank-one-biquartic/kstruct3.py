"""Kernel of E1 (target A) inside a structured subspace, modulo J, rational generators via several primes.
Subspace: D = P7(X) + L*M*beta3(X) + c*M*(extra)  (blocks given below); d=7, Y-even."""
import sys, json, time
import lib
from lib import *
from flint import nmod_mat
_argv = sys.argv; sys.argv = [sys.argv[0]]
from kstruct2 import kernel_coords, canon_complement, xmons, inter_dim
sys.argv = _argv
NPR = int(sys.argv[1]) if len(sys.argv) > 1 else 6
SUB = sys.argv[2] if len(sys.argv) > 2 else 'S7+LMS3'
d = 7
blocks = []
if SUB == 'S7+LMS3':
    blocks = [('S7', [{m: Fraction(1)} for m in xmons(7)]), ('LM.S3', [core.mul(core.mul(L, M), {m: Fraction(1)}) for m in xmons(3)])]
elif SUB == 'S7+MS5':
    blocks = [('S7', [{m: Fraction(1)} for m in xmons(7)]), ('M.S5', [core.mul(M, {m: Fraction(1)}) for m in xmons(5)])]
elif SUB == 'X03S4+LMS3':
    blocks = [('X0^3.S4', [core.mul(P('X0**3'), {m: Fraction(1)}) for m in xmons(4)]), ('LM.S3', [core.mul(core.mul(L, M), {m: Fraction(1)}) for m in xmons(3)])]
elif SUB == 'LS5+LMS3':
    blocks = [('L.S5', [core.mul(L, {m: Fraction(1)}) for m in xmons(5)]), ('LM.S3', [core.mul(core.mul(L, M), {m: Fraction(1)}) for m in xmons(3)])]
labels = []; Bpolys = []
for lab, polys in blocks:
    for pp_ in polys: labels.append(lab); Bpolys.append(pp_)
nB = len(Bpolys)
print(f'SUB={SUB}: dim {nB}, blocks {[(l, len(p_)) for l, p_ in blocks]}', flush=True)
t0 = time.time()
kers = {}; primes = core.PRIMES[:NPR]
for p in primes:
    E1, rkE, Jr, rkJ = kernel_coords(A, Bpolys, p)
    ker = canon_complement(E1, rkE, Jr, rkJ, nB, p)
    print(f'  p={p}: dim(E1 cap B)={rkE}, dim(J cap B)={rkJ}, kernel mod J = {len(ker)} ({time.time()-t0:.0f}s)', flush=True)
    if p == primes[0]:
        for lab, _ in blocks:
            keep = [j for j, l in enumerate(labels) if l == lab]
            print(f'    block {lab} alone: kernel mod J = {inter_dim(E1, rkE, keep, nB, p) - inter_dim(Jr, rkJ, keep, nB, p)}')
    kers[p] = ker
    if len(ker) == 0: break
    # try reconstruction with primes so far
    gens = []; ok_all = True
    for idx in range(len(ker)):
        poly = {}; coords = []; ok = True
        for j in range(nB):
            a, mm = kers[primes[0]][idx][j], primes[0]
            for pp in primes[1:primes.index(p) + 1]:
                a, mm = core.crt_pair(a, mm, kers[pp][idx][j], pp)
            r = core.ratrecon(a, mm)
            if r is None: ok = False; break
            if r: coords.append((labels[j], j, r)); poly = core.add(poly, core.scale(Bpolys[j], r))
        ok_all &= ok
        if ok: gens.append((poly, coords))
    if ok_all and len(gens) == len(ker):
        print(f'  reconstruction succeeded with {primes.index(p) + 1} primes', flush=True)
        out = []
        for poly, coords in gens:
            # verify E1 exactly mod another prime via e1_solvable
            okE, info = fixed_d.e1_solvable(F, A, q, poly, core.PRIMES[-1])
            hts = max(max(abs(r.numerator), r.denominator) for _, _, r in coords)
            print('   generator: nterms', len(coords), 'max height', hts, 'E1 check (last prime):', okE, flush=True)
            print('     coords:', [(lab, core.to_str({next(iter(Bpolys[j].keys())): 1}) if lab == 'S7' else lab + '#' + str(j), str(r)) for lab, j, r in coords], flush=True)
            out.append(core.to_str(poly))
        json.dump({'target': 'X0**3*Y0**3', 'd': 7, 'sub': SUB, 'generators': out}, open(f'{WORK}/kstruct3_{SUB}.json', 'w'))
        break
print('done', time.time() - t0)
