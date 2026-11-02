"""Degree-d cones with a line vertex (polynomials in 4 variables) that are E1-admissible: reconstruct a rational basis from k primes
and run the rigorous fixed-D test on nmax basis elements (+ a few random rational combinations).
Usage: python3 conescreen.py <d> <absent vars e.g. Y0,Y1> <chi ab|all> <nprimes> <nmax> [nrand]"""
import sys, json, time, random
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
import rigorfix
from hodge import core
from hodge.core import monomials, shift
from fractions import Fraction
d = int(sys.argv[1]); absent = sys.argv[2].split(','); chiarg = sys.argv[3]; k = int(sys.argv[4]); nmax = int(sys.argv[5])
nrand = int(sys.argv[6]) if len(sys.argv) > 6 else 3
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']; ab = [V.index(v) for v in absent]
chis = L.CHARS if chiarg == 'all' else [tuple(int(x) for x in chiarg)]
tag = f'cone-d{d}-no{"".join(absent)}'
for c in chis:
    Um = [m for m in L.mons_chi(d, c) if all(m[i] == 0 for i in ab)]
    if not Um: continue
    per_prime = []; dims = []
    for p in L.PRIMES[2:2 + k]:
        Ap = core.reduce_mod(L.A, p)
        RT = L.CharReducer(d + 6, L.cadd(c, L.CHI_A), p, verbose=False)
        K = RT.cls_many([shift(Ap, m, 1, p) for m in Um])
        X, nul = K.transpose().nullspace()
        kers = [{Um[j]: int(X[j, kk]) for j in range(len(Um)) if int(X[j, kk])} for kk in range(nul)]
        # canonical form: rref of the kernel vectors over Um
        from flint import nmod_mat
        M = nmod_mat(len(kers), len(Um), [v.get(m, 0) for v in kers for m in Um], p) if kers else None
        reps = []
        if M is not None:
            Mr, rk = M.rref()
            for i in range(rk):
                reps.append({Um[j]: int(Mr[i, j]) for j in range(len(Um)) if int(Mr[i, j])})
        per_prime.append((p, reps)); dims.append(len(reps))
        print(f'chi={c} p={p}: E1-admissible cones dim {len(reps)}', flush=True)
    if not dims or dims[0] == 0: continue
    assert len(set(dims)) == 1, dims
    rat = []
    for idx in range(dims[0]):
        ok, Pq = L.ratrecon_polys([(p, reps[idx]) for p, reps in per_prime])
        rat.append({e: v for e, v in Pq.items()} if ok else None)
    nok = sum(1 for r in rat if r)
    L.receipt('conescreen', exp='cone-family', d=d, chi=c, absent=absent, dim_admissible=dims[0], reconstructed=nok, primes=L.PRIMES[2:2 + k])
    good = [r for r in rat if r]
    random.seed(d)
    tests = good[:nmax]
    for t in range(nrand):
        D = {}
        for r in good:
            cf = random.randint(-2, 2)
            if cf: D = core.add(D, core.scale(r, Fraction(cf)))
        if D: tests.append(D)
    for i, D in enumerate(tests):
        try:
            rigorfix.run(f'{tag}-chi{c[0]}{c[1]}-{i}', core.to_str(D), full=False)
        except Exception as ex:
            L.receipt('rigorfix', exp='rigorfix', family=f'{tag}-chi{c[0]}{c[1]}-{i}', error=repr(ex))
print('CONESCREEN DONE', flush=True)
