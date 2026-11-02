"""Step 9: (a) Fermat-deformation denominators (f6 = X0^6+X1^6+X2^6 serves A on the Fermat sextic; test f6, g6, f6+g6, Qx times
monomials/blocks on off-paired); (b) the cubic-cover discriminant Delta_t = 4 Qx^3 + 27 f6^2 g6 (degree 18, G-invariant)
and its irreducible factors.  All candidates: E1 (reduced) first, then the full reduced fixed-D system (dense for small,
compressed-light for large)."""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/off-paired')
import oplib
from hodge import core, fixed_d
from hodge.core import mul, add, sub, scale, power, monomials, shift, pderiv
import numpy as np
from flint import nmod_mat, fmpq_mpoly_ctx, fmpq, Ordering

F, A, q = oplib.load()
p = oplib.P0
V = core.VARS
def mstr(e): return '*'.join(f'{V[i]}**{k}' if k > 1 else V[i] for i, k in enumerate(e) if k) or '1'

ctx = fmpq_mpoly_ctx.get(tuple(V), Ordering.lex)
def topoly(P): return ctx.from_dict({e: fmpq(c.numerator, c.denominator) for e, c in P.items()})
def frompoly(P):
    out = {}
    for e, c in P.to_dict().items():
        out[tuple(int(x) for x in e)] = core.Fraction(int(c.p), int(c.q)) if hasattr(c, 'p') else core.Fraction(str(c))
    return out
def factors(P):
    c, fs = topoly(P).factor()
    return [(frompoly(f), m) for f, m in fs]

def compressed_light(S, k_extra=40, seed=1, verbose=True):
    """Like fixed_d.solve_compressed but builds the nmod_mat from a generator (no giant Python list)."""
    t0 = time.time()
    M = fixed_d._sparse_matrix(S)
    C = S.ncols
    k = min(S.nrows, C + 1 + k_extra)
    Cm = fixed_d._compress(M, S.p, k, seed)
    del M
    t1 = time.time()
    def gen():
        for i in range(k):
            row = Cm[i].tolist()
            for v in row: yield v
    Mn = nmod_mat(k, C + 1, gen(), S.p)
    del Cm
    t2 = time.time()
    Mr, rank = Mn.rref()
    t3 = time.time()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        pivots.append(j)
    consistent = not (pivots and pivots[-1] == C)
    rank0 = rank - (0 if consistent else 1)
    if verbose: print(f'[compressed-light] k={k} C={C} consistent={consistent} rank={rank0} nullity={C-rank0} (compress {t1-t0:.0f}s, build {t2-t1:.0f}s, rref {t3-t2:.0f}s)', file=sys.stderr, flush=True)
    return {'consistent': consistent, 'rank': rank0, 'nullity': C - rank0, 'pivots': pivots, 'rows': S.nrows, 'cols': C}

def test(name, fam, D, do_solve=True):
    d = core.degree(D)
    c = oplib.char_of(D)
    if c is None:
        print(f'{name}: not an eigenvector, skip', file=sys.stderr, flush=True); return
    t0 = time.time()
    ok, info = oplib.e1_reduced(F, A, q, D, p)
    if not ok:
        dt = time.time() - t0
        oplib.log(core.to_str(D) if len(D) <= 40 else name, d, fam, False, None, info['rank'], 0, dt, f'{name}: E1 fails (char {c}; rank {info["rank"]}/{info["rank_aug"]}, rows {info["rows"]} cols {info["cols"]})')
        print(f'{name:45s} d={d} char={c} E1=False rank {info["rank"]}/{info["rank_aug"]} {dt:.1f}s', file=sys.stderr, flush=True); return
    if not do_solve:
        print(f'{name:45s} d={d} char={c} E1=True (solve skipped)', file=sys.stderr, flush=True); return
    S = oplib.ReducedSystem(F, A, q, D, p, verbose=True)
    if S.ncols <= 6000:
        r = S.solve()
    else:
        r = compressed_light(S)
    dt = time.time() - t0
    oplib.log(core.to_str(D) if len(D) <= 40 else name, d, fam, True, r['consistent'], r['rank'], r['nullity'], dt, f'{name}: char={c} rows={S.nrows} cols={S.ncols}')
    print(f'{name:45s} d={d} char={c} E1=True consistent={r["consistent"]} rank={r["rank"]} nullity={r["nullity"]} {dt:.1f}s', file=sys.stderr, flush=True)
    if r['consistent']:
        print('!!! CONSISTENT:', name, file=sys.stderr, flush=True)
        st, data = oplib.lift_solution(F, A, q, D)
        print('lift:', st, file=sys.stderr, flush=True)
        if st == 'solved':
            json.dump({'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}, open('/tmp/hodge/results/off-paired/witness-candidate.json', 'w'))
            print('!!! WITNESS SAVED', file=sys.stderr, flush=True)

f6 = core.parse('X0**6+X1**6+X2**6'); g6 = core.parse('Y0**6+Y1**6+Y2**6')
Qx = sub(sub(F, f6), g6)
H = core.parse('X0**6+X1**6+Y0**6+Y1**6+X0**3*X1*Y0*Y1/2'); C = core.parse('X1**3*Y1/2+X0**2*X1*Y0/3')
Fi = [pderiv(F, i) for i in range(6)]
D0 = {}
for e, cc in Fi[0].items():
    e2 = list(e); e2[0] -= 1; D0[tuple(e2)] = cc / 6
which = sys.argv[1] if len(sys.argv) > 1 else 'all'

if which in ('all', 'fermat'):
    blocks = {'f6': f6, 'g6': g6, 'f6+g6': add(f6, g6), 'Qx': Qx}
    for bn, B in blocks.items():
        for m in monomials(1) + monomials(2):
            test(f'{bn}*{mstr(m)}', 'fermat-deformation*monomial', shift(B, m))
    for (an, Pa), (bn, Pb) in itertools.combinations_with_replacement(list(blocks.items()) + [('H', H), ('C', C), ('D0', D0)], 2):
        if an in blocks or bn in blocks:
            if an == bn: continue   # same support as the single block (degree <= 6 excluded)
            test(f'{an}*{bn}', 'fermat-deformation*block', mul(Pa, Pb))

if which in ('all', 'disc'):
    Dt = add(scale(power(Qx, 3), 4), scale(mul(mul(f6, f6), g6), 27))
    fs = factors(Dt)
    print('Delta_t factors:', [(core.degree(f), m, len(f)) for f, m in fs], file=sys.stderr, flush=True)
    with open('/tmp/hodge/results/off-paired/delta_t_factors.txt', 'w') as fh:
        for f, m in fs: fh.write(f'mult {m} deg {core.degree(f)}: {core.to_str(f)}\n')
    if len(fs) > 1 or fs[0][1] > 1:
        for f, m in fs:
            if core.degree(f) >= 7: test(f'factor(deg {core.degree(f)}) of Delta_t', 'cubic-cover-discriminant-factor', f)
    test('Delta_t = 4Qx^3+27f6^2g6', 'cubic-cover-discriminant', Dt)
print('done', file=sys.stderr, flush=True)
