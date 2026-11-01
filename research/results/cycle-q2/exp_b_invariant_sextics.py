"""Experiment 1/3: invariant sextic supports.  D = subset sums of the 9 monomials of F (with F's
coefficients), the pure-power part P, the coupling sum M, plus the three invariant 'square' monomials.
Each support is tested with D (d=6) and D^2 (d=12, the definitive test at q=2)."""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, numpy as np
from hodge import core
from fractions import Fraction

F, A, q, rec = core.load_instance('cycle-q2')
mons = sorted(F.keys())
names = [core.to_str({m: 1}) for m in mons]
print('F monomials:', names, flush=True)
# complementary subsets give the same support (D + (F-D) = F); keep the representative containing X0**6 or not: just dedupe by support via frozenset pairs
seen = set()
subsets = []
for r in range(1, 9):
    for S in itertools.combinations(range(9), r):
        comp = tuple(sorted(set(range(9)) - set(S)))
        key = min(S, comp)
        if key in seen: continue
        seen.add(key)
        subsets.append(S)
print('distinct supports from F-subsets:', len(subsets), flush=True)

def Dsum(S, coeffs=True):
    D = {}
    for i in S:
        D[mons[i]] = F[mons[i]] if coeffs else Fraction(1)
    return D

results = []
t0 = time.time()
def run(D, fam, note):
    r6 = eq.test_D(F, A, q, D, fam, note=note + ' [D]', verbose=False)
    D2 = core.mul(D, D)
    r12 = eq.test_D(F, A, q, D2, fam, note=note + ' [D^2]', verbose=False, Dstr='(' + core.to_str(D) + ')**2')
    s = f'{fam:14s} {core.to_str(D)[:70]:70s} | d=6: e1={r6.get("e1")} cons={r6.get("consistent")} rk={r6.get("rank")} nul={r6.get("nullity")} | d=12: e1={r12.get("e1")} cons={r12.get("consistent")} rk={r12.get("rank")} nul={r12.get("nullity")} cols={r12.get("cols")} [{time.time()-t0:.0f}s]'
    print(s, flush=True)
    results.append({'D': core.to_str(D), 'fam': fam, 'r6': {k: v for k, v in r6.items() if not k.startswith('_')}, 'r12': {k: v for k, v in r12.items() if not k.startswith('_')}})
    if r6.get('consistent') or r12.get('consistent'):
        print('!!!!! CONSISTENT SYSTEM FOUND !!!!!', core.to_str(D), flush=True)

# special named ones first
P = core.parse('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')
M = core.parse('X0**3*X1*Y0*Y1 + X1**3*X2*Y1*Y2 + X0*X2**3*Y0*Y2')
run(P, 'pure-power', 'P = sum of sixth powers (same support as M)')
run(M, 'coupling-sum', 'M = sum of the three couplings')
run(core.parse('X0**6 + X1**6 + X2**6'), 'X-block', 'sum X_i^6')
run(core.parse('Y0**6 + Y1**6 + Y2**6'), 'Y-block', 'sum Y_i^6')
for i in range(3):
    run(core.parse(f'X{i}**6 + Y{i}**6'), 'pair-fermat', f'X{i}^6 + Y{i}^6')
# edge pieces F restricted to a pair of pairs
run(core.parse('X0**6 + X1**6 + Y0**6 + Y1**6 + X0**3*X1*Y0*Y1/2'), 'edge', 'F|_{X2=Y2=0}')
run(core.parse('X1**6 + X2**6 + Y1**6 + Y2**6 + X1**3*X2*Y1*Y2/2'), 'edge', 'F|_{X0=Y0=0}')
run(core.parse('X0**6 + X2**6 + Y0**6 + Y2**6 + X0*X2**3*Y0*Y2/2'), 'edge', 'F|_{X1=Y1=0}')
# F minus one coupling (same support as that coupling alone)
for c in ['X0**3*X1*Y0*Y1', 'X1**3*X2*Y1*Y2', 'X0*X2**3*Y0*Y2']:
    run(core.parse(c), 'single-coupling', f'{c} (support = F - coupling/2)')
# the three invariant squares and sums
sq = ['X1**2*Y0**2*Y1**2', 'X2**2*Y1**2*Y2**2', 'X0**2*Y0**2*Y2**2']
run(core.parse(' + '.join(sq)), 'inv-squares', 'sum of three invariant squares')
for s_ in sq:
    run(core.add(P, core.parse(s_)), 'P+square', f'P + {s_}')
    run(core.add(M, core.parse(s_)), 'M+square', f'M + {s_}')
run(core.add(core.add(P, core.parse(' + '.join(sq))), {}), 'P+squares', 'P + all three squares')
run(core.add(M, core.parse(' + '.join(sq))), 'M+squares', 'M + all three squares')
# P + t M for a few t (all same support except t=1/2 which is F) -- t=1/2 excluded; t = -1 and 2 as checks of gauge invariance of the test
run(core.add(P, core.scale(M, Fraction(2))), 'pencil-check', 'P + 2M (should match P)')
# all subsets with F coefficients
for S in subsets:
    D = Dsum(S)
    if len(S) == 9: continue
    run(D, 'F-subset', 'subset ' + ','.join(names[i] for i in S))
# subsets with unit coefficients on couplings (different supports when mixed)
for S in subsets:
    if any(F[mons[i]] != 1 for i in S) and any(F[mons[i]] == 1 for i in S):
        D = Dsum(S, coeffs=False)
        run(D, 'F-subset-unit', 'unit-coeff subset ' + ','.join(names[i] for i in S))
json.dump(results, open('/tmp/hodge/results/cycle-q2/exp_b_invariant_sextics.json', 'w'), indent=1, default=str)
print('DONE', flush=True)
