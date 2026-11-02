"""Exact fixed-D tests (Koszul normalization) of structured mu6-homogeneous degree-7 elements of J for the K3 target
(and the full target): D = sum_i P_i a_i with a_i quadrics of prescribed character. Random small-integer coefficients
=> generic members of the character-b slice of J_7 (irreducible for generic choices). Receipts fixedD-koszul-struct-d7.jsonl."""
import sys, json, random
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder/wave2-excl')
from fixedD_koszul import run, z, mon
from fractions import Fraction as Fr
from hodge.core import monomials, add, scale
random.seed(2026)
def vdeg(m): return m[3] + m[4] + m[5]
def rand_quadric(vd):
    q = {}
    for m in monomials(2):
        if vdeg(m) == vd:
            c = random.randrange(-3, 4)
            if c: q = add(q, {tuple(m): Fr(c)})
    return q
out = 'fixedD-koszul-struct-d7.jsonl'
# character b: a_i (i<3) needs vdeg == b, a_i (i>=3) needs vdeg == b+1 (mod 6); quadrics have vdeg in {0,1,2}
cases = []
for b in (0, 1, 2):
    for t in range(2):
        a = []
        for i in range(6):
            want = b if i < 3 else (b + 1) % 6
            a.append(rand_quadric(want) if want <= 2 else {})
        if all(not x for x in a): continue
        cases.append((f'rand-J7-char{b}-{t}', a))
for label, a in cases:
    for kind in ('k3', 'full'):
        try:
            res = run(label + '-' + kind, 7, a, kind, 65521, False, verbose=False)
            print(json.dumps({k: res[k] for k in ('label', 'target', 'consistent_Q', 'consistent_p', 'seconds')}), [(b['char'], b['rows'], b['cols'], b['rank_Q'], b['aug_rank_Q']) for b in res['blocks']], flush=True)
            with open(out, 'a') as fh: fh.write(json.dumps(res) + '\n')
        except Exception as ex:
            print('ERROR', label, kind, repr(ex), flush=True)
