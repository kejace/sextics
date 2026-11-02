"""Exact Koszul fixed-D screen of random mu6-homogeneous members of J_d (degree d), all characters, both targets.
usage: struct_dn.py d [ntrials] [--k3-only]"""
import sys, json, random
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder/wave2-excl')
from fixedD_koszul import run
from fractions import Fraction as Fr
from hodge.core import monomials, add
d = int(sys.argv[1]); ntr = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 1
kinds = ('k3',) if '--k3-only' in sys.argv else ('k3', 'full')
random.seed(100 + d)
def vdeg(m): return m[3] + m[4] + m[5]
def rand_form(deg, vd):
    q = {}
    for m in monomials(deg):
        if vdeg(m) == vd:
            c = random.randrange(-3, 4)
            if c: q = add(q, {tuple(m): Fr(c)})
    return q
out = f'fixedD-koszul-struct-d{d}.jsonl'
for b in range(6):
    for t in range(ntr):
        a = []
        for i in range(6):
            want = b if i < 3 else (b + 1) % 6
            a.append(rand_form(d - 5, want) if want <= d - 5 else {})
        if all(not x for x in a): continue
        label = f'rand-J{d}-char{b}-{t}'
        for kind in kinds:
            try:
                res = run(label + '-' + kind, d, a, kind, 65521, False, verbose=False)
                print(json.dumps({k: res[k] for k in ('label', 'target', 'consistent_Q', 'consistent_p', 'seconds')}), [(bl['char'], bl['rows'], bl['cols'], bl['rank_Q'], bl['aug_rank_Q']) for bl in res['blocks']], flush=True)
                with open(out, 'a') as fh: fh.write(json.dumps(res) + '\n')
            except Exception as ex:
                print('ERROR', label, kind, repr(ex), flush=True)
