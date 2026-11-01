"""Fixed-D consistency tests (engine v2, compressed rows, mod p) for kernel members.
Usage: python3 ktest.py gens <jsonfile>       -- rational generators from kstruct2 json (exact lift attempted if consistent)
       python3 ktest.py random <d> <par> <n>  -- n random mod-p members of the saved E1 basis at degree d / parity
       python3 ktest.py expr '<D expr>' '<label>' '<family>'"""
import sys, json, time, random
import lib
from lib import *
p = core.PRIMES[0]

def test_modp(Dp, label, family, d, note=''):
    t0 = time.time()
    S = fixed_d.FixedDSystem(F, A, q, Dp, p)
    r = fixed_d.solve_compressed(S)
    sec = round(time.time() - t0, 1)
    print(f'{label}: consistent={r["consistent"]} rank={r["rank"]} nullity={r["nullity"]} rows={S.nrows} cols={S.ncols} ({sec}s)', flush=True)
    log({'D': label, 'd': d, 'family': family, 'e1': True, 'consistent': bool(r['consistent']), 'rank': int(r['rank']),
         'nullity': int(r['nullity']), 'rows': S.nrows, 'cols': S.ncols, 'seconds': sec, 'note': note, 'engine': 2})
    return r

mode = sys.argv[1]
if mode == 'gens':
    J = json.load(open(sys.argv[2]))
    for i, s in enumerate(J['generators']):
        D = core.parse(s)
        label = f'K{J["d"]} generator {i} (rational, base-type; target {J["target"]})'
        r = run_D(D, label, f'e1-kernel-d{J["d"]}-gen', note='rational generator of annihilator mod J; engine v2', skip_e1=False)
        show(r)
elif mode == 'random':
    d, par, n = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    E = [core.parse(s, mod=p) for s in json.load(open(f'{WORK}/e1kernel_d{d}_par{par}.json'))['basis_modp']]
    random.seed(int(sys.argv[5]) if len(sys.argv) > 5 else 1)
    for k in range(n):
        D = {}
        for e in E: D = core.add(D, core.scale(e, random.randint(1, p - 1), p), p)
        test_modp(D, f'random member #{k} of E1 family d={d} par={par} ({len(E)}-dim, mod p)', f'e1-family-d{d}-random', d,
                  note='mod-p generic member of the affine E1-admissible family (kernel + J)')
elif mode == 'expr':
    D = core.parse(sys.argv[2])
    r = run_D(D, sys.argv[3], sys.argv[4], note='engine v2')
    show(r)
