"""(c) Sector-restricted full E1+E2 solves (character-twisted reduced system, charsolve.CharSystem) of every
non-Jacobian E1-kernel element in the S3xS3 character sectors (triv, sgnX, sgnY, sgnsgn), both Y-parities,
degrees d = 8..dmax, plus random rational combinations of the kernel classes and random points of the full
E1-admissible space (kernel + J^chi_d).  All results are mod p (first prime) unless consistent, then CRT-lifted.

Usage: python3 c_sectors.py dmin dmax [types...]
"""
import sys, time, json, pickle, os, random
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual/wave2-excl')
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from xlib import *
from charsolve import Group, e1_kernel_char, solve_char, J_char_basis, CharSystem
from fractions import Fraction

F, A, q = load()
dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
types = sys.argv[3:] or ['triv', 'sgnX', 'sgnY', 'sgnsgn']
groups = {'triv': Group([[0,1,2],[3,4,5]], [1,1]), 'sgnX': Group([[0,1,2],[3,4,5]], [-1,1]),
          'sgnY': Group([[0,1,2],[3,4,5]], [1,-1]), 'sgnsgn': Group([[0,1,2],[3,4,5]], [-1,-1]),
          'S22_++': Group([[1,2],[4,5]], [1,1]), 'S22_-+': Group([[1,2],[4,5]], [-1,1]),
          'S22_+-': Group([[1,2],[4,5]], [1,-1]), 'S22_--': Group([[1,2],[4,5]], [-1,-1]),
          'S3S2_++': Group([[0,1,2],[4,5]], [1,1]), 'S3S2_-+': Group([[0,1,2],[4,5]], [-1,1]),
          'S3S2_+-': Group([[0,1,2],[4,5]], [1,-1]), 'S3S2_--': Group([[0,1,2],[4,5]], [-1,-1]),
          'S2S3_++': Group([[1,2],[3,4,5]], [1,1]), 'S2S3_-+': Group([[1,2],[3,4,5]], [-1,1]),
          'S2S3_+-': Group([[1,2],[3,4,5]], [1,-1]), 'S2S3_--': Group([[1,2],[3,4,5]], [-1,-1])}
t00 = time.time()
def say(*a): print(f'[{time.time()-t00:7.1f}s]', *a, flush=True)
rng = random.Random(20261102)
def rand_comb(polys, k=None):
    out = {}
    for P in polys:
        c = rng.randint(-5, 5)
        if c: out = core.add(out, core.scale(P, Fraction(c)))
    return out
def to_frac(P):
    return {e: Fraction(c) for e, c in P.items()}
def run(D, label, d, pD, typ, note=''):
    t0 = time.time()
    try:
        status, data = solve_char(F, A, q, D, groups[typ], max_primes=4)
    except Exception as ex:
        say(label, 'ERROR', repr(ex)[:200]); log(experiment='c:sector-solve', D=label, d=d, pD=pD, type=typ, error=repr(ex)[:200]); return
    rec = dict(experiment='c:sector-solve', D=label, d=d, pD=pD, type=typ, status=status, consistent=(status not in ('inconsistent',)),
               rank=data.get('rank'), nullity=data.get('nullity'), rows=data.get('rows'), cols=data.get('cols'), prime=data.get('prime'),
               seconds=round(time.time() - t0, 1), note=note)
    log(**rec)
    say(label, status, f"rank={data.get('rank')} null={data.get('nullity')} rows={data.get('rows')} cols={data.get('cols')} ({rec['seconds']}s)")
    if status == 'solved':
        wit = dict(d=d, D=core.to_str(D), N=[core.to_str(x) for x in data['N']], Q=core.to_str(data['Q']))
        with open(f'{PK}/WITNESS_{label}.json', 'w') as f: json.dump(wit, f)
        say('!!!!!!!! WITNESS FOUND', label)
    elif status not in ('inconsistent',):
        with open(f'{PK}/CONSISTENT_{label}.pkl', 'wb') as f: pickle.dump((D, data), f)
        say('!!!!!!!! consistent mod p but not lifted', label)
for d in range(dmin, dmax + 1):
    for typ in types:
        G_ = groups[typ]
        for pD in (0, 1):
            pk = f'{RES}/kchar_{typ}_d{d}_p{pD}.pkl'
            if os.path.exists(pk):
                info, polys = pickle.load(open(pk, 'rb'))
            else:
                t0 = time.time(); info, polys = e1_kernel_char(F, A, q, d, pD, G_); info['seconds'] = round(time.time() - t0, 1); info['type'] = typ
                pickle.dump((info, polys), open(pk, 'wb'))
            say(f'd={d} {typ} pD={pD}: kernel beyond J = {info["ker"]} (dimJ={info["dimJ"]})')
            log(experiment='c:sector-kernel', d=d, pD=pD, type=typ, **{k: v for k, v in info.items() if k not in ('d', 'pD', 'type')})
            if not polys: continue
            polys = [to_frac(P) for P in polys]
            for i, P in enumerate(polys):
                run(P, f'K{typ}_d{d}_p{pD}_rep{i}', d, pD, typ, note='kernel basis element')
            nrand = 4 if len(polys) > 1 else 0
            for r in range(nrand):
                run(rand_comb(polys), f'K{typ}_d{d}_p{pD}_rand{r}', d, pD, typ, note='random small-integer combination of kernel basis')
            # random points of the full E1-admissible space (kernel + J^chi_d)
            try:
                Jb = [to_frac({e: sym_int(c, P1) for e, c in Pj.items()}) for Pj in J_char_basis(core.reduce_mod(F, P1), d, pD, G_, P1)]
            except Exception as ex:
                Jb = []
            for r in range(2):
                D = core.add(rand_comb(polys), rand_comb(Jb[:40]))
                if D: run(D, f'KJ{typ}_d{d}_p{pD}_rand{r}', d, pD, typ, note='random point of kernel + J^chi_d (first 40 J generators)')
say('DONE')
