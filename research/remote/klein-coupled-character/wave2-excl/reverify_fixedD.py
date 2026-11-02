"""(a) re-verify wave-1 fixed-D exclusions with the engine-v2 character-restricted system (independent of the
wave-1 Koszul formulation) at a NEW prime, and upgrade each to an exact rational Farkas certificate.
Families: polars F_i (6), Z_j*F_i (36, both parities), F_i*F_j (21), xyztrs*F_i (6, d=11), Jhat (d=9).
Certificates -> certs/<name>.json; one log line per D."""
import sys, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
import kx
from hodge import core
from hodge.core import mul, shift, parse
p = int(sys.argv[1]); which = sys.argv[2].split(',') if len(sys.argv) > 2 else ['polar', 'ZF', 'FF', 'xyztrsF', 'jhat']
F, A, q = kx.load()
Fi = kx.polars(F)
names = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
tests = []
if 'jhat' in which:
    tests.append(('Jhat', kx.jhat(), 'jhat'))
if 'polar' in which:
    for i in range(6): tests.append((f'F_{i}', Fi[i], 'polar'))
if 'ZF' in which:
    for j in range(6):
        for i in range(6):
            tests.append((f'{names[j]}*F_{i}', mul(parse(names[j]), Fi[i]), 'Z*polar'))
if 'FF' in which:
    for i in range(6):
        for j in range(i, 6):
            tests.append((f'F_{i}*F_{j}', mul(Fi[i], Fi[j]), 'polar*polar'))
if 'xyztrsF' in which:
    xyztrs = parse('X0*X1*X2*Y0*Y1*Y2')
    for i in range(6): tests.append((f'xyztrs*F_{i}', mul(xyztrs, Fi[i]), 'xyztrs*polar'))
summary = []
for name, D, fam in tests:
    t0 = time.time()
    try:
        res = kx.farkas(F, A, q, D, p, save=f'certs/{name.replace("*","_")}_p{p}.json', Dstr=name)
    except Exception as ex:
        res = {'status': 'error', 'error': repr(ex), 'D': name}
    res['family'] = fam; res['name'] = name
    kx.log(exp='fixedD-farkas', **{k: v for k, v in res.items() if k not in ('pivot_cols',)})
    print(name, res.get('status'), 'd', res.get('d'), 'rows', res.get('rows'), 'cols', res.get('cols'), 'rankM', res.get('rank_M'),
          'nullity', res.get('nullity'), 'farkas_rows', res.get('farkas_rows'), 'bits', res.get('bits_max'), round(time.time() - t0, 1), 's', flush=True)
    summary.append({k: res.get(k) for k in ('name', 'family', 'd', 'status', 'rows', 'cols', 'rank_M', 'nullity', 'farkas_rows', 'bits_max', 'sha256_y', 'p', 'charD', 'parD', 'stab')})
json.dump(summary, open(f'fixedD_summary_p{p}_{"-".join(which)}.json', 'w'), indent=1)
print('DONE')
