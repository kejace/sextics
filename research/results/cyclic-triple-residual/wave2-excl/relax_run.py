"""Hessian-relaxation driver.  usage: python3 relax_run.py d p threads [tag]
Builds the relaxation quadrics in degree d (mod p), writes msolve input, computes the reduced GB (mod p), then for
every linear form l vanishing on the illegal locus runs msolve on I + (l - 1) (empty <=> V(I) inside {l = 0}).
Outputs: relax_d{d}_p{p}.json receipt, relax_d{d}_p{p}.ms (ideal), relax_d{d}_p{p}_gb.txt, per-form outputs in a subdir."""
import sys, os, time, json, subprocess
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from w2lib import *
from multiprocessing import Pool

d = int(sys.argv[1]); p = int(sys.argv[2]); threads = int(sys.argv[3]); tag = sys.argv[4] if len(sys.argv) > 4 else ''
S = Setup()
P = D8 = None
if d == 8:
    dat = json.load(open(PK + '/e1_d8_exact.json'))
    D8 = {tuple(k_): v for k_, v in ((tuple(map(int, k.split(','))), v) for k, v in dat['D8z'].items())}
    P = [{tuple(map(int, k.split(','))): v for k, v in Pi.items()} for Pi in dat['Pz']]
t0 = time.time()
th, quads, info = relaxation_quadrics(S, d, p, P=P, D8=D8)
names = [f't{k}' for k in range(len(th))]
base = f'{PK}/relax_d{d}_p{p}{tag}'
polys = [quad_str(q, names, p) for q in quads]
write_msolve(base + '.ms', len(th), p, polys, names)
json.dump({'d': d, 'p': p, 'theta': [list(t) if t[0] != 'c' else ['c'] for t in th], 'n_quadrics': len(quads), 'sectors': info},
          open(base + '_system.json', 'w'))
print(f'built {len(quads)} quadrics in {len(th)} vars ({time.time()-t0:.1f}s)', flush=True)
# 1. reduced Groebner basis of I (mod p)
rc, dt, so, se = run_msolve(base + '.ms', base + '_gb.txt', threads=threads, gb=2)
gb = parse_msolve_gb(base + '_gb.txt') if rc == 0 else []
nlin = sum(1 for g in gb if g and '^' not in g and '*' not in g.replace('*t', '#') ) if False else None
# count linear elements: no monomial of degree >= 2 (crude: no '^' and no product of two variables)
def is_linear(s):
    s2 = s.replace(' ', '')
    import re
    return not re.search(r't\d+\^', s2) and not re.search(r't\d+\*t\d+', s2)
nlin = sum(1 for g in gb if is_linear(g))
print(f'GB: rc={rc} {len(gb)} elements, {nlin} linear, {dt:.1f}s; msolve: {so.strip()[-300:]}', flush=True)
log({'exp': f'hessian-relaxation-gb', 'family': f'E1-admissible degree {d}', 'd': d, 'prime': p, 'result': f'GB has {len(gb)} elements, {nlin} linear',
     'rank': None, 'seconds': round(dt, 1), 'job': f'excl-ctr-relax{d}{tag}', 'note': f'{len(quads)} quadrics in {len(th)} unknowns; msolve -g 2 mod p; rc={rc}'})
# 2. for each illegal-locus linear form l: I + (l - 1) empty?
forms = illegal_forms(d, th)
os.makedirs(base + '_forms', exist_ok=True)
def one(k):
    f = forms[k]
    inp = f'{base}_forms/f{k}.ms'; out = f'{base}_forms/f{k}.out'
    write_msolve(inp, len(th), p, polys + [lin_str(f, names) + '-1'], names)
    rc, dt, so, se = run_msolve(inp, out, threads=1, gb=0)
    txt = open(out).read().strip() if os.path.exists(out) else ''
    empty = txt.startswith('[-1]')
    return {'k': k, 'form': lin_str(f, names), 'rc': rc, 'seconds': round(dt, 1), 'empty': empty, 'out': txt[:200]}
t1 = time.time()
with Pool(min(len(forms), threads)) as pool:
    res = pool.map(one, range(len(forms)))
n_empty = sum(1 for r in res if r['empty'])
json.dump({'d': d, 'p': p, 'n_theta': len(th), 'n_quadrics': len(quads), 'gb_size': len(gb), 'gb_linear': nlin,
           'n_forms': len(forms), 'n_empty': n_empty, 'forms': res, 'seconds_forms': round(time.time() - t1, 1)}, open(base + '.json', 'w'), indent=1)
print(f'forms: {n_empty}/{len(forms)} give empty V(I + (l-1)) ({time.time()-t1:.1f}s)', flush=True)
log({'exp': 'hessian-relaxation-illegal-locus', 'family': f'E1-admissible degree {d}', 'd': d, 'prime': p,
     'result': f'{n_empty}/{len(forms)} linear forms l vanishing on the illegal locus have V(I+(l-1)) empty mod p'
               + (' => V(I) subset of illegal locus (mod p)' if n_empty == len(forms) else ' => extra components!'),
     'rank': None, 'seconds': round(time.time() - t1, 1), 'job': f'excl-ctr-relax{d}{tag}', 'note': 'msolve zero-dim emptiness test mod p per form; modular evidence, exact certificates separate'})
