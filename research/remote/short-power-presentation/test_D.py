"""Test a denominator: E1, full fixed-D consistency (pure parity: engine fast solver; mixed: MixedSystem),
optional target_image.  Usage: python3 test_D.py <json-file-with-list-of-{D,...}> [--image] [--family name]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hlib import *
from mixed_solver import mixed_consistent

def test(F, A, q, D, family, note='', image=False):
    d = core.degree(D); par = core.parity_of(D)
    t0 = time.time()
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    t1 = time.time()
    print(f'[E1] d={d} par={par} {family}: {ok} {info} ({t1-t0:.0f}s)', flush=True)
    if not ok:
        log(core.to_str(D), d, family, False, False, info['rank'], info['cols'] - info['rank'], t1 - t0, note + ' E1 fails')
        return None
    if par is None:
        out = mixed_consistent(F, A, q, D)
        status = 'consistent-modp' if out['consistent'] else 'inconsistent'
    else:
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=True)
        status = 'solved' if out.get('exact') else ('inconsistent' if not out['consistent'] else 'lift-failed')
    t2 = time.time()
    print(f'[E2] {family}: {status} rank={out["rank"]} nullity={out["nullity"]} rows={out["rows"]} cols={out["cols"]} ({t2-t1:.0f}s)', flush=True)
    extra = ''
    if image and par is not None:
        dim, basis, S, dt = served(F, q, D)
        extra = f' served-dim={dim} ({dt:.0f}s)'
        print(f'[image] {family}: served dim {dim} ({dt:.0f}s)', flush=True)
    log(core.to_str(D), d, family, True, out['consistent'], out['rank'], out['nullity'], time.time() - t0,
        note + f' status={status} rows={out["rows"]} cols={out["cols"]}' + extra)
    if status == 'solved':
        with open(f'{RES}/witness-candidate.json', 'w') as f:
            json.dump({'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q'])}, f)
        print('!!! WITNESS FOUND — saved witness-candidate.json', flush=True)
    return out

if __name__ == '__main__':
    F, A, q, rec = load()
    path = sys.argv[1]; image = '--image' in sys.argv
    fam = 'E1-kernel'
    if '--family' in sys.argv: fam = sys.argv[sys.argv.index('--family') + 1]
    items = json.load(open(path))
    for it in items:
        D = core.parse(it['D'])
        test(F, A, q, D, it.get('family', fam), note=it.get('note', ''), image=image)
