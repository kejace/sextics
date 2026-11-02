"""Job C (remote): full E1+E2 fixed-D tests at d = 9, 10 for every (parity, sigma-character) sector, prime 1000003:
every E1-kernel representative (canonical rref basis mod p of (J:A)_d modulo J_d), random kernel combinations and
random kernel + J_d combinations; then exact rational dual certificates (rigorous exclusions over Q) for Jhat,
the six polars F_i, the 36 Z_j F_i, the rational d=9 kernel representatives of wave 1 and the 21 polar products.
Outputs: jobC_*.json, jobC_duals/*.json, log lines in wave2-excl.log.jsonl."""
import sys, os, json, time, random
import multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
from w2lib import *
P = int(sys.argv[1]) if len(sys.argv) > 1 else 1000003
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 24
F, A, q = load()
OUT = PACKET + f'/jobC_p{P}'; os.makedirs(OUT, exist_ok=True); os.makedirs(OUT + '/duals', exist_ok=True)
def say(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)

def task_fixedD(args):
    name, Dser, family, note = args
    D = {tuple(json.loads(k)): v for k, v in Dser.items()}
    try:
        e = test_D_restricted(F, A, q, D, P, name, family, note=note, extra={'job': 'jobC'})
        return {k: v for k, v in e.items() if not k.startswith('_')}
    except Exception as ex:
        log({'exp': 'fixedD', 'D': name, 'family': family, 'prime': P, 'result': 'error', 'note': repr(ex), 'job': 'jobC'})
        return {'D': name, 'error': repr(ex)}

def ser(D): return {json.dumps(list(m)): int(v) for m, v in D.items()}

if __name__ == '__main__':
    ctx = Ctx(F, A, P)
    rng = random.Random(20261101)
    tasks = []
    kernel_info = {}
    for d in (9, 10):
        for par in (0, 1):
            for k in range(7):
                t = time.time()
                cs = colon_sector(ctx, d, par, k, want_reps=True)
                kernel_info[f'{d},{par},{k}'] = {'dimS': cs['dimS'], 'dimJ': cs['dimJ'], 'dimColon': cs['dimColon'], 'kernel': cs['kernel']}
                say(f'd={d} par={par} k={k}: kernel {cs["kernel"]} J {cs["dimJ"]} [{time.time()-t:.0f}s]')
                reps = cs['reps']
                for i, D in enumerate(reps):
                    tasks.append((f'E1-kernel rep d={d} par={par} k={k} #{i} (p={P})', ser(D), 'e1-kernel-rep', 'canonical rref kernel basis element mod p'))
                if len(reps) >= 2:
                    for tr in range(3):
                        tasks.append((f'E1-kernel random combo d={d} par={par} k={k} trial {tr} (p={P})', ser(rand_combo(reps, P, rng)), 'e1-kernel-random', 'random GF(p) combination of the kernel reps'))
                if reps:
                    for tr in range(3):
                        D = add(rand_combo(reps, P, rng), rand_combo(cs['jpolys'], P, rng), P)
                        tasks.append((f'E1-kernel+J_d random d={d} par={par} k={k} trial {tr} (p={P})', ser(D), 'e1-kernel-plus-J-random', 'random kernel combo plus random element of J_d in the sector'))
                if d == 10 and not reps and cs['dimJ']:
                    for tr in range(2):
                        tasks.append((f'J_d random d={d} par={par} k={k} trial {tr} (p={P})', ser(rand_combo(cs['jpolys'], P, rng)), 'J_d-random', 'random element of J_d in the sector (kernel empty)'))
    json.dump(kernel_info, open(f'{OUT}/kernel_info.json', 'w'), indent=1)
    say(f'{len(tasks)} fixed-D tasks; running with {NW} workers')
    results = []
    with mp.Pool(NW) as pool:
        for r in pool.imap_unordered(task_fixedD, tasks):
            results.append(r)
            say('done', r.get('D'), r.get('result'), 'rank', r.get('rank'), 'rows', r.get('rows'), 'cols', r.get('cols'), r.get('seconds'), 's')
            json.dump(results, open(f'{OUT}/fixedD_results.json', 'w'), indent=1, default=str)
    say('SECTOR_TESTS_DONE')

    # exact rational dual certificates
    Fstr = open('/home/kejace/devel/sextics/data/klein-coupled-character/instance.rec').read().split('F: ')[1].split('\n')[0]
    Jhat = core.parse(f'7*X0*X2**2*(9*(Y0**6 - Y1**6 + Y2**6) - 9*({Fstr})) + 15*X1*Y2*Y0*(2*Y0**6 - 4*Y1**6 + 8*Y2**6 - 2*({Fstr}))')
    Fi = [core.pderiv(F, i) for i in range(6)]
    Z = [core.parse(v) for v in core.VARS]
    cands = [('Jhat (fqlx, degree 9, char 2)', Jhat, 'companion')]
    cands += [(f'F_{i}', Fi[i], 'polar') for i in range(6)]
    cands += [(f'{core.VARS[j]}*F_{i}', core.mul(Z[j], Fi[i]), 'polar-times-coordinate') for i in range(6) for j in range(6)]
    kr = json.load(open(LANE + '/kernel_d9_rational.json'))
    for key, lst in kr.items():
        for i, s in enumerate(lst):
            cands.append((f'rational E1-kernel rep d=9 (par,k)={key} #{i}', core.parse(s), 'e1-kernel-rep-rational'))
    cands += [(f'F_{i}*F_{j}', core.mul(Fi[i], Fi[j]), 'polar-product') for i in range(6) for j in range(i, 6)]
    def task_dual(args):
        name, Dser, family = args
        D = {tuple(json.loads(k)): Fraction(v) for k, v in Dser.items()}
        t = time.time()
        try:
            r = exact_dual(F, A, q, D, P)
            y = r.pop('y', None)
            fn = f'{OUT}/duals/{sha(name)}.json'
            json.dump({'D': name, 'D_poly': core.to_str(D), 'prime_for_pivots': P, 'certificate': r, 'y': y}, open(fn, 'w'))
            r['file'] = fn
            log({'exp': 'exact-dual', 'D': name, 'd': core.degree(D), 'k': char_of(D), 'par': core.parity_of(D), 'prime': P, 'family': family,
                 'result': 'RIGOROUS-INCONSISTENT (exact rational dual y: y^T M = 0, y^T b = 1 verified over Q)' if r['ok'] else 'dual-failed',
                 'rank': r.get('rank_mod_p'), 'seconds': round(time.time() - t, 1), 'job': 'jobC', 'note': json.dumps({k: v for k, v in r.items() if k != 'file'})[:600]})
            return name, r
        except Exception as ex:
            log({'exp': 'exact-dual', 'D': name, 'prime': P, 'family': family, 'result': 'error', 'note': repr(ex), 'job': 'jobC'})
            return name, {'error': repr(ex)}
    dtasks = [(n, {json.dumps(list(m)): str(v) for m, v in D.items()}, fam) for n, D, fam in cands]
    duals = {}
    with mp.Pool(min(NW, 16)) as pool:
        for name, r in pool.imap_unordered(task_dual, dtasks):
            duals[name] = r
            say('dual', name, r.get('ok'), 'support', r.get('support'), 'bits', r.get('max_bits'), 't', r.get('t_solve'), r.get('t_verify'), r.get('error', ''))
            json.dump(duals, open(f'{OUT}/duals_summary.json', 'w'), indent=1, default=str)
    say('JOBC_DONE')
