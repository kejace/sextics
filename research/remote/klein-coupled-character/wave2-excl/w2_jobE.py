"""Job E (remote): d = 11, 12 sector probes (prime 1000003): E1-kernel dimensions per sector, 2 random kernel combos and
2 random kernel + J_d combos per sector through the full restricted E1+E2 system (dense or compressed rref).
Modular screening only (no exact duals at this size)."""
import sys, os, json, time, random
import multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
from w2lib import *
P = int(sys.argv[1]) if len(sys.argv) > 1 else 1000003
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 16
DEGS = [int(a) for a in sys.argv[3].split(',')] if len(sys.argv) > 3 else [11, 12]
F, A, q = load()
OUT = PACKET + f'/jobE_p{P}'; os.makedirs(OUT, exist_ok=True)
def say(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)
def task(args):
    name, Dser, family, note = args
    D = {tuple(json.loads(k)): v for k, v in Dser.items()}
    try:
        e = test_D_restricted(F, A, q, D, P, name, family, note=note, compressed=True, extra={'job': 'jobE'})
        return {k: v for k, v in e.items() if not k.startswith('_')}
    except Exception as ex:
        log({'exp': 'fixedD', 'D': name, 'family': family, 'prime': P, 'result': 'error', 'note': repr(ex), 'job': 'jobE'})
        return {'D': name, 'error': repr(ex)}
def ser(D): return {json.dumps(list(m)): int(v) for m, v in D.items()}
if __name__ == '__main__':
    ctx = Ctx(F, A, P); rng = random.Random(7)
    tasks = []; info = {}
    for d in DEGS:
        for par in (0, 1):
            for k in range(7):
                t = time.time()
                cs = colon_sector(ctx, d, par, k, want_reps=True)
                info[f'{d},{par},{k}'] = {'dimS': cs['dimS'], 'dimJ': cs['dimJ'], 'dimColon': cs['dimColon'], 'kernel': cs['kernel'], 'seconds': round(time.time() - t, 1)}
                say(f'd={d} par={par} k={k}: kernel {cs["kernel"]} J {cs["dimJ"]} colon {cs["dimColon"]} [{time.time()-t:.0f}s]')
                log({'exp': 'colon-kernel-sector', 'd': d, 'par': par, 'k': k, 'prime': P, 'result': info[f'{d},{par},{k}'], 'job': 'jobE'})
                reps = cs['reps']
                if reps:
                    for tr in range(2):
                        tasks.append((f'E1-kernel random combo d={d} par={par} k={k} trial {tr} (p={P})', ser(rand_combo(reps, P, rng)), 'e1-kernel-random', f'random GF(p) combination of the {len(reps)} kernel reps'))
                        tasks.append((f'E1-kernel+J_d random d={d} par={par} k={k} trial {tr} (p={P})', ser(add(rand_combo(reps, P, rng), rand_combo(cs['jpolys'], P, rng), P)), 'e1-kernel-plus-J-random', 'random kernel combo plus random element of J_d'))
                json.dump(info, open(f'{OUT}/kernel_info.json', 'w'), indent=1)
    say(len(tasks), 'tasks')
    results = []
    with mp.Pool(NW) as pool:
        for r in pool.imap_unordered(task, tasks):
            results.append(r); say('done', r.get('D'), r.get('result'), 'rank', r.get('rank'), 'rows', r.get('rows'), 'cols', r.get('cols'), r.get('seconds'), 's')
            json.dump(results, open(f'{OUT}/results.json', 'w'), indent=1, default=str)
    say('JOBE_DONE')
