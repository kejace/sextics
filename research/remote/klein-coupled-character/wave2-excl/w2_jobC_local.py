"""Local reduced version of job C (remote saturated): d = 9 all sectors: every kernel rep + 2 random kernel combos +
2 random kernel+J combos; d = 10: 2 random kernel combos + 1 kernel+J per sector.  Prime 1000003, 4 workers, nice."""
import sys, os, json, time, random
import multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
from w2lib import *
P = 1000003; NW = 4
F, A, q = load()
OUT = PACKET + f'/jobC_local_p{P}'; os.makedirs(OUT, exist_ok=True)
def say(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)
def task(args):
    name, Dser, family, note = args
    D = {tuple(json.loads(k)): v for k, v in Dser.items()}
    try:
        e = test_D_restricted(F, A, q, D, P, name, family, note=note, extra={'job': 'jobC-local'})
        return {k: v for k, v in e.items() if not k.startswith('_')}
    except Exception as ex:
        log({'exp': 'fixedD', 'D': name, 'family': family, 'prime': P, 'result': 'error', 'note': repr(ex), 'job': 'jobC-local'})
        return {'D': name, 'error': repr(ex)}
def ser(D): return {json.dumps(list(m)): int(v) for m, v in D.items()}
if __name__ == '__main__':
    ctx = Ctx(F, A, P); rng = random.Random(20261101)
    tasks = []; info = {}
    for d in (9, 10):
        for par in (0, 1):
            for k in range(7):
                t = time.time()
                cs = colon_sector(ctx, d, par, k, want_reps=True)
                info[f'{d},{par},{k}'] = {'dimS': cs['dimS'], 'dimJ': cs['dimJ'], 'dimColon': cs['dimColon'], 'kernel': cs['kernel'], 'HF_R': hilbert_R_sector(d, par, k), 'HF_ok': cs['dimS'] - cs['dimJ'] == hilbert_R_sector(d, par, k)}
                say(f'd={d} par={par} k={k}: kernel {cs["kernel"]} J {cs["dimJ"]} HF_ok {info[f"{d},{par},{k}"]["HF_ok"]} [{time.time()-t:.0f}s]')
                log({'exp': 'colon-kernel-sector', 'd': d, 'par': par, 'k': k, 'prime': P, 'result': info[f'{d},{par},{k}'], 'seconds': round(time.time()-t, 1), 'job': 'jobC-local'})
                reps = cs['reps']
                if d == 9:
                    for i, D in enumerate(reps):
                        tasks.append((f'E1-kernel rep d=9 par={par} k={k} #{i} (p={P})', ser(D), 'e1-kernel-rep', 'canonical rref kernel basis element mod p'))
                nr = 2 if d == 9 else 2
                if len(reps) >= 2:
                    for tr in range(nr):
                        tasks.append((f'E1-kernel random combo d={d} par={par} k={k} trial {tr} (p={P})', ser(rand_combo(reps, P, rng)), 'e1-kernel-random', 'random GF(p) combination of the kernel reps'))
                if reps:
                    for tr in range(2 if d == 9 else 1):
                        tasks.append((f'E1-kernel+J_d random d={d} par={par} k={k} trial {tr} (p={P})', ser(add(rand_combo(reps, P, rng), rand_combo(cs['jpolys'], P, rng), P)), 'e1-kernel-plus-J-random', 'random kernel combo plus random element of J_d'))
                json.dump(info, open(f'{OUT}/kernel_info.json', 'w'), indent=1)
    say(len(tasks), 'tasks')
    results = []
    with mp.Pool(NW) as pool:
        for r in pool.imap_unordered(task, tasks):
            results.append(r); say('done', r.get('D'), r.get('result'), 'rank', r.get('rank'), 'rows', r.get('rows'), 'cols', r.get('cols'), r.get('seconds'), 's')
            json.dump(results, open(f'{OUT}/fixedD_results.json', 'w'), indent=1, default=str)
    say('JOBC_LOCAL_DONE')
