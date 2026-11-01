import sys, time, random, json
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from kcc2 import *
F, A, q = load(); p = P0
ctx = Ctx(F, A, p)
random.seed(2026)
def rnd(): return random.randint(1, p - 1)
def combo(polys, coeffs):
    acc = {}
    for P, c in zip(polys, coeffs):
        if c: acc = add(acc, scale(P, c, p), p)
    return acc
import pickle
store = {}
for d in [9, 10]:
    for parD in [0, 1]:
        for k in range(7):
            t = time.time()
            reps, (jk, jp) = e1_kernel(ctx, d, k, parD, verbose=True)
            store[(d, parD, k)] = (reps, jk)
            print(f'd={d} parD={parD} k={k}: nontrivial kernel {len(reps)}, J_d block {len(jp)} gens  [{time.time()-t:.1f}s]', flush=True)
            if not reps: continue
            # (i) each kernel representative
            for idx, D in enumerate(reps):
                e = full_test(F, A, D, f'E1-kernel rep d={d} parD={parD} k={k} #{idx}', 'e1-kernel-rep',
                              note=f'canonical kernel representative mod p (mod J_d), {len(D)} terms', verbose=False)
                print(f'   rep #{idx}: e1={e["e1"]} consistent={e["consistent"]} rank={e["rank"]} nullity={e["nullity"]} rows={e.get("rows")} cols={e.get("cols")} {e["seconds"]}s', flush=True)
            # (ii) random combinations of kernel reps (pure kernel part)
            for trial in range(3 if len(reps) > 1 else 0):
                D = combo(reps, [rnd() for _ in reps])
                e = full_test(F, A, D, f'E1-kernel random combo d={d} parD={parD} k={k} trial {trial}', 'e1-kernel-random',
                              note='random F_p combination of the nontrivial kernel reps (no J_d part)')
                print(f'   random kernel combo {trial}: consistent={e["consistent"]} rank={e["rank"]} nullity={e["nullity"]} {e["seconds"]}s', flush=True)
            # (iii) random kernel + random J_d part
            for trial in range(3):
                D = combo(reps + jp, [rnd() for _ in reps] + [rnd() for _ in jp])
                e = full_test(F, A, D, f'E1-kernel+J_d random d={d} parD={parD} k={k} trial {trial}', 'e1-kernel-plus-J-random',
                              note=f'random F_p point of (kernel reps + J_d block, {len(reps)}+{len(jp)} gens)')
                print(f'   random kernel+J {trial}: consistent={e["consistent"]} rank={e["rank"]} nullity={e["nullity"]} {e["seconds"]}s', flush=True)
    pickle.dump({str(k): v for k, v in store.items()}, open(f'kernel_store_d{d}.pkl', 'wb'))
