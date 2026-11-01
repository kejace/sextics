import sys, time, random, pickle
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from kcc2 import *
F, A, q = load(); p = P0
ctx = Ctx(F, A, p)
random.seed(1010)
def rnd(): return random.randint(1, p - 1)
def combo(polys, coeffs):
    acc = {}
    for P, c in zip(polys, coeffs):
        if c: acc = add(acc, scale(P, c, p), p)
    return acc
store = {}
d = 10
for parD in [0, 1]:
    for k in range(7):
        t = time.time()
        reps, (jk, jp) = e1_kernel(ctx, d, k, parD, verbose=True)
        store[(d, parD, k)] = reps
        print(f'd={d} parD={parD} k={k}: nontrivial kernel {len(reps)}, J_d block {len(jp)} gens  [{time.time()-t:.1f}s]', flush=True)
        if not reps: continue
        for idx, D in enumerate(reps[:2]):
            e = full_test(F, A, D, f'E1-kernel rep d={d} parD={parD} k={k} #{idx}', 'e1-kernel-rep', note=f'canonical kernel representative mod p (mod J_d); kernel dim {len(reps)}')
            print(f'   rep #{idx}: consistent={e["consistent"]} rank={e["rank"]} nullity={e["nullity"]} rows={e.get("rows")} cols={e.get("cols")} {e["seconds"]}s', flush=True)
        for trial in range(2):
            D = combo(reps, [rnd() for _ in reps])
            e = full_test(F, A, D, f'E1-kernel random combo d={d} parD={parD} k={k} trial {trial}', 'e1-kernel-random', note=f'random F_p combination of the {len(reps)} kernel reps (no J_d part)')
            print(f'   random kernel combo {trial}: consistent={e["consistent"]} rank={e["rank"]} nullity={e["nullity"]} {e["seconds"]}s', flush=True)
        for trial in range(2):
            D = combo(reps + jp, [rnd() for _ in reps] + [rnd() for _ in jp])
            e = full_test(F, A, D, f'E1-kernel+J_d random d={d} parD={parD} k={k} trial {trial}', 'e1-kernel-plus-J-random', note=f'random F_p point of kernel({len(reps)}) + J_d({len(jp)})')
            print(f'   random kernel+J {trial}: consistent={e["consistent"]} rank={e["rank"]} nullity={e["nullity"]} {e["seconds"]}s', flush=True)
pickle.dump({str(kk): v for kk, v in store.items()}, open('kernel_store_d10.pkl', 'wb'))
