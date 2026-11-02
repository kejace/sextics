"""Part (1): E1-admissible denominator spaces per degree d = 4..12, both Y-parities, all sigma-weights k.
Stores dims + small-height rational kernel representatives (mod J_d) + exact E1 lifts in w2_sectors.pkl,
appends one JSON line per sector to wave2-search.log.jsonl, writes w2_kernels.out."""
import os, sys, time, json, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import w2lib as L
from hodge import core
from multiprocessing import Pool
OUT = '/tmp/hodge/results/klein-coupled-character/w2_kernels.out'
PKL = '/tmp/hodge/results/klein-coupled-character/w2_sectors.pkl'
F, A, q = L.load()

def work(task):
    d, parD, k = task
    t = time.time()
    info = L.sector_space(F, A, d, parD, k, nprimes=12, verbose=False, want_rational=True)
    lifts = []
    for D in info['reps']:
        N = L.rational_e1_lift(F, A, D, nprimes=12, restrict=True)
        lifts.append(N)
    info['lifts'] = lifts
    info['lift_fail'] = sum(1 for N in lifts if N is None)
    if info['nontrivial'] > 0:
        try:
            reps_f, lifts_f, finfo = L.float_kernel_reps(F, A, d, parD, k)
            info['reps_float'] = reps_f; info['lifts_float'] = lifts_f; info['float_info'] = finfo
        except Exception as ex:
            info['float_error'] = repr(ex)
    info['seconds_total'] = round(time.time() - t, 1)
    return info

if __name__ == '__main__':
    degs = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else list(range(4, 13))
    procs = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    tasks = [(d, parD, k) for d in degs for parD in (0, 1) for k in range(7)]
    store = {}
    if os.path.exists(PKL):
        store = pickle.load(open(PKL, 'rb'))
    tasks = [t for t in tasks if str(t) not in store]
    with Pool(procs) as pool, open(OUT, 'a') as out:
        for info in pool.imap_unordered(work, tasks):
            key = str((info['d'], info['parD'], info['k']))
            store[key] = info
            pickle.dump(store, open(PKL + '.tmp', 'wb')); os.replace(PKL + '.tmp', PKL)
            line = (f"d={info['d']} parD={info['parD']} k={info['k']}: nD={info['nD']} (R/J)_(d+6)={info['ring_dim_target']} "
                    f"kernel={info['kernel_total']} J_d={info['J_d']} nontrivial={info['nontrivial']} reps={len(info['reps'])} "
                    f"heights={info.get('rep_heights')} lift_ok={info.get('lift_ok')} lift_fail={info['lift_fail']} [{info['seconds_total']}s]")
            print(line, flush=True); out.write(line + '\n'); out.flush()
            L.wlog({'family': 'sector-dims', 'd': info['d'], 'parD': info['parD'], 'k': info['k'], 'nD': info['nD'],
                    'ring_dim_target': info['ring_dim_target'], 'kernel_total': info['kernel_total'], 'J_d': info['J_d'],
                    'nontrivial': info['nontrivial'], 'reps': len(info['reps']), 'rep_heights': info.get('rep_heights'),
                    'lift_ok': info.get('lift_ok'), 'lift_fail': info['lift_fail'], 'seconds': info['seconds_total'], 'p': info['p'],
                    'note': 'E1-admissible sector space: kernel of [A] on (R/J)_d^{parD}[k] mod p=2147483647 (upper bound over Q), rational reps via CRT(12 primes)+LLL'})
    print('KERNELS_DONE', flush=True)
