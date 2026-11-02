"""Compute one sector space locally (light) into w2_sectors.pkl: python3 w2_sector_local.py d parD k [nprimes]"""
import os, sys, time, pickle
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import w2lib as L
F, A, q = L.load()
d, parD, k = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
nprimes = int(sys.argv[4]) if len(sys.argv) > 4 else 8
PKL = L.LANE + '/w2_sectors.pkl'
t = time.time()
info = L.sector_space(F, A, d, parD, k, nprimes=nprimes, verbose=False, want_rational=True)
lifts = [L.rational_e1_lift(F, A, D, nprimes=nprimes, restrict=True) for D in info['reps']]
info['lifts'] = lifts; info['lift_fail'] = sum(1 for N in lifts if N is None)
if info['nontrivial'] > 0:
    reps_f, lifts_f, finfo = L.float_kernel_reps(F, A, d, parD, k)
    info['reps_float'] = reps_f; info['lifts_float'] = lifts_f; info['float_info'] = finfo
info['seconds_total'] = round(time.time() - t, 1)
store = pickle.load(open(PKL, 'rb')) if os.path.exists(PKL) else {}
store[str((d, parD, k))] = info
pickle.dump(store, open(PKL + '.tmp', 'wb')); os.replace(PKL + '.tmp', PKL)
line = (f"d={d} parD={parD} k={k}: nD={info['nD']} (R/J)_(d+6)={info['ring_dim_target']} kernel={info['kernel_total']} J_d={info['J_d']} "
        f"nontrivial={info['nontrivial']} reps={len(info['reps'])} heights={info.get('rep_heights')} lift_ok={info.get('lift_ok')} lift_fail={info['lift_fail']} [{info['seconds_total']}s]")
print(line, flush=True)
L.wlog({'family': 'sector-dims', 'd': d, 'parD': parD, 'k': k, 'nD': info['nD'], 'ring_dim_target': info['ring_dim_target'],
        'kernel_total': info['kernel_total'], 'J_d': info['J_d'], 'nontrivial': info['nontrivial'], 'reps': len(info['reps']),
        'rep_heights': info.get('rep_heights'), 'lift_ok': info.get('lift_ok'), 'lift_fail': info['lift_fail'], 'seconds': info['seconds_total'],
        'p': info['p'], 'nprimes': nprimes, 'note': 'E1-admissible sector space (local run): kernel of [A] on (R/J)_d^{parD}[k] mod p (upper bound over Q), rational reps via CRT+LLL'})
