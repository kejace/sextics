"""Sector dims (1 prime) + float kernel reps/lifts into w2_sectors.pkl: python3 w2_sector_float.py d parD k"""
import os, sys, time, pickle
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import w2lib as L
F, A, q = L.load()
d, parD, k = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
t = time.time()
info = L.sector_space(F, A, d, parD, k, nprimes=1, verbose=False, want_rational=False)
print('dims', {kk: info[kk] for kk in ('nD', 'ring_dim_target', 'kernel_total', 'J_d', 'nontrivial')}, f'[{time.time()-t:.1f}s]', flush=True)
info.update({'reps': [], 'lifts': [], 'lift_fail': 0, 'lift_ok': False})
if info['nontrivial'] > 0:
    reps_f, lifts_f, finfo = L.float_kernel_reps(F, A, d, parD, k, verbose=True)
    info.update({'reps_float': reps_f, 'lifts_float': lifts_f, 'float_info': finfo})
info['seconds_total'] = round(time.time() - t, 1)
PKL = L.LANE + '/w2_sectors.pkl'
store = pickle.load(open(PKL, 'rb')) if os.path.exists(PKL) else {}
store[str((d, parD, k))] = info
pickle.dump(store, open(PKL + '.tmp', 'wb')); os.replace(PKL + '.tmp', PKL)
L.wlog({'family': 'sector-dims', 'd': d, 'parD': parD, 'k': k, 'nD': info['nD'], 'ring_dim_target': info['ring_dim_target'], 'kernel_total': info['kernel_total'],
        'J_d': info['J_d'], 'nontrivial': info['nontrivial'], 'reps': 0, 'rep_heights': None, 'lift_ok': False, 'lift_fail': 0, 'seconds': info['seconds_total'], 'p': info['p'],
        'note': 'E1-admissible sector space (local, 1 prime): dims mod p; floating-point orthonormal kernel reps stored'})
print('stored', f'[{time.time()-t:.1f}s]', flush=True)
