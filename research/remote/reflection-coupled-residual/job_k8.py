"""K8 (unique non-Jacobian E1 class at d=8, Y-even, Klein-invariant): full fixed-D certificate test.
1. sector-restricted compressed solve (symmetry.solve_fixed_D_symmetric), 2. unrestricted compressed solve as cross-check."""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d, symmetry, jacobian
LOG = 'wave2-search.log.jsonl'
def log(rec):
    rec['lane'] = 'reflection-coupled-residual'; rec['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S'); rec['engine'] = fixed_d.ENGINE_VERSION
    open(LOG, 'a').write(json.dumps(rec, default=str) + '\n')
F, A, q, rec = core.load_instance('reflection-coupled-residual')
K8 = core.parse(json.load(open('kernel_d8_X03Y03.json'))[0])
t0 = time.time()
ok, info = fixed_d.e1_solvable(F, A, q, K8)
print('e1', ok, info, round(time.time() - t0, 1), flush=True)
log({'event': 'k8-e1', 'e1': bool(ok), 'info': info, 'seconds': round(time.time() - t0, 1)})
t0 = time.time()
status, data = symmetry.solve_fixed_D_symmetric(F, A, q, K8, verbose=True, compressed=True)
print('symmetric:', status, {k: v for k, v in data.items() if k not in ('N', 'Q')}, round(time.time() - t0, 1), flush=True)
log({'event': 'k8-fixedD-symmetric', 'status': status, 'seconds': round(time.time() - t0, 1),
     **{k: v for k, v in data.items() if k not in ('N', 'Q')}})
if status == 'solved':
    ok2, info2 = core.verify_witness(F, A, q, K8, data['N'], data['Q'])
    w = {'instance': 'reflection-coupled-residual', 'd': 8, 'D': core.to_str(K8), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}
    json.dump(w, open('witness-candidate.json', 'w'), indent=1)
    print('!!!!!!!! WITNESS FOUND (K8) verified =', ok2, flush=True)
    log({'event': 'WITNESS', 'D': 'K8', 'verified': bool(ok2)})
t0 = time.time()
out = fixed_d.solve_fixed_D_fast(F, A, q, K8, verbose=True, want_exact=(status == 'solved'))
print('full compressed:', {k: v for k, v in out.items() if k not in ('N', 'Q')}, round(time.time() - t0, 1), flush=True)
log({'event': 'k8-fixedD-full', 'seconds': round(time.time() - t0, 1), **{k: v for k, v in out.items() if k not in ('N', 'Q')}})
