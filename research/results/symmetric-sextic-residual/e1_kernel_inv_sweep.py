import sys, time, json, pickle
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from invsolve import *
F, A, q = load()
dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 16
for d in range(1, dmax + 1):
    for pD in (0, 1):
        t0 = time.time()
        info, polys = e1_kernel_inv(F, A, q, d, pD)
        info['seconds'] = round(time.time() - t0, 1)
        print(json.dumps(info), flush=True)
        pickle.dump((info, polys), open(f'{RES}/e1inv_d{d}_p{pD}.pkl', 'wb'))
        for Pk in polys[:12]:
            s = core.to_str(Pk)
            print('   inv ker elt (', len(Pk), 'terms):', s[:300], flush=True)
