import sys, json, pickle
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from e1kernel import *
d, parD = int(sys.argv[1]), int(sys.argv[2])
info, basis = e1_kernel(d, parD, p=core.PRIMES[1]); info['p'] = core.PRIMES[1]; info['engine'] = 2
print(json.dumps(info), flush=True)
with open(f'{WD}/e1kernel_remote.jsonl', 'a') as f: f.write(json.dumps(info) + '\n')
if basis:
    with open(f'{WD}/e1kernel_basis_d{d}_par{parD}.pkl', 'wb') as f: pickle.dump(basis, f)
