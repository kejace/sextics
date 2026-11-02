import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from e1kernel import *
jobs = [(8,1),(9,0),(9,1)]
for d, parD in jobs:
    t0=time.time()
    info, basis = e1_kernel(d, parD)
    info['engine']=2
    print(json.dumps(info), flush=True)
    with open(f'{WD}/e1kernel.jsonl','a') as f: f.write(json.dumps(info)+'\n')
    if basis:
        with open(f'{WD}/e1kernel_basis_d{d}_par{parD}.txt','w') as f:
            for b in basis:
                f.write(core.to_str({k:(v if v<2**30 else v-core.PRIMES[0]) for k,v in b.items()})+'\n')
