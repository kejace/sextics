import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from e1kernel import *
import pickle
import pickle
for d, parD in [(10,0),(10,1)]:
    info, basis = e1_kernel(d, parD, p=core.PRIMES[1]); info["p"]=core.PRIMES[1]
    info['engine']=2
    print(json.dumps(info), flush=True)
    with open(f'{WD}/e1kernel.jsonl','a') as f: f.write(json.dumps(info)+'\n')
    if basis:
        with open(f'{WD}/e1kernel_basis_d{d}_par{parD}.pkl','wb') as f:
            pickle.dump(basis, f)
