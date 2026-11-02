"""E1-kernel ladder at degree d, parity par (mod PRIMES[1]); saves kernel basis pickle like k10one.py."""
import sys, json, pickle, time
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from e1kernel import *
d, parD = int(sys.argv[1]), int(sys.argv[2])
t0 = time.time()
info, basis = e1_kernel(d, parD, p=core.PRIMES[1]); info['p'] = core.PRIMES[1]; info['engine'] = 2; info['wave'] = 2
print(json.dumps(info), flush=True)
with open(f'{WD}/e1kernel_remote.jsonl', 'a') as f: f.write(json.dumps(info) + '\n')
with open(f'{WD}/wave2-search.log.jsonl', 'a') as f:
    f.write(json.dumps({'batch': 'e1-ladder', 'd': d, 'parity': parD, 'p': core.PRIMES[1], 'kerdim': info['kerdim'], 'dimJ': info['dimJ_d'], 'kermodJ': info['kermodJ'], 'rows': info['rows'], 'cols': info['cols'], 'seconds': round(time.time()-t0,1), 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}) + '\n')
if basis:
    with open(f'{WD}/e1kernel_basis_d{d}_par{parD}.pkl', 'wb') as f: pickle.dump(basis, f)
