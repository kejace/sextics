cd /tmp/hodge/results/eight-plane-configuration
cat > k10one.py <<'PY'
import sys, json, pickle
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from e1kernel import *
d, parD = int(sys.argv[1]), int(sys.argv[2])
info, basis = e1_kernel(d, parD, p=core.PRIMES[1]); info['p'] = core.PRIMES[1]; info['engine'] = 2
print(json.dumps(info), flush=True)
with open(f'{WD}/e1kernel_remote.jsonl', 'a') as f: f.write(json.dumps(info) + '\n')
if basis:
    with open(f'{WD}/e1kernel_basis_d{d}_par{parD}.pkl', 'wb') as f: pickle.dump(basis, f)
PY
nohup python3 k10one.py 10 0 > k10_par0.out 2>&1 &
nohup python3 k10one.py 10 1 > k10_par1.out 2>&1 &
nohup python3 k10one.py 9 0 > k9_par0.out 2>&1 &
nohup python3 k10one.py 9 1 > k9_par1.out 2>&1 &
sleep 1; ps aux | grep k10one | grep -v grep | wc -l
