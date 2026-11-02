import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
LOG='/tmp/hodge/results/polarized-k3-remainder/wave4/log.jsonl'
def log(**kw):
    kw['ts']=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()); open(LOG,'a').write(json.dumps(kw, default=str)+'\n')
F, A, q, rec = core.load_instance('polarized-k3-remainder')
Amix=core.parse(json.load(open('/tmp/hodge/results/polarized-k3-remainder/wave4/amix_pieces.json'))['A_mix'])
p=core.PRIMES[0]
for d in [8,9]:
    for par in [0,1]:
        t0=time.time()
        km=jacobian.mult_kernel(F,Amix,d,p,parity=par)
        row={'d':d,'parity':par,'Amix_kernel_in_RJ':km['dim'],'Amix_admissible':km['kernel_dim'],'J_d':km['jacobian_dim'],'RJ_d':km['ring_dim'],'sec':round(time.time()-t0,1)}
        print(row, flush=True); log(what='e1-mult-kernel-Amix', prime=p, d=d, parity=par, result=row)
