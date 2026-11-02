import sys, time
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl-referee')
from reflib import *
LOG = '/tmp/hodge/results/one-reflection-skew/wave2-excl-referee/referee.log.jsonl'
t0 = time.time()
for d, p in ((8, None), (9, 1000000007)):
    K = 0; J = 0; per = {}
    for ch in all_sectors():
        k, rJ, nr = e1_kernel_dim(d, ch, p=p); jd = j_dim(d, ch, p=p); K += k; J += jd
        if k: per[str(ch)] = (k, jd)
    print(f"6var d={d} ({'exact' if p is None else 'mod '+str(p)}): dim E1-kernel = {K}, dim J_d = {J}, extra = {K-J}  ({time.time()-t0:.1f}s)", flush=True)
    log(LOG, exp='e1-kernel-6var', d=d, prime=p, dimK=K, dimJ=J, extra=K-J, per_sector=per, method=('fmpz_mat.rank exact' if p is None else 'nmod_mat rank (kernel dim is an upper bound... see note)'), note='mod-p kernel dimension >= Q-kernel dimension in general; equality with the exact table is consistency evidence')
for d in (8, 9):
    K = 0; J = 0
    for ch in all_sectors():
        k, rJ, nr = e1_kernel_dim(d, ch, nv=5, p=1000000007); K += k; J += j_dim(d, ch, nv=5, p=1000000007)
    print(f"5var d={d} (mod 1000000007): dim (J_f:A)_d = {K}, dim (J_f)_d = {J}, extra = {K-J}  ({time.time()-t0:.1f}s)", flush=True)
    log(LOG, exp='e1-kernel-5var', d=d, prime=1000000007, dimK=K, dimJ=J, extra=K-J)
