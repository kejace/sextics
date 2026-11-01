import json, sys, time
from sympy import primerange
from ffcount import twisted_sums, f_B, g_E
BAD = {2, 3, 7, 11, 2957}
pmin, pmax = int(sys.argv[1]), int(sys.argv[2])
out = open(f'sums_{pmin}_{pmax}.jsonl', 'w')
for p in primerange(pmin, pmax + 1):
    if p in BAD: continue
    t0 = time.time()
    rec = {'p': p}
    ks = [1, 2] if p % 3 == 1 else [2]
    if p % 3 == 1 and p <= 110: ks.append(3)
    if p % 3 == 2 and p <= 30: ks.append(4)
    for k in ks:
        rb = twisted_sums(p, k, f_B); re_ = twisted_sums(p, k, g_E)
        rec[f'TB_k{k}'] = {kk: rb[kk] for kk in ('S', 'npts', 'n0', 'omega')}
        rec[f'E_k{k}'] = {kk: re_[kk] for kk in ('S', 'npts', 'n0', 'omega')}
    rec['seconds'] = round(time.time() - t0, 1)
    out.write(json.dumps(rec) + '\n'); out.flush()
    print(p, rec['seconds'], flush=True)
