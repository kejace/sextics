"""D8 exclusion re-verification (lighter than d8full.py): (ii) sector system, modular ranks at two primes + exact rational
null vectors (CRT/ratrecon) or exact fmpz ranks; (iii) exact fmpq rref consistency of the augmented sector system."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from w2lib import *
import sect
S = Setup(); F, A = S.F, S.A
D8 = core.parse(json.load(open(WD + '/kernel_d8_0_4_0.json'))['D'])
t0 = time.time()
out = sealed_exclusion(F, A, 1, D8, exact_fallback=True)
print('(ii)', out, flush=True)
log({'exp': 'D8-sector-sealed-reverify', 'D': 'D8', 'd': 8, 'prime': ','.join(map(str, out['primes'])),
     'result': ('SEALED INCONSISTENT over Q' if out['sealed_inconsistent'] else 'not sealed') + f" (method {out.get('method')})", 'rank': out['rank'][0],
     'seconds': out['seconds'], 'job': 'excl-ctr-chainA', 'note': json.dumps({k: v for k, v in out.items() if k not in ('primes',)})[:600]})
t1 = time.time()
Sy = sect.Sys(F, A, 1, D8)
Maug = fmpq_mat(Sy.dense(None, aug=True)); Mr, rk = Maug.rref()
piv = []
for i in range(rk):
    j = 0
    while j <= Sy.ncols and Mr[i, j] == 0: j += 1
    piv.append(j)
cons = not (piv and piv[-1] == Sy.ncols)
print(f'(iii) exact fmpq rref of augmented sector system: rank_aug={rk} consistent={cons} ({time.time()-t1:.1f}s)', flush=True)
log({'exp': 'D8-sector-exact-fmpq-rref', 'D': 'D8', 'd': 8, 'prime': None, 'result': f'rank_aug over Q = {rk}, consistent={cons}', 'rank': rk,
     'seconds': round(time.time() - t1, 1), 'job': 'excl-ctr-chainA', 'note': f'rows={Sy.nrows} cols={Sy.ncols}; wave-1: fmpz rank 1155 / 1156'})
