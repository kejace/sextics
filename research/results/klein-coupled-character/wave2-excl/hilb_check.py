"""Regular-sequence check mod p: dim_Fp J(F)_k per (parity, sigma-character) block summed must equal
dim S_k - coeff_k ((1-t^5)/(1-t))^6 for every k <= KMAX.  Also records the per-block dims (used later)."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
import kx
from hodge import core
p = int(sys.argv[1]); KMAX = int(sys.argv[2])
F, A, q = kx.load()
out = {'p': p, 'blocks': {}, 'ok': True}
t0 = time.time()
for k in range(5, KMAX + 1):
    tot = 0
    for par in (0, 1):
        for ch in range(7):
            keys, polys = kx.jac_block(F, k, p, par, ch)
            mons = kx.wmonomials(k, par, ch)
            rows = {m: i for i, m in enumerate(mons)}
            rk = kx.dense(polys, rows, p).rank() if polys else 0
            out['blocks'][f'{k},{par},{ch}'] = {'dimS': len(mons), 'dimJ': rk, 'dimR': len(mons) - rk}
            tot += rk
    exp = kx.dimS(k) - kx.hilbert_ci(k)
    ok = (tot == exp)
    out['ok'] = out['ok'] and ok
    print(k, 'dimJ_p', tot, 'expected', exp, 'OK' if ok else 'MISMATCH', round(time.time() - t0, 1), flush=True)
json.dump(out, open(f'hilbert_blocks_p{p}.json', 'w'))
kx.log(exp='hilbert-regular-sequence', p=p, kmax=KMAX, ok=out['ok'], seconds=round(time.time() - t0, 1),
       note='dim_Fp J_k (sum over 14 parity/character blocks) equals the complete-intersection Hilbert function for 5<=k<=KMAX')
print('ALL OK' if out['ok'] else 'MISMATCH FOUND')
