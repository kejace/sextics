"""V1: re-verify the E1-kernel profile (both parities, d=1..9) at a NEW prime with a different code path
(hodge.jacobian.mult_kernel), check the J-block rank equals the regular-sequence value (needed for the char-0 lift),
and check K_9 = K_7*S_2 + K_8*S_1 (mod J) at that prime.  Output: verify_e1kernel_<p>.json + log lines."""
import sys, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000007
OUT = f'{PK}/verify_e1kernel_{p}.json'
LOGF = f'{PK}/verify_e1kernel_{p}.log.jsonl'
def fermat_hf(k, par):
    """dim (S/(z_i^5))_k with Y-parity par (the regular-sequence value for any smooth sextic's Jacobian ring)."""
    n = 0
    for e in monomials(k, par):
        if max(e) <= 4: n += 1
    return n
results = {}
kers = {}
for d in range(1, 10):
    for par in (0, 1):
        t0 = time.time()
        out = jacobian.mult_kernel(F, A, d, p, parity=par)
        kers[(d, par)] = out['basis']
        # J-block rank in degree d+6, parity (1+par)%2
        tpar = (1 + par) % 2
        Jm, mons = jacobian.jacobian_matrix(F, d + 6, p, tpar)
        rJ = Jm.rank()
        expJ = len(mons) - fermat_hf(d + 6, tpar)
        rec_ = dict(exp='reverify-e1-kernel', d=d, parity=par, prime=p, kernel_mod_J=out['dim'], kernel_dim=out['kernel_dim'],
                    J_d=out['jacobian_dim'], ring_d=out['ring_dim'], J_d_expected=len(monomials(d, par)) - fermat_hf(d, par),
                    rankJ_target=rJ, rankJ_target_expected=expJ, rankJ_ok=(rJ == expJ), seconds=round(time.time() - t0, 1),
                    note='mult_kernel (hodge.jacobian) at a new prime; rankJ_ok=True means the J block has the regular-sequence rank so dim_p kernel >= dim_Q kernel')
        print(json.dumps(rec_), flush=True)
        results[f'{d},{par}'] = rec_
        log(rec_, path=LOGF)
        json.dump(results, open(OUT, 'w'), indent=1)
# generation check: K_9 (both parities) vs products K_7 * S_2 + K_8 * S_1 modulo J_9
t0 = time.time()
Fp = core.reduce_mod(F, p)
for par in (0, 1):
    prods = []
    for b in kers[(7, 0)]:
        for m in monomials(2, par):
            prods.append(mul(b, {m: 1}, p))
    for b in kers[(8, 0)]:
        for m in monomials(1, par):
            prods.append(mul(b, {m: 1}, p))
    for b in kers[(8, 1)]:
        for m in monomials(1, (par + 1) % 2):
            prods.append(mul(b, {m: 1}, p))
    red = jacobian._Reducer(F, 9, p, par)
    nf = [red.normal_form(x) for x in prods]
    sm = red.standard; sidx = {m: i for i, m in enumerate(sm)}
    flat = [0] * (len(nf) * len(sm))
    for r_, v in enumerate(nf):
        for e, c in v.items(): flat[r_ * len(sm) + sidx[e]] = c
    rk = nmod_mat(len(nf), len(sm), flat, p).rank() if nf else 0
    rec_ = dict(exp='reverify-K9-generation', d=9, parity=par, prime=p, span_products_mod_J=rk, K9_mod_J=len(kers[(9, par)]),
                generated=(rk == len(kers[(9, par)])), seconds=round(time.time() - t0, 1),
                note='rank mod J_9 of K_7*S_2 + K_8*S_1 versus dim K_9 (dossier 4.6: no new (J:A) generators in degree 9)')
    print(json.dumps(rec_), flush=True); log(rec_, path=LOGF)
    results[f'gen9,{par}'] = rec_
json.dump(results, open(OUT, 'w'), indent=1)
print('done')
