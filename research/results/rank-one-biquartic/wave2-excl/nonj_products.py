"""nonj_products: build the full non-Jacobian E1-kernel generator list for degree d from exact lower-degree generators:
products (m * D_k, m * N_k) for monomials m of degree d - deg D_k (exact: (m N_k)(F) = m A D_k), plus new generators
(nonj_d<d>.json from k7exact.py), then prune to a set independent modulo J_d (mod p) so that the family basis is a
Q-basis of V_d (dimension check against the modular upper bound).
Usage: python3 nonj_products.py <d> <out.json> --in nonj_d7.json [--in nonj_d8.json ...] [--prime P]"""
import sys, json, argparse
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
ap = argparse.ArgumentParser(); ap.add_argument('d', type=int); ap.add_argument('out'); ap.add_argument('--in', dest='inp', action='append', default=[]); ap.add_argument('--prime', type=int, default=1000000007)
a = ap.parse_args(); d, p = a.d, a.prime
cands = []
for fn in a.inp:
    for rec_ in json.load(open(fn)):
        D = core.parse(rec_['D']); N = [core.parse(s) for s in rec_['N']]
        dd = core.degree(D)
        if dd > d: continue
        for m in monomials(d - dd):
            cands.append((f"{rec_.get('label','K')}*{core.to_str(mono(m))}", mul(D, mono(m)), [mul(Nk, mono(m)) for Nk in N]))
# prune modulo J_d (per parity) mod p, keeping a maximal independent set
out = []
for par in (0, 1):
    red = jacobian._Reducer(F, d, p, par)
    sm = red.standard; sidx = {m: i for i, m in enumerate(sm)}
    cur = []; cur_rank = 0
    for lab, D, N in cands:
        if core.parity_of(D) != par: continue
        nf = red.normal_form(core.reduce_mod(D, p))
        if not nf: continue
        row = [0] * len(sm)
        for e, v in nf.items(): row[sidx[e]] = v
        rk = nmod_mat(len(cur) + 1, len(sm), [x for r_ in cur for x in r_] + row, p).rank()
        if rk > cur_rank:
            cur.append(row); cur_rank = rk
            assert sub(apply_field(N, F), mul(A, D)) == {}
            out.append(dict(label=lab, D=core.to_str(D), N=[core.to_str(Nk) for Nk in N]))
    print(f'parity {par}: {cur_rank} independent non-J directions mod J_{d}', flush=True)
json.dump(out, open(a.out, 'w'))
log(dict(exp='nonJ-products', d=d, prime=p, n=len(out), inputs=a.inp, file=a.out, note='exact non-J E1-kernel directions of degree d (products of lower-degree exact generators + new ones), pruned to a J-independent set mod p'), path=f'{PK}/nonj_products_d{d}.log.jsonl')
print('wrote', a.out, len(out))
