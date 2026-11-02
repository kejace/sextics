"""Hessian-relaxation family exclusion: python3 hess_run.py <d> <prime> [--cert] [--nonj FILE] [--tag TAG] [--blocks 0,1]
Per parity block: E' = [J_{2d} | P_kl+P_lk]; coordinates u on V_d/(F*S_{d-6}); checks u^2 in rowspace mod p and (with --cert)
produces exact rational certificates lambda^T E' = u^2 (lambda vanishes on J_{2d}, hence on all Koszul/Q columns)."""
import sys, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
ap = argparse.ArgumentParser()
ap.add_argument('d', type=int); ap.add_argument('prime', type=int)
ap.add_argument('--cert', action='store_true'); ap.add_argument('--tag', default=''); ap.add_argument('--nonj', default='')
ap.add_argument('--blocks', default='0,1'); ap.add_argument('--parity', default='all')
a = ap.parse_args()
d, p = a.d, a.prime
par = None if a.parity == 'all' else int(a.parity)
tag = a.tag or f'hess_d{d}_p{p}_{a.parity}'
LOGF = f'{PK}/{tag}.log.jsonl'
t0 = time.time()
basis = polar_basis(d, par)
nonj = []
if a.nonj:
    for k, rec_ in enumerate(json.load(open(a.nonj))):
        D = core.parse(rec_['D']); N = [core.parse(s) for s in rec_['N']]
        if par is not None and core.parity_of(D) != par: continue
        basis.append((('K', rec_.get('label', k)), D, N)); nonj.append(k)
K = len(basis); labidx = {b[0]: k for k, b in enumerate(basis)}
print(f'd={d} p={p}: basis {K} (polar {K-len(nonj)}, non-J {len(nonj)})', flush=True)
HF = HessianFamily(d, basis)
sp_groups = []
if d >= 6:
    for m in monomials(d - 6):
        grp = [labidx[('J', i, tuple(m[j] + (1 if j == i else 0) for j in range(NV)))] for i in range(NV) if ('J', i, tuple(m[j] + (1 if j == i else 0) for j in range(NV))) in labidx]
        if len(grp) == NV: sp_groups.append(grp)
in_group = {k for grp in sp_groups for k in grp}
coords = [(f'a[{basis[k][0]}]', {k: Fraction(1)}) for k in range(K) if k not in in_group]
for grp in sp_groups:
    for i in range(1, NV):
        coords.append((f'a[{basis[grp[0]][0]}]-a[{basis[grp[i]][0]}]', {grp[0]: Fraction(1), grp[i]: Fraction(-1)}))
def square_vec(lin, blk):
    out = {}; ks = sorted(lin); cidx = HF.blocks[blk]['colidx']
    for x in range(len(ks)):
        for y in range(x, len(ks)):
            k, l = ks[x], ks[y]
            j = cidx[('aa', k, l)]
            out[j] = out.get(j, 0) + lin[k] * lin[l] * (1 if k == l else 2)
    return out
summary = dict(exp='hessian-relaxation-weight2', d=d, prime=p, K=K, nonJ=len(nonj), n_spurious_dirs=len(sp_groups), n_coords=len(coords), blocks={})
allok = True
for blk in [int(x) for x in a.blocks.split(',')]:
    info, Mr, piv = HF.analyze(blk, p)
    C = HF.blocks[blk]['ncols']
    tg = [(name, lin) for name, lin in coords if (1 + 2 * HF.par[next(iter(lin))]) % 2 == blk]
    okm = [reduce_by_rref(Mr, piv, square_vec(lin, blk), p, C) for name, lin in tg]
    n_ok = sum(okm)
    binfo = dict(info, n_coords=len(tg), coord_squares_in_rowspace_modp=n_ok, failing=[name for (name, _), ok in zip(tg, okm) if not ok][:30])
    print('[block]', json.dumps(binfo), flush=True)
    summary['blocks'][blk] = binfo
    log(dict(exp='hessian-relaxation-weight2-block', d=d, prime=p, **binfo, job=tag, seconds=round(time.time() - t0)), path=LOGF)
    if n_ok < len(tg): allok = False
    if a.cert and n_ok:
        targets = [(name, square_vec(lin, blk)) for (name, lin), ok in zip(tg, okm) if ok]
        certs = hess_certificates(HF, blk, Mr, piv, targets, p)
        nok = sum(1 for v in (certs or {}).values() if v)
        if certs:
            json.dump({name: ({str(r): str(v) for r, v in lam.items()} if lam else None) for name, lam in certs.items()}, open(f'{PK}/{tag}.blk{blk}.certs.json', 'w'))
        binfo['certified'] = nok
        log(dict(exp='hessian-relaxation-exact-certificates', d=d, prime=p, blk=blk, certified=nok, of=len(tg), job=tag, seconds=round(time.time() - t0),
                 result=('all coordinates of this block certified over Q' if nok == len(tg) else 'incomplete')), path=LOGF)
        if nok < len(tg): allok = False
    del Mr
summary['result'] = 'PROVED (mod p rowspace; exact if --cert and all certified): no certificate of degree %d' % d if allok else 'incomplete: some coordinate squares not in the J-relaxed weight-2 ideal'
summary['seconds'] = round(time.time() - t0)
print(json.dumps(summary, default=str), flush=True)
json.dump(dict(summary=summary, coords=[(n, {str(k): str(v) for k, v in lin.items()}) for n, lin in coords], basis_labels=[str(b[0]) for b in basis]), open(f'{PK}/{tag}.json', 'w'), default=str)
log(dict(summary, job=tag), path=LOGF)
print('done', round(time.time() - t0), 's')
