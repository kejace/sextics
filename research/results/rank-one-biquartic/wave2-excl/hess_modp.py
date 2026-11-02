"""hess_modp: memory-light MODULAR Hessian-relaxation screening for larger d (8, 9): per prime and parity block,
rref of the J_{2d} block once, normal forms of the sparse quadratic columns (Hessian representatives) in the standard
monomial basis, row space of the resulting (standard monomials x aa-columns) matrix = quadratic forms in the J-relaxed
weight-2 ideal (mod p); then coordinate-square membership and the residual Singular radical test.  No exact step.
Usage: python3 hess_modp.py <d> <prime> [--nonj FILE] [--blocks 1] [--tag TAG]"""
import sys, json, time, argparse, subprocess
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
ap = argparse.ArgumentParser(); ap.add_argument('d', type=int); ap.add_argument('prime', type=int); ap.add_argument('--nonj', default=''); ap.add_argument('--blocks', default='1'); ap.add_argument('--tag', default='')
a = ap.parse_args(); d, p = a.d, a.prime
tag = a.tag or f'hmodp_d{d}_p{p}'
LOGF = f'{PK}/{tag}.log.jsonl'
T0 = time.time(); el = lambda: round(time.time() - T0)
basis = polar_basis(d)
if a.nonj:
    for k, rec_ in enumerate(json.load(open(a.nonj))):
        basis.append((('K', rec_.get('label', k)), core.parse(rec_['D']), [core.parse(s) for s in rec_['N']]))
K = len(basis); labidx = {b[0]: k for k, b in enumerate(basis)}
par = [core.parity_of(b[1]) for b in basis]
sp_groups = []
if d >= 6:
    for m in monomials(d - 6):
        sp_groups.append([labidx[('J', i, tuple(m[j] + (1 if j == i else 0) for j in range(NV)))] for i in range(NV)])
in_group = {k for grp in sp_groups for k in grp}
coords = [(f'a{k}', {k: 1}) for k in range(K) if k not in in_group]
for grp in sp_groups:
    for i in range(1, NV): coords.append((f'a{grp[0]}-a{grp[i]}', {grp[0]: 1, grp[i]: -1}))
centers = [grp[0] for grp in sp_groups]
print(f'd={d} p={p} K={K} coords={len(coords)}', flush=True)
results = {}
for blk in [int(x) for x in a.blocks.split(',')]:
    t1 = time.time()
    red = jacobian._Reducer(F, 2 * d, p, blk)
    sm = red.standard; sidx = {m: i for i, m in enumerate(sm)}
    print(f'[blk{blk}] J_{2*d} rref: rank {red.rank}, standard monomials {len(sm)} ({el()}s)', flush=True)
    pairs = [(k, l) for k in range(K) for l in range(k, K) if (1 + par[k] + par[l]) % 2 == blk]
    # normal forms of the sparse Hessian representatives -> matrix NF (len(sm) x len(pairs))
    flat = [0] * (len(sm) * len(pairs))
    npairs = len(pairs)
    for j, (k, l) in enumerate(pairs):
        col = hess_column(basis[k], basis[l])
        nf = red.normal_form(col)
        for e, v in nf.items(): flat[sidx[e] * npairs + j] = v
        if j % 5000 == 0: print(f'   nf {j}/{npairs} ({el()}s)', flush=True)
    NF = nmod_mat(len(sm), npairs, flat, p)
    del flat
    Mr, rank = NF.rref(); piv = rref_pivots(Mr, rank, npairs)
    pidx = {pr: j for j, pr in enumerate(pairs)}
    print(f'[blk{blk}] forms space dim {rank} (of {len(sm)} functionals, {npairs} columns) ({el()}s)', flush=True)
    rows_ = []
    for i, j in enumerate(piv):
        rows_.append({pairs[k_]: int(Mr[i, k_]) for k_ in range(j, npairs) if int(Mr[i, k_])})
    sq_ok = {}
    if blk == 1:
        for name, lin in coords:
            tv = {}
            ks = sorted(lin)
            for x in range(len(ks)):
                for y in range(x, len(ks)):
                    k, l = ks[x], ks[y]; tv[pidx[(k, l)]] = tv.get(pidx[(k, l)], 0) + lin[k] * lin[l] * (1 if k == l else 2)
            sq_ok[name] = reduce_by_rref(Mr, piv, {j: Fraction(v) for j, v in tv.items()}, p, npairs)
        cert = [n for n, _ in coords if sq_ok[n]]; fail = [n for n, _ in coords if not sq_ok[n]]
        print(f'[blk{blk}] coordinate squares in forms space: {len(cert)} of {len(coords)}; failing {fail[:15]}', flush=True)
        results['cert'] = cert; results['fail'] = fail
    results[f'forms{blk}'] = rows_
    results[f'rank{blk}'] = rank
    log(dict(exp='hmodp-forms', d=d, prime=p, blk=blk, K=K, n_forms=rank, n_functionals=len(sm), n_pairs=npairs, squares_in=len(results.get('cert', [])), n_coords=len(coords),
             failing=results.get('fail', [])[:40], job=tag, seconds=el()), path=LOGF)
    del Mr, NF
# residual Singular test (mod p) using all forms from the computed blocks
fail = results.get('fail', [n for n, _ in coords]); cert = results.get('cert', [])
res_vars = list(fail) + [f'a{g}' for g in centers]; rv_idx = {n: i for i, n in enumerate(res_vars)}
subst = {}
for k in range(K):
    if k in in_group:
        grp = next(gr for gr in sp_groups if k in gr); g0 = grp[0]
        expr = {rv_idx[f'a{g0}']: 1}
        if k != g0 and f'a{g0}-a{k}' in rv_idx: expr[rv_idx[f'a{g0}-a{k}']] = -1
        subst[k] = expr
    else:
        subst[k] = {rv_idx[f'a{k}']: 1} if f'a{k}' in rv_idx else {}
lines = [f'ring r = {p}, (' + ','.join(f'v{i}' for i in range(len(res_vars))) + '), dp;', 'ideal I;']
nnz = 0
for blk in [int(x) for x in a.blocks.split(',')]:
    for fm in results[f'forms{blk}']:
        out = {}
        for (k, l), v in fm.items():
            for ik, ck in subst[k].items():
                for il, cl in subst[l].items():
                    key = (min(ik, il), max(ik, il)); out[key] = (out.get(key, 0) + v * ck * cl) % p
        out = {kk: vv for kk, vv in out.items() if vv}
        if not out: continue
        nnz += 1
        lines.append('I[size(I)+1] = ' + '+'.join((f'{v}*v{i}*v{j}' if i != j else f'{v}*v{i}^2') for (i, j), v in out.items()) + ';')
lines += ['ideal G = std(I);', '"dim V(I):", dim(G), " #forms:", size(I);', 'int bad = 0; ideal T;']
for nm in fail:
    lines.append(f'T = G, v{rv_idx[nm]} - 1; if (dim(std(T)) >= 0) {{ bad++; "residual coordinate {nm} NOT in radical"; }}')
lines += [f'"residual coordinates not in radical:", bad, "of", {len(fail)};', 'quit;']
open(f'{PK}/{tag}.sing', 'w').write('\n'.join(lines) + '\n')
r = subprocess.run(['Singular', '-q', f'{PK}/{tag}.sing'], capture_output=True, text=True, timeout=7200)
open(f'{PK}/{tag}.sing.out', 'w').write(r.stdout + r.stderr)
bad = next((int(l.split(':')[1].split('of')[0]) for l in r.stdout.splitlines() if l.startswith('residual coordinates not in radical:')), None)
dimline = next((l for l in r.stdout.splitlines() if l.startswith('dim V(I):')), '')
result = (f'SCREENED mod {p}: V(J-relaxed forms) subset Sp => no certificate of degree {d} modulo p' if bad == 0 else f'inconclusive mod p: {bad} residual coordinates not in radical')
print('[residual]', dimline, '|', result, f'({el()}s)', flush=True)
log(dict(exp='hmodp-residual-groebner', d=d, prime=p, n_residual_vars=len(res_vars), n_restricted_forms=nnz, result=result, singular=dimline, job=tag, seconds=el()), path=LOGF)
json.dump(dict(cert=cert, fail=fail, res_vars=res_vars, ranks={b: results.get(f'rank{b}') for b in (0, 1)}), open(f'{PK}/{tag}.json', 'w'))
print('done', el(), 's')
