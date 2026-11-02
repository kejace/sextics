"""hess_pipeline: complete Hessian-relaxation exclusion for degree d (all parities, mixed allowed).
  A (mod p_0..p_{n-1}): per parity block, rref of E' = [J_{2d} | quadratic columns]; forms = rref rows with aa-pivots.
     Block 1 at p_0: which coordinate squares u^2 (u on V_d/Sp) lie in the row space.
  B (mod p_0): residual variables = failing coordinates + spurious centers; restrict all forms; Singular radical test;
     greedy selection of full rows whose restrictions span the restricted span.
  C (exact): CRT-lift selected rows; exact certificates for coordinate squares and selected forms (lambda^T E' = (0|phi)).
  D (exact): Singular over Q on (certified linear coordinates) + (exact forms): dim(std(I + (u-1))) = -1 for all residual u.
Usage: python3 hess_pipeline.py <d> <nprimes> [--nonj FILE] [--sparse] [--tag TAG] [--skip-exact]"""
import sys, json, time, argparse, subprocess, random
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
from sympy import prevprime
from collections import Counter
ap = argparse.ArgumentParser(); ap.add_argument('d', type=int); ap.add_argument('nprimes', type=int); ap.add_argument('--nonj', default=''); ap.add_argument('--sparse', action='store_true')
ap.add_argument('--tag', default=''); ap.add_argument('--skip-exact', action='store_true'); ap.add_argument('--p0', type=int, default=1000000007)
a = ap.parse_args(); d = a.d
tag = a.tag or f'hpipe_d{d}'
LOGF = f'{PK}/{tag}.log.jsonl'
T0 = time.time()
def el(): return round(time.time() - T0)
basis = polar_basis(d)
if a.nonj:
    for k, rec_ in enumerate(json.load(open(a.nonj))):
        basis.append((('K', rec_.get('label', k)), core.parse(rec_['D']), [core.parse(s) for s in rec_['N']]))
K = len(basis); labidx = {b[0]: k for k, b in enumerate(basis)}
HF = HessianFamily(d, basis, sparse=a.sparse)
# coordinates
sp_groups = []
if d >= 6:
    for m in monomials(d - 6):
        grp = [labidx[('J', i, tuple(m[j] + (1 if j == i else 0) for j in range(NV)))] for i in range(NV)]
        sp_groups.append(grp)
in_group = {k for grp in sp_groups for k in grp}
coords = [(f'a{k}', {k: Fraction(1)}) for k in range(K) if k not in in_group]
for grp in sp_groups:
    for i in range(1, NV): coords.append((f'a{grp[0]}-a{grp[i]}', {grp[0]: Fraction(1), grp[i]: Fraction(-1)}))
centers = [grp[0] for grp in sp_groups]
def square_vec(lin, blk):
    out = {}; ks = sorted(lin); cidx = HF.blocks[blk]['colidx']
    for x in range(len(ks)):
        for y in range(x, len(ks)):
            k, l = ks[x], ks[y]; j = cidx[('aa', k, l)]
            out[j] = out.get(j, 0) + lin[k] * lin[l] * (1 if k == l else 2)
    return out
primes = [a.p0]; x = a.p0
while len(primes) < a.nprimes:
    x = prevprime(x); primes.append(int(x))
# ---------------- Step A ----------------
per = {0: [], 1: []}
sq_ok = {}
for blk in (0, 1):
    B = HF.blocks[blk]; C = B['ncols']; nJ = HF.nJ[blk]
    for p in primes:
        info, Mr, piv = HF.analyze(blk, p, verbose=False)
        rows_ = []
        for i, j in enumerate(piv):
            if j < nJ: continue
            rows_.append({k: int(Mr[i, k]) for k in range(j, C) if int(Mr[i, k])})
        per[blk].append(dict(p=p, rank=info['rank'], rank_J=info['rank_J'], rows=rows_, pivs=tuple(min(r_) for r_ in rows_)))
        print(f'[A blk{blk} p={p}] rank {info["rank"]} rank_J {info["rank_J"]} forms {len(rows_)} ({el()}s)', flush=True)
        if blk == 1 and p == primes[0]:
            for name, lin in coords:
                sq_ok[name] = reduce_by_rref(Mr, piv, square_vec(lin, 1), p, C)
        del Mr
    log(dict(exp='hpipe-A-forms', d=d, blk=blk, primes=[x_['p'] for x_ in per[blk]], rank=per[blk][0]['rank'], rank_J=per[blk][0]['rank_J'], n_forms=len(per[blk][0]['rows']),
             nrows=B['nrows'], ncols=C, job=tag, seconds=el()), path=LOGF)
cert_names = [n for n, _ in coords if sq_ok[n]]; fail_names = [n for n, _ in coords if not sq_ok[n]]
print(f'[A] coordinate squares in row space mod p: {len(cert_names)} of {len(coords)}; failing: {fail_names[:20]}', flush=True)
log(dict(exp='hpipe-A-squares', d=d, prime=primes[0], K=K, n_coords=len(coords), squares_in_I2=len(cert_names), failing=fail_names, job=tag, seconds=el()), path=LOGF)
# ---------------- Step B (mod p0) ----------------
p0 = primes[0]
# residual variables: failing coordinates + centers.  Substitution: certified non-group a_k -> 0; certified difference a_g0 - a_gi -> a_gi = a_g0.
coord_lin = dict(coords)
subst = {}   # k -> linear combination over residual variable indices (dict var->coef)
res_vars = []  # list of (name, kind)
for name in fail_names:
    lin = coord_lin[name]
    res_vars.append(name)
for g in centers: res_vars.append(f'a{g}')
rv_idx = {n: i for i, n in enumerate(res_vars)}
# express each a_k in terms of residual variables (mod p): a_k = u_k if non-group (u_k = 0 if certified); group: a_gi = a_g0 - u_{g0,gi}
for k in range(K):
    if k in in_group:
        grp = next(gr for gr in sp_groups if k in gr); g0 = grp[0]
        expr = {rv_idx[f'a{g0}']: 1}
        if k != g0:
            nm = f'a{g0}-a{k}'
            if nm in rv_idx: expr[rv_idx[nm]] = expr.get(rv_idx[nm], 0) - 1
        subst[k] = expr
    else:
        nm = f'a{k}'
        subst[k] = {rv_idx[nm]: 1} if nm in rv_idx else {}
def restrict(form, blk, p):
    """quadratic form (dict colidx->val) in a -> dict (i<=j in residual vars) -> val mod p"""
    out = {}
    cols = HF.blocks[blk]['cols']
    for j, v in form.items():
        _, k, l = cols[j]
        ek, elc = subst[k], subst[l]
        for ik, ck in ek.items():
            for il, cl in elc.items():
                key = (min(ik, il), max(ik, il)); out[key] = (out.get(key, 0) + v * ck * cl) % p
    return {kk: vv for kk, vv in out.items() if vv}
restricted = {blk: [restrict(r_, blk, p0) for r_ in per[blk][0]['rows']] for blk in (0, 1)}
nres = len(res_vars)
# Singular mod p radical test on all restricted forms
def sing_script(forms_by_blk, char, resnames, test_names, exact=False):
    lines = [f'ring r = {char}, (' + ','.join(f'v{i}' for i in range(len(resnames))) + '), dp;', 'ideal I;']
    for blk, fl in forms_by_blk.items():
        for fm in fl:
            if not fm: continue
            terms = []
            for (i, j), v in fm.items():
                if exact:
                    v = Fraction(v); cs = f'({v.numerator}/{v.denominator})'
                else: cs = str(v)
                terms.append(f'{cs}*v{i}*v{j}' if i != j else f'{cs}*v{i}^2')
            lines.append('I[size(I)+1] = ' + '+'.join(terms) + ';')
    lines += ['int tt = timer;', 'ideal G = std(I);', '"dim V(I):", dim(G), " #forms:", size(I), " time:", timer - tt;', 'int bad = 0; ideal T;']
    for nm in test_names:
        lines.append(f'T = G, v{rv_idx[nm]} - 1; if (dim(std(T)) >= 0) {{ bad++; "residual coordinate {nm} NOT in radical"; }}')
    lines += [f'"residual coordinates not in radical:", bad, "of", {len(test_names)};', 'quit;']
    return '\n'.join(lines) + '\n'
open(f'{PK}/{tag}.B.sing', 'w').write(sing_script(restricted, p0, res_vars, fail_names))
r = subprocess.run(['Singular', '-q', f'{PK}/{tag}.B.sing'], capture_output=True, text=True, timeout=3600)
open(f'{PK}/{tag}.B.out', 'w').write(r.stdout + r.stderr)
badB = next((int(l.split(':')[1].split('of')[0]) for l in r.stdout.splitlines() if l.startswith('residual coordinates not in radical:')), None)
dimB = next((l for l in r.stdout.splitlines() if l.startswith('dim V(I):')), '')
print(f'[B] residual vars {nres} ({len(fail_names)} failing + {len(centers)} centers); Singular mod p: {dimB}; not in radical: {badB} ({el()}s)', flush=True)
log(dict(exp='hpipe-B-groebner-modp', d=d, prime=p0, n_residual_vars=nres, failing=len(fail_names), result=('V(forms) subset Sp mod p' if badB == 0 else f'{badB} residual coordinates not in radical'), singular=dimB, job=tag, seconds=el()), path=LOGF)
# greedy selection of full rows whose restrictions span the restricted span (per block)
selected = {}
for blk in (0, 1):
    keys = sorted({kk for fm in restricted[blk] for kk in fm}); kidx = {kk: i for i, kk in enumerate(keys)}
    chosen = []; cur = []; cur_rank = 0
    for i, fm in enumerate(restricted[blk]):
        if not fm: continue
        row = [0] * max(1, len(keys))
        for kk, vv in fm.items(): row[kidx[kk]] = vv
        rk = nmod_mat(len(cur) + 1, len(row), [x_ for r_ in cur for x_ in r_] + row, p0).rank()
        if rk > cur_rank: cur.append(row); cur_rank = rk; chosen.append(i)
    selected[blk] = chosen
    print(f'[B] blk{blk}: restricted span rank {cur_rank}, selected {len(chosen)} full rows', flush=True)
json.dump(dict(coords=[n for n, _ in coords], certified_modp=cert_names, failing=fail_names, res_vars=res_vars, selected={b: selected[b] for b in (0, 1)}), open(f'{PK}/{tag}.B.json', 'w'))
if a.skip_exact or badB != 0:
    print('stopping before exact step', flush=True); sys.exit(0)
# ---------------- Step C (exact) ----------------
def crt_rows(blk, idxs):
    good = [x_ for x_ in per[blk] if x_['pivs'] == per[blk][0]['pivs'] and x_['rank'] == per[blk][0]['rank']]
    out = {}
    for i in idxs:
        keys = set()
        for x_ in good: keys |= set(x_['rows'][i])
        form = {}; ok = True
        for k in sorted(keys):
            a_, m_ = good[0]['rows'][i].get(k, 0), good[0]['p']
            for x_ in good[1:]:
                a_, m_ = core.crt_pair(a_, m_, x_['rows'][i].get(k, 0), x_['p'])
            rr = core.ratrecon(a_, m_)
            if rr is None: ok = False; break
            if rr: form[k] = rr
        out[i] = form if ok else None
    return out, len(good)
exact_forms = {0: {}, 1: {}}; cert_ok = {}
for blk in (0, 1):
    lifted, ng = crt_rows(blk, selected[blk])
    nfail = sum(1 for v in lifted.values() if v is None)
    targets = [(f'form{blk}_{i}', fm) for i, fm in lifted.items() if fm is not None]
    if blk == 1:
        targets += [(f'sq:{name}', square_vec(coord_lin[name], 1)) for name in cert_names]
    print(f'[C blk{blk}] lifted {len(lifted)-nfail}/{len(lifted)} selected forms with {ng} primes; certifying {len(targets)} targets ({el()}s)', flush=True)
    if not targets: continue
    info, Mr, piv = HF.analyze(blk, p0, verbose=False)
    certs = hess_certificates(HF, blk, Mr, piv, targets, p0, verbose=False)
    del Mr
    nok = 0
    for name, tv in targets:
        ok = bool(certs and certs.get(name))
        nok += ok
        if name.startswith('sq:'): cert_ok[name[3:]] = ok
        elif ok: exact_forms[blk][name] = tv
    hts = max((max(abs(x_.numerator), x_.denominator) for lam in (certs or {}).values() if lam for x_ in lam.values()), default=1)
    print(f'[C blk{blk}] exact certificates verified {nok} of {len(targets)}; max certificate height digits {len(str(hts))} ({el()}s)', flush=True)
    if certs:
        json.dump({name: ({str(r_): str(v) for r_, v in lam.items()} if lam else None) for name, lam in certs.items()}, open(f'{PK}/{tag}.blk{blk}.certs.json', 'w'))
    log(dict(exp='hpipe-C-exact-certificates', d=d, blk=blk, prime=p0, targets=len(targets), verified=nok, max_height_digits=len(str(hts)), job=tag, seconds=el()), path=LOGF)
json.dump(dict(certified_coords=[n for n, ok in cert_ok.items() if ok], exact_forms={b: {n: {str(HF.blocks[b]['cols'][j][1:]): str(v) for j, v in fm.items()} for n, fm in exact_forms[b].items()} for b in (0, 1)}),
          open(f'{PK}/{tag}.C.json', 'w'))
# ---------------- Step D (exact, Singular over Q) ----------------
# variables: all K coordinates a0..a{K-1}; ideal = certified linear coordinates + exact forms (full, unrestricted)
lines = ['ring r = 0, (' + ','.join(f'a{k}' for k in range(K)) + '), dp;', 'ideal I;']
for name, ok in cert_ok.items():
    if ok:
        lin = coord_lin[name]
        lines.append('I[size(I)+1] = ' + '+'.join(f'({v})*a{k}' for k, v in lin.items()) + ';')
for blk in (0, 1):
    cols = HF.blocks[blk]['cols']
    for name, fm in exact_forms[blk].items():
        terms = []
        for j, v in fm.items():
            _, k, l = cols[j]; v = Fraction(v)
            terms.append(f'({v.numerator}/{v.denominator})*a{k}*a{l}')
        lines.append('I[size(I)+1] = ' + '+'.join(terms) + ';')
lines += ['int tt = timer;', 'ideal G = std(I);', '"dim V(I) over Q:", dim(G), " #gens:", size(I), " time:", timer - tt;', 'int bad = 0; ideal T;']
for name, lin in coords:
    u = '+'.join(f'({v})*a{k}' for k, v in lin.items())
    lines.append(f'T = G, {u} - 1; if (dim(std(T)) >= 0) {{ bad++; "coordinate {name} NOT in radical over Q"; }}')
lines += [f'"coordinates not in radical over Q:", bad, "of", {len(coords)};', 'quit;']
open(f'{PK}/{tag}.D.sing', 'w').write('\n'.join(lines) + '\n')
r = subprocess.run(['Singular', '-q', f'{PK}/{tag}.D.sing'], capture_output=True, text=True, timeout=7200)
open(f'{PK}/{tag}.D.out', 'w').write(r.stdout + r.stderr)
badD = next((int(l.split(':')[1].split('of')[0]) for l in r.stdout.splitlines() if l.startswith('coordinates not in radical over Q:')), None)
dimD = next((l for l in r.stdout.splitlines() if l.startswith('dim V(I) over Q:')), '')
result = ('PROVED over Q: no certificate of degree %d (all parities)' % d) if badD == 0 else f'incomplete over Q ({badD} coordinates not in radical)'
print('[D]', dimD, '|', result, f'({el()}s)', flush=True)
log(dict(exp='hpipe-D-groebner-Q', d=d, result=result, singular=dimD, n_certified_coords=sum(cert_ok.values()), n_exact_forms=sum(len(v) for v in exact_forms.values()), job=tag, seconds=el(),
         note='Singular char 0: ideal generated by exactly certified linear coordinates u (u^2 in I) and exactly certified quadratic forms in I; every coordinate of V_d/Sp is in the radical'), path=LOGF)
print('done', el(), 's')
