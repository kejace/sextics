"""hess_lift: exact rational quadratic forms in a spanning (a lower bound for) the J-relaxed weight-2 ideal, via
multi-prime CRT of the rref rows with aa-pivots, exact certification (lambda^T E' = (0|phi) over Q), and a Singular-over-Q
script testing V(forms) subset Sp.  Usage: python3 hess_lift.py <d> <nprimes> [--nonj FILE] [--blocks 0,1] [--tag TAG]"""
import sys, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
from sympy import prevprime
from collections import Counter
ap = argparse.ArgumentParser(); ap.add_argument('d', type=int); ap.add_argument('nprimes', type=int); ap.add_argument('--nonj', default=''); ap.add_argument('--blocks', default='0,1'); ap.add_argument('--tag', default='')
ap.add_argument('--singular', action='store_true', help='also run Singular over Q on the exact forms'); ap.add_argument('--sparse', action='store_true', help='use mod-J sparse column representatives (Hessian identity)')
a = ap.parse_args(); d = a.d
tag = a.tag or f'hesslift_d{d}'
LOGF = f'{PK}/{tag}.log.jsonl'
t0 = time.time()
basis = polar_basis(d)
if a.nonj:
    for k, rec_ in enumerate(json.load(open(a.nonj))):
        basis.append((('K', rec_.get('label', k)), core.parse(rec_['D']), [core.parse(s) for s in rec_['N']]))
K = len(basis); labidx = {b[0]: k for k, b in enumerate(basis)}
HF = HessianFamily(d, basis, sparse=a.sparse)
primes = []; x = 1000000007 + 1
while len(primes) < a.nprimes:
    x = prevprime(x); primes.append(int(x))
exact_forms = {}
for blk in [int(x) for x in a.blocks.split(',')]:
    B = HF.blocks[blk]; C = B['ncols']; nJ = HF.nJ[blk]
    per = []
    for p in primes:
        info, Mr, piv = HF.analyze(blk, p, verbose=False)
        rows_ = []
        for i, j in enumerate(piv):
            if j < nJ: continue
            rows_.append({k: int(Mr[i, k]) for k in range(j, C) if int(Mr[i, k])})
        per.append(dict(p=p, rank=info['rank'], rank_J=info['rank_J'], rows=rows_, pivs=tuple(min(r_) for r_ in rows_)))
        print(f'[blk{blk} p={p}] rank {info["rank"]} rank_J {info["rank_J"]} forms {len(rows_)} ({time.time()-t0:.0f}s)', flush=True)
        del Mr
    piv0 = Counter(x['pivs'] for x in per).most_common(1)[0][0]
    good = [x for x in per if x['pivs'] == piv0]
    nforms = len(good[0]['rows'])
    lifted = []
    fail = 0
    for i in range(nforms):
        form = {}
        keys = set()
        for x in good: keys |= set(x['rows'][i])
        ok = True
        for k in sorted(keys):
            a_, m_ = good[0]['rows'][i].get(k, 0), good[0]['p']
            for x in good[1:]:
                a_, m_ = core.crt_pair(a_, m_, x['rows'][i].get(k, 0), x['p'])
            rr = core.ratrecon(a_, m_)
            if rr is None: ok = False; break
            if rr: form[k] = rr
        if ok: lifted.append((i, form))
        else: fail += 1
    print(f'[blk{blk}] lifted {len(lifted)} of {nforms} forms with {len(good)} primes; ratrecon failures {fail}', flush=True)
    # exact certification at the first prime's pivot data
    info, Mr, piv = HF.analyze(blk, primes[0], verbose=False)
    targets = [(f'form{blk}_{i}', form) for i, form in lifted]
    certs = hess_certificates(HF, blk, Mr, piv, targets, primes[0], verbose=False) if targets else {}
    del Mr
    nok = sum(1 for v in (certs or {}).values() if v)
    hts = max((max(abs(v.numerator), v.denominator) for _, form in lifted for v in form.values()), default=1)
    print(f'[blk{blk}] exact certificates verified: {nok} of {len(lifted)}; max height digits {len(str(hts))} ({time.time()-t0:.0f}s)', flush=True)
    exact_forms[blk] = [{str(B['cols'][k][1:]): str(v) for k, v in form.items()} for (i, form), name in zip(lifted, [n for n, _ in targets]) if certs and certs[name]]
    log(dict(exp='hessian-relaxation-exact-forms', d=d, blk=blk, primes=len(good), forms_modp=nforms, lifted=len(lifted), certified=nok, max_height_digits=len(str(hts)),
             rank=good[0]['rank'], rank_J=good[0]['rank_J'], job=tag, seconds=round(time.time() - t0),
             note='rref rows with aa-pivots CRT-lifted and certified exactly: each is a quadratic form in a lying in the E2 ideal over Q'), path=LOGF)
json.dump(dict(d=d, K=K, basis_labels=[str(b[0]) for b in basis], forms=exact_forms), open(f'{PK}/{tag}.forms.json', 'w'))
# Singular over Q script
sp_groups = []
if d >= 6:
    for m in monomials(d - 6):
        sp_groups.append([labidx[('J', i, tuple(m[j] + (1 if j == i else 0) for j in range(NV)))] for i in range(NV)])
in_group = {k for grp in sp_groups for k in grp}
coords = [(f'a{k}', {k: 1}) for k in range(K) if k not in in_group]
for grp in sp_groups:
    for i in range(1, NV): coords.append((f'a{grp[0]}-a{grp[i]}', {grp[0]: 1, grp[i]: -1}))
lines = ['ring r = 0, (' + ','.join(f'a{k}' for k in range(K)) + '), dp;', 'ideal I;']
for blk, forms in exact_forms.items():
    for fm in forms:
        terms = []
        for key, v in fm.items():
            k, l = eval(key); v = Fraction(v)
            terms.append(f'({v.numerator}/{v.denominator})*a{k}*a{l}')
        lines.append('I[size(I)+1] = ' + '+'.join(terms) + ';')
lines += ['int tt = timer;', 'ideal G = std(I);', '"dim V(I) (affine, over Q):", dim(G), " #forms:", size(I), " time:", timer - tt;', 'int bad = 0; ideal T;']
for name, lin in coords:
    u = '+'.join(f'({v})*a{k}' for k, v in lin.items())
    lines.append(f'T = G, {u} - 1; if (dim(std(T)) >= 0) {{ bad++; "coordinate {name} NOT in radical"; }}')
lines += [f'"coordinates not in radical:", bad, "of", {len(coords)};', 'quit;']
open(f'{PK}/{tag}.Q.sing', 'w').write('\n'.join(lines) + '\n')
print('wrote', f'{PK}/{tag}.Q.sing', flush=True)
if a.singular:
    import subprocess
    r = subprocess.run(['Singular', '-q', f'{PK}/{tag}.Q.sing'], capture_output=True, text=True)
    open(f'{PK}/{tag}.Q.out', 'w').write(r.stdout + r.stderr)
    print(r.stdout[-2000:], flush=True)
    bad = None
    for line in r.stdout.splitlines():
        if line.startswith('coordinates not in radical:'):
            bad = int(line.split(':')[1].split('of')[0])
    log(dict(exp='hessian-relaxation-groebner-Q', d=d, result=('PROVED over Q: V(exact forms) subset Sp => no certificate of degree %d' % d) if bad == 0 else f'incomplete ({bad} coordinates not in radical)',
             n_forms=sum(len(v) for v in exact_forms.values()), job=tag, seconds=round(time.time() - t0), singular_out=r.stdout[-500:]), path=LOGF)
print('done', round(time.time() - t0), 's')
