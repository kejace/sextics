"""Sector table for a degree from the canonical log (lens excl). usage: python3 wave2-excl/mktable.py --deg 7"""
import json, argparse, collections, glob, os
ap = argparse.ArgumentParser(); ap.add_argument('--deg', type=int, required=True); ap.add_argument('--prime', type=int, default=1000000007)
a = ap.parse_args(); d = a.deg; p = a.prime
rows = [json.loads(l) for l in open('/tmp/hodge/results/cycle-q2/wave2-excl.log.jsonl')]
secs = {}
for r in rows:
    if r.get('exp') == f'sector-d{d}' and r.get('prime') == p and r.get('engine') == 2:
        secs[tuple(r['sector'])] = r
els = collections.defaultdict(list)
for r in rows:
    if r.get('exp') == f'sector-d{d}-element' and r.get('engine') == 2 and r.get('sector') is not None:
        els[tuple(r['sector'])].append(r)
certs = collections.Counter()
for f in glob.glob(f'/tmp/hodge/results/cycle-q2/wave2-excl/certs/d{d}_p{p}_k*.json'):
    try:
        c = json.load(open(f))
        if c.get('verified_exact'): certs[tuple(int(ch) for ch in os.path.basename(f).split('_k')[1][:3])] += 1
    except Exception: pass
print(f'# Degree-{d} sector table (p = {p}); {len(secs)} of 72 sectors analysed so far\n')
print('| sector | dim S | dim K | dim J | non-J | relaxation quadrics | relaxation affine dim | Farkas certificates | status |')
print('|---|---|---|---|---|---|---|---|---|')
totK = totJ = 0
for k, r in sorted(secs.items()):
    K, J = r['dimK'], r['dimJ']; totK += K; totJ += J
    if K == 0: continue
    nc = certs.get(k, 0)
    st = ('single element: Farkas-excluded' if K == 1 and nc else ('single element: pending' if K == 1 else
          ('family EXCLUDED (relaxation empty)' if r.get('relax_dim_affine') == 0 else f"family not decided ({K}-dim; {nc} elements/points excluded)")))
    print(f"| {k} | {r['dimS']} | {K} | {J} | {K-J} | {r.get('relax_nquadrics','-')} | {r.get('relax_dim_affine','-')} | {nc} | {st} |")
print(f'\nTotals over analysed sectors: dim K = {totK}, dim J = {totJ}, non-Jacobian = {totK-totJ}; certificates: {sum(certs.values())}')
