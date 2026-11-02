"""Summarise wave-2 logs (local mirror of the remote lane dir) into tables for the report.

Usage: python3 w2_report.py [--mirror /tmp/hodge/remote/one-reflection-skew] [--merge]
  --merge  append new lines of the mirror's wave2-search.log.jsonl into the local lane log (dedup by line).
Prints: search-end/batch summary per job, dims table, struct summary, candidates, witnesses.
"""
import sys, os, json, argparse, collections
ap = argparse.ArgumentParser()
ap.add_argument('--mirror', default='/tmp/hodge/remote/one-reflection-skew')
ap.add_argument('--local', default='/tmp/hodge/results/one-reflection-skew')
ap.add_argument('--merge', action='store_true')
a = ap.parse_args()

def read_jsonl(path):
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except Exception:
            pass
    return out

if a.merge:
    rl, ll = f'{a.mirror}/wave2-search.log.jsonl', f'{a.local}/wave2-search.log.jsonl'
    have = set(open(ll).read().splitlines()) if os.path.exists(ll) else set()
    added = 0
    if os.path.exists(rl):
        with open(ll, 'a') as fh:
            for line in open(rl).read().splitlines():
                if line and line not in have:
                    fh.write(line + '\n'); have.add(line); added += 1
    print(f'merged {added} new log lines into {ll}')

rows = read_jsonl(f'{a.local}/wave2-search.log.jsonl')
print(f'{len(rows)} log lines')
# dims
dims5 = [r for r in rows if r.get('event') == 'e1-dims-5var']
if dims5:
    print('\n## E1 dims (5-var): d, (rho,sigma): nD, kerD, J_d, ker_mod_J, fmult, rank[A]')
    seen = set()
    for r in sorted(dims5, key=lambda r: (r['d'], r.get('pD', 0), r.get('sD') or 0)):
        key = (r['d'], r.get('pD'), r.get('sD'))
        if key in seen: continue
        seen.add(key)
        print(f"  d={r['d']} ({r.get('pD')},{r.get('sD')}): nD={r.get('nD')} ker={r.get('kerD')} J_d={r.get('dimJ_d')} "
              f"ker/J={r.get('ker_mod_J')} fmult={r.get('dim_fmult')} rank[A]={r.get('rank_AD_mod_J')} J_d+6={r.get('dimJ_d6')} {r.get('error','')}")
dims6 = [r for r in rows if r.get('event') == 'e1-dims-6var']
if dims6:
    print('\n## E1 dims (6-var, Y-even D, engine mult_kernel):')
    for r in dims6:
        print(f"  d={r['d']}: kernel={r.get('kernel_dim')} J_d={r.get('jacobian_dim')} ker_mod_J={r.get('ker_mod_J')} ring={r.get('ring_dim')} {r.get('error','')}")
# searches
starts = {r['job']: r for r in rows if r.get('event') == 'search-start'}
ends = {r['job']: r for r in rows if r.get('event') == 'search-end'}
batches = collections.defaultdict(list)
for r in rows:
    if r.get('event') == 'batch':
        batches[r['job']].append(r)
print('\n## Searches: job | inst | d (rho,sigma) family target | K rows cols | method | restarts | best rel | median | hist(<1e-8,<1e-6,<1e-4,<1e-2)')
for job, s in sorted(starts.items(), key=lambda kv: (kv[1].get('d', 0), kv[0])):
    e = ends.get(job)
    bl = batches.get(job, [])
    if e:
        h = e.get('hist', {})
        hs = f"{h.get('<1e-08', 0)}/{h.get('<1e-06', 0)}/{h.get('<0.0001', 0)}/{h.get('<0.01', 0)}"
        print(f"  {job:22s} {s.get('inst','')[:16]:16s} d={s['d']} ({s['rho']},{s['sigma']}) {s['family'][:14]:14s} {s['target']:5s} K={s['K']:<4d} {s['rows']}x{s['cols']} {s['method']:5s} "
              f"n={e['restarts']:<5d} best={e['best_rel'] if e['best_rel'] is None else format(e['best_rel'], '.2e')} med={e['median_rel'] if e['median_rel'] is None else format(e['median_rel'], '.2e')} hist={hs} {e['elapsed']}s cands={e.get('n_cands')} conf={e.get('confirmed')}")
    elif bl:
        b = bl[-1]
        print(f"  {job:22s} {s.get('inst','')[:16]:16s} d={s['d']} ({s['rho']},{s['sigma']}) {s['family'][:14]:14s} {s['target']:5s} K={s['K']:<4d} {s['rows']}x{s['cols']} {s['method']:5s} "
              f"RUNNING n={b['done_total']} best={format(b['best_rel_total'], '.2e') if b['best_rel_total'] is not None else None} cands={b.get('n_cands')}")
    else:
        print(f"  {job:22s} {s.get('inst','')[:16]:16s} d={s['d']} ({s['rho']},{s['sigma']}) {s['family'][:14]:14s} K={s['K']} STARTED (no batch yet)")
# structured
sc = [r for r in rows if r.get('event') == 'struct-candidate']
if sc:
    print(f'\n## Structured candidates: {len(sc)}')
    fam = collections.defaultdict(lambda: [0, 0, 0, 0])
    for r in sc:
        k = (r['family'].split(':')[0], r['d'])
        fam[k][0] += 1
        if r.get('e1'): fam[k][1] += 1
        if r.get('served_mod_J'): fam[k][2] += 1
        if r.get('serves_A_mod_J'): fam[k][3] += 1
    for k in sorted(fam):
        v = fam[k]
        print(f'  {k[0]:28s} d={k[1]:<3d} n={v[0]:<4d} e1-pass={v[1]:<4d} served_mod_J>0: {v[2]}  serves A: {v[3]}')
    for r in sc:
        if r.get('served_mod_J') or r.get('serves_A_mod_J') or r.get('confirm'):
            print('  !!', r['name'], r['d'], r.get('served_dim'), r.get('served_mod_J'), r.get('serves_A_mod_J'), r.get('confirm'))
# candidates / witnesses
for ev in ('candidate', 'exact-confirm', 'WITNESS'):
    L = [r for r in rows if r.get('event') == ev]
    if L:
        print(f'\n## {ev}: {len(L)}')
        for r in L[-20:]:
            print('  ', json.dumps(r)[:300])
