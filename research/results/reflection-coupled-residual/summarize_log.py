"""Summarise wave2-search.log.jsonl (lane reflection-coupled-residual) as markdown tables."""
import json, sys, collections
path = sys.argv[1] if len(sys.argv) > 1 else '/tmp/hodge/results/reflection-coupled-residual/wave2-search.log.jsonl'
recs = [json.loads(l) for l in open(path) if l.strip()]
recs = [r for r in recs if 'control_instance' not in r]
by = collections.defaultdict(list)
for r in recs:
    by[r.get('event')].append(r)
print('events:', {k: len(v) for k, v in by.items()})
print('\n### E1 spaces (exact mod p)')
print('| d | sec | dim S | dim J | dim R | E1-admissible | non-J |')
seen = set()
for r in sorted(by['e1space'], key=lambda r: (r['d'], r['sec'])):
    if (r['d'], r['sec']) in seen: continue
    seen.add((r['d'], r['sec']))
    print(f"| {r['d']} | {r['sec']} | {r['dimS']} | {r['dimJ']} | {r['dimR']} | {r['kernel_dim']} | {r['nonJ']} |")
print('\n### Hessian relaxation, exact (msolve)')
print('| d | sectors | family | K | polys | empty |')
for r in by['msolve'] + by['msolve-family']:
    print(f"| {r['d']} | {r['secs']} | {r.get('family','J')} | {r['K']} | {r.get('nquadrics', r.get('npolys'))} | {r.get('empty_projective')} {'(timeout)' if r.get('timed_out') else ''} |")
print('\n### Relaxation LM landscapes')
print('| d | sectors | family/affine | K | search dim | m | restarts | best | median |')
for r in by['relax'] + by['famrelax']:
    print(f"| {r['d']} | {r['secs']} | {r.get('family', r.get('affine'))} | {r['K']} | {r.get('dimP')} | {r['m']} | {r['restarts']} | {r['best_rel']:.2e} | {r['median_rel']:.2e} |")
print('\n### Koszul-reduced full-E2 varpro')
print('| d | sectors | family | K | search dim | rows | n_a | restarts | best rel | median | time |')
for r in by['kvp']:
    print(f"| {r['d']} | {r['secs']} | {r['family']} | {r['K']} | {r.get('dimP')} | {r.get('rows')} | {r.get('n_a')} | {r['restarts']} | {r['best_rel']:.2e} | {r['median_rel']:.2e} | {r['seconds']}s |")
print('\n### Exact confirmations')
for r in by['confirm']:
    print(f"- d={r['d']} {r['status']} rank={r.get('rank')} nullity={r.get('nullity')} note={r.get('note')} D={r['D'][:70]}")
for ev in ('k8-e1', 'k8-fixedD-symmetric', 'k8-fixedD-full', 'served', 'screen', 'kvp-candidate', 'relax-reconstruct'):
    for r in by[ev]:
        print(f"- {ev}: " + json.dumps({k: v for k, v in r.items() if k not in ('event', 'lane', 'ts', 'engine', 'D')})[:300])
