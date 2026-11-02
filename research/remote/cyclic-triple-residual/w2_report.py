"""Render markdown tables from wave2-search.log.jsonl (stdout)."""
import json, sys, math, collections
LOG = '/tmp/hodge/results/cyclic-triple-residual/wave2-search.log.jsonl'
rows = [json.loads(l) for l in open(LOG) if l.strip()]
def fmt(x):
    if x is None: return '-'
    if isinstance(x, float):
        return f'{x:.2e}' if (abs(x) < 1e-3 or abs(x) >= 1e4) and x != 0 else f'{x:.4g}'
    return str(x)
# E1 spaces
e1 = [r for r in rows if r.get('batch') == 'e1space' and 'totals' in r]
if e1:
    print('### E1-admissible spaces by degree (G-sectors; dims mod p = rigorous upper bounds, bases exact over Q)\n')
    print('| d | dim S_d | dim J_d | kerA in S_d | ker mod J | dimension-forced | sporadic | sectors with ker mod J > 0 | time |')
    print('|---|---|---|---|---|---|---|---|---|')
    seen = set()
    for r in sorted(e1, key=lambda r: r['d']):
        if r['d'] in seen: continue
        seen.add(r['d'])
        t = r['totals']
        secs = ', '.join(f"{k}:{v['ker_mod_J']}" for k, v in sorted(r['sectors'].items()) if v.get('ker_mod_J'))
        print(f"| {r['d']} | {t['dim_S']} | {t['dim_J']} | {t['kerA_in_S']} | {t['ker_mod_J']} | {t['forced']} | {t['sporadic']} | {secs or '-'} | {r['seconds']}s |")
    print()

mb = [r for r in rows if r.get('batch') == 'macaulay-blocks']
if mb:
    print('### Rigorous relaxation exclusions (block Macaulay ranks mod p; HF_k = 0 is rigorous over Q)\n')
    print('| d | Y-parity | restricted to complement of F-gauge locus | #w vars | quadric span rank | k | HF_k | rigorous exclusion | time |')
    print('|---|---|---|---|---|---|---|---|---|')
    for r in mb:
        print(f"| {r['d']} | {r['parity']} | {r.get('restricted_to_complement_of_spurious', False)} | {r['nW']} | {r.get('quadric_span_rank', '-')} | {r['k']} | {r['HF_k_total']} | {'YES' if r['rigorous_exclusion'] else 'no'} | {r['seconds']}s |")
    print()
rp = [r for r in rows if r.get('batch') in ('relaxation-point-d6', 'relax-d6-even-local')]
for r in rp:
    print('- ' + r.get('batch') + ': ' + json.dumps({k: v for k, v in r.items() if k not in ('best_a', 'best_D_float', 'ts', 'engine', 'lane', 'batch', 'rels_sorted')}, default=str)[:900])
if rp: print()
# searches
sb = [r for r in rows if r.get('batch') == 'search']
if sb:
    print('### Joint searches (variable projection; rel = f/|b|^2 at the best restart; F-multiples projected out)\n')
    print('| label | mode | group | K (dimP) | rows x cols (c,Q) | restarts | best rel | median rel | #<1e-6 | #<1e-10 | spurious | mean iters | time |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for r in sb:
        mode = 'mod-J relax' if r.get('relax') else 'full E2'
        print(f"| {r['label']} | {mode} | {r.get('group', r.get('group_order'))} | {r['K']} ({r['dimP']}) | {r['rows']} x {r['cols']} ({r['nC']},{r['nQ']}) | {r['restarts']} | {fmt(r['best_rel'])} | {fmt(r['median_rel'])} | {r['n_rel_below_1e-6']} | {r['n_rel_below_1e-10']} | {r['n_spurious']} | {fmt(r['mean_iters'])} | {r['seconds']}s |")
    print()
rf = [r for r in rows if r.get('batch') == 'refine']
if rf:
    print('### Relaxation minima refined by the full search\n')
    for r in rf:
        print(f"- {r['label']}: " + '; '.join(f"{x['from_relax_rel']:.1e} -> {x['full_rel']:.1e}" for x in r['results']))
    print()
pr = [r for r in rows if r.get('batch') == 'promote']
if pr:
    print('### Rational reconstructions tested exactly\n')
    print('| label | d | descr | recon err | num rel | E1 | consistent | rank/aug | witness |')
    print('|---|---|---|---|---|---|---|---|---|')
    for r in pr:
        print(f"| {r['tag']} | {r['d']} | {r['descr']} | {fmt(r.get('recon_err'))} | {fmt(r.get('num_rel'))} | {r.get('e1')} | {r.get('consistent')} | {r.get('rank')}/{r.get('rank_aug')} | {'YES' if r.get('WITNESS') else 'no'} |")
    print()
sw = [r for r in rows if r.get('batch') == 'sweep']
if sw:
    print('### Structured candidates (sweep; exact sector solver for d <= 9, two primes above)\n')
    print('| family | candidate | d | terms | F|D | E1 | consistent | rank/aug (cols) | exact | time |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    for r in sorted(sw, key=lambda r: (r['d'], r['family'])):
        print(f"| {r['family']} | {r['name']} | {r['d']} | {r['terms']} | {r.get('F_divides_D')} | {r.get('e1')} | {r.get('consistent')} | {r.get('rank')}/{r.get('rank_aug')} ({r.get('cols')}) | {r.get('exact_rank')} | {r.get('seconds')}s |" + (f" error {r['error']}" if r.get('error') else ''))
    print()
others = [r for r in rows if r.get('batch') not in ('e1space', 'search', 'refine', 'promote', 'sweep')]
if others:
    print('### Other log entries\n')
    for r in others:
        print('-', r.get('batch'), r.get('note', '')[:400])
