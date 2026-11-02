"""Summarize the merged wave-2 search log into markdown tables (per family/degree/sector: K, dimP, system size,
restarts, best/median residual, residual histogram, candidates) and control hit rates.  Writes report_tables.md."""
import json, collections, math, sys
L = '/tmp/hodge/results/four-plus-two-residual/wave2-search.log.jsonl'
recs = [json.loads(l) for l in open(L)]
best = {}
for r in recs:
    if r.get('kind') != 'batch' or not r.get('restarts'):
        continue
    key = (r.get('formulation', 'full'), r['family'], r['d'], tuple(r['sector']))
    b = best.setdefault(key, {'K': r.get('K'), 'dimP': r.get('dimP'), 'rows': r.get('rows'), 'cols': r.get('cols'), 'mode': set(),
                              'restarts': 0, 'best': 1.0, 'cands': 0, 'secs': 0.0, 'hist': collections.Counter(), 'medians': []})
    b['mode'].add(r.get('mode')); b['restarts'] += r['restarts']; b['best'] = min(b['best'], r.get('best_rel', 1.0))
    b['cands'] += r.get('n_cands', 0); b['secs'] += r.get('seconds', 0) or 0
    for k, v in (r.get('hist') or {}).items():
        b['hist'][k] += v
    if r.get('median_rel') is not None:
        b['medians'].append(r['median_rel'])
lines = ['| formulation | family | d | sector (h,r2,s) | K | dimP | rows x cols | mode | restarts | best rel | median rel | hist (log10 bins) | cands<1e-8 | CPU s |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
tot_restarts = 0
for key in sorted(best, key=lambda k: (k[0], k[1], k[2], k[3])):
    b = best[key]; tot_restarts += b['restarts']
    med = sorted(b['medians'])[len(b['medians']) // 2] if b['medians'] else float('nan')
    hist = ' '.join(f'{k}:{v}' for k, v in sorted(b['hist'].items(), key=lambda kv: float(kv[0].replace('<', '').replace('>=', '')) if kv[0] != '>=1' else 9))
    lines.append(f"| {key[0]} | {key[1]} | {key[2]} | {key[3]} | {b['K']} | {b['dimP']} | {b['rows']}x{b['cols']} | {'/'.join(sorted(m for m in b['mode'] if m))} | {b['restarts']} | {b['best']:.2e} | {med:.1e} | {hist} | {b['cands']} | {b['secs']:.0f} |")
search_tab = '\n'.join(lines)
# controls
ctl = collections.defaultdict(lambda: {'restarts': 0, 'hits': 0, 'solved': 0, 'K': None, 'dimP': None, 'sub': None})
for r in recs:
    if r.get('kind') == 'control-batch':
        key = (r.get('instance'), r.get('formulation', 'full'), r.get('method', '?'), r.get('subgroup_size'))
        c = ctl[key]; c['restarts'] += r['restarts']; c['hits'] += r['n_rel_lt_1e-10']; c['solved'] += r['n_solved']
        c['K'] = r.get('K', c['K']); c['dimP'] = r.get('dimP', c['dimP'])
    elif r.get('kind') == 'control':
        key = (r.get('instance'), r.get('formulation', 'full'), r.get('method', '?'), r.get('subgroup_size'))
        ctl[key]['K'] = r.get('K'); ctl[key]['dimP'] = r.get('dimP'); ctl[key]['sys'] = f"{r.get('rows')}x{r.get('cols')}"
clines = ['| instance | formulation | method | subgroup | K | dimP | system | restarts | hits (rel<1e-10) | exactly solved | hit rate |', '|---|---|---|---|---|---|---|---|---|---|---|']
for key in sorted(ctl, key=str):
    c = ctl[key]
    if not c['restarts']:
        continue
    clines.append(f"| {key[0]} | {key[1]} | {key[2]} | {key[3]} | {c['K']} | {c['dimP']} | {c.get('sys','')} | {c['restarts']} | {c['hits']} | {c['solved']} | {100.0*c['hits']/c['restarts']:.1f}% |")
ctl_tab = '\n'.join(clines)
# candidates / reconstructions
cands = [r for r in recs if r.get('kind') in ('cand', 'reconstruct')]
cl = [f"- {r.get('kind')} {r.get('tag', r.get('family'))} d={r.get('d')} rel={r.get('rel', r.get('polished_rel', ''))} status={r.get('status', r.get('solved'))}" for r in cands][:50]
contact = [r for r in recs if r.get('kind', '').startswith('contact')]
open('/tmp/hodge/results/four-plus-two-residual/wave2-search/report_tables.md', 'w').write(
    f'## search ({tot_restarts} restarts total)\n' + search_tab + '\n\n## controls\n' + ctl_tab + '\n\n## candidates\n' + ('\n'.join(cl) or 'none') +
    f'\n\n## contact records: {len(contact)}\n')
print(f'{tot_restarts} restarts in {len(best)} (formulation,family,d,sector) cells; {len(cands)} candidate records; controls: {len(ctl)}')
print(ctl_tab)
