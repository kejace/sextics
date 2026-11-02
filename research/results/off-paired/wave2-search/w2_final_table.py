#!/usr/bin/env python3
"""Local (light): build the consolidated search table from all raw restart records (wave2-<job>.w*.jsonl in the remote mirror),
grouped by job; includes jobs without a summary (killed/running).  Prints markdown.  python3 w2_final_table.py [--md out.md]"""
import os, sys, json, glob, re, collections
MIR = '/tmp/hodge/remote/off-paired/wave2-search'
jobs = collections.defaultdict(list)
for fn in glob.glob(f'{MIR}/wave2-*.w*.jsonl'):
    m = re.match(r'.*/wave2-(.+)\.w\d+\.jsonl$', fn)
    if not m: continue
    for line in open(fn):
        try: jobs[m.group(1)].append(json.loads(line))
        except Exception: pass
summ = {}
for fn in glob.glob(f'{MIR}/wave2-*.summary.json'):
    try:
        s = json.load(open(fn)); summ[s['job']] = s
    except Exception: pass
rows = []
for job, recs in sorted(jobs.items()):
    rel = sorted(r['relres'] for r in recs if 'relres' in r)
    if not rel: continue
    s = summ.get(job, {})
    hist = collections.Counter(f'{x:.3e}' for x in rel)
    d = s.get('d'); sec = s.get('sector'); fam = s.get('family', '?'); grp = s.get('group', 'G'); dimP = s.get('dimP', '?')
    if d is None:
        m = re.match(r'(?:ft|hess-|t01|fe)?d(\d+)', job)
        d = int(m.group(1)) if m else '?'
    if sec is None:
        m = re.search(r'c(\d)(\d)', job)
        sec = [int(m.group(1)), int(m.group(2))] if m else None
    if fam == '?' and job.startswith(('d', 'ft')):
        fam = 'colon'
    step = s.get('step', 'varpro')
    rows.append({'job': job, 'step': step, 'd': d, 'sector': sec, 'group': grp, 'family': fam, 'dimP': dimP, 'restarts': len(rel),
                 'best': rel[0], 'median': rel[len(rel)//2], 'hits_1e-9': sum(1 for x in rel if x < 1e-9), 'distinct_minima': len(hist),
                 'top_min_count': hist.most_common(1)[0][1], 'mean_time': round(sum(r.get('time', 0) for r in recs) / len(recs), 1),
                 'done': 'yes' if job in summ else 'partial'})
rows.sort(key=lambda r: (str(r['d']), str(r['sector']), r['job']))
lines = ['| job | step | d | sector | grp | family | dimP | restarts | best relres | median | hits<1e-9 | distinct minima | s/restart | status |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
for r in rows:
    lines.append(f"| {r['job']} | {r['step']} | {r['d']} | {r['sector']} | {r['group']} | {r['family']} | {r['dimP']} | {r['restarts']} | {r['best']:.3e} | {r['median']:.2e} | {r['hits_1e-9']} | {r['distinct_minima']} | {r['mean_time']} | {r['done']} |")
table = '\n'.join(lines)
print(table)
tot = sum(r['restarts'] for r in rows); print(f'\ntotal restarts: {tot}; jobs: {len(rows)}; numerical hits (<1e-9): {sum(r["hits_1e-9"] for r in rows)}')
if '--md' in sys.argv:
    open(sys.argv[sys.argv.index('--md') + 1], 'w').write(table + f'\n\nTotal restarts {tot} over {len(rows)} jobs; numerical hits (relres < 1e-9): {sum(r["hits_1e-9"] for r in rows)}.\n')
