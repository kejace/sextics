"""Regenerate the auto sections of the wave-2 report from the merged wave2 log (search table, structured summary)."""
import json, os, collections, subprocess, sys
REP = '/tmp/hodge/notes/wave2-one-reflection-skew-search.md'
LOG = '/tmp/hodge/results/one-reflection-skew/wave2-search.log.jsonl'
rows = [json.loads(l) for l in open(LOG) if l.strip()]
MIRROR = '/tmp/hodge/remote/one-reflection-skew'
import glob
# primary source: per-job status files in the mirror (complete even when log lines were lost)
jobs = {}
for p in glob.glob(f'{MIRROR}/*.status.json') + glob.glob('/tmp/hodge/results/one-reflection-skew/*.status.json'):
    try:
        s = json.load(open(p))
    except Exception:
        continue
    j = s['job']
    if j.startswith('smoke') or j.startswith('test-'):
        continue
    if j not in jobs or s.get('done', 0) >= jobs[j].get('done', 0):
        jobs[j] = s
starts = {r['job']: r for r in rows if r.get('event') == 'search-start'}
def inst_of(job, s):
    if job in starts:
        return starts[job].get('inst', 'one-reflection-skew')
    if job.startswith('calib-hesse'):
        return 'family-hesse-web'
    if job.startswith('calib-xtype'):
        return 'xtype-strict-rational-113'
    if job.startswith('ctrl-onerefl'):
        return 'one-reflection'
    return 'one-reflection-skew'
lines = ['| job | instance | d | sector | family | target | K | system rows×cols | method | restarts | best rel | median rel | #<1e-8 / <1e-6 / <1e-4 | wall |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
n_restarts = 0
for job, s in sorted(jobs.items(), key=lambda kv: (inst_of(kv[0], kv[1]) != 'one-reflection-skew', kv[1].get('d', 0), kv[0])):
    inst = inst_of(job, s)
    st = starts.get(job, {})
    sysz = f"{st.get('rows', '?')}×{st.get('cols', '?')}"
    h = s.get('hist', {}) or {}
    br, mr = s.get('best_rel'), s.get('median_rel')
    if inst == 'one-reflection-skew':
        n_restarts += s.get('done', 0)
    state = '' if s.get('final') else ' (running)'
    lines.append(f"| {job} | {inst} | {s['d']} | ({s['rho']},{s['sigma']}) | {s['family']} | {s['target']} | {s['K']} | {sysz} | {s['method']} | {s.get('done', 0)}{state} | "
                 f"{'—' if br is None else format(br, '.2e')} | {'—' if mr is None else format(mr, '.2e')} | {h.get('<1e-08', 0)} / {h.get('<1e-06', 0)} / {h.get('<0.0001', 0)} | {s.get('elapsed', '')} s |")
# empty-family searches only appear in the log
ends = {r['job']: r for r in rows if r.get('event') == 'search-end'}
for job, e in sorted(ends.items()):
    if job in jobs or job.startswith('smoke') or job.startswith('test-'):
        continue
    st = starts.get(job, {})
    lines.append(f"| {job} | {st.get('inst', '')} | {st.get('d', '')} | ({st.get('rho', '')},{st.get('sigma', '')}) | {st.get('family', '')} | {st.get('target', '')} | {st.get('K', 0)} | | | 0 | — | — | | {e.get('note', '')} |")
table = '\n'.join(lines)
# structured summary
sc = [r for r in rows if r.get('event') == 'struct-candidate']
fam = collections.defaultdict(lambda: [0, 0, 0, 0])
for r in sc:
    k = (r['family'].split(':')[0], r['d']); fam[k][0] += 1
    if r.get('e1'): fam[k][1] += 1
    if r.get('served_mod_J'): fam[k][2] += 1
    if r.get('serves_A_mod_J'): fam[k][3] += 1
sl = ['| family | d | candidates | pass E1 | served_mod_J > 0 | serves A |', '|---|---|---|---|---|---|']
for k in sorted(fam):
    v = fam[k]; sl.append(f'| {k[0]} | {k[1]} | {v[0]} | {v[1]} | {v[2]} | {v[3]} |')
hits = [r for r in sc if r.get('served_mod_J') or r.get('serves_A_mod_J')]
stable = '\n'.join(sl) + ('\n\nCandidates serving anything beyond J: ' + '; '.join(f"{r['name']} (d={r['d']}, served_mod_J={r.get('served_mod_J')}, A={r.get('serves_A_mod_J')})" for r in hits) if hits else '\n\nNo structured candidate serves any target beyond the trivially served Jacobian part.') if sc else '(no structured candidates processed yet)'
cands = [r for r in rows if r.get('event') in ('candidate', 'exact-confirm', 'WITNESS', 'unconfirmed-candidate') and r.get('inst', 'one-reflection-skew') == 'one-reflection-skew']
ctext = '\n'.join('* ' + json.dumps({k: v for k, v in r.items() if k not in ('ts', 'lane', 'engine', 'J_svals_tail')})[:400] for r in cands[-30:]) or '(none: no skew-lane restart reached the 1e-7 confirmation threshold)'
rep = open(REP).read()
def repl(rep, tag, body):
    a, b = f'<!-- AUTO-{tag} -->', f'<!-- /AUTO-{tag} -->'
    if a not in rep:
        return rep
    i, j = rep.index(a) + len(a), rep.index(b)
    return rep[:i] + '\n' + body + '\n' + rep[j:]
rep = repl(rep, 'SEARCH', table + f'\n\nTotal restarts on the skew lane so far: {n_restarts}.')
rep = repl(rep, 'STRUCT', stable)
rep = repl(rep, 'CANDS', ctext)
open(REP, 'w').write(rep)
print('report refreshed; skew-lane restarts:', n_restarts, 'jobs:', len(starts), 'struct candidates:', len(sc))

# exact kernel table (from exact-kernel events and cache files)
import glob as _g
kern = {}
for r in rows:
    if r.get('event') == 'exact-kernel' and 'error' not in r:
        kern[(r['d'], r['rho'], r['sigma'], r['target'])] = r
for p in _g.glob('/tmp/hodge/remote/one-reflection-skew/w2cache/kernel-one-reflection-skew-*.json'):
    try:
        cj = json.load(open(p)); inf = cj['info']
        tgt = 'A' if p.endswith('da2df3f3.json') else ('Aplus' if p.endswith('094329ba.json') else '?')
        key = (inf['d'], inf['pD'], inf['sD'], tgt)
        kern.setdefault(key, {'d': inf['d'], 'rho': inf['pD'], 'sigma': inf['sD'], 'target': tgt, 'kerD': inf['kerD'],
                              'primes_used': inf.get('primes_used'), 'max_height_bits': inf.get('max_height_bits'), 'secs': inf.get('secs')})
    except Exception:
        pass
kl = ['| d | sector | target | kernel dim | primes used | max height (bits) | seconds |', '|---|---|---|---|---|---|---|']
for key in sorted(kern):
    r = kern[key]
    kl.append(f"| {r['d']} | ({r['rho']},{r['sigma']}) | {r['target']} | {r['kerD']} | {r.get('primes_used')} | {r.get('max_height_bits')} | {r.get('secs')} |")
rep = open(REP).read()
rep = repl(rep, 'KERNELS', '\n'.join(kl) if kern else '(none yet)')
open(REP, 'w').write(rep)
print('kernel table entries:', len(kern))
