"""Regenerate the wave-2 report from the log + receipts.  usage: python3 w2_report.py [extra-notes.md]"""
import json, os, sys, glob, time, collections
LANE = '/tmp/hodge/results/off-paired'; PK = LANE + '/wave2-excl'
REPORT = '/tmp/hodge/notes/wave2-off-paired-excl.md'
rows = []
for path in [LANE + '/wave2-excl.log.jsonl'] + glob.glob('/tmp/hodge/remote/off-paired/wave2-excl.log.jsonl'):
    if os.path.exists(path):
        for l in open(path):
            try: rows.append(json.loads(l))
            except Exception: pass
# dedupe
seen = set(); uniq = []
for r in rows:
    k = json.dumps({a: b for a, b in r.items() if a not in ('ts',)}, sort_keys=True, default=str)
    if k not in seen: seen.add(k); uniq.append(r)
rows = uniq
witness = any(r.get('result') in ('WITNESS',) or (r.get('exp') == 'fixedD' and r.get('consistent') is True) for r in rows)
status = 'WITNESS FOUND (see log)' if witness else 'no witness'
out = []
out.append('# Wave 2 — lane `off-paired` — lens EXCLUSION / RIGOR\n')
out.append(f'STATUS: {status}; last update {time.strftime("%Y-%m-%d %H:%M")}; {len(rows)} logged experiments.\n')
out.append('Packet: `/tmp/hodge/results/off-paired/wave2-excl/` (README.md there explains every receipt and how to replay). '
           'Log: `/tmp/hodge/results/off-paired/wave2-excl.log.jsonl`. Wave-2 primes 2147483629 / 2147483587 (fixed-D), 2147483647 (msolve).\n')
out.append('## Rigor conventions\n* modular rank <= rational rank; compressed rank <= modular rank.\n'
           '* Fixed-D exclusion is PROVEN when rank_p([M|b]) = rank_p(M)+1 and rank_Q(M) <= rank_p(M) is certified by the explicit rational kernel '
           '(Euler gauge L*Z and F-multiples (F W, B_D(W)), dimension e computed exactly) or by exact verification of a rationally reconstructed basis of the whole modular kernel.\n'
           '* E1 failure is PROVEN when rank_p(J_{d+6}[chi]) equals the exact Koszul dimension of J(F)_{d+6}[chi] and the augmented rank is one more.\n'
           '* Bilinear sector results (msolve, charts) are modular SCREENS unless a rational certificate is recorded.\n')
# fixed-D summary by family
fam = collections.defaultdict(lambda: collections.Counter())
for r in rows:
    if r.get('exp') in ('fixedD', 'e1'):
        fam[(r.get('family') or '?', r.get('exp'))][r.get('result')] += 1
out.append('## Fixed-denominator receipts by family (counts of (label, prime) records)\n')
out.append('| family | exp | results |\n|---|---|---|')
for (f, e), cnt in sorted(fam.items()):
    out.append(f'| {f} | {e} | ' + ', '.join(f'{k}: {v}' for k, v in sorted(cnt.items(), key=lambda x: str(x[0]))) + ' |')
# per-label rigorous status
lab = {}
for r in rows:
    if r.get('exp') in ('fixedD', 'e1'):
        L = r.get('label') or r.get('D')
        s = lab.setdefault(L, {'d': r.get('d'), 'char': r.get('char'), 'e1': set(), 'full': set(), 'family': r.get('family')})
        (s['e1'] if r['exp'] == 'e1' else s['full']).add(str(r.get('result')))
nprov = sum(1 for L, s in lab.items() if any(x.startswith('PROVEN') for x in s['e1'] | s['full']))
nmod = sum(1 for L, s in lab.items() if not any(x.startswith('PROVEN') for x in s['e1'] | s['full']) and any(('inconsistent' in x or 'fail' in x) for x in s['e1'] | s['full']))
nsolv = [L for L, s in lab.items() if any('solvable' in x or 'WITNESS' in x for x in s['full'])]
out.append(f'\nDistinct denominators tested: {len(lab)}; PROVEN excluded (rational rank sandwich): {nprov}; excluded modulo p only: {nmod}; solvable mod p: {nsolv}\n')
# bilinear
bil = [r for r in rows if r.get('exp') == 'bilinear']
if bil:
    out.append('## Bilinear sector systems (D free in a character sector; msolve on charts)\n')
    out.append('| d | char | mode | tsub | nt | ns | neq | charts run | no-solution | other | max s |\n|---|---|---|---|---|---|---|---|---|---|---|')
    grp = collections.defaultdict(list)
    for r in bil: grp[(r['d'], str(r.get('char')), r.get('mode'), r.get('tsub'))].append(r)
    for k, rs in sorted(grp.items(), key=lambda x: (x[0][0], x[0][1])):
        ok = sum(1 for r in rs if str(r.get('result', '')).startswith('no-solution'))
        other = [f"{r.get('chart')}:{r.get('result')}" for r in rs if not str(r.get('result', '')).startswith('no-solution') and r.get('result') != 'written']
        out.append(f"| {k[0]} | {k[1]} | {k[2]} | {k[3]} | {rs[0].get('nt')} | {rs[0].get('ns')} | {rs[0].get('neq_nonzero')} | {len(rs)} | {ok} | {' '.join(other)[:80]} | {max(float(r.get('seconds') or 0) for r in rs):.0f} |")
# kernel dims
kd = [r for r in rows if r.get('exp') == 'kernel-dims']
if kd:
    out.append('\n## Colon / non-Jacobian E1-kernel dimensions (colon mod p = upper bound for Q; dim J exact by Koszul)\n')
    out.append('| d | colon | dim J_d | non-Jacobian kernel | gauge | prime | nonjac per sector (nonzero) |\n|---|---|---|---|---|---|---|')
    seen_d = set()
    for r in sorted(kd, key=lambda r: r['d']):
        if r['d'] in seen_d: continue
        seen_d.add(r['d']); t = r['total']
        out.append(f"| {r['d']} | {t['colon']} | {t['J']} | {t['nonjac']} | {t['gauge']} | {r['prime']} | " + ' '.join(f"{k}:{v}" for k, v in r['nonjac_per_sector'].items() if v) + ' |')
# notable single results
out.append('\n## Notable single receipts\n')
for r in rows:
    if r.get('family') in ('discriminant', 'colon-quartic', 'fermat-cone', 'genus2-discriminant-component') or r.get('d', 0) >= 15:
        out.append('* ' + json.dumps({k: r.get(k) for k in ('label', 'D', 'd', 'char', 'exp', 'prime', 'rows', 'cols', 'rank_M', 'rank_aug', 'rank_J', 'dimJ_exact', 'explicit_kernel', 'kernel_verified', 'kernel_total', 'result', 'rigor', 'seconds') if r.get(k) is not None}, default=str))
if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
    out.append('\n' + open(sys.argv[1]).read())
open(REPORT, 'w').write('\n'.join(out) + '\n')
print('report written', REPORT, len(rows), 'rows')
