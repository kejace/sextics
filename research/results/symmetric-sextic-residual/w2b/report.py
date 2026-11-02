"""Rewrite /tmp/hodge/notes/wave2-symmetric-sextic-residual-search.md from the merged log + static narrative."""
import json, os, time, collections, math
LANE = 'symmetric-sextic-residual'; RES = f'/tmp/hodge/results/{LANE}'; W2B = f'{RES}/w2b'
REPORT = f'/tmp/hodge/notes/wave2-{LANE}-search.md'
recs = []
for line in open(f'{RES}/wave2-search.log.jsonl'):
    try: recs.append(json.loads(line))
    except Exception: pass
wit = os.path.exists(f'{RES}/witness-candidate.json')
witness_lines = [r for r in recs if r.get('WITNESS') or r.get('event') == 'WITNESS']
out = []
out.append(f'# Wave 2 search report: `{LANE}`\n')
if wit and witness_lines:
    out.append(f'STATUS: **WITNESS FOUND** — see `{RES}/witness-candidate.json` and the WITNESS log lines ({len(witness_lines)}).\n')
else:
    out.append(f'STATUS: NO WITNESS FOUND (wave 2 started 2026-11-01 ~19:10 EST; report regenerated {time.strftime("%Y-%m-%d %H:%M")} EST from {len(recs)} log records).\n')
out.append(f'Log: `{RES}/wave2-search.log.jsonl` (one JSON line per batch; engine tag `w2lib-v1`). Scripts: `{W2B}/` (`w2lib.py`, `w2_varpro.py`, `w2_kernels.py`, `w2_single*.py`, `w2_struct.py`, `scheduler.py`, `report.py`).\n')
out.append(open(f'{W2B}/report_static.md').read())
# ---- kernel dims
K = [r for r in recs if r.get('family') == 'e1-kernel-dims' and r.get('job') == 'kernels']
if K:
    tab = collections.defaultdict(dict)
    for r in K: tab[(r['chi'], r['pD'])][r['d']] = (r['dimK'], r['dimJ'] or 0, r['ker_beyond_J'])
    out.append('\n## 5. E1-admissible denominator spaces (job `kernels`, mod p = 2^31−1)\n')
    out.append('`dimK / dimJ / beyond` = dim{D of this type : A·D ∈ J(F)} / dim J(F)_d inside the type / kernel classes modulo J.\n')
    out.append('| sector (χ, pD) | ' + ' | '.join(f'd={d}' for d in range(4, 15)) + ' |')
    out.append('|---|' + '---|' * 11)
    for key in [('triv', 0), ('triv', 1), ('sgnX', 0), ('sgnX', 1), ('sgnY', 0), ('sgnY', 1), ('sgnsgn', 0), ('sgnsgn', 1)]:
        row = tab.get(key, {}); out.append(f'| {key[0]}, pD={key[1]} | ' + ' | '.join(f'{row[d][0]}/{row[d][1]}/{row[d][2]}' if d in row else '' for d in range(4, 15)) + ' |')
    subs = [(k, v) for k, v in sorted(tab.items()) if k[0].startswith('S') and 10 in v]
    if subs:
        out.append('\nSubgroup types at d = 10 (`dimK/dimJ/beyond`): ' + '; '.join(f'{k[0]} pD={k[1]}: {v[10][0]}/{v[10][1]}/{v[10][2]}' for k, v in subs) + '.')
    S = [r for r in K if 'served_random_points' in r]
    nb = sum(1 for r in S for s in r['served_random_points'] if s.get('beyond_J', 0) > 0)
    out.append(f'\n* Served-target check at 2 random admissible points per sector (d ≤ 10): {len(S)} sectors, points serving something beyond `J(F)_6^G`: **{nb}**.')
    KF = [r for r in recs if r.get('family') == 'e1-kernel-dims-unrestricted']
    if KF:
        out.append('* Unrestricted kernels (job `kernfull`): ' + '; '.join(f"d={r['d']} pD={r['pD']}: beyond J = {r.get('ker_beyond_J')} (kernel {r.get('kernel_dim')}, J {r.get('jacobian_dim')}, R/J {r.get('ring_dim')})" for r in sorted(KF, key=lambda r: (r['d'], r['pD'])) if 'error' not in r) + '.')
Fm = [r for r in recs if r.get('family') == 'fermat-deformation-kernel']
if Fm:
    tab2 = collections.defaultdict(dict)
    for r in Fm: tab2[(r['chi'], r['d'], r['pD'])][r['t']] = r['ker_beyond_J']
    same = all(len(set(v[t] for t in v if t != '0')) == 1 for v in tab2.values())
    out.append(f'* **Fermat deformation** `F_t = p6X + p6Y + t·e2X e2Y (e2X + 2e2Y)`, t ∈ {{0, 5/7, 1, 2, 1/3, −1, 7/5, 3/11, −5/7}}: beyond-J kernels identical for all t ≠ 0: **{same}**; they jump only at t = 0. '
               'Examples (t=0 vs t≠0): ' + '; '.join(f'{k}: {v.get("0")} vs {v.get("5/7")}' for k, v in sorted(tab2.items()) if k[0] == 'triv') + '. So 5/7 is not special for E1.')
# ---- varpro landscape
V = [r for r in recs if r.get('family') == 'varpro-sector' and 'best_rel' in r]
if V:
    out.append('\n## 7. Joint (varpro) search landscape — all batches so far\n')
    out.append('Per sector: admissible dimension `dimP` (= dimK − F-multiples), total restarts (all jobs/seeds), best and median relative residual, number of restarts below 1e-6 / 1e-10, and the distinct minima values (rounded, count). A solution has rel ≤ 1e-16.\n')
    out.append('| sector | nD | dimK (J+ker) | dimP | rows×cols | restarts | best rel | median | <1e-6 | <1e-10 | minima (value: count) | jobs |')
    out.append('|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---|---|')
    g = collections.defaultdict(list)
    for r in V: g[(r['chi'], r['d'], r['pD'])].append(r)
    def skey(k):
        order = {'jac': 0, 'jacJ': 0, 'triv': 1, 'sgnX': 2, 'sgnY': 3, 'sgnsgn': 4}
        return (order.get(k[0], 9), k[1], k[2], k[0])
    for key in sorted(g, key=skey):
        rs = g[key]; n = sum(r.get('n', 0) for r in rs); best = min(r['best_rel'] for r in rs)
        med = sorted(r.get('median_rel', r['best_rel']) for r in rs)[len(rs) // 2]
        b6 = sum(r.get('n_below_1e-6', 0) for r in rs); b10 = sum(r.get('n_below_1e-10', 0) for r in rs)
        mins = collections.Counter()
        for r in rs:
            for val, cnt in r.get('top_minima', []): mins[val] += cnt
        mtxt = ', '.join(f'{v}: {c}' for v, c in sorted(mins.items(), key=lambda kv: float(kv[0]))[:5])
        r0 = rs[0]
        out.append(f"| {key[0]} d={key[1]} pD={key[2]} | {r0.get('nD')} | {r0.get('dimK')} ({r0.get('dimJ')}+{r0.get('ker_beyond_J')}) | {r0.get('dimP')} | {r0.get('rows')}×{(r0.get('nN') or 0) + (r0.get('nQ') or 0)} | {n} | {best:.2e} | {med:.1e} | {b6} | {b10} | {mtxt} | {','.join(sorted(set(r['job'] for r in rs)))} |")
    conf = [r for r in V if r.get('confirmations')]
    out.append(f'\nCandidates with rel < 1e-10 submitted to exact confirmation: {sum(len(r["confirmations"]) for r in conf)}; verified witnesses: {sum(1 for r in conf for c in r["confirmations"] if c.get("ok"))}.')
# ---- partial (in-progress) sectors from worker files
import glob, re
done_keys = set((r['chi'], r['d'], r['pD'], r['job']) for r in V) if V else set()
wf = collections.defaultdict(list)
for pat in (f'{W2B}/*_w*.jsonl', f'/tmp/hodge/remote/{LANE}/w2b/*_w*.jsonl'):
    for fn in glob.glob(pat):
        m = re.match(r'.*/(\w+)_([A-Za-z0-9+\-]+)-(\d+)-(\d)_w\d+\.jsonl$', fn)
        if not m: continue
        job, chi, d, pD = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        if job.startswith('local'): continue
        for line in open(fn):
            try: wf[(chi, d, pD, job)].append(json.loads(line))
            except Exception: pass
partial = {k: v for k, v in wf.items() if k not in done_keys and v}
if partial:
    out.append('\n### 7a. Sectors still running (partial worker output, deduplicated by file)\n')
    out.append('| sector | job | restarts so far | best rel | median | <1e-6 | <1e-10 | minima (value: count) |')
    out.append('|---|---|---:|---:|---:|---:|---:|---|')
    for k in sorted(partial, key=lambda k: (k[0], k[1], k[2])):
        rs = partial[k]; rels = sorted(r['rel'] for r in rs if not r.get('spurious')) or sorted(r['rel'] for r in rs)
        mins = collections.Counter(f'{r:.3g}' for r in rels)
        mtxt = ', '.join(f'{v}: {c}' for v, c in sorted(mins.items(), key=lambda kv: float(kv[0]))[:5])
        out.append(f"| {k[0]} d={k[1]} pD={k[2]} | {k[3]} | {len(rs)} | {rels[0]:.2e} | {rels[len(rels)//2]:.1e} | {sum(1 for r in rels if r < 1e-6)} | {sum(1 for r in rels if r < 1e-10)} | {mtxt} |")
# ---- single solves
Sg = [r for r in recs if r.get('family') in ('kernel-class', 'kernel-random-combo', 'structured-single')]
if Sg:
    out.append('\n## 8. Single-denominator exact solves (character-reduced systems mod p, CRT lift on consistency)\n')
    c = collections.Counter((r['family'], r.get('e1'), r.get('consistent')) for r in Sg)
    out.append('| family | E1 pass & consistent | E1 pass & inconsistent | E1 fails | errors |')
    out.append('|---|---:|---:|---:|---:|')
    for fam in ('kernel-class', 'kernel-random-combo', 'structured-single'):
        rr = [r for r in Sg if r['family'] == fam]
        if not rr: continue
        out.append(f"| {fam} | {sum(1 for r in rr if r.get('consistent'))} | {sum(1 for r in rr if r.get('e1') and r.get('consistent') is False)} | {sum(1 for r in rr if r.get('e1') is False)} | {sum(1 for r in rr if 'error' in r)} |")
    cons = [r for r in Sg if r.get('consistent')]
    if cons:
        out.append('\n**Consistent single denominators:** ' + '; '.join(f"{r.get('label')} (d={r.get('d')}, {r.get('chi')}, lift={r.get('lift_status')}, verified={r.get('verified')})" for r in cons))
    st = [r for r in Sg if r['family'] == 'structured-single']
    if st:
        out.append('\nStructured candidates: ' + '; '.join(f"`{r['label']}` d={r.get('d')} {r.get('chi')}: " + ('E1 fails' if r.get('e1') is False else ('consistent' if r.get('consistent') else f"inconsistent (rank {r.get('rank')}, null {r.get('nullity')}, {r.get('rows')}×{r.get('cols')})")) for r in sorted(st, key=lambda r: (r.get('d') or 0, r['label']))))
    kc = [r for r in Sg if r['family'] in ('kernel-class', 'kernel-random-combo')]
    if kc:
        bys = collections.defaultdict(lambda: [0, 0, 0])
        for r in kc:
            k = (r.get('chi'), r.get('d'), r.get('pD')); bys[k][0] += 1
            if r.get('consistent'): bys[k][1] += 1
            if 'error' in r: bys[k][2] += 1
        out.append('\nKernel classes / combos per sector (tested, consistent, errors): ' + '; '.join(f'{k[0]} d={k[1]} pD={k[2]}: {v[0]}/{v[1]}/{v[2]}' for k, v in sorted(bys.items(), key=lambda kv: (str(kv[0][0]), kv[0][1] or 0, kv[0][2] or 0))))
# ---- structured families
St = [r for r in recs if r.get('family') == 'structured-family']
if St:
    out.append('\n## 9. Structured vertex families (cones with G-stable vertex; restricted E1 kernel + varpro)\n')
    out.append('| family | χ | d | pD | monomials | nD | dimK | dimP | restarts | best rel | single solves (consistent/total) | note |')
    out.append('|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|')
    for r in sorted(St, key=lambda r: (r.get('fam'), r.get('d'), r.get('chi'), r.get('pD'))):
        ss = r.get('single_solves') or []
        out.append(f"| {r.get('fam')} | {r.get('chi')} | {r.get('d')} | {r.get('pD')} | {r.get('n_monomials')} | {r.get('nD')} | {r.get('dimK')} | {r.get('dimP')} | {r.get('restarts', '')} | {r['best_rel']:.2e} | {sum(1 for s in ss if s.get('consistent'))}/{len(ss)} | {r.get('note', r.get('error', ''))[:60]} |" if 'best_rel' in r else
                   f"| {r.get('fam')} | {r.get('chi')} | {r.get('d')} | {r.get('pD')} | {r.get('n_monomials')} | {r.get('nD')} | {r.get('dimK')} | {r.get('dimP')} | | | | {r.get('note', r.get('error', ''))[:60]} |")
# ---- jobs
J = collections.defaultdict(lambda: [0, 0.0])
for r in recs:
    if r.get('job'): J[r['job']][0] += 1; J[r['job']][1] += float(r.get('seconds') or 0)
out.append('\n## 10. Batches logged per job\n')
out.append(', '.join(f'`{j}`: {v[0]} records, {v[1]/3600:.1f} h' for j, v in sorted(J.items())))
open(REPORT, 'w').write('\n'.join(out) + '\n')
print('report written', REPORT, len(recs), 'records; witness:', wit)
