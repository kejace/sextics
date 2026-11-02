"""Summarize all pointwise-search batches from the mirror (/tmp/hodge/remote/diagonal-permutation/w2-pt-d*.jsonl):
per (d, tag): restarts, best/median relative residual, best absolute-scaled residual, iteration stats, time/restart, histogram."""
import glob, json, math, os, statistics
MIR = '/tmp/hodge/remote/diagonal-permutation'
rows = []
for fn in sorted(glob.glob(os.path.join(MIR, 'w2-pt-d*.jsonl'))):
    recs = []
    for l in open(fn):
        try: recs.append(json.loads(l))
        except Exception: pass
    if not recs: continue
    rels = [r['rel'] for r in recs]; ra = [r.get('rel_abs', float('nan')) for r in recs]
    its = [r['iters'] for r in recs]; secs = [r['seconds'] for r in recs]
    fam = recs[0]['family']
    edges = [-300, -14, -12, -10, -8, -6, -4, -3, -2, -1, 0, 10]; hist = [0] * 11
    for v in rels:
        lv = math.log10(max(v, 1e-300))
        for i in range(11):
            if edges[i] <= lv < edges[i + 1]: hist[i] += 1; break
    rows.append({'file': os.path.basename(fn), 'd': recs[0]['d'], 'tag': recs[0]['tag'], 'family': fam, 'restarts': len(recs),
                 'best_rel': min(rels), 'median_rel': statistics.median(rels), 'best_rel_abs': min(x for x in ra if x == x) if any(x == x for x in ra) else None,
                 'iters_median': statistics.median(its), 'iters_max': max(its), 'sec_per_restart': round(statistics.mean(secs), 1), 'hist_log10': hist})
for r in rows:
    print(f"d={r['d']} {r['tag']:6s} n={r['restarts']:5d} best={r['best_rel']:.3e} med={r['median_rel']:.3e} best_abs={r['best_rel_abs']} iters~{r['iters_median']} (max {r['iters_max']}) {r['sec_per_restart']}s/restart hist={r['hist_log10']}  {r['family']}")
json.dump(rows, open('/tmp/hodge/results/diagonal-permutation/w2-batch-stats.json', 'w'), indent=1)
