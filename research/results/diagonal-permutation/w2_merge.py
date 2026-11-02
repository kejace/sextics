"""Merge remote wave-2 JSONL outputs (mirror /tmp/hodge/remote/diagonal-permutation/w2-*.jsonl) into the local
wave2-search.log.jsonl (dedup by line) and print per-file counts.  Per-restart lines (w2-pt-d*.jsonl) are NOT merged
line by line (too many); instead one summary line per (d, tag) with counts/best/histogram is appended (dedup by key)."""
import os, json, glob, math
MIR = '/tmp/hodge/remote/diagonal-permutation'
LOCAL = '/tmp/hodge/results/diagonal-permutation/wave2-search.log.jsonl'
have = set()
if os.path.exists(LOCAL):
    have = set(l.rstrip('\n') for l in open(LOCAL) if l.strip())
keys = set()
for l in have:
    try:
        o = json.loads(l)
        if o.get('exp') == 'pt-batch': keys.add((o['d'], o['tag']))
    except Exception: pass
added = 0
with open(LOCAL, 'a') as out:
    for fn in sorted(glob.glob(os.path.join(MIR, 'w2-*.jsonl'))):
        base = os.path.basename(fn)
        if base.startswith('w2-near-'):
            continue
        if base.startswith('w2-pt-d'):
            rels = []; its = []; secs = []; fams = set(); d = None; tag = None
            for l in open(fn):
                try: o = json.loads(l)
                except Exception: continue
                rels.append(o['rel']); its.append(o['iters']); secs.append(o['seconds']); fams.add(o['family']); d = o['d']; tag = o['tag']
            if not rels: continue
            hist = [0] * 11
            edges = [-300, -14, -12, -10, -8, -6, -4, -3, -2, -1, 0, 10]
            for r in rels:
                v = math.log10(max(r, 1e-300))
                for i in range(11):
                    if edges[i] <= v < edges[i + 1]: hist[i] += 1; break
            line = {'exp': 'pt-batch', 'd': d, 'tag': tag, 'family': sorted(fams), 'restarts': len(rels), 'best_rel': min(rels),
                    'median_rel': sorted(rels)[len(rels) // 2], 'mean_iters': round(sum(its) / len(its), 1), 'mean_seconds': round(sum(secs) / len(secs), 1),
                    'hist_log10': hist, 'hist_bins': '[<-14,-12,-10,-8,-6,-4,-3,-2,-1,0,>0]', 'source': base}
            s = json.dumps(line)
            # replace older summary for the same key: we simply append (latest wins when reading)
            out.write(s + '\n'); added += 1
            print(f'{base}: {len(rels)} restarts, best {min(rels):.3e}, median {line["median_rel"]:.3e}, hist {hist}')
        else:
            n = 0
            for l in open(fn):
                l = l.rstrip('\n')
                if not l.strip(): continue
                try: o = json.loads(l)
                except Exception: continue
                o['source'] = base
                s = json.dumps(o)
                if s in have: continue
                out.write(s + '\n'); have.add(s); added += 1; n += 1
            print(f'{base}: {n} new lines')
print('added', added)
