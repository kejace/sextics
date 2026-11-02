"""Step 1: E1-admissible denominator spaces per degree and character sector (exact over Q).
Writes wave2-dims.jsonl (one line per (d, chi)) and one summary batch line to the lane log."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/off-paired')
import w2lib as W
F, A, q = W.load()
d0, d1 = int(sys.argv[1]), int(sys.argv[2])
out = open('/tmp/hodge/results/off-paired/wave2-dims.jsonl', 'a')
for d in range(d0, d1 + 1):
    t0 = time.time()
    tab = {}
    for par in (0, 1):
        for c in W.all_chars_parity(par):
            S = W.SectorSpaces(F, A, q, d, c)
            s = S.summary(); s['parity'] = par
            tab[str(c)] = s
            out.write(json.dumps(s) + '\n'); out.flush()
            print(f"d={d} chi={c} par={par} R={s['dimR']} E1={s['dimE1']} gauge={s['gauge']} jac={s['jac_mod_gauge']} ker={s['ker']} ({s['seconds']}s)", flush=True)
    totE1 = sum(v['dimE1'] for v in tab.values()); totK = sum(v['ker'] for v in tab.values())
    totJ = sum(v['jac_mod_gauge'] for v in tab.values()); totG = sum(v['gauge'] for v in tab.values())
    W.log_batch({'step': 'dims', 'd': d, 'sectors': 36, 'dimE1_total': totE1, 'gauge_total': totG,
                 'jac_mod_gauge_total': totJ, 'ker_total': totK,
                 'nonzero_sectors': sum(1 for v in tab.values() if v['dimE1'] - v['gauge'] > 0),
                 'seconds': round(time.time() - t0, 1), 'detail_file': 'wave2-dims.jsonl'})
