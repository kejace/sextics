"""R1: E1 kernel table Ann([A])_d mod J_d, all 18 sectors, exact over Q (d <= 8) and mod a fresh prime (d <= 11)."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-referee')
import ref_lib as L
P = 999999937
LOG = '/tmp/hodge/results/cyclic-triple-residual/wave2-referee/referee.log.jsonl'
def log(rec):
    rec['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S'); rec['script'] = 'ref_e1.py'
    with open(LOG, 'a') as f: f.write(json.dumps(rec) + '\n')
L.check_setup()
table = {}
for d in range(0, 12):
    for mode in (['exact', P] if d <= 8 else [P]):
        t0 = time.time()
        r = L.e1_kernel(d, p=None if mode == 'exact' else mode)
        r['seconds'] = round(time.time() - t0, 1)
        nz = {k: v['ann'] for k, v in r['sectors'].items() if v['ann']}
        print(f"d={d} {r['mode']}: Ann mod J = {r['ann_total']} (kerA {r['kerA_total']}, dimJ_d {r['dimJ_d']}={r['ci_dimJ_d']}, dimJ_{d+6} {r['dimJ_d6']}={r['ci_dimJ_d6']}, ci_ok {r['ci_ok']}) sectors {nz} {r['seconds']}s", flush=True)
        log({'check': f'e1_kernel_d{d}', 'mode': r['mode'], 'ann_mod_J': r['ann_total'], 'kerA': r['kerA_total'], 'dimJ_d': r['dimJ_d'],
             'ci_dimJ_d': r['ci_dimJ_d'], 'dimJ_d6': r['dimJ_d6'], 'ci_dimJ_d6': r['ci_dimJ_d6'], 'ci_ok': r['ci_ok'], 'nonzero_sectors': nz, 'seconds': r['seconds']})
        table[f'{d}:{r["mode"]}'] = r
        json.dump(table, open('/tmp/hodge/results/cyclic-triple-residual/wave2-referee/e1_table.json', 'w'), indent=1)
print('done', flush=True)
