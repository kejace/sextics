#!/usr/bin/env python3
"""Cone / coordinate-subspace families: for every (variable subset, sector) pair whose E1 kernel exceeds the
gauge (table e1_sectors_wave2.json written by w2_e1profile.py --cones), run the full-E2 Koszul varpro search
over that (small) family.   python3 w2_cones.py [--dmin 6] [--dmax 12] [--restarts 60]
"""
import sys, os, json, argparse
os.environ['OPENBLAS_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/quadric-contact-residual')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import w2lib as W
import w2_varpro as V

ap = argparse.ArgumentParser()
ap.add_argument('--dmin', type=int, default=6); ap.add_argument('--dmax', type=int, default=12)
ap.add_argument('--restarts', type=int, default=60); ap.add_argument('--procs', type=int, default=8)
ap.add_argument('--seed', type=int, default=9000); ap.add_argument('--tag', default='cone')
ap.add_argument('--table', default=W.RES + '/e1_sectors_wave2.json')
a = ap.parse_args()
T = json.load(open(a.table))
n = 0
for d in range(a.dmin, a.dmax + 1):
    cones = T.get(str(d), {}).get('cones', {})
    for key, r in sorted(cones.items()):
        if r['kernel'] <= 0:   # masked families contain no F-multiples (F involves all six variables)
            continue
        name, cs = key.split(':')
        c = W.parse_chi(cs)
        vars_mask = [v in [name[i:i + 2] for i in range(0, len(name), 2)] for v in W.core.VARS]
        n += 1
        print(f'=== cone family {name} sector {cs} d={d}: kernel {r["kernel"]} gauge {r["gauge"]}', flush=True)
        try:
            V.run_sector(d, c, None, False, 'basis', vars_mask, a.restarts, a.procs, a.seed + 13 * n, 300, a.tag, 1e-10, 10, 1200, 8, 12)
        except Exception as ex:
            import traceback; traceback.print_exc()
            W.log_line({'task': 'varpro', 'd': d, 'chi': cs, 'vars': name, 'tag': a.tag, 'status': 'family-error', 'error': str(ex)[:300]})
print('cones done', n)
