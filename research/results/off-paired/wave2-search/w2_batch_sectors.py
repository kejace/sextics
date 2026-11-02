#!/usr/bin/env python3
"""Run w2_varpro.py over many (d, sector, kind) specs sequentially inside one remote job.
  python3 w2_batch_sectors.py --d 9 --sectors all|h0|LIST --kind colon --restarts 400 --workers 8 --seed0 0 --maxiter 120 --job d9all [--per-time 600]
  sectors: 'all' = 36 chars ordered by eff (largest first); 'h0' = (a,0); 'even' = Y-even (a even); explicit 'a,b;a,b;...'
"""
import os, sys, json, subprocess, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True); ap.add_argument('--sectors', default='all'); ap.add_argument('--kind', default='colon')
ap.add_argument('--restarts', type=int, default=400); ap.add_argument('--workers', type=int, default=8); ap.add_argument('--seed0', type=int, default=0)
ap.add_argument('--maxiter', type=int, default=120); ap.add_argument('--job', required=True); ap.add_argument('--per-time', type=float, default=900)
ap.add_argument('--opt', default='lm'); ap.add_argument('--F', default='offpaired'); ap.add_argument('--min-eff', type=int, default=1)
ap.add_argument('--free-target', action='store_true'); ap.add_argument('--group', default='G')
args = ap.parse_args()
dims = json.load(open('dims.json'))
if args.sectors == 'all':
    secs = sorted([tuple(v['c']) for v in dims.values() if v['d'] == args.d], key=lambda c: -dims[f'{args.d}:{c[0]},{c[1]}']['eff'])
elif args.sectors == 'h0':
    secs = [(a, 0) for a in range(6)]
elif args.sectors == 'even':
    secs = sorted([tuple(v['c']) for v in dims.values() if v['d'] == args.d and v['c'][0] % 2 == 0], key=lambda c: -dims[f'{args.d}:{c[0]},{c[1]}']['eff'])
else:
    secs = [tuple(int(x) for x in s.split(',')) for s in args.sectors.split(';')]
T0 = time.time()
for c in secs:
    eff = dims.get(f'{args.d}:{c[0]},{c[1]}', {}).get('eff', None) if args.group == 'G' else None
    if eff is not None and eff < args.min_eff:
        continue
    job = f'{args.job}-c{c[0]}{c[1]}'
    if os.path.exists(f'wave2-{job}.summary.json'):
        print(f'skip {job} (done)', flush=True); continue
    cmd = [sys.executable, 'w2_varpro.py', '--d', str(args.d), '--c', f'{c[0]},{c[1]}', '--kind', args.kind, '--restarts', str(args.restarts),
           '--workers', str(args.workers), '--seed0', str(args.seed0), '--maxiter', str(args.maxiter), '--time-limit', str(args.per_time),
           '--job', job, '--opt', args.opt, '--F', args.F] + (['--free-target'] if args.free_target else []) + ['--group', args.group]
    print(f'[{time.strftime("%H:%M:%S")}] start {job} (eff {eff}) elapsed {time.time()-T0:.0f}s', flush=True)
    with open(f'{job}.out', 'w') as fo:
        subprocess.run(cmd, stdout=fo, stderr=subprocess.STDOUT)
    try:
        s = json.load(open(f'wave2-{job}.summary.json'))
        print(f'   -> best {s["best_relres"]:.3e} median {s["median_relres"]:.3e} restarts {s["restarts"]} dimP {s["dimP"]} t/restart {s["mean_time_per_restart"]} witness {s["witness"]}', flush=True)
        if s['witness']:
            print('!!!!!!!! WITNESS FOUND !!!!!!!!', flush=True)
    except Exception as ex:
        print('   -> no summary:', ex, flush=True)
print('batch done', flush=True)
