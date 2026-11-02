#!/usr/bin/env python3
"""Local (light): pull the remote lane dir, merge remote batch-log lines into wave2-search.log.jsonl, print a results table.
  python3 w2_report.py [--no-pull]"""
import os, sys, json, subprocess, glob, time
HERE = os.path.dirname(os.path.abspath(__file__)); LANE = os.path.dirname(HERE)
LOCAL_LOG = os.path.join(LANE, 'wave2-search.log.jsonl')
MIRROR = '/tmp/hodge/remote/off-paired/wave2-search'
if '--no-pull' not in sys.argv:
    subprocess.run([sys.executable, '-m', 'hodge.rjob', 'pull', '--lane', 'off-paired'], cwd='/home/kejace/devel/sextics-work',
                   capture_output=True, text=True, timeout=600)
# merge remote batch log
have = set(l.strip() for l in open(LOCAL_LOG)) if os.path.exists(LOCAL_LOG) else set()
added = 0
rl = os.path.join(MIRROR, w2lib.REMOTE_LOG)
if os.path.exists(rl):
    with open(LOCAL_LOG, 'a') as f:
        for line in open(rl):
            line = line.strip()
            if line and line not in have:
                f.write(line + '\n'); have.add(line); added += 1
print(f'merged {added} new remote batch lines into {LOCAL_LOG}')
# summaries
rows = []
for fn in sorted(glob.glob(os.path.join(MIRROR, 'wave2-*.summary.json'))):
    try:
        s = json.load(open(fn))
    except Exception:
        continue
    rows.append(s)
rows.sort(key=lambda s: (s.get('d', 0), str(s.get('sector')), s.get('family', '')))
print(f"{'job':22s} {'step':18s} {'d':>2s} {'sector':8s} {'grp':3s} {'fam':12s} {'K':>4s} {'dimP':>4s} {'rows':>6s} {'rst':>5s} {'best':>9s} {'median':>9s} {'t/rst':>6s} {'hits<1e-9':>9s} wit")
for s in rows:
    print(f"{s.get('job',''):22s} {s.get('step',''):18s} {s.get('d',''):>2} {str(tuple(s.get('sector',[]))):8s} {s.get('group','G'):3s} {s.get('family',''):12s} "
          f"{s.get('K', s.get('K_w','')):>4} {s.get('dimP',''):>4} {s.get('rows',''):>6} {s.get('restarts',''):>5} "
          f"{(s.get('best_relres') if s.get('best_relres') is not None else float('nan')):9.2e} {(s.get('median_relres') if s.get('median_relres') is not None else float('nan')):9.2e} "
          f"{str(s.get('mean_time_per_restart', s.get('mean_time',''))):>6s} {s.get('n_numerical_hits_1e-9',''):>9} {s.get('witness', s.get('served_hit',''))}")
# other steps from the merged log
for line in open(LOCAL_LOG):
    try: r = json.loads(line)
    except Exception: continue
    if r.get('step') in ('control-fermat', 'homotopy', 'structured', 'hess'):
        print({k: v for k, v in r.items() if k not in ('ts',)} if r.get('step') != 'homotopy' else {'step': 'homotopy', 'job': r.get('job'), 'd': r.get('d'), 'path': r.get('relres_path')})
