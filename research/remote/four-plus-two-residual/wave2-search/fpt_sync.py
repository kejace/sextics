"""Pull the remote lane dir and merge all wave-2 search records (mine: wave2-search/*.jsonl on the remote mirror;
plus any pre-existing lines of the shared wave2-search.log.jsonl) into the lane's wave2-search.log.jsonl (deduplicated).
usage: python3 fpt_sync.py [--no-pull]"""
import os, sys, json, subprocess, glob
LANE = '/tmp/hodge/results/four-plus-two-residual'
MIRROR = '/tmp/hodge/remote/four-plus-two-residual'
LOG = f'{LANE}/wave2-search.log.jsonl'
if '--no-pull' not in sys.argv:
    r = subprocess.run(['python3', '-m', 'hodge.rjob', 'pull', '--lane', 'four-plus-two-residual'], cwd='/home/kejace/devel/sextics-work',
                       capture_output=True, text=True, timeout=600)
    print(r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-300:])
have = set()
if os.path.exists(LOG):
    have = set(l.rstrip('\n') for l in open(LOG) if l.strip())
new = []
sources = [f'{MIRROR}/wave2-search.log.jsonl'] + sorted(glob.glob(f'{MIRROR}/wave2-search/*.jsonl')) + sorted(glob.glob(f'{LANE}/wave2-search/*.jsonl')) + sorted(glob.glob(f'{MIRROR}/wave2-ctl*.jsonl')) + sorted(glob.glob(f'{MIRROR}/wave2-kJ*.jsonl')) + sorted(glob.glob(f'{MIRROR}/wave2-dims*.jsonl'))
for src in sources:
    if not os.path.exists(src):
        continue
    tag = os.path.basename(src)
    for line in open(src):
        line = line.rstrip('\n')
        if not line.strip():
            continue
        if line in have:
            continue
        try:
            o = json.loads(line)
        except Exception:
            continue
        if 'src' not in o:
            o['src'] = tag
        line2 = json.dumps(o, default=str)
        if line2 in have:
            continue
        have.add(line2); have.add(line)
        new.append(line2)
with open(LOG, 'a') as f:
    for l in new:
        f.write(l + '\n')
print(f'merged {len(new)} new lines into {LOG} (total {sum(1 for _ in open(LOG))})')
