"""Pull the remote lane dir and merge wave2 batch lines into the lane's wave2-search.log.jsonl (dedup by content)."""
import json, os, subprocess, sys, glob
LANE = '/tmp/hodge/results/rank-one-biquartic'
MIR = '/tmp/hodge/remote/rank-one-biquartic'
LOG = f'{LANE}/wave2-search.log.jsonl'
if '--no-pull' not in sys.argv:
    subprocess.run(['python3', '-m', 'hodge.rjob', 'pull', '--lane', 'rank-one-biquartic'], cwd='/home/kejace/devel/sextics-work',
                   capture_output=True, text=True, timeout=600)
have = set()
if os.path.exists(LOG):
    for line in open(LOG):
        try:
            o = json.loads(line); o.pop('merged_from', None); have.add(json.dumps(o, sort_keys=True))
        except Exception:
            pass
added = 0
for path in sorted(glob.glob(f'{MIR}/wave2/out/*.jsonl')):
    name = os.path.basename(path)
    if name.endswith('.restarts.jsonl') or name.endswith('.cands.jsonl') or name.startswith('smoke') or name.startswith('planted'):
        continue
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        try:
            o = json.loads(line)
        except Exception:
            continue
        if str(o.get('tag', '')).startswith(('planted', 'smoke')):
            continue
        key = json.dumps(o, sort_keys=True)
        if key in have:
            continue
        have.add(key)
        o['merged_from'] = name
        with open(LOG, 'a') as fh:
            fh.write(json.dumps(o) + '\n')
        added += 1
print(f'merged {added} new lines into {LOG}')
