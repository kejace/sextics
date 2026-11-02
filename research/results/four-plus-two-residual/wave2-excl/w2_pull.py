"""Local helper: pull the remote lane mirror, merge all wave2-excl/remote/*.jsonl lines into the local lane log
/tmp/hodge/results/four-plus-two-residual/wave2-excl.log.jsonl (deduplicated), copy receipts into the packet, print a summary."""
import subprocess, json, os, glob, sys, collections
LANE = '/tmp/hodge/results/four-plus-two-residual'
MIR = '/tmp/hodge/remote/four-plus-two-residual/wave2-excl/remote'
LOG = f'{LANE}/wave2-excl.log.jsonl'
if '--nopull' not in sys.argv:
    r = subprocess.run(['python3', '-m', 'hodge.rjob', 'pull', '--lane', 'four-plus-two-residual'], cwd='/home/kejace/devel/sextics-work',
                       capture_output=True, text=True, timeout=600)
    print(r.stdout.strip()[-300:], r.stderr.strip()[-300:] if 'engine v2' not in r.stderr else '')
seen = set()
if os.path.exists(LOG):
    for line in open(LOG):
        seen.add(line.strip())
new = 0
os.makedirs(f'{LANE}/wave2-excl/receipts', exist_ok=True)
with open(LOG, 'a') as out:
    for f in sorted(glob.glob(f'{MIR}/*.jsonl')):
        for line in open(f):
            s = line.strip()
            if s and s not in seen:
                out.write(s + '\n'); seen.add(s); new += 1
    for f in glob.glob(f'{MIR}/*.json') + glob.glob(f'{MIR}/*.sing') + glob.glob(f'{MIR}/*.jsonl'):
        subprocess.run(['cp', f, f'{LANE}/wave2-excl/receipts/'])
print('merged new lines:', new, 'total:', len(seen))
rows = [json.loads(l) for l in open(LOG) if l.strip()]
cnt = collections.Counter((r.get('job'), r.get('exp')) for r in rows)
for k, v in sorted(cnt.items(), key=str): print(f'{v:5d}  {k}')
res = collections.Counter((r.get('job'), str(r.get('result'))[:60]) for r in rows if 'result' in r)
for k, v in sorted(res.items(), key=str): print(f'{v:5d}  {k}')
