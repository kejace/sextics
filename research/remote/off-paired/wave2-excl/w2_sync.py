"""Pull the remote lane dir and merge the remote wave2-excl.log.jsonl into the local one (set union), so that a later
`rjob submit` (which rsyncs local -> remote without --delete) never clobbers lines written on the remote by jobs or other agents."""
import subprocess, os, sys
sys.path.insert(0, '/home/kejace/devel/sextics-work')
LANE = '/tmp/hodge/results/off-paired'; MIR = '/tmp/hodge/remote/off-paired'
r = subprocess.run(['python3', '-m', 'hodge.rjob', 'pull', '--lane', 'off-paired'], cwd='/home/kejace/devel/sextics-work', capture_output=True, text=True, timeout=600)
print(r.stdout.strip()[-300:], r.stderr.strip()[-300:])
import glob
ll = f'{LANE}/wave2-excl.log.jsonl'
for rl in sorted(glob.glob(f'{MIR}/wave2-excl*.log.jsonl')):
    name = os.path.basename(rl)
    if os.path.exists(rl):
        have = set(open(ll).read().splitlines()) if os.path.exists(ll) else set()
        new = [l for l in open(rl).read().splitlines() if l and l not in have]
        with open(ll, 'a') as f:
            for l in new: f.write(l + '\n')
        print(f'merged {len(new)} new lines of {name} into local ({len(have) + len(new)} total)')
print('receipts stay in the mirror dir', MIR + '/wave2-excl/')
