"""Pull the remote lane dir and merge all w2-*.jsonl lines into the local wave2-search.log.jsonl (dedup)."""
import subprocess, glob, json, os, sys
LANE = 'polarized-k3-remainder'
subprocess.run([sys.executable, '-m', 'hodge.rjob', 'pull', '--lane', LANE], cwd='/home/kejace/devel/sextics-work', capture_output=True, text=True, timeout=600)
local = f'/tmp/hodge/results/{LANE}/wave2-search.log.jsonl'
have = set()
if os.path.exists(local):
    have = set(l.rstrip('\n') for l in open(local) if l.strip())
added = 0
with open(local, 'a') as out:
    for fn in sorted(glob.glob(f'/tmp/hodge/remote/{LANE}/w2-*.jsonl') + glob.glob(f'/tmp/hodge/results/{LANE}/w2-*.jsonl')):
        for line in open(fn):
            line = line.rstrip('\n')
            if line and line not in have:
                out.write(line + '\n'); have.add(line); added += 1
print('merged', added, 'new lines; total', len(have))
