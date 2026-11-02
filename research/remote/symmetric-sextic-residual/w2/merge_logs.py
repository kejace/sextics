"""Merge per-job logs (local w2/*.log.jsonl and the pulled remote mirror) into wave2-search.log.jsonl (dedupe)."""
import glob, os, json
LANE = 'symmetric-sextic-residual'
RES = f'/tmp/hodge/results/{LANE}'; MIR = f'/tmp/hodge/remote/{LANE}'
main = f'{RES}/wave2-search.log.jsonl'
have = set(open(main).read().splitlines()) if os.path.exists(main) else set()
added = 0
with open(main, 'a') as out:
    for pat in (f'{RES}/w2/*.log.jsonl', f'{MIR}/w2/*.log.jsonl'):
        for fn in sorted(glob.glob(pat)):
            for line in open(fn):
                line = line.rstrip('\n')
                if line and line not in have:
                    out.write(line + '\n'); have.add(line); added += 1
print(f'merged {added} new lines; total {len(have)}')
