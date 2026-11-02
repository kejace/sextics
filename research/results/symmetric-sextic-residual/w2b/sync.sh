#!/bin/bash
# pull remote lane dir, merge job logs into wave2-search.log.jsonl, regenerate the report
cd /home/kejace/devel/sextics-work && timeout 600 python3 -m hodge.rjob pull --lane symmetric-sextic-residual 2>&1 | grep -v "hodge engine" | tail -1
cd /tmp/hodge/results/symmetric-sextic-residual && python3 - <<'PY'
import glob, os
LANE='symmetric-sextic-residual'; RES=f'/tmp/hodge/results/{LANE}'; MIR=f'/tmp/hodge/remote/{LANE}'
main=f'{RES}/wave2-search.log.jsonl'
have=set(open(main).read().splitlines()) if os.path.exists(main) else set()
added=0
with open(main,'a') as out:
    for pat in (f'{RES}/w2/*.log.jsonl', f'{RES}/w2b/*.log.jsonl', f'{MIR}/w2/*.log.jsonl', f'{MIR}/w2b/*.log.jsonl'):
        for fn in sorted(glob.glob(pat)):
            for line in open(fn):
                line=line.rstrip('\n')
                if line and 'w2lib-v1' in line and line not in have:
                    out.write(line+'\n'); have.add(line); added+=1
print('merged', added, 'new lines; total', len(have))
# witness candidate from remote?
if os.path.exists(f'{MIR}/witness-candidate.json') and not os.path.exists(f'{RES}/witness-candidate.json'):
    import shutil; shutil.copy(f'{MIR}/witness-candidate.json', f'{RES}/witness-candidate.json'); print('*** witness-candidate.json pulled from remote ***')
PY
python3 w2b/report.py
