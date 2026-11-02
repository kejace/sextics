#!/bin/bash
# pull the remote lane dir and merge remote log lines (both the legacy remote copy of wave2-search.log.jsonl and the remote-only file)
cd /home/kejace/devel/sextics-work && timeout 600 python3 -m hodge.rjob pull --lane dense-rational > /tmp/hodge/results/dense-rational/pull.log 2>&1
python3 - <<'PY'
import os
loc='/tmp/hodge/results/dense-rational/wave2-search.log.jsonl'
have=set(open(loc).read().splitlines()) if os.path.exists(loc) else set()
added=0
for fn in ['/tmp/hodge/remote/dense-rational/wave2-search.log.jsonl','/tmp/hodge/remote/dense-rational/wave2-search.remote.log.jsonl']:
    if os.path.exists(fn):
        with open(loc,'a') as out:
            for line in open(fn).read().splitlines():
                if line and line not in have:
                    out.write(line+'\n'); have.add(line); added+=1
print('merged', added, 'new log lines')
PY
