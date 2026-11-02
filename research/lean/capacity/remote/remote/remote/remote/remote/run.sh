#!/usr/bin/env bash
# Usage: run.sh <Module.Name> [label]   -- source at src/<Module/Name>.lean
# Mimics scripts/proposals/sandbox.py: lean -j 2 -M 14336, then leanchecker replay.
set -u
export PATH="$(cat /tmp/hodge/lean/capacity/prefix.txt)/bin:$PATH"
CAP=/tmp/hodge/lean/capacity
MOD="$1"; LABEL="${2:-$1}"
REL="${MOD//.//}.lean"
export LEAN_PATH="$(cat $CAP/leanpath.txt):$CAP/out:$CAP/src"
mkdir -p "$CAP/out/$(dirname "$REL")" "$CAP/logs"
cd "$CAP/src"
LOG="$CAP/logs/${MOD}.compile.log"; KLOG="$CAP/logs/${MOD}.kernel.log"
/usr/bin/time -f "WALL=%e USER=%U SYS=%S MAXRSS_KB=%M" -o "$LOG.time" \
  lean -j 2 -M 14336 -o "$CAP/out/${REL%.lean}.olean" "$REL" > "$LOG" 2>&1
CRC=$?
/usr/bin/time -f "WALL=%e USER=%U SYS=%S MAXRSS_KB=%M" -o "$KLOG.time" \
  leanchecker --verbose "$MOD" > "$KLOG" 2>&1
KRC=$?
python3 - "$MOD" "$LABEL" "$CRC" "$KRC" "$LOG" "$KLOG" <<'PY'
import sys,json,re,time
mod,label,crc,krc,log,klog=sys.argv[1:]
def parse(p):
    try: t=open(p+'.time').read()
    except: return {}
    m=re.search(r'WALL=([\d.]+) USER=([\d.]+) SYS=([\d.]+) MAXRSS_KB=(\d+)',t)
    return {'wall':float(m[1]),'user':float(m[2]),'sys':float(m[3]),'maxrss_gb':round(int(m[4])/1e6,2)} if m else {}
errs=[l for l in open(log,errors='replace') if 'error' in l][:3]
rec={'module':mod,'label':label,'compile_rc':int(crc),'compile':parse(log),'kernel_rc':int(krc),'kernel':parse(klog),'errors':[e.strip()[:200] for e in errs],'ts':time.strftime('%H:%M:%S')}
open('/tmp/hodge/lean/capacity/timings.jsonl','a').write(json.dumps(rec)+'\n')
print(json.dumps(rec))
PY
