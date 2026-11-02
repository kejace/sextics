#!/usr/bin/env bash
# Usage: run_module.sh <Module.Name> <label> [cpuset]
# Mimics scripts/proposals/sandbox.py: `lean -j 2 -M 14336 -o <olean> <file>` pinned to 2 cores, then
# `leanchecker --verbose <Module>` pinned to the same cores. Source at $CAP/src/<Module/Name>.lean.
# Appends one JSON record to $CAP/results/timings.jsonl.
set -u
LOAD0=$(cut -d" " -f1 /proc/loadavg)
CAP=${CAP:-/tmp/hodge/lean/capacity}
ROOT=${ROOT:-$HOME/devel/sextics/formal/lean}
export PATH="$HOME/.elan/bin:$PATH"
MOD="$1"; LABEL="${2:-$1}"; CPUSET="${3:-}"
REL="${MOD//.//}.lean"
cd "$ROOT"
PREFIX="$(lake env lean --print-prefix)"
LP="$(lake env printenv LEAN_PATH)"
export PATH="$PREFIX/bin:$PATH"
export LEAN_PATH="$LP:$CAP/out"
mkdir -p "$CAP/out/$(dirname "$REL")" "$CAP/logs" "$CAP/results"
cd "$CAP/src"
PIN=""
if [ -n "$CPUSET" ]; then PIN="taskset -c $CPUSET"; fi
LOG="$CAP/logs/${LABEL}.compile.log"; KLOG="$CAP/logs/${LABEL}.kernel.log"
T0=$(date +%s.%N)
/usr/bin/time -f "WALL=%e USER=%U SYS=%S MAXRSS_KB=%M" -o "$LOG.time" \
  $PIN lean -j 2 -M 14336 -o "$CAP/out/${REL%.lean}.olean" "$REL" > "$LOG" 2>&1
CRC=$?
KRC=-1
if [ $CRC -eq 0 ]; then
  /usr/bin/time -f "WALL=%e USER=%U SYS=%S MAXRSS_KB=%M" -o "$KLOG.time" \
    $PIN leanchecker --verbose "$MOD" > "$KLOG" 2>&1
  KRC=$?
fi
LOAD1=$(cut -d" " -f1 /proc/loadavg)
python3 - "$MOD" "$LABEL" "$CRC" "$KRC" "$LOG" "$KLOG" "$CAP" "$LOAD0" "$LOAD1" <<'PY'
import sys,json,re,time,os
mod,label,crc,krc,log,klog,cap,load0,load1=sys.argv[1:]
def parse(p):
    try: t=open(p+'.time').read()
    except Exception: return {}
    m=re.search(r'WALL=([\d.]+) USER=([\d.]+) SYS=([\d.]+) MAXRSS_KB=(\d+)',t)
    return {'wall':float(m[1]),'user':float(m[2]),'sys':float(m[3]),'maxrss_gb':round(int(m[4])/1e6,2)} if m else {}
errs=[l for l in open(log,errors='replace') if 'error' in l][:4]
meta={}
mp=os.path.join(cap,'src',mod.replace('.','/')+'.json')
if os.path.exists(mp):
    try: meta=json.load(open(mp))
    except Exception: pass
rec={'module':mod,'label':label,'compile_rc':int(crc),'compile':parse(log),'kernel_rc':int(krc),'kernel':parse(klog),
     'errors':[e.strip()[:300] for e in errs],'meta':meta,'host':os.uname().nodename,'load0':float(load0),'load1':float(load1),'ts':time.strftime('%Y-%m-%dT%H:%M:%S')}
open(os.path.join(cap,'results','timings.jsonl'),'a').write(json.dumps(rec)+'\n')
print(json.dumps(rec))
PY
