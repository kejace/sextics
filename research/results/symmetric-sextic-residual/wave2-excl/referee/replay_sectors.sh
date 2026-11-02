#!/bin/bash
# Replay msolve (-g 1, 1 thread) on every stored nonempty sector system at d=7,8 and compare pure-power sets with the stored .lead
P=/tmp/hodge/results/symmetric-sextic-residual/wave2-excl
for f in $P/relax_sector_{triv,sgnX,sgnY}_d{7,8}_par{0,1}_p2147483647.ms; do
  b=$(basename $f .ms); [ -s $f ] || { echo "$b: empty .ms (empty family)"; continue; }
  s=$(date +%s.%N); timeout 900 nice -n 15 msolve -f $f -o /tmp/hodge/referee-ssr/$b.lead -g 1 -t 1 >/tmp/hodge/referee-ssr/$b.msolve.log 2>&1; rc=$?; e=$(date +%s.%N)
  python3 -I - "$f" "$P/$b.lead" "/tmp/hodge/referee-ssr/$b.lead" "$b" "$rc" "$(echo "$e - $s" | bc)" <<'PY'
import re,sys
ms,lead0,lead1,b,rc,secs=sys.argv[1:7]
vars_=open(ms).readline().strip().split(',')
def pure(path):
    txt=open(path).read(); found={}
    for tok in re.split(r'[,\n\[\]]',txt):
        tok=tok.strip()
        m=re.fullmatch(r'(\d+\*)?([A-Za-z_]\w*)(\^(\d+))?',tok)
        if m: v=m.group(2); k=int(m.group(4)) if m.group(4) else 1; found[v]=min(found.get(v,99),k)
    return found
def allmons(path):
    txt=open(path).read().split('\n'); body=[l for l in txt if not l.startswith('#')]
    return sorted(set(t.strip() for t in re.split(r'[,\[\]\n]', '\n'.join(body)) if t.strip()))
p0,p1=pure(lead0),pure(lead1)
same=allmons(lead0)==allmons(lead1)
print(f"{b}: rc={rc} {float(secs):.1f}s vars={len(vars_)} replay pure powers {sum(v in p1 for v in vars_)}/{len(vars_)} (stored {sum(v in p0 for v in vars_)}/{len(vars_)}) identical leading ideal={same} n_lead={len(allmons(lead1))}")
PY
done
