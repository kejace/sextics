#!/usr/bin/env bash
# runs ON THE REMOTE: launch Kronecker scaling jobs (after the noncomputable fix); one nohup per job
CAP=/tmp/hodge/lean/capacity
CPU=60
launch() { # name cpuset args...
  local name=$1; local cs=$2; shift 2
  nohup setsid bash -c "bash $CAP/kjob.sh $name $cs $*" > /tmp/hodge/jobs/cap-$name.log 2>&1 < /dev/null &
  echo $! > /tmp/hodge/jobs/cap-$name.pid; echo launched $name pid $!
}
launch K2_m3_d10_k300  "$CPU,$((CPU+1))"  --F m3 --d 10 --kD 300 --kH 2 --kT 2 --parts e1,e2,homog; CPU=$((CPU+2))
launch K2_dn_d6_k40    "$CPU,$((CPU+1))"  --F dense --d 6 --kD 40 --kH 2 --parts e1,e2,homog; CPU=$((CPU+2))
launch K2_dn_d8_k100   "$CPU,$((CPU+1))"  --F dense --d 8 --kD 100 --kH 2 --parts e1,e2,homog; CPU=$((CPU+2))
launch K2_m3_d10_k1000 "$CPU,$((CPU+1))"  --F m3 --d 10 --kD 1000 --kH 2 --kT 3 --parts e1,e2,homog; CPU=$((CPU+2))
launch K2_dn_d10_k300  "$CPU,$((CPU+1))"  --F dense --d 10 --kD 300 --kH 2 --kT 2 --parts e1,e2,homog; CPU=$((CPU+2))
launch K2_m3_d10_dense "$CPU,$((CPU+1))"  --F m3 --d 10 --denseD --kH 1 --kT 1 --parts e1,e2,homog; CPU=$((CPU+2))
launch K2_m3_q2_d8_k100 "$CPU,$((CPU+1))" --F m3 --q 2 --d 8 --kD 100 --kH 1 --parts e1,e2,homog; CPU=$((CPU+2))
