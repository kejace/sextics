#!/usr/bin/env bash
CAP=/tmp/hodge/lean/capacity
CPU=80
launch() { local name=$1; local cs=$2; shift 2
  nohup setsid bash -c "bash $CAP/kjob.sh $name $cs $*" > /tmp/hodge/jobs/cap-$name.log 2>&1 < /dev/null &
  echo $! > /tmp/hodge/jobs/cap-$name.pid; echo launched $name pid $!; }
launch K3_dn_d6_k40   "$CPU,$((CPU+1))" --F dense --d 6 --kD 40 --kH 2 --parts e1,e2,homog; CPU=$((CPU+2))
launch K3_dn_d8_k100  "$CPU,$((CPU+1))" --F dense --d 8 --kD 100 --kH 2 --parts e1,e2,homog; CPU=$((CPU+2))
launch K3_dn_d10_k300 "$CPU,$((CPU+1))" --F dense --d 10 --kD 300 --kH 2 --kT 2 --parts e1,e2,homog; CPU=$((CPU+2))
launch K3_dn_d4_k10   "$CPU,$((CPU+1))" --F dense --d 4 --kD 10 --kH 1 --parts e1,e2,homog; CPU=$((CPU+2))
