#!/usr/bin/env bash
# runs ON THE REMOTE: recompile Kron (slicing + Option degrees), then sliced Kronecker jobs at d=10
CAP=/tmp/hodge/lean/capacity
export PATH="$HOME/.elan/bin:$PATH"
bash $CAP/run_module.sh Reflect.Kron kron-lib-v3 100,101 > /tmp/hodge/jobs/kron-lib-v3.log 2>&1
CPU=112
launchS() { local name=$1; local cs=$2; shift 2
  nohup setsid bash -c "cd $CAP/gen && python3 kron_slice_gen.py $* --name $name && bash $CAP/run_module.sh Cap.$name $name $cs" > /tmp/hodge/jobs/cap-$name.log 2>&1 < /dev/null &
  echo $! > /tmp/hodge/jobs/cap-$name.pid; echo launched $name pid $!; }
launchS KS_m3_d10_k300  "$CPU,$((CPU+1))" --F m3 --d 10 --kD 300 --kH 2 --kT 2 --parts e1,e2,homog; CPU=$((CPU+2))
launchS KS_m3_d10_k1000 "$CPU,$((CPU+1))" --F m3 --d 10 --kD 1000 --kH 2 --kT 3 --parts e1,e2,homog; CPU=$((CPU+2))
launchS KS_dn_d10_k300  "$CPU,$((CPU+1))" --F dense --d 10 --kD 300 --kH 2 --kT 2 --parts e1,e2,homog; CPU=$((CPU+2))
launchS KS_m3_d8_k150   "$CPU,$((CPU+1))" --F m3 --d 8 --kD 150 --kH 2 --parts e1,e2,homog; CPU=$((CPU+2))
