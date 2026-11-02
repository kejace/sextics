#!/usr/bin/env bash
# runs ON THE REMOTE (after batch5's kron-lib-v3): dense-F d=10 cases, unsliced (psiN) and sliced
CAP=/tmp/hodge/lean/capacity
export PATH="$HOME/.elan/bin:$PATH"
for i in $(seq 1 200); do [ -f /tmp/hodge/jobs/kron-lib-v3.log ] && grep -q '"compile_rc": 0' /tmp/hodge/jobs/kron-lib-v3.log && break; sleep 15; done
CPU=120
launchG() { local gen=$1; local name=$2; local cs=$3; shift 3
  nohup setsid bash -c "cd $CAP/gen && python3 $gen $* --name $name && bash $CAP/run_module.sh Cap.$name $name $cs" > /tmp/hodge/jobs/cap-$name.log 2>&1 < /dev/null &
  echo $! > /tmp/hodge/jobs/cap-$name.pid; echo launched $name pid $!; }
launchG kron_slice_gen.py KS_dn_d10_k20 "$CPU,$((CPU+1))" --F dense --d 10 --kD 20 --kH 1 --kT 1 --parts e1,e2,homog; CPU=$((CPU+2))
launchG kron_gen.py       K5_dn_d10_k20 "$CPU,$((CPU+1))" --F dense --d 10 --kD 20 --kH 1 --kT 1 --parts e1,e2,homog; CPU=$((CPU+2))
launchG kron_slice_gen.py KS_dn_d10_k300 "$CPU,$((CPU+1))" --F dense --d 10 --kD 300 --kH 2 --kT 2 --parts e1,e2,homog; CPU=$((CPU+2))
