#!/usr/bin/env bash
# runs ON THE REMOTE: psiN-based Kronecker modules at d=10 + sliced template at k150 (needs fresh Kron olean)
CAP=/tmp/hodge/lean/capacity
export PATH="$HOME/.elan/bin:$PATH"
bash $CAP/run_module.sh Reflect.Kron kron-lib-v2 100,101 > /tmp/hodge/jobs/kron-lib-v2.log 2>&1
CPU=102
launch() { local name=$1; local cs=$2; shift 2
  nohup setsid bash -c "bash $CAP/kjob.sh $name $cs $*" > /tmp/hodge/jobs/cap-$name.log 2>&1 < /dev/null &
  echo $! > /tmp/hodge/jobs/cap-$name.pid; echo launched $name pid $!; }
launch K4_m3_d10_k300  "$CPU,$((CPU+1))" --F m3 --d 10 --kD 300 --kH 2 --kT 2 --parts e1,e2,homog; CPU=$((CPU+2))
launch K4_m3_d10_k1000 "$CPU,$((CPU+1))" --F m3 --d 10 --kD 1000 --kH 2 --kT 3 --parts e1,e2,homog; CPU=$((CPU+2))
launch K4_dn_d10_k300  "$CPU,$((CPU+1))" --F dense --d 10 --kD 300 --kH 2 --kT 2 --parts e1,e2,homog; CPU=$((CPU+2))
launch K4_m3_q2_d8_k100 "$CPU,$((CPU+1))" --F m3 --q 2 --d 8 --kD 100 --kH 1 --parts e1,e2,homog; CPU=$((CPU+2))
# sliced template at the size where the plain template failed
nohup setsid bash -c "cd $CAP/gen && python3 slice_gen.py --F m3 --d 10 --kD 150 --kH 1 --parts e2 --name SL_e2m3k150 && bash $CAP/run_module.sh Cap.SL_e2m3k150 SL_e2m3k150 $CPU,$((CPU+1))" > /tmp/hodge/jobs/cap-SL_e2m3k150.log 2>&1 < /dev/null &
echo $! > /tmp/hodge/jobs/cap-SL_e2m3k150.pid; echo launched SL_e2m3k150
