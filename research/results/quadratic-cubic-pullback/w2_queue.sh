#!/bin/bash
# next submissions for lane quadratic-cubic-pullback (run one at a time when a remote slot frees)
cd /home/kejace/devel/sextics-work
case "$1" in
  d89odd)
    python3 -m hodge.rjob submit --force --lane quadratic-cubic-pullback --name search-quadratic-cubic-pullback-jvp-d89 --cpus 8 -- "python3 w2_jvp.py --d 8 --parity 1 --stage 12 --restarts 1600 --restarts2 96 --procs 8 --tag r1 --seed 81 --time-limit 2400 --thresh1 1e-4 > w2_jvp_d8odd_r1.out 2>&1; python3 w2_jvp.py --d 9 --parity 0 --stage 12 --restarts 240 --restarts2 160 --procs 8 --tag r1 --seed 9 --time-limit 5400 --thresh1 1e-6 --maxiter2 200 > w2_jvp_d9even_r1.out 2>&1; python3 w2_jvp.py --d 9 --parity 1 --stage 12 --restarts 240 --restarts2 160 --procs 8 --tag r1 --seed 91 --time-limit 5400 --thresh1 1e-6 --maxiter2 200 > w2_jvp_d9odd_r1.out 2>&1" ;;
  struct1112)
    python3 -m hodge.rjob submit --force --lane quadratic-cubic-pullback --name search-quadratic-cubic-pullback-struct1112 --cpus 1 -- "python3 w2_struct.py --d 10 --mirror /tmp/hodge/results/quadratic-cubic-pullback > w2_struct_d10b.out 2>&1; python3 w2_struct.py --d 12 --mirror /tmp/hodge/results/quadratic-cubic-pullback > w2_struct_d12.out 2>&1; python3 w2_struct.py --d 11 --mirror /tmp/hodge/results/quadratic-cubic-pullback > w2_struct_d11.out 2>&1" ;;
  *) echo "usage: $0 {d89odd|struct1112}" ;;
esac
