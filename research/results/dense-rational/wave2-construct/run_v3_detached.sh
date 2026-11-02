#!/bin/bash
# detached duplicate of the V3 (reflection locus) certificate run; survives the agent session
cd /tmp/hodge/results/dense-rational/wave2-construct
f=V3_reflection_p1073741789
/usr/bin/time -f "%e s %M KB" nice -n 19 msolve -g 1 -t 1 -f $f.ms -o ${f}_dup.g1.out 2> ${f}_dup.time.log
python3 ../wave2_construct_ltcheck.py $f.ms ${f}_dup.g1.out construct-V3-detached > ${f}_dup.ltcheck.txt 2>&1
touch ${f}_dup.done
