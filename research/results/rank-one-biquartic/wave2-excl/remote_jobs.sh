#!/bin/bash
# Remote job commands (submit when the 128-core box has capacity; each is independent).  Lane: rank-one-biquartic.
cd /home/kejace/devel/sextics-work
L=rank-one-biquartic
# degree-8 non-Jacobian generators (exact) then the degree-8 pipeline (both need nonj_d7.json from k7exact.py 7)
python3 -m hodge.rjob submit --wait 7200 --lane $L --name excl-$L-k8exact --cpus 4 -- 'python3 wave2-excl/k7exact.py 8 --known wave2-excl/nonj_d7.json > wave2-excl/k7exact_d8.out 2>&1'
python3 -m hodge.rjob submit --wait 7200 --lane $L --name excl-$L-hpipe8 --cpus 8 -- 'python3 wave2-excl/hess_pipeline.py 8 2 --sparse --nonj wave2-excl/nonj_d8_all.json --tag hpipe_d8 > wave2-excl/hpipe_d8.out 2>&1'
# modular screening of degree 9 (no exact step)
python3 -m hodge.rjob submit --wait 7200 --lane $L --name excl-$L-hmodp9 --cpus 8 -- 'python3 wave2-excl/hess_modp.py 9 1000000007 --nonj wave2-excl/nonj_d9_all.json --tag hmodp_d9 > wave2-excl/hmodp_d9.out 2>&1'
# re-verification: dense fixed-D at a new prime; Singular quotient at 32003
python3 -m hodge.rjob submit --wait 7200 --lane $L --name excl-$L-V3fixedD --cpus 12 -- 'python3 wave2-excl/verify_fixedD.py 2147483423 12 > wave2-excl/verify_fixedD.out 2>&1'
python3 -m hodge.rjob submit --wait 7200 --lane $L --name excl-$L-V2quot32003 --cpus 1 -- 'cd wave2-excl && Singular -q quotJA_p32003.sing > quotJA_p32003.out 2>&1'
