#!/bin/bash
cd /tmp/hodge/results/reflection-coupled-residual/wave2-excl
bash job_companion.sh
export W2LIB_SUFFIX=.r4
python3 e1recon.py 9 00 8 > e1recon_d9_00.out 2>&1
python3 e1recon.py 9 10 8 > e1recon_d9_10.out 2>&1
echo MISC DONE >> job_misc.log
