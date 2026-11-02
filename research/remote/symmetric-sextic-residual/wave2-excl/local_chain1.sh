#!/bin/bash
# sequential, single process, niced: small enough for the shared box while the remote is saturated
cd /tmp/hodge/results/symmetric-sextic-residual
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLINT_NUM_THREADS=1
for chi in triv sgnX sgnY sgnsgn; do for par in 0 1; do
  nice -n 15 python3 w2x_relax_job.py --kind sector --chi $chi --d 8 --par $par --prime 2147483647 --threads 2 --timeout 1800 --singular 2>&1 | grep -v "engine v2"
done; done
nice -n 15 python3 w2x_relax_job.py --kind full --d 7 --par 0 --prime 2147483647 --threads 2 --timeout 3600 2>&1 | grep -v "engine v2"
nice -n 15 python3 w2x_relax_job.py --kind full --d 7 --par 1 --prime 2147483647 --threads 2 --timeout 3600 2>&1 | grep -v "engine v2"
nice -n 15 python3 w2x_e1kernel_job.py --prime 2147483629 --degrees 1,2,3,4,5,6,7,8 --tag e1kernel_P2 2>&1 | grep -v "engine v2"
echo CHAIN1 DONE
