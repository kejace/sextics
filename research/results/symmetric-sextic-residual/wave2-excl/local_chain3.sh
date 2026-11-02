#!/bin/bash
cd /tmp/hodge/results/symmetric-sextic-residual
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLINT_NUM_THREADS=1
for chi in sgnsgn sgnX sgnY triv; do for par in 0 1; do
  nice -n 15 python3 w2x_relax_job.py --kind sector --chi $chi --d 9 --par $par --prime 2147483647 --threads 2 --timeout 1800 2>&1 | grep -v "engine v2"
done; done

nice -n 15 python3 w2x_exactkernel_small.py 1,2,3 2>&1 | grep -v "engine v2"
nice -n 15 python3 w2x_relax_job.py --kind full --d 6 --par mixed --prime 2147483629 --threads 2 --timeout 1800 2>&1 | grep -v "engine v2"
echo CHAIN3b DONE
