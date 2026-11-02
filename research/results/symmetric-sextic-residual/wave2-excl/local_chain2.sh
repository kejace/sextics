#!/bin/bash
cd /tmp/hodge/results/symmetric-sextic-residual
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLINT_NUM_THREADS=1
nice -n 15 python3 w2x_kernel_structure.py 2147483647 8,9 2>&1 | grep -v "engine v2"
nice -n 15 python3 w2x_farkas_job.py D7 D7Y EX EY DX*DY 2>&1 | grep -v "engine v2"
nice -n 15 python3 w2x_sector_screen_job.py --prime 2147483647 --degrees 8,9,10 --npoints 4 --seed 1 2>&1 | grep -v "engine v2"
echo CHAIN2 DONE
