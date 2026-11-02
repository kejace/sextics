#!/bin/bash
# sector relaxations: d=7 (consistency with the full run), d=8, d=9; 4 characters x 2 parities; two primes for d=8
cd /tmp/hodge/results/symmetric-sextic-residual
for d in 7 8; do for chi in triv sgnX sgnY sgnsgn; do for par in 0 1; do
  python3 w2x_relax_job.py --kind sector --chi $chi --d $d --par $par --prime 2147483647 --threads 8 --timeout 7200 --singular
done; done; done
for chi in triv sgnX sgnY sgnsgn; do for par in 0 1; do
  python3 w2x_relax_job.py --kind sector --chi $chi --d 9 --par $par --prime 2147483647 --threads 16 --timeout 14400 --singular
done; done
for chi in triv sgnX sgnY sgnsgn; do for par in 0 1; do
  python3 w2x_relax_job.py --kind sector --chi $chi --d 8 --par $par --prime 2147483629 --threads 8 --timeout 7200
done; done
