#!/bin/bash
# served-space profiling of the 121 structured denominators (degree <= 8), chunked to keep each process short
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
cd /tmp/hodge/results/reflection-coupled-residual
split -l 12 -d cands_served.txt cands_served_chunk_
for f in cands_served_chunk_*; do timeout 600 python3 job_served.py $f 8 2>&1 | grep -v "engine v2" | grep "served\]\|error" | cut -c1-330; done
echo LOCAL_SERVED_DONE
