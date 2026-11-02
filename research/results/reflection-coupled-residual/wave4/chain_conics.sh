#!/bin/bash
# w4-conics: weighted conics on W through a random point (two seeds, two primes) and inside a random hyperplane section
cd /tmp/hodge/results/reflection-coupled-residual/wave4
T=${HODGE_JOB_CPUS:-16}
echo "=== conics start $(date -u)"
python3 conics_pt.py rpt_s1_p31 2147483647 1 $T pt > rpt_s1_p31.out 2>&1; tail -1 rpt_s1_p31.out
python3 conics_pt.py rpt_s2_p31 2147483647 2 $T pt > rpt_s2_p31.out 2>&1; tail -1 rpt_s2_p31.out
python3 conics_pt.py rhyp_s1_p31 2147483647 1 $T hyp > rhyp_s1_p31.out 2>&1; tail -1 rhyp_s1_p31.out
echo "=== conics done $(date -u)"
