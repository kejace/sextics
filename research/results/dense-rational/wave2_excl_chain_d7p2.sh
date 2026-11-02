#!/bin/bash
# second-prime cross-check of the degree-7 pure-sector certificates (p = 1073741827, the wave-1 prime), sequential
cd /tmp/hodge/results/dense-rational
export PYTHONPATH=/home/kejace/devel/sextics-work
for sec in even odd; do bash wave2_excl_run_quadsys.sh 7 $sec 1073741827 16 excl-dense-rational-d7p2-$sec; done
echo done > wave2-excl/chain_d7p2.done
