#!/bin/bash
# Job excl-misc: wave-1 replays at P2 (sector solves of D8*, Z_i D8*; engine cross-check), served spaces of the degree-6/8 blocks,
# then the degree-10 extras screen (all sectors, 8 workers).
cd /tmp/hodge/results/quadric-contact-residual/wave2-excl
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-excl/excl-misc.log.jsonl
python3 fixedD_batch.py 1 2 2>&1 | grep -v "engine v2" | tail -3
ONLY=rad0,rad1,rad2,X1^6+Y1^6,X2^6+Y2^6,D8* python3 served_local.py 8 2>&1 | grep -v "engine v2" | tail -3
python3 extras_screen.py 10 8 2>&1 | grep -v "engine v2" | tail -3
echo MISC DONE
