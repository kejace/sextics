#!/bin/bash
# Job excl-lattice7: (1) bilinear zero-pair tests for all 28 parity-class pairs at degree 7 (resumable; skips pairs already
# logged in the master log copy), (2) the 8 diagonal affine planes (60 unknowns) by msolve -g 1 (16 threads, 2 h cap each).
cd /tmp/hodge/results/quadric-contact-residual/wave2-excl
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-excl/excl-lattice7.log.jsonl
# seed the job log with pairs already certified locally (so they are skipped)
grep '"bilinear-pair"' /tmp/hodge/results/quadric-contact-residual/wave2-excl.log.jsonl 2>/dev/null | grep "NO ZERO PAIR" >> $HODGE_W2LOG
python3 bilinear_pairs2.py d7mixed.quads.local.json 2 2 all 2>&1 | grep -v "engine v2" | tail -30
bash planes7.sh
echo LATTICE7 DONE
