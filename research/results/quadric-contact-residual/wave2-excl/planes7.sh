#!/bin/bash
# Job excl-planes7: the 8 "diagonal" affine planes of the parity cube at degree 7 (4 classes each, ~60 unknowns; the 6 coordinate
# half-cubes are in excl-sub7).  msolve -g 1, 16 threads, 2 h cap each.
cd /tmp/hodge/results/quadric-contact-residual/wave2-excl
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-excl/excl-planes7.log.jsonl
J=d7mixed.quads.local.json
run() { MSOLVE_TIMEOUT=7200 python3 sub_closure.py $J "$1" "$2" 16 2>&1 | grep -v "engine v2" | tail -1; }
run d7_pl_even 'b0=0&b2=0&y0|b0=0&b2=1&y1|b0=1&b2=0&y1|b0=1&b2=1&y0'
run d7_pl_odd  'b0=0&b2=0&y1|b0=0&b2=1&y0|b0=1&b2=0&y0|b0=1&b2=1&y1'
run d7_pl_b0b2_0 'b0=0&b2=0&y0|b0=0&b2=0&y1|b0=1&b2=1&y0|b0=1&b2=1&y1'
run d7_pl_b0b2_1 'b0=0&b2=1&y0|b0=0&b2=1&y1|b0=1&b2=0&y0|b0=1&b2=0&y1'
run d7_pl_b0y_0 'b0=0&b2=0&y0|b0=0&b2=1&y0|b0=1&b2=0&y1|b0=1&b2=1&y1'
run d7_pl_b0y_1 'b0=0&b2=0&y1|b0=0&b2=1&y1|b0=1&b2=0&y0|b0=1&b2=1&y0'
run d7_pl_b2y_0 'b0=0&b2=0&y0|b0=1&b2=0&y0|b0=0&b2=1&y1|b0=1&b2=1&y1'
run d7_pl_b2y_1 'b0=0&b2=0&y1|b0=1&b2=0&y1|b0=0&b2=1&y0|b0=1&b2=1&y0'
echo PLANES7 DONE
