#!/bin/bash
# Exact mod-p decision (msolve) of the joint certificate problem for every nonempty sector at d=6, then d=7 (small K first).
cd /tmp/hodge/results/four-plus-two-residual
python3 wave2-search/fpt_msolve.py four-plus-two-residual 6 auto 1000003 4
python3 wave2-search/fpt_msolve.py four-plus-two-residual 6 auto 2147483629 4
python3 wave2-search/fpt_msolve.py four-plus-two-residual 7 auto 1000003 4
echo "[msolve-small] done"
