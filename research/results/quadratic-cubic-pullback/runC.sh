#!/bin/bash
ulimit -v 14000000
cd /tmp/hodge/results/quadratic-cubic-pullback
exec python3 run_list2.py ours listC_ours.txt
