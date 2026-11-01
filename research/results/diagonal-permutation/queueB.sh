#!/bin/sh
cd /tmp/hodge/results/diagonal-permutation
while pgrep -f "batch2b.py" >/dev/null; do sleep 10; done
python3 batch3.py > batch3.out 2>&1
python3 batch5.py > batch5.out 2>&1
python3 batch4.py > batch4.out 2>&1
