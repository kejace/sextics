#!/bin/bash
cd /tmp/hodge/results/dense-rational
while pgrep -f "e1kernel.py direct 9" > /dev/null; do sleep 5; done
python3 deform.py > deform.log 2>&1
