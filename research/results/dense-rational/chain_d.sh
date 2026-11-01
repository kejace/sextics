#!/bin/bash
cd /tmp/hodge/results/dense-rational
while pgrep -f "chain_b.sh" > /dev/null; do sleep 10; done
python3 timage.py > timage.log 2>&1
