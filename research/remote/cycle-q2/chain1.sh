#!/bin/bash
cd /tmp/hodge/results/cycle-q2
while kill -0 4027624 2>/dev/null; do sleep 10; done
python3 exp_g_deg5.py > exp_g_deg5.out 2>&1
