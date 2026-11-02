#!/bin/bash
cd /tmp/hodge/results/cycle-q2
while kill -0 3931661 2>/dev/null; do sleep 10; done
python3 exp_h_served12.py a > exp_h_served12.out 2>&1
python3 exp_e_families.py subsets > exp_e_subsets.out 2>&1
