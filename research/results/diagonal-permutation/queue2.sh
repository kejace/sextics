#!/bin/sh
cd /tmp/hodge/results/diagonal-permutation
python3 /tmp/hodge/results/diagonal-permutation/batch3.py > batch3.out 2>&1
python3 /tmp/hodge/results/diagonal-permutation/batch5.py > batch5.out 2>&1
python3 /tmp/hodge/results/diagonal-permutation/batch2b.py > batch2b.out 2>&1
python3 /tmp/hodge/results/diagonal-permutation/batch4.py > batch4.out 2>&1
