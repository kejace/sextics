#!/bin/bash
cd /tmp/hodge/results/quadric-contact-residual/wave2-excl
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-excl/excl-triples7.log.jsonl
python3 triples7.py 8 3600 2>&1 | grep -v "engine v2" | tail -60
