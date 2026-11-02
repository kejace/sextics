#!/bin/sh
# follow-up: degrees 11,12 of the restricted-kernel families (heavy: reducer nullspaces up to ~17k x 26k mod p); checkpoint per (d,par) via the jsonl
cd /tmp/hodge/results/short-power-presentation
FLINT_THREADS=16 python3 w2c_vert.py construct-spp-vert-hi.out.jsonl 11 12 > construct-spp-vert-hi.log 2>&1
