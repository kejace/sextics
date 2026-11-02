#!/bin/bash
# Pipeline E: EXACT (mod p) Hessian-congruence closures for character-mixed pure-parity families (Singular std).
cd /tmp/hodge/results/quadric-contact-residual
export PYTHONPATH=/home/kejace/devel/sextics-work
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-search.pipeE.log.jsonl
run() { echo "=== $(date '+%F %T') $*"; python3 "$@" 2>&1 | grep -v 'hodge engine v2'; echo "=== exit $? $(date '+%F %T')"; }
run w2_hessgb.py --d 7 --parity 0 --timeout 5400
run w2_hessgb.py --d 7 --parity 1 --timeout 5400
run w2_hessgb.py --d 8 --parity 0 --timeout 7200
run w2_hessgb.py --d 8 --parity 1 --timeout 7200
run w2_hessgb.py --d 9 --parity 0 --timeout 7200
echo "PIPELINE E DONE $(date)"
