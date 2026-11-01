#!/bin/bash
cd /tmp/hodge/results/dense-rational
while pgrep -f "e1kernel.py direct 9" > /dev/null; do sleep 15; done
timeout 3600 python3 deform.py > deform.log 2>&1
timeout 3600 python3 timage.py > timage.log 2>&1
timeout 5400 python3 families2.py deg6b > fam_deg6b.log 2>&1
timeout 7200 python3 families2.py deg7x > fam_deg7x.log 2>&1
timeout 9000 python3 families2.py deg8x > fam_deg8x.log 2>&1
echo done > chain_f.done
