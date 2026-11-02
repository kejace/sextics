#!/bin/bash
# poll until chain D2 has ended (cooperative STOP files), then submit chain E into the freed slot
cd /home/kejace/devel/sextics-work
until timeout 90 python3 -m hodge.rjob status --name search-rank-one-biquartic-chainD2 2>/dev/null | grep -q "^DONE"; do sleep 180; done
echo "chain D2 ended at $(date)"
timeout 300 python3 -m hodge.rjob submit --lane rank-one-biquartic --name search-rank-one-biquartic-chainE --cpus 8 --force -- 'bash wave2/chainE.sh > wave2/out/chainE.log 2>&1'
echo "chain E submitted at $(date)"
