#!/bin/sh
# Construction-lens remote jobs for short-power-presentation (queued but not started on 2026-11-01 because the remote stayed at load ~1000).
# Submit from brunos (each job: --cpus 8, ~1-2 h):
cd /home/kejace/devel/sextics-work
python3 -m hodge.rjob submit --wait 7200 --lane short-power-presentation --name construct-spp-vert  --cpus 8 -- "FLINT_THREADS=8 python3 w2c_vert.py construct-spp-vert.out.jsonl 7 8 9 10 > construct-spp-vert.log 2>&1"
python3 -m hodge.rjob submit --wait 7200 --lane short-power-presentation --name construct-spp-k7var --cpus 8 -- "FLINT_THREADS=8 python3 w2c_k7var.py construct-spp-k7var.out.jsonl > construct-spp-k7var.log 2>&1"
# follow-ups once the above are in: degrees 11-12 (w2c_vert_hi.sh, --cpus 16) and, for any E1-passing vertical candidate, the two-prime reconstruction + fixed-D solve:
#   python3 w2c_vert2.py construct-spp-vert2.out.jsonl vert 10 0   (family d parity)
