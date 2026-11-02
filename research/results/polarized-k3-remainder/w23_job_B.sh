#!/bin/bash
# wave-2 construction lens (repurposed): extra F_5 seeds (other kx) + random F_11 / F_13 seed search on W_{23,f}
cd /tmp/hodge/results/polarized-k3-remainder
for kx in 1 2 3; do python3 w23_seed_fastA.py --p 5 --kx $kx --out construct-w23seedA-p5-kx$kx.jsonl > construct-w23seedA-p5-kx$kx.out 2>&1; done
python3 w23_seed_randomA.py --p 11 --kx 1 --out construct-w23seedA-p11-kx1.jsonl --trials 300000000 --seed 11 > construct-w23seedA-p11-kx1.out 2>&1 &
python3 w23_seed_randomA.py --p 13 --kx 1 --out construct-w23seedA-p13-kx1.jsonl --trials 600000000 --seed 13 > construct-w23seedA-p13-kx1.out 2>&1 &
wait
echo ALLDONE > construct-w23seed.done
