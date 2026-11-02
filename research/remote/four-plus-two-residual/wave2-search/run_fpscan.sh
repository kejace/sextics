#!/bin/bash
# Exact finite-field exhaustive scans of small E1-admissible sector families (Koszul-form consistency mod p).
cd /tmp/hodge/results/four-plus-two-residual
# first: polish-and-reconstruct every sub-1e-4 point of the core-only degree-8 family (best 7.5e-5 seen in kstruct)
python3 wave2-search/fpt_search.py --procs 8 --threads 1 --formulation koszul --method lm --maxiter 300 --job core8 --family vars:X0,X1,Y0,Y1 --degrees 8 --restarts 25 --batches 12 --seed0 18000 --cand-rel 1e-4
S="python3 wave2-search/fpt_fpscan.py"
for d in 5 6 7 8; do $S $d auto 1019 3 8; done          # K <= 3: all of P^(K-1)(F_1019)
for d in 7 8; do $S $d auto 61 4 8; $S $d auto 67 4 8; done   # K = 4: two primes
for d in 9 10; do $S $d auto 1019 2 8; done            # pencils at higher degree
echo "[fpscan] done"
# rank / nullity statistics per sector (full and Koszul-form systems)
python3 wave2-search/fpt_ranks.py 5 8
# harder positive controls for the LM pipeline: known witnesses of degree 12 (K=68) and degree 10 (K=235)
python3 wave2-search/fpt_control.py kctl-sylv-lm family-sylvester-norm 800 8 lm 300 0 koszul
python3 wave2-search/fpt_control.py kctl-m3-lm m3-coupled 800 8 lm 300 0 koszul
echo "[fpscan+misc] done"
