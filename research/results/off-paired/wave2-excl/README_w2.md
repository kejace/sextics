# Wave-2 exclusion/rigor packet for `off-paired` (q = 1, A = X0^3 Y0^3) — rank-sandwich receipts and chart-complete bilinear sweeps

Instance: `data/off-paired/instance.rec`; F = X0^6+X1^6+X2^6+Y0^6+Y1^6+Y2^6 + X0^3X1Y0Y1/2 + X1^3X2Y1Y2/2 + X0^2X1X2Y0Y2/3.
Diagonal group G = Z/6 x Z/6 with generators g = diag(1,z^2,1,z,z^3,z^3), h = diag(1,1,z,1,1,z^5) (z = e^{2 pi i/6});
characters written (g-weight, h-weight) mod 6; chi(A) = (3,0); the residue of A is G-invariant.

(Another wave-2 agent keeps its own independent scripts in this directory: `exlib.py`, `stepA1_*.py`, `cert_*.json`
— rational Farkas certificates.  Everything below is produced by the `w2*.py` scripts only.)

## Files
* `w2lib.py` — library: characters, exact equivariant Koszul Hilbert function of J(F), explicit rational kernel of the fixed-D
  system, rigorous receipts (`e1_receipt`, `fixedD_receipt_certified`), kernel-lift certification, logging (env `W2LOG`).
* `w2_cands.py` — candidate lists: `cands_wave1.jsonl` (621 wave-1 structured denominators), `cands_hessian.jsonl` (213 2x2
  Hessian minors, E1 only), `cands_blocks.jsonl` (680 squarefree semi-invariant block products, degrees 7–14, the dossier's E2 list).
* `w2_batch.py` — multiprocess batch of receipts (two primes). `w2_deltat.py` — degree-18 discriminant full solve (chunked compression).
* `w2_bilinear.py` (+ `w2_fast.py`) — bilinear sector systems: D free in a character sector of the colon space modulo the Euler gauge
  (`--tsub all`) or in the non-Jacobian part (`--tsub nonjac`); `--mode sector` (lossless E1+E2) or `--mode relaxJ` (necessary
  pure-t quadrics modulo J(F)); charts t_b = 0 (b < a), t_a = 1 (`--charts all|rev`) or a random affine normalization (`--charts random`);
  msolve mod 2147483647.  `--fast` uses flint nmod_mpoly products and matrix normal forms (verified identical to the slow path).
* `w2_kernels.py` — colon / non-Jacobian kernel dimensions per sector and degree (`kernel_dims_p*.json`).
* `w2_sync.py` (pull + merge logs), `w2_report.py` (report generator).
* Receipts: `receipts_<job>.jsonl` (one JSON line per (denominator, prime) or per chart). Log: `../wave2-excl*.log.jsonl`.

## Exact statements certified by a receipt
1. `{"exp":"e1", "result":"PROVEN-E1-fail"}`: rank_p of the Jacobian generators of degree d+6 in the sector chi(A)+chi(D) equals the
   exact dimension dim_Q J(F)_{d+6}[chi] (Koszul complex of the smooth sextic: F is smooth by the coefficient-load bound), and the
   augmented rank is one more.  Since rank_p <= rank_Q: A·D is not in J(F) over Q, so no certificate with denominator D exists
   (E1 is necessary for any certificate, after Euler normalization N(F) = A D).
2. `{"exp":"fixedD", "result":"PROVEN-inconsistent", "rigor":"explicit-kernel"}`: in the character-restricted fixed-D system
   M x = b (lossless for an eigenvector D), rank_p([M|b]) = rank_p(M)+1 and rank_p(M) = C − e with e the exact dimension of the
   explicit rational kernel {(L Z, 0)} + {(F W, B_D(W))} (q = 1).  Hence rank_Q(M) <= C − e = rank_p(M) < rank_Q([M|b]): the
   rational system is inconsistent.  By the prescribed-support theorem (q = 1) every denominator with the same reduced support is excluded.
3. `{"rigor":"kernel-certified"}`: as in 2, but rank_Q(M) <= rank_p(M) is certified by a rationally reconstructed basis of the
   full modular kernel, each vector verified EXACTLY to satisfy −D N(F) + F B_D(N) − F^2 Q = 0 (`kernel_verified == kernel_total == nullity`).
4. `{"exp":"bilinear", "result":"no-solution(modp)"}` on chart a: the bilinear system (t, s) with t_b = 0 (b < a), t_a = 1 has no
   solution over the algebraic closure of F_p.  All charts together: no denominator D != 0 of degree d in the sector modulo p.
   This is a SCREEN (two primes together with wave 1 / deanodotdotdot); a rational certificate is needed for a theorem.

## Replay
```
python3 w2_batch.py --cands cands_blocks.jsonl --out receipts_X.jsonl --job X --workers 8
python3 w2_bilinear.py --d 8 --chars 00 --mode sector --charts all --threads 8 --job Y --addJ [--fast]
python3 w2_deltat.py Z
python3 w2_kernels.py 4 12
```
