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

## Exact (theorem-level) sector exclusions: `w2_elim2.py` (+ `w2_local.py`, `w2_elim.py`)
Statement proved for a sector (d, c) when the receipt says `PROVEN-no-denominator-in-sector(exact: ...)`:
*no legal denominator of degree d that is a (g,h)-eigenvector of character c exists.*  Proof structure (all over Q):
1. D_a = rational basis of (J:A)_{d,c} modulo the Euler gauge F·S_{d-6} (reconstructed from a modular basis and VERIFIED: N_a(F) = A·D_a
   exactly; independence modulo p implies independence over Q; the modular colon dimension bounds the rational one, so the D_a span).
   S_k = rational Koszul kernel basis (verified S_k(F) = 0).  B-products and normal forms modulo F·S_{2d-6}[sector] computed exactly
   (fmpq_mpoly, rational rref), giving the exact bilinear equations eq_r(t, s) of the Euler-normalized certificate problem.
2. Degree-2 relaxation: every lambda in the left kernel of the s-coefficient blocks (exact integer kernel, LLL-reduced) gives a pure-t
   quadric lambda^T q(t) satisfied by every certificate.  Singular over Q: if the ideal J2 of these quadrics has `dim(std(J2)) = 0`
   (affine cone = {0}) the sector is excluded.
3. Otherwise `minAssGTZ(J2)` over Q gives the components; each component lies in its linear span L (basis from the linear generators of
   its prime).  On L, exact relaxations of degree 2, 3, 4 (`w2_local.py`: cofactors mu_r(x) of degree k-2 whose combination has no
   s-terms) give pure-x consequences; Singular over Q checks that they generate an m-primary ideal in Q[x].  A component of span
   dimension 1 is a single rational D and is excluded by the fixed-D rank-sandwich receipt instead.
4. Empty colon sectors (`PROVEN-empty-sector`) are excluded because the modular colon dimension (0) bounds the rational one.
No modular computation enters the conclusion; modular arithmetic is only used to FIND the rational bases and the candidate structure.
Receipts: `receipts_excl-off-paired-elim2Q-d7.jsonl`, `receipts_excl-off-paired-chain-rigor-d{7,8}.jsonl`; Singular scripts in `elim2_*/`.
