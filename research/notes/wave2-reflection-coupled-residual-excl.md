# Wave 2 — `reflection-coupled-residual` — lens EXCLUSION / RIGOR

STATUS: in progress — no witness; RIGOROUS: deg D ≥ 7 for every denominator (d ≤ 4: E1; d = 5, 6: E1 + Hessian relaxation; E1 re-verified at a third prime in all four characters); degree 7 excluded for three of the four pure Klein characters so far; K8 (the unique non-Jacobian degree-8 E1 class) certified not a denominator

Started 2026-11-01 ~19:15 local. Work dir `/tmp/hodge/results/reflection-coupled-residual/wave2-excl/`,
log `/tmp/hodge/results/reflection-coupled-residual/wave2-excl.log.jsonl` (one JSON line per experiment).
Remote jobs named `excl-reflection-coupled-residual-*`.

## Operational note
Remote receipts are merged with `wave2-excl/merge_remote.py` after every `rjob pull`; `rjob submit` rsyncs the local lane dir
over the remote, which clobbered the remote `e1char.jsonl` once (degree-9 receipts for three characters lost; rerun queued with
W2LIB_SUFFIX so remote and local receipt files never share a name again).

## Plan
(a) re-verify wave-1/dossier rigorous claims with a different prime (p = 2147483587 / 2147483579; wave 1 used 2147483647, 2147483629; dossier used 32003), with explicit good-reduction (Hilbert function) checks;
(b) rigorous lower bound on deg D via E1 by Klein character + Hessian relaxation of E2 on the Jacobian family;
(c) non-Jacobian E1-kernel classes and structured families, degrees 7-12, character sectors, full E1+E2 (symmetric solver);
(d) packet; (e) draft post.

## Results so far

### R1. Degree 5 (Hessian relaxation, local, p = 2147483587) — DONE
For D = Σ c_i F_i (c_i constants; this is the whole E1-admissible space at d = 5 once job e1char confirms
"E1 ⇒ D ∈ J_5" in all four characters), the E2-mod-J obstruction is ob(c) = −[A Σ c_i c_j F_ij] ∈ (S/J)_10.
The Hessian has exactly one vanishing off-diagonal entry, F_{X2 Y2} = 0; the other 20 classes [A F_ij] are
linearly independent in (S/J)_10 (rank 20 in each of the four target-character blocks: sizes 6/10/2/2;
good reduction verified: rank J_10^ψ = dim S_10^ψ − Fermat count for all ψ). Hence c_i c_j = 0 for all 20 pairs,
in particular c_i² = 0, so c = 0. **No denominator of degree 5 (any character, mixed parity included).**
Receipts: `wave2-excl/hess.jsonl` (exp phi-kernel/phi-verdict, d=5).

### R2. Degree 6 (Hessian relaxation, local, p = 2147483587) — DONE: deg-6 excluded
The linear criterion fails: ker Φ_6 has dimension 75 (per target character 17/26/13/19; trivial part only 1): the
Euler relations Σ_j z_j F_ij = 5F_i put every symmetrised tuple z⊗a + a⊗z (a ∈ S_1^6, 36 dims) into ker Φ. So the
genuine quadratic system is needed: on W = V/T (V = S_1^6 = 36, T = C·z, 35 unknowns, all four characters jointly,
so mixed-parity D included) the obstruction gives 963 nonzero quadrics of which 345 are independent (< 630 = dim S_2(W)),
so the degree-2 part alone does not decide. msolve (grevlex GB mod p, 709 elements, 1.8 s) + Singular: **dim 0, mult 1060,
minimal associated primes = {(w0,…,w34)}**, Hilbert function HF(0..7) = 1, 35, 285, 477, 225, 36, 1, **0**.
Hence I_p ⊇ S_7(W); since HF_Q(k) ≤ HF_p(k) (rank argument, good reduction in degree 12 verified), I_Q ⊇ S_7(W),
so V(ob) = T over ℂ: every g with ob(g) = 0 has g(F) ∈ (F). **No denominator of degree 6** (any character, mixed
parity included), given E1 ⇒ D ∈ J_6 in all characters (dossier mod 32003; re-checked by job e1char).
Files `wave2-excl/hess_d6_W00011011_p2147483587.{ms,json,lt,gb}`, `hess_d6_singular_minass.txt`.
Proof that B(T,·)=0 exactly: B(zG, g) = [A G Σ_j g_j Σ_i z_i F_ij] = 5[A G g(F)] = 0 since g(F) ∈ J.


### R3. Degree 7, pure Klein characters (Hessian relaxation + msolve/Singular, local, p = 2147483587) — 3 of 4 DONE
D ∈ J_7 (E1) of a single character χ ⇒ g ∈ V_χ, W_χ = V_χ/T_χ; the quadrics live in (S/J)_14^{(0,1)} (dim 426).
| χ | dim V | dim T | dim W | indep. quadrics | GB size | dim | mult | HF | verdict |
|---|---|---|---|---|---|---|---|---|---|
| (1,1) | 18 | 0 | 18 | 117 | 135 | 0 | 85 | 1,18,54,12,0 | V = {0} |
| (1,0) | 23 | 1 | 22 | 132 | 301 | 0 | 239 | 1,22,121,83,12,0 | V = {0} |
| (0,1) | 40 | 2 | 38 | 281 | 2987 | 0 | 1109 | 1,38,460,556,54,0 | V = {0} |
| (0,0) | 45 | 3 | 42 | 274 | (running) | | | | |
Hence (modulo the (0,0) block) **no degree-7 denominator with a single Klein character**. Mixed-character degree-7
denominators: the joint systems (pairs 40–80 unknowns, triples 78–102, full 120) are being exported and run
(chain `wave2-excl/chain_d7mixed.sh`; receipts `gbcheck.jsonl`).

### R4. The degree-8 E1-kernel class K8 (remote job excl-rcr-k8job) — DONE: K8 rigorously excluded
* K8 ∉ J(F): verified mod 2147483587 and 2147483579 (rank 168 vs 169); a rational polynomial outside J mod p is outside J over ℚ.
* A·K8 ∈ J(F) **exactly over ℚ**: the character-restricted E1 system (3276 × 3388, rank 2850) was solved mod two primes,
  CRT-lifted (62 bits, max denominator 18 bits) and the identity Σ N_i ∂_iF = A·K8 verified exactly (`wave2-excl/K8_E1_exact.json`,
  N-term counts 28/65/42/56/23/0). So the wave-1 class is a genuine rational E1-admissible non-Jacobian denominator.
* Fixed-D certificate test (rigorfix protocol): Klein-invariant-sector system 22016 × 4256 (of 7497 unrestricted columns,
  stabiliser order 4): **inconsistent mod p**, rank_p M = 3780, nullity 476; the 482 a-priori rational kernel vectors
  (zG, 0) and (F m e_i, D ∂_i m − m ∂_i D) have rank 476 = nullity, hence rank_ℚ M = rank_p M and the inconsistency holds
  over ℚ (and ℂ). Second prime agrees; the full unrestricted system (40301 × 7497) is also inconsistent mod p (rank 6703).
  **K8 itself is not a certificate denominator.** (427 s.) The affine family K8 + J_8 is treated in R6.

### R5. E1 by Klein character, independent re-verification (remote job excl-rcr-e1char, p = 2147483587; reps lifted with 2147483579)
Kernel of [A]: (S/J)_d^χ → (S/J)_{d+6}^{χ+(0,1)}; "ker" = dim{D ∈ S_d^χ : AD ∈ J}, "mod J" = its dimension modulo J_d^χ.
Good reduction (rank J_n^χ = dim S_n^χ − Fermat count) holds in every row, so each modular kernel dimension bounds the
rational one from above and the zeros are theorems.
| d | χ=(0,0): dim S, (S/J), ker, mod J | (0,1) | (1,0) | (1,1) |
|---|---|---|---|---|
| 1–4 | ker = 0 in all characters (dims 3/2/1/0, 10/6/3/2, 22/18/10/6, 48/38/22/18) | | | |
| 5 | 88, 85, 3, **0** | 78, 76, 2, **0** | 48, 47, 1, **0** | 38, 38, 0, **0** |
| 6 | 158, 144, 14, **0** | 138, 126, 12, **0** | 88, 82, 6, **0** | 78, 74, 4, **0** |
| 7 | 258, 213, 45, **0** | 238, 198, 40, **0** | 158, 135, 23, **0** | 138, 120, 18, **0** |
| 8 | 413, 301, 113, **1** (= K8, two primes, rational lift agrees with wave 1) | 378, 274, 104, **0** | 258, 194, 64, **0** | 238, 182, 56, **0** |
| 9 | 623, 381, 263, **21** | (receipts lost to an rsync clobber; rerun queued; reps files show kernel 0 for (0,1),(1,1) and > 0 for (1,0)) | | |
| 10 | —, 460, 650, **186** (dimension-forced: 460 → 274) | (running) | | |
Agreement with wave 1 (Y-even sector, primes 2^31−1, 2^31−19) and the dossier (all characters mod 32003, d ≤ 6): complete.
New: the Y-odd characters at d = 7, 8 have no non-Jacobian E1-admissible class either, so at degree 8 the whole
E1-admissible space is J_8 ⊕ ℂ·K8 (mixed characters allowed).

### R7. Structured-family E1 screen (remote job excl-rcr-struct, p = 2147483587, d = 7..12) — d = 7 done
48 monomial-defined shapes (cones with a line or plane vertex, multiples of coordinates/squares, polynomials in squares,
X-only/Y-only, bounded X/Y-degree, the Klein-invariant subring, (X0,Y0)-ideal powers, …), each intersected with the
E1-admissible space by character. Degree 7: only `multiple-of-X0` (36), `multiple-of-X0^2` (6), `multiple-of-X1` (36),
`multiple-of-X1^2` (6) contain E1-admissible elements, all inside J_7 (consistent with R5); no cone, binary/ternary form,
polynomial in squares, u = Y2²-polynomial, … of degree 7 passes E1. Degree 8: only `cone-line: independent of X0,Y2` (1 element) and `multiple-of-X0` (126) meet the E1-admissible space, both
inside J_8: **no structured shape contains a member of the K8 family** (extends wave-1's 55-shape test). Degrees 9–12: running
(survivors outside J are saved as `struct_d<d>_<shape>_chi<ab>.json` for fixed-D tests).

### Jobs submitted (remote)
* `excl-rcr-e1char`: E1 kernels by Klein character, d = 1..10, p = 2147483587 (+ 2147483579 for reconstruction).
* `excl-rcr-hess`: Φ-kernel at d = 6 (two primes); Macaulay-3 rank of the quadratic system per character at d = 7;
  msolve exports of the pair/full systems at d = 7.
* `excl-rcr-k8job`: rigorous fixed-D test of K8 (restricted + full), K8 ∉ J, exact E1 lift, e2lin.

