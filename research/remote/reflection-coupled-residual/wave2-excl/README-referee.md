# Referee replay guide — `reflection-coupled-residual`, wave-2 exclusion packet

Written by a referee pass (Claude, Anthropic, for `kejace`; no human review). Everything below was re-derived with code that
does **not** import the lane engine (`referee/reflib.py`: dict polynomials over `fractions.Fraction`, python-flint `nmod_mat` /
`fmpq_mat` for linear algebra), at primes nobody else used: `p = 1000000021` (python-flint), `p = 1000003` (Singular).
Pinned engines: python-flint 0.9.0, msolve 0.6.5, Singular 4.3.2. Receipts: `referee/ref_*.jsonl` (one JSON line each).

Instance: `data/reflection-coupled-residual/instance.rec` (q = 1, A = X0^3 Y0^3, F the coupled sextic with coupling 5/7),
parsed directly from the record by `reflib.load_instance()`. Variables X0,X1,X2,Y0,Y1,Y2; Klein character of a monomial
chi(m) = (deg_Y2 m mod 2, (deg_Y0 m + deg_Y1 m) mod 2); chi(F) = (0,0), chi(A) = (0,1).

## What the statements mean (contract, `formal/README.md`)
A certificate is (D, N, Q) with E := A D^2 - D sum_i N_i F_i + F (D sum_i d_i N_i - sum_i N_i d_i D) = F^2 Q (q = 1),
D coprime to F, deg N_i = deg D + 1. Since F is irreducible and F does not divide D, reducing E mod F gives
**E1: A D in J(F)** (J = (F_0,...,F_5)); J is Klein-graded, so E1 holds character by character.
For D = g(F) = sum g_i F_i the natural E1 solution is N = A g, any other E1 solution differs by Koszul terms and multiples of z,
and E = F^2 Q forces **ob(g) := -[A g^T Hess(F) g] = 0 in (S/J)_{2d}** ("Hessian relaxation of E2"); ob(g + zG) = ob(g).

## Replay commands (run from this directory; each one appends a receipt)
```sh
cd referee
python3 ref_e1.py 1 10 1000000021          # E1 kernels by character; ~15 min for d = 9, 10 (a few seconds for d <= 8)
python3 ref_k8.py 1000000021                # exact A*K8 = sum N_i F_i over Q (fractions); K8 not in J mod p; A*K8 in J mod p
python3 ref_hess.py rank20 1000000021       # degree 5: the 20 classes [A F_ij] (F_{X2 Y2} = 0) independent in (S/J)_10
python3 ref_hess.py export 6 1000000021 00 01 10 11      # degree 6, all characters jointly (35 unknowns)
python3 ref_hess.py export 7 1000000021 11               # degree 7 single characters: 11, 10, 01, 00; pair: "10 11"
python3 ref_gb.py ref_hess_d6_W00011011_p1000000021.ms 2 600 8     # msolve GB + Singular dim/vdim/Hilbert function
python3 ref_gb.py ref_hess_d7_W11_p1000000021.ms 2 600 8           # likewise W10 (1 s), W01 (~150 s), W00 (~1-2 h), W1011 (~3 min)
python3 ref_fixedD.py K8 ../K8_E1_exact.json:D 1000000021 5200     # fixed-D test of K8, certified over Q (see below)
python3 ref_fixedD.py d9b0 ../rigorfix.r2.jsonl:0/D 1000000021 8400 # one degree-9 kernel representative
python3 ref_k8shapes.py                                             # exact (fmpq_mat): which shapes meet K8 + J_8
Singular -q ref_smooth.sing                                         # smoothness of X and of X cap {z_i = 0}, p = 1000003
```

## Why a modular result is a proof here (transfer arguments, checked)
1. **Ranks.** For an integer matrix M, rank_p M <= rank_Q M. "Good reduction" in degree n and character chi means
   rank_p J_n^chi = dim S_n^chi - #{monomials of degree n, all exponents <= 4} (the Fermat count = dim (S/J)_n for any smooth
   sextic). It is recorded in every receipt (`good_d`, `good_target`, `good_reduction`). With it, the saturated lattice
   J_Q cap Z_(p)[z] reduces onto J_p in that degree, so the modular E1 kernel dimension bounds the rational one from above:
   `kernel_mod_J = 0` mod p is a theorem over Q and C; `kernel_mod_J = 1` is an upper bound (K8's membership is then exact, `ref_k8.py`).
2. **Quadratic systems.** With good reduction in degree 2d, the mod-p quadrics are reductions of a Z_(p)-basis of the rational
   coordinate quadrics of ob on W = V/T, so HF_Q(k) <= HF_p(k) (Macaulay-matrix ranks). A finite quotient dimension mod p
   (msolve "Dimension of quotient", Singular `dim = 0`) gives HF_p(k0) = 0 for k0 = 1 + max standard degree, hence
   HF_Q(k0) = 0, hence V(ob) = T over C: every D = g(F) with ob(g) = 0 lies in (F) and is not coprime to F.
3. **Fixed-D inconsistency.** (zG, 0) and (F m e_i, D d_i m - m d_i D) are exact rational kernel vectors of the homogeneous
   fixed-D system (checked symbolically and mod p). If their rank r_ap satisfies rank_p(M_sub) = ncols - r_ap for some row subset,
   then rank_Q M = ncols - r_ap, and rank_p([M|b]_sub) = rank_Q M + 1 proves inconsistency over Q and C.
4. **Smoothness.** `dim = 0` of the Jacobian ideal of X cap {z_i = 0} mod p means the ideal is m-primary mod p; by the
   Hilbert-function inequality of item 1 it is m-primary over Q, so the section is smooth over C.
5. **Reduction lemma (served(z_i G) = served(G)).** Not a computation: it needs the prescribed-support theorem of the book,
   Chapter 1 (`docs/01-certificates.md`, statements 1 <=> 3 with denominator D_red^q, q = 1) and the Gysin description of
   N^1 H^4 in Chapter 3 §3.1 (a smooth hyperplane section Z_i has H^{2,0} = 0 and H^2(Z_i) = Q h, and h^2 vanishes off any
   hypersurface section). Those book theorems rest on Hodge-module strictness (Mustata-Popa) and are not Lean-verified.

## Referee results (all agree with the lane's receipts unless stated)
* E1 kernels mod J (d; chi 00/01/10/11): d <= 7: 0/0/0/0; d = 8: 1/0/0/0; d = 9: 21/0/21/0; d = 10: 186/126/136/108 (upper bounds); good reduction in every row.
* A*K8 = sum N_i F_i holds exactly over Q (37-term K8, N_i with 28/65/42/56/23/0 terms); K8 not in J_8 (mod p, good reduction).
* Degree 5: 20 nonzero Hessian pairs, rank 20 mod J_10.  Degree 6: W = 35, 345 independent quadrics, GB 709, quotient 1060,
  HF 1,35,285,477,225,36,1,0.  Degree 7: (1,1) 18/117 -> quotient 85, HF 1,18,54,12,0; (1,0) 22/132 -> 239, HF 1,22,121,83,12,0;
  (0,1) 38/281 -> 1109, HF 1,38,460,556,54,0; pair (1,0)+(1,1) 40/300 -> quotient 4399, HF 1,40,520,1543,1695,567,33,0 (64 s);
  (0,0) 42/274: the lane's run gives quotient 3028, HF 1,42,629,1986,359,11,0 (4224 s at 2147483587); the referee run at 10^9+21 is
  recorded in `ref_gb.jsonl` when finished (it was still running when the post was written).
* K8 fixed-D: restricted system 22016 x 4256, 482 a-priori kernel vectors of rank 476, rank_p(M_sub) = 3780 = 4256 - 476,
  rank_p([M|b]_sub) = 3781: inconsistent over Q and C.
* Smoothness (Singular, p = 1000003): X: dim 0, vdim 15625 = 5^6; all six sections: dim 0, vdim 3125 = 5^5.
* Shapes meeting the affine family K8 + J_8, exact over Q: **the family contains a member with all Y-exponents even** (the lane's
  prose "no structured shape contains a member of the K8 family" is wrong for that shape and for the Klein-invariant subring,
  which contains K8 itself; the lane's own `structscreen.jsonl` rows for d = 8 say the same). No member is a polynomial in squares,
  a cone independent of (X0, Y2), a multiple of X0, or of X- or Y-degree <= 4.
