# Wave 2 — `short-power-presentation` — lens EXCLUSION / RIGOR

STATUS: in progress — no witness. Rigorous so far (this session, independent prime 1073741789, two Gröbner engines):
**no denominator of degree ≤ 5 of any parity; no Y-pure denominator of degree 6 (other than multiples of F).**
Mixed parity d = 6: excluded mod p by the annihilator-functional argument (§0 R4); exact certificate queued (remote).

Lane: `short-power-presentation` (q=1, A = X0^3 Y0^3, F = eX + eY + c((L+M)^6+(L-M)^6), c = 5/7^7,
L = X0+2X1+3X2, M = Y0+3Y1+5Y2). Log: `/tmp/hodge/results/short-power-presentation/wave2-excl.log.jsonl`.
Scripts (prefix `excl_`, lane dir is shared with the search lens): `excl_lib.py`, `excl_e1kernel.py`, `excl_relax.py`, `excl_slice.py`.
Packet: `/tmp/hodge/results/short-power-presentation/wave2-excl/` (assembled at the end).

## 0. Rigor framework used (all statements below are over Q-bar unless marked "mod p")

* (R1) Hilbert-series check: F is smooth, so dim_Q J(F)_e^{±} is known exactly
  (H_R(t) = (1+t+…+t^4)^6, trace of τ: (1+t^2+…+t^8)^3; dim R_e^± = (H_e ± Θ_e)/2). Every modular rank of a
  Jacobian generator matrix is compared with this value; equality means the prime is "good" in that degree and the
  mod-p quotient S_e/J_e is the reduction of the Q-quotient.
* (R2) rank_p ≤ rank_Q. E1 failure (A·D ∉ J) mod a good prime proves it over Q; dim_Q(K_d/J_d) ≤ dim_p(K_d/J_d) when
  the prime is good in degrees d and d+6.
* (R3) Relaxation of E2 (derivation): for a certificate with D = Σ λ_a D_a + w(F) (D_a kernel representatives with
  particular N_a, N_a(F) = A D_a; w a vector field of degree d-5), N can be normalised to N = Σ λ_a N_a + A w + (G/6)E + K
  (E Euler, K Koszul), and the identity forces
  Ξ(λ,w) := Σ λ_aλ_b B_{D_a}(N_b) + Σ_a λ_a (D_a w(A) − A w(D_a) − N_a^T Hess(F) w) − A·w^T Hess(F) w ∈ J(F)_{2d}
  (B_D(N) = D div N − N(D)). Ξ is a homogeneous quadratic map; the condition is invariant under w ↦ w + G·E (which
  changes D by 6GF, a legal change of denominator with the same N), so we work in the slice {w_{0,m} = 0 : X0 | m}.
  If the mod-p ideal I_p of the slice system is m-primary (Hilbert function 0 in some degree) then, since
  h_{S/I_p}(t) ≥ h_{S/I_Q}(t) for all t (rank_p ≤ rank_Q on the Macaulay matrices), the Q-variety is empty: every
  denominator of degree d with that parity structure is a multiple of F (not coprime). Nothing else is assumed.

* (R4) Mixed parity. For D = D^+ + D^- the odd part of the relaxation reads Ξ^+(x) + Ξ^-(y) ∈ J_{2d}^- with NO cross
  terms (x = (λ^+,w^+), y = (λ^-,w^-)); the cross terms sit in the even part. Let L^± ⊂ R_{2d}^- be the spans of the
  coordinate polynomials of Ξ^±. Every linear functional φ on S_{2d}^- vanishing on J + L^+ kills Ξ^-(y) for every
  mixed solution, giving dim L^- − dim(L^+ ∩ L^-) quadrics in y alone (and symmetrically in x). If that ideal is
  m-primary mod p then y = 0 for every mixed solution over Q-bar (Hilbert-function inequality), and then Ξ^+(x) ∈ J
  forces x = 0 by the pure result. Rigor needs the functionals to be reductions of rational ones, i.e.
  dim_Q L^+ = dim_p L^+: certified by exact rational dependency certificates (c, h) with Σ c_k P_k = Σ h_i F_i,
  verified with fmpq_mpoly (`excl_lib.exact_dependencies`).

* (R5) Block elimination (generalises R4; `excl_blocks.py`). Because Hess(F) = 30 Σ_k l_k^4 ∇l_k ∇l_k^T has rank ≤ 8,
  Φ(w) = 30 Σ_k l_k^4 g_k(w)^2 couples the six component blocks of w only through the two small cross spans
  (l_7^4 ± l_8^4)·S. For a block B, functionals vanishing on J + span{pair polynomials not inside B} applied to Ξ give
  quadrics in the B-unknowns only; if m-primary in the unknowns that occur, those unknowns vanish on every solution and
  are dropped; the leftover ("invisible", typically the diagonal fields Z_i∂_i) directions go to a small residual
  system. Validated on d = 6 even: six blocks of 2–3 unknowns, residual 3 unknowns, verdict proj-empty (matches §1).

## 1. Results so far (prime 1073741789 unless stated; all Hilbert checks passed)

| exp | d | parity | unknowns (slice) | indep. quadrics | engine | result | status |
|---|---|---|---|---|---|---|---|
| relaxation | 5 | even | 3 | — | msolve + Singular | h = [1,3,2,0]: proj-empty | **rigorous** |
| relaxation | 6 | even | 17 | 39 | msolve + Singular agree | h = [1,17,114,358,501,225,0]: proj-empty (Singular vdim 1216) | **rigorous** |
| relaxation | 6 | odd | 18 | — | Singular (msolve count stopped at t=5) | dim 0, vdim 540 (wave 1: 540 with p=536870909), h(6)=0 | **rigorous** |
| relaxation | 6 | mixed | 35 | 132 | msolve: >100 s locally; slices P^k, k ≤ 26, all empty | full 35-variable GB not finished | screened |
| (R4) mixed d=6, y-side | 6 | mixed | 18 (all y visible) | 36 | msolve + Singular: dim 0, vdim 4096 | y = 0 for every mixed solution ⇒ x = 0 by pure even | **mod p; exact certificate (114 dependencies) queued** |
| (R4) mixed d=6, x-side | 6 | mixed | 12 of 17 visible | 18 | dim 0, vdim 729 | x ∈ Π_x = span of the 5 diagonal fields Z_i∂_i; restricted system on Π_x: h=[1,5,0] | consistent |
| L^+ ∩ L^- | 6 | mixed | dim L^+ = 39, dim L^- = 57, intersection 21 (mod p) | | | naive span argument insufficient | info |
| bilinear block | 6 | mixed | 17×18 | rank 57/306 | — | cross-term shortcut fails (Hessian of 8 sixth powers has rank ≤ 8) | info |
| E1 kernel d=1..3 | | both | | | own code, p=1073741789 | 0 (agrees with wave 1) | rigorous |

Combined with the E1-kernel results (K_d^±/J_d^± = 0 for d ≤ 6, wave 1 with p = 2^31−1, re-verified by the sibling job
excl-spp-e1k with p = 2147483629 for d ≤ 6 and by this session for d ≤ 3 with p = 1073741789):
**every certificate denominator that is Y-pure has degree ≥ 7, and every denominator has degree ≥ 6
(≥ 7 once the mixed d = 6 exact certificate is in).** (Degrees 1–4: J_d = 0 and K_d = 0, nothing to test;
degree 5: w constant, Φ(w) = 30 Σ_k l_k^4 (w·∇l_k)^2 ∈ J_4 = 0 forces w = 0 since the eight quartics l_k^4 are
independent — also confirmed computationally.)

## 2. Wave-1 claims: re-verification status
1. E1 kernel dims (p = 2^31−1): d ≤ 3 confirmed with p = 1073741789; d = 4..10 queued (job excl-spp-e1k).
2. d = 6 Hessian relaxation (p = 536870909, Singular): confirmed with p = 1073741789 by both msolve and Singular;
   odd vdim 540 reproduced.
3. Wave-1 d = 7 relaxation imposed Φ(w) ∈ J_8 — invalid since K_8^+ ≠ 0; replaced by the degree-14 system with the
   kernel element as an extra unknown (queued: excl-spp-relax-d7).
4. D7^+ full E2: the search lens ran it (inconsistent mod p, rank 4035, nullity 462); my trivial-nullity count is
   456, so 6 extra homogeneous solutions exist and that single test is modular evidence only. The d = 7 relaxation
   covers D7 + J_7 entirely and will supersede it.

## 2b. Other agents on this lane (for the record, not re-derived here)
* search lens: K_8^+/J and K_8^-/J both 6-dimensional (primes 2^31−1 and 2147483629); full fixed-D test of D7^+ inconsistent mod 2^31−1.
* verify lens: jobs verify-spp-e1k, verify-spp-relax6, verify-spp-d7exact (running).
* sibling excl instance: excl-spp-e1k (d ≤ 6 agree with wave 1 and with this session), excl-spp-quad6/7/8 (running).

## 3. Open / queued
* exclB-spp-relax-mixed (queued): exact D7^± + exact dependency certificates for mixed d=6; second prime; mixed d=7 by (R4).
* exclB-spp-kernel (queued): exact D7^±, K_8/K_9 generation by K_7, λ-only sub-families at d=8,9.
* exclB-spp-relax-d7 (queued): d=7 even/odd pure relaxation (61 unknowns + λ) with p = 1073741789, msolve + Singular.
* The remote is saturated (load ≈ 800 on 128 cores, ≈ 90 jobs); submissions wait for admission.
* d=8 pure: ~162 unknowns / 468 quadrics — semi-regular estimate puts the Gröbner regularity near 10 (infeasible);
  will be attempted only if the d=7 systems turn out much easier than generic.
