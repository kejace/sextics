# Wave 2 — `short-power-presentation` — lens EXCLUSION / RIGOR

STATUS: in progress — no witness. **RIGOROUS: every certificate denominator has degree ≥ 7 (any Y-parity).**
(E1 kernels zero through degree 6 with four primes + Hilbert check; E2 relaxations m-primary for d ≤ 6 pure; mixed d = 6
by annihilator functionals with 114 + 114 exact rational dependency certificates; prime 1073741789, msolve + Singular.)
Degree 7, Y-even: λ = 0 (no D7^+ component) and the Y-coefficient blocks vanish; a 42-unknown residual is under Gröbner.

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
| (R4) mixed d=6, y-side | 6 | mixed | 18 (all y visible) | 36 | msolve + Singular: dim 0, vdim 4096 | y = 0 for every mixed solution ⇒ x = 0 by pure even | **rigorous: 114/114 exact dependencies verified (dim_Q L^+ = 39), remote job exclB-spp-relax-mixed** |
| (R4) mixed d=6, x-side | 6 | mixed | 12 of 17 visible | 18 | dim 0, vdim 729 | exactness: 114/114 dependencies (dim_Q L^- = 57) | rigorous (redundant with the y-side) |
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

### 1a. Degree 7, Y-even, block elimination (R5), prime 1073741789, pass 1 (local run, ~7 min)
61 unknowns in the slice (λ for D7^+ plus 60 field coefficients), 1891 pairs, span of all pair polynomials in R_14^- = 734 (of 744).
* Y-component blocks w3, w4, w5 (9 unknowns each): 18 quadrics in 6 visible unknowns, m-primary → eliminated; 3 invisible each.
* λ: the λ^2 coefficient B_{D7}(N7) survives all functionals → **λ = 0: no degree-7 Y-even denominator involves D7^+
  (every solution of the relaxation is a pure Jacobian denominator w(F))** (mod p; exact certificate pending).
* X-component blocks w0, w1, w2: 29 quadrics in 9 visible unknowns each, NOT m-primary in pass 1.
* Residual: 42 unknowns, 231 independent quadrics. A second pass and the single-unknown scheme ('fine') eliminate exactly
  the same set: the 18 Y-coefficients w_{k, X_i Y_j} (j ≠ k) and λ; the 9 "diagonal-type" Y-coefficients w_{k, X_i Y_k}
  and all 33 X-coefficients resist (their pair polynomials modulo J lie in the span of the other pairs). The residual family is
  w = (arbitrary quadratic X-components, 33 coefficients) + Σ_k a_k(X)·Y_k∂_{Y_k} (a_k linear in X, 9 coefficients):
  fields that are "diagonal in Y with X-linear coefficients". Random linear
  slices of the residual up to codimension 11 (P^30 ⊂ P^41) are all empty (exploratory, not a proof); msolve's F4 on the
  full residual: degree 3 (9680 × 13222), degree 4 (139915 × 116404, 7581 new, 897 s with 16 threads on the remote),
  degree 5 (959426 × 685669, in progress). The local attempt was OOM-killed (brunos has no free memory). Singular `std`
  on the same residual runs on the remote as a backup. The remote job exclB-spp-blocks repeated the d = 7 even chain (same
  block steps) but its residual msolve died under memory pressure; the standalone 16-thread msolve on the residual is
  the live attempt. Processing λ first (others = all w-pairs, span 304) also eliminates it, so its certificate only needs
  dependencies among the w-pair polynomials (moderate height). Exact certificates for the eliminating steps (λ, w3, w4,
  w5) are being produced by the replay job exclB-spp-certify7e (λ first, up to 12 primes for the lifts) (`excl_certify.py`: dim_Q span(other pairs mod J_14) = dim_p, rational dependency certificates).

### 1b. Degree 7, Y-odd, block elimination (remote job exclB-spp-blocks, prime 1073741789)
Mirror image of the even case: λ (the D7^- coefficient) is eliminated, and in the second pass the three X-component
blocks w0, w1, w2 are eliminated; the Y-component blocks resist. So **no degree-7 Y-odd denominator involves D7^-**
(mod p), and the residual lives in the Y-components (plus invisible X-directions).

### 1c. λ-only sub-families at d = 8 (remote job exclB-spp-kernel2, prime 1073741789)
For D in the span of the K_8^±/J bases (w = 0, 6 unknowns each) the relaxed E2 defect ideal is m-primary (Singular dim 0,
vdim 7 in both parities): no denominator in the kernel span itself works. K_9^+/J = 18 with S_2·K_7 spanning 15.

## 2. Wave-1 claims: re-verification status
1. E1 kernel dims (p = 2^31−1): d ≤ 3 confirmed with p = 1073741789; d = 4..10 queued (job excl-spp-e1k).
2. d = 6 Hessian relaxation (p = 536870909, Singular): confirmed with p = 1073741789 by both msolve and Singular;
   odd vdim 540 reproduced.
3. Wave-1 d = 7 relaxation imposed Φ(w) ∈ J_8 — invalid since K_8^+ ≠ 0; replaced by the degree-14 system with the
   kernel element as an extra unknown (queued: excl-spp-relax-d7).
4. D7^+ full E2: the search lens ran it (inconsistent mod p, rank 4035, nullity 462); my trivial-nullity count is
   456, so 6 extra homogeneous solutions exist and that single test is modular evidence only. The d = 7 relaxation
   covers D7 + J_7 entirely and will supersede it.

## 2b. Other agents on this lane (for the record, read from their logs; not re-derived here)
* E1-kernel dimensions now agree across FOUR primes: 2^31−1 (wave 1), 2147483629 (sibling excl job, d ≤ 7 even),
  1000000007 (verify lens, d ≤ 7 both parities), 1073741789 (this session: d ≤ 3 both parities, d = 7 even).
* verify lens: Singular re-run of the wave-1 d = 6 relaxation (no Euler slice): even dim 1 / mult 1216 (the Euler line),
  odd dim 0 / mult 540 — matches this session's sliced computation (vdim 1216 / 540).
* search lens (`wave2-search.log.jsonl`): K_8^-/J = S_1·K_7 exactly (Y_i D7^+, X_i D7^-); K_8^+/J = 6 of which S_1·K_7 gives
  only 4 because **L·D7^+ ∈ J_8 and M·D7^- ∈ J_8** (D7^± are residuals of degree-8 polar divisors containing the
  hyperplane sections L = 0 resp. M = 0), plus 2 new even elements; K_9^+/J = 18; Gorenstein-duality forced values
  dim Ann_10 = (300, 267), Ann_11 = (511, 511), Ann_12 = (677, 648). Full fixed-D test of D7^+ inconsistent mod 2^31−1
  (served space = J_6^- only). Structured points (D7 alone, polar-residual fields, pencils) fail the d = 7 relaxation.
* construct lens: exact identities (e.g. eX·M^6 − eY·L^6 = (L^2+M^2)·G10 − L^6·F), Fano scheme of lines, rationality of the
  τ-quotient discriminant curve, a non-domination lemma for X/τ.
* sibling excl instance: excl-spp-quad6/7/8 (Gröbner on the unsliced relaxations; no output yet after >1 h).

## 3. Open / queued
* exclB-spp-relax-mixed (queued): exact D7^± + exact dependency certificates for mixed d=6; second prime; mixed d=7 by (R4).
* exclB-spp-kernel (queued): exact D7^±, K_8/K_9 generation by K_7, λ-only sub-families at d=8,9.
* exclB-spp-relax-d7 (queued): d=7 even/odd pure relaxation (61 unknowns + λ) with p = 1073741789, msolve + Singular.
* The remote is saturated (load ≈ 800 on 128 cores, ≈ 90 jobs); submissions wait for admission.
* d=8 pure: ~162 unknowns / 468 quadrics — semi-regular estimate puts the Gröbner regularity near 10 (infeasible);
  will be attempted only if the d=7 systems turn out much easier than generic.


## 4. Summary for the orchestrator (updated as results arrive)
1. Best rigorous bound: **deg D ≥ 7 for every certificate denominator, any Y-parity** (E1 kernels zero through d = 6
   with four primes + Hilbert check; E2 relaxations m-primary for d ≤ 6 Y-pure; mixed d = 6 by exact-certified
   annihilator functionals).
2. Degree 7: Y-even: λ = 0 (no D7^+ component) and 18 Y-coefficients vanish; Y-odd: λ = 0 and the X-blocks vanish
   (mod 1073741789); the even residual (42 unknowns / 231 quadrics) is at F4 degree 5 on the remote (16 threads);
   exact certificates for the eliminating steps in progress (λ first).
3. E1 kernel: dim K_7^±/J = 1 exactly (exact rational D7^± with exact numerator certificates N7^± over Q, both parities,
   remote job exclB-spp-kernel2; the odd rref representative has 324 terms and is not itself a cone, the class has a cone
   representative in Q[Y,L] per the search lens); K_8^+/J = 6 of which S_1·K_7 spans 4, and K_8^-/J = 6 = S_1·K_7 exactly (this
   session, p = 1073741789, job exclB-spp-kernel2, agreeing with the search lens);
   D7^+ ∈ Q[X0,X1,X2,M] (cone over the line {X = M = 0} = pullback from the base of the plane-sextic fibration),
   D7^- ∈ Q[Y0,Y1,Y2,L]; K_8^- = S_1·K_7, K_8^+ = 6 with L·D7^+ ∈ J_8, M·D7^- ∈ J_8 (search lens).
4. Methods: (R1) Hilbert-series prime check, (R2) modular rank inequality, (R3) E2 relaxation with the Euler slice,
   (R4) parity annihilator functionals, (R5) block elimination exploiting rank(Hess F) ≤ 8, exact dependency
   certificates with fmpq_mpoly; msolve `u*u` parsing pitfall avoided (`u^2`).
5. Obstacles: d = 8 pure needs ≈160 unknowns; the relaxation systems are far from semi-regular (regularity ≥ 5–6),
   so only block/functional reductions make degree ≥ 7 feasible.
