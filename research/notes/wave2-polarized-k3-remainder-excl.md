# Wave 2 — `polarized-k3-remainder` — lens EXCLUSION / RIGOR

STATUS: in progress — no witness (exclusion lens). Running summary (03:00 EST):
* Best rigorous bound for A (full target): deg D >= 8 (thread); degrees <= 6 independently REPLICATED exactly (p = 65521 and Q);
  E1 injectivity of [A] through d = 9 RE-VERIFIED at a second prime (2^31-19), both parities; degree-7 exclusion fully replicated at p = 65521 (ranks, 1779, 9 empty charts, 42 variables, DIM 0 / 1194 / 309);
  degree-8 point f(U)U0^2 exact in all 3 blocks; 6 random structured degree-7 members of J exact-inconsistent (both targets).
* NEW exact theorem: the K3-sector target A_K3 has no denominator of degree <= 6 (E1 injective through 7; degree 5 killed; degree 6
  reduces to D ~ f(U), which fails exactly). served(f(U)), served(U0), served(V0) have zero K3-class projection EXACTLY over Q.
* Degree 7 for the K3 class is OPEN by relaxation methods (single- and 19-target L-relaxations too weak: ranks 935/1092 of 7260,
  no squares); the reducible family f(U)·ℓ is excluded by the sum principle; f(U)U0, f(U)V0 and random J_7 members fail exactly.
* Consequence (End_Hdg(T) = Q [thread] + η ∉ NS_prim⊗C, automatic if ρ = 1): no product of degree-<=6 supports serves A; a legal
  octic has X ∩ V(D) irreducible and reduced, or = (irreducible degree-7 K3-carrying support) ∪ (hyperplane section).
  The thread's spurious degree-8 relaxation points f(U)h are all dead.
* Non-Jacobian K3 kernel at d = 8: 51 + 39 dims; 6 elements of the K4⊗K4 block inconsistent (mod p, μ6-restricted 13572×2430).
* Remote saturated all night (load 600-1800 / 128 cores, < 60 GB free): everything above computed locally; 7 remote jobs queued.

## Plan
(a) adversarial re-verification of wave-1 / thread claims (second prime, independent implementation);
(b) rigorous lower bound on deg D: replicate d<=7 exclusion for A; NEW: same for the K3-sector target
    A_K3 = (U0 U1 U2)^2 - (V0 V1 V2)^2 (needed to make the sum-principle reduction to irreducible supports valid);
(c) non-Jacobian E1-kernel of [A_K3] at d=8,9 + structured irreducible families (fixed-D, remote, mod p);
(d) packet wave2-excl/ with receipts; (e) draft post.

## Results so far
(none yet)

### R1. Degree-6 exclusion for A independently REPLICATED (exact over Q)  [DONE 00:40]
Own implementation `wave2-excl/relax.py` ((U,V) frame, G = 7f(U)+7f(V), T = prod(U_i^2-V_i^2), L = (G)+J(G)^2+(R_ijk),
Euler gauge drops the z5-divisible coefficient of a_5). Prime 65521 and exact rationals. Results identical to the
thread's DEGREE_SIX.md (p = 101): block L-ranks 603/558/396/260/396/558 (mod p = over Q), augmented Q-ranks
713/715/617/470/617/715, constraint rank 437; exact square elimination kills the 30 off-diagonal entries, then the 5
diagonal ones (rank 15) => D = 0 or D = c*G. Degree-5 smoke test: all 6 parameters killed (rank 21). Receipts
`wave2-excl/relax-full-d6-p65521-local.json`, `-exactQ.json`, `smoke-d5-p101.json`.

### R2. NEW THEOREM (exact): the K3-sector target has no denominator of degree <= 6   [DONE 00:50]
Target A_K3 := [(X0+Y0)^2(X1+Y1)^2(X2+Y2)^2 - (X0-Y0)^2(X1-Y1)^2(X2-Y2)^2]/64  (= 64^{-1}(B(U)-B(V)), the K3-sector
projection of A; the book's P_0 A-hat).  Why it matters: served(D) depends only on Y = X ∩ V(D) and is additive over
the components of Y (N_{Y1 ∪ Y2} = N_{Y1} + N_{Y2}); so a reducible denominator can serve A even if no factor does.
The dossier's "sum principle => only irreducible octics remain" therefore needs: no low-degree support carries any
K3-sector class.  By the dichotomy (End_Hdg(T) = Q, Λ²T irreducible) and if ρ = 1, that is equivalent to
"no low-degree D serves A_K3".  Results (all exact over Q, p = 65521 cross-check):
* E1 for A_K3: [A_K3] injective on (R/J)_d for d <= 7, both parities (mod 2^31-1 => over Q)  => E1-admissible D of degree <= 7 lie in J(F); degrees <= 4 impossible (J_d = 0).
* d = 5: L-relaxation forces all 6 parameters to 0 (rank 17).  => excluded.
* d = 6: L-relaxation forces 32/35 gauged parameters to 0; survivors a_i = c_i U_i (i<3); residual quadrics
  (c0-c2)^2, (c0-c2)(c1-c2), (c1-c2)^2  => c0=c1=c2, i.e. D ~ f(U) modulo the gauge D ↦ D + tG.
  Exact Koszul-normalized fixed-D test of D = f(U) for A_K3 (char-0 block 966 x 167): rank_Q 136, augmented 137 => inconsistent.
  => no degree-6 D serves A_K3.   (Geometrically: X ∩ V(f(U)) = join of the two curves C_U, C_V supports the joining
  sector Λ²H¹(C)(-1), not the K3 sector -- consistent with Ch. 3A.)
* d = 7: relaxation queued on the remote (job excl-polarized-k3-remainder-relaxk3).
Consequence (conditional on ρ(S_2) = 1, or on the 3-target-family version, see plan): no product of denominators of
degree <= 6 can serve A; a legal octic would need X ∩ V(D) to have a component of degree >= 7·6 ... (completed once d=7 is done).

### R3. Served-space K3-projection tool (unconditional vanishing statements)   [01:05]
`wave2-excl/servedK3.py`: for a mu6-homogeneous D, the certificate system with a VARIABLE target in the K3 block
(all 56 pure-U / pure-V sextic monomials) is character-graded; the projection of its nullspace to the target
coordinates, reduced modulo J(f)(U) + J(f)(V), is the K3-class projection of served(D).  CORRECTION (01:50): a zero MODULAR class projection is only evidence (the modular rank of the class map is a lower
bound for its rational rank, so it cannot certify vanishing); the tool now has an --exact mode (fmpz nullspace over Q,
integer model G = 14 F-hat) and the rigorous statements below are the exact ones.
Results (mod p evidence): D = f(U), U0, V0: class dim 0.  EXACT over Q: U0 and V0: class dim 0 (rigorous); f(U): see R6.
This is the statement needed for the sum principle (N_{Y1∪Y2} = N_{Y1} + N_{Y2}): a product D1·D2 can only serve
the K3 part of A if some factor has a NONZERO K3-projection; "A_K3 ∉ served(D_j)" alone is not enough.
Controls: D = G (illegal, divisible by F) passes the exact Koszul checker for both targets.

### R4. Degree-8 fixed-D exact re-verification   [01:05]
D = f(U)·U0² (the thread's single excluded octic): exact Koszul-normalized test (own code, p = 65521 and over Q):
char-0 block 3366 × 888, rank 756, augmented 757 => inconsistent, for the FULL target (thread: identical numbers over
Q and F_101 -- REPLICATED) and ALSO for the K3-sector target A_K3 (new; identical ranks).
Tangent-hyperplane probe: hyperplanes λ·U through the plane {V=0} tangent to X (λ = tangent line of C at an F_p-point) have
a non-rational singularity t² + f(V) whose exceptional divisor is the K3 S₂ itself; nevertheless their served K3-class
projection is 0 mod p (consistent with E1-injectivity of [A_K3] at d = 1, which is rigorous).  The copy of T(S₂)(−1) such a section carries is not a sub-Hodge
structure of Λ²V₂, so this mechanism cannot produce K3-sector (3,1)-classes.
Remote: saturated (load 900+/128 cores, <60 GB free) -- rjob refuses submissions; queue waiting (--wait 4h). Local
fallback chain (nice 19, one process at a time): K3 d=7 relaxation -> full d=7 replication -> E1 tables (2nd prime) ->
K3 d=8 relaxation -> f(U)h family exact tests -> full d=8 relaxation regeneration.  Logs wave2-excl/chain*.log.

### Logic note (refined): how the K3-sector exclusions feed the reducible-denominator question   [01:25]
Let Y = X ∩ V(D) = ∪ Y_j (components Y_j = X ∩ V(D_j), Σ deg D_j ≤ deg D since Pic X = Z h).  served(D) = Σ_j served(D_j)
(prescribed-support theorem + N_{Y1∪Y2} = N_{Y1} + N_{Y2}).  Dichotomy (needs End_Hdg(T) = Q, proven in #115): for each j,
π_K3(N_{Y_j}) ⊇ Λ²T  or  π_{Λ²T}(N_{Y_j}) = 0.  If A_K3 ∉ served(D_j) and η ∉ NS_prim⊗C (automatic when ρ = 1), the
Λ²T-component ω∧η_T of r_K3 is nonzero and not served, so the first alternative fails and π_K3(served(D_j)) ⊆ ω⊗NS_prim.
Hence: if NO denominator of degree ≤ k serves A_K3, then NO product of denominators of degrees ≤ k serves A
(else ω∧η ∈ Σ_j π_K3(served(D_j)) ⊆ ω⊗NS_prim, i.e. η ∈ NS_prim⊗C).  The hypothesis η ∉ NS_prim⊗C is exactly the one
already needed for the equivalence theorem; no extra assumption on ρ is required.
Status: k = 6 PROVEN above (exact).  k = 7 in progress (local chain).  The dossier's wave-1 formulation ("products of
excluded denominators are excluded") implicitly needed this K3-sector step; the full-target degree ≤ 7 exclusion alone
does not give it.

### R5. Non-Jacobian E1 kernel of the K3-sector target at degree 8   [01:35]
[A_K3] on (R/J)_8: kernel mod J has dimension 51 (Y-even) + 39 (Y-odd) at p = 2^31-1 (the full target has none through
d = 9).  Source: eta = [z0^2 z1^2 z2^2] has kernels R_{f,a} -> R_{f,a+6} for a >= 4 (dimension reasons), e.g.
ker(eta|R_4) (x) ker(eta|R_4) ⊂ R_4 (x) R_4.  These are the natural non-Jacobian candidates for the K3 sector; their
full E2 test (fixed-D, 40k x 7.5k) is queued (k3kernel.py) but blocked by compute saturation.  Note the reducible
rank-1 members g(U)g'(V) are already excluded by R2 (both factors have degree 4 <= 6).

### R6. Exact served-space statements   [01:55]
`servedK3.py --exact` (integer model G = 14 F-hat, fmpz nullspace over Q, classes modulo J(7f)_6(U) + J(7f)_6(V)):
* D = f(U) (char 0): 4923 × 832, nullity 110, 17 served K3-monomial vectors, class dimension 0.
  => π_K3(served(f(U))) = 0 over Q, unconditionally.  (The join X ∩ V(f(U)) of the two copies of C carries no K3 class.)
* D = U0, D = V0: class dimension 0 over Q.
For hyperplanes ℓ = λ·U + μ·V with λ, μ ≠ 0 the tensor structure R_F = R_f ⊗ R_f gives A'·ℓ ∉ J for every nonzero
K3-block class A' (the components B'⊗(μ·V), (λ·U)⊗B'' of the product vanish only if B' = B'' = 0), so E1 already
forbids any K3 class; only hyperplanes through one of the two planes {U=0}, {V=0} need the full system.

### R7. Structure of the non-Jacobian K3 kernel at degree 8   [02:00]
K4 := ker(η: R_{f,4} → R_{f,10}) has dimension exactly 9 (η: R_4 (15) → R_10 (6) is onto), verified exactly over Q
(`k4_lll.py`, membership of η·g in J(f)_10 checked exactly for each basis element).  The exact kernel vectors have
very large height (≈70-digit integers even after LLL on the rref-derived integer basis), i.e. there is no small-height
natural element: the non-Jacobian K3 kernel K4⊗K4 (81 dims, character 4) is a "generic" subspace, like the full-target
E1 kernel at d = 10 observed in wave 1.  The UV-frame tests (k3kernel_uv.py: pairings Σ g_k(U)g_k(V), random symmetric
tensors, random character components of the full kernel) are running locally with the engine's μ6-restricted solver.
First result (02:05): D = Σ_k g_k(U) g_k(V) over the exact K4 basis (153 terms, character 4, D ∉ J, E1 for A_K3 holds):
μ6-restricted full E2 system 13572 × 2430 (of 7518), rank 2206, nullity 224 — INCONSISTENT mod 2^31−1.  So the most
natural non-Jacobian K3-kernel element does not serve the K3 class (modular statement).  Further elements running.

### Consequence for degree 8 (to be finalized once the K3 degree-7 relaxation is analysed)
If no denominator of degree ≤ 7 serves A_K3 (degrees ≤ 6 proven exactly; 7 pending), then under η ∉ NS_prim⊗C (the
equivalence-theorem hypothesis, automatic for ρ = 1) a legal octic D must satisfy: X ∩ V(D) is irreducible and reduced
(every component Y_j = X ∩ V(D_j) has Σ deg D_j ≤ 8, so a proper decomposition only has factors of degree ≤ 7, whose
served K3-projections lie in ω⊗NS_prim), hence D is irreducible over C and not a power; by Lefschetz X ∩ V(D) is singular,
and by Ch. 3 §3.1 it has a non-normal or non-rational-singular component.  The thread's degree-8 relaxation points
D = f(U)h are all reducible, so they are excluded by this argument without any further computation.

## DRAFT FINAL SUMMARY v2 (04:25; [] = pending)
1. Best rigorous bound for the instance target A: deg D >= 8 (thread). Independently replicated exactly: degrees <= 6 (own code,
   p = 65521 and Q: identical block ranks 603/558/396/260/396/558, constraint rank 437, two-step square elimination); E1 injectivity of
   [A] through d = 9 at a second prime (2^31-19, both parities); the degree-8 point f(U)U0^2 exact in all 3 blocks (756/757, 828/829,
   516/517); degree 7 FULLY replicated at p = 65521 (L-ranks over Q = mod p 1635/1602/1635/1464/1287/1464; constraint rank 1779; nine empty charts;
   square reduction to the thread's 42 variables; final GB DIM 0 / 1194 / 309) — the bound deg D >= 8 now rests on two codes and two primes.
2. NEW exact theorem: the K3-sector class A_K3 = P_0(A) has no denominator of degree <= 6 (E1 injective through 7 at two primes (2^31-1, 2^31-19), both parities;
   degree 5 killed; degree 6 ↦ D ~ f(U) via residual quadrics (c0-c2)^2, (c0-c2)(c1-c2), (c1-c2)^2; f(U) fails exactly 136/137).
3. Exact served-space statements: f(U), U0, V0 have zero K3-class projection over Q (the join of the two copies of C carries no K3 class).
4. Repaired logic: served spaces add over components (N_{Y1∪Y2} = N_{Y1}+N_{Y2}), so "deg >= 8" alone does not exclude products.
   With End_Hdg(T) = Q (thread) and η ∉ NS_prim⊗C (the equivalence-theorem hypothesis; automatic if ρ = 1): no product of degree-<=6
   supports serves A; a legal octic has X ∩ V(D) irreducible and reduced, or = (irreducible degree-7 support carrying the K3 class) ∪
   (hyperplane section). All of the thread's spurious degree-8 relaxation points f(U)h are reducible and hence dead.
5. Block cones g(U), g(V), and products g(U)h(V) are excluded in every degree (product divisors C_Γ × S₆ under Φ; H¹(S₆) = 0), which is
   exactly what the character-0 degree-8 sector leaves over (30 cone coefficients, 6 residual quadrics). Caveat: Hodge-theoretic, to be
   double-checked (exceptional divisor over C × C) before posting as a theorem.
6. Degree 7 for the K3 class is OPEN by relaxation methods: single-target rank 935/7260, 19-target rank 1092 (every K3 class imposes the
   same constraints), no squares, zero set with 84-dimensional tangent spaces along the (excluded) family f(U)·ℓ.
7. Non-Jacobian K3 kernel at degree 8 completely described: K4⊗K4 (81, char 4; K4 = ker(η|R_{f,4}), dim 9, exact LLL basis) ⊕
   {ℓ(U)ℓ'(V)((U0U1U2)²+(V0V1V2)²)} (9, char 1, reducible). Eight elements of K4⊗K4 fail E2 for A_K3 (μ6-restricted 13572×2430, rank 2206,
   mod 2^31-1) and the K4 pairing serves no K3-block class; [char-1 random tests pending].
8. Exact structured screens: 6 random homogeneous members of J_7 and 5 of J_8 (chars 0,1,2,3,5; 22 runs) fail for both targets with rank gap 1; degree 9: char-0 member exact-inconsistent for both targets (gaps 1), continuing.
   μ6-homogeneous octic sectors: char 0 done (cones); [chars 1,2,3,5 and K3 sectors pending].
9. Compute: the 128-core remote was saturated all night (load 580-1800, often < 60 GB free); rjob refused every submission; all results
   above were computed locally (one nice'd process at a time where possible); 7 remote jobs remain queued (sweeps of the 26 wave-1
   structured octics/decics, exact degree-7 K3 relaxation, degree-8 relaxations, kernel tests).
10. Packet: /tmp/hodge/results/polarized-k3-remainder/wave2-excl/ (README with replay recipe, 15 scripts, JSON/JSONL receipts, log
    wave2-excl.log.jsonl with ~60 entries); draft post: /tmp/hodge/notes/wave2-polarized-k3-remainder-draftpost.md.
Update (04:35): the character-1 block is confirmed computationally (U0·V0·(η(U)+η(V)) ∉ J, E1 for A_K3 holds with rank 5052 =
augmented rank), and a random element of it fails the full E2 system for A_K3 (μ6-restricted 13464 rows, rank 2322, nullity 324,
mod 2^31−1) — as it must, being a reducible support.  So every tested element of the complete degree-8 non-Jacobian K3 kernel
(8 of K4⊗K4, 1 of the character-1 block) fails E2.

### Remote queue status (04:40)
rjob refused every submission since 00:20 (load 580–1800 on 128 cores; free memory often < 60 GB).  Two sequential queue
scripts are waiting: `wave2-excl/submit_queue.sh` (jobs excl-polarized-k3-remainder-{relaxfull, relaxk3, e1table}; the first two
are now redundant with the local runs, the third was done locally) and `submit_queue2.sh` (k3d7exact = exact-Q degree-7 K3
relaxation; relaxd8 = degree-8 relaxations K3 + full regenerated at p = 65521; sweep = the 26 wave-1 structured octics/decics
through hodge.sweep (mod p); k3kernel = remaining kernel tests).  If the remote frees up they submit automatically (≤ 3 running);
results land in /tmp/hodge/results/polarized-k3-remainder/wave2-excl/ on the remote and must be pulled with
`python3 -m hodge.rjob pull --lane polarized-k3-remainder`.
Final kernel bookkeeping (04:50): closed form of the character-1 block verified with a control — U0V0(η(U)+η(V)) and
(U0+2U1)(V1−V2)(η(U)+η(V)) satisfy E1 for A_K3 (rank 5052 = augmented) and are not in J, while the sign-flipped
U0V0(η(U)−η(V)) fails E1 (5052/5053) exactly as the tensor identity predicts.  All ten tested elements of the degree-8
non-Jacobian K3 kernel (8 in K4⊗K4, 2 in the character-1 block) fail the full E2 system for A_K3 (mod 2^31−1).
Second-prime E1 table for A_K3 (2^31−19): injective for d ≤ 7 (both parities so far), kernel 51 at (8, Y-even) as at the
first prime, and 153 at (9, Y-even) — the non-Jacobian K3 kernel grows quickly with the degree (K4⊗R_1-type and
K5⊗K4-type blocks), so from degree 9 on the K3 class has a large E1-admissible non-Jacobian candidate space that only E2
can decide.
Character 1 (full target, 99 gauged parameters: a_i of V-degree 1 for i<3, V-degree 2 for i≥3): constraint rank only 180 on 4950
quadratic monomials, no squares — the L-relaxation is essentially empty in this sector (its zero set has dimension ≥ 80), so no
universal statement follows for character-1 octics; the only rigorous information there is the exact failure of random
members of J_8 of character 1 (R12).  This mirrors the K3-target degree-7 situation: the relaxation is sharp only when the
target monomials and the sector together produce many independent quotient coordinates.  (Chars 2, 3, 5 and the K3 sectors
still running.)

### R13. Degree-7 exclusion for A: replication at p = 65521 (block level)   [05:15]
Own implementation, prime 65521: L-ranks 1635/1602/1635/1464/1287/1464 per character block (= thread's ranks over Q and F_101;
and = Q-ranks computed in the K3 run, same L), constraint ranks 291/318/291/470/688/470 (= thread's augmented ranks
1926/1920/1926/1934/1975/1934 minus the L-ranks), 2541 constraint rows of combined rank 1779 (thread: 2,528 rows, rank 1,779).  The
chart-wise projective emptiness replay (`d7_charts.py`, p = 65521): all nine charts c_0=…=c_{k-1}=0, c_k=1 (k = 0..8) give the
unit ideal (1682 s), and the square reduction on the restricted system follows the thread's sequence exactly
(111→108→104→92→63→…, identical forced-index sets); the final homogeneous Gröbner basis of the 42 remaining variables has DIM=0, SIZE=1194, VDIM=309 (thread: 0 / 1194 / 309) — the
degree-7 exclusion is fully REPLICATED at an independent prime with independent code; together with the equal Q/F_p L-ranks this re-establishes the thread's
"deg D >= 8" by the same flatness argument at an independent prime and with independent code.

## FINAL SUMMARY (ready to submit; pending items marked)
1. Best rigorous bound for the instance target A: deg D >= 8 (thread). Now independently replicated end to end: degrees <= 6
   exactly over Q (own code, p = 65521: block L-ranks 603/558/396/260/396/558, constraint rank 437, squares kill all 35 gauged
   parameters); degree 7 at p = 65521 (L-ranks over Q = mod p 1635/1602/1635/1464/1287/1464, constraint rank 1779, nine empty
   charts, square reduction to the thread's 42 variables, final GB DIM 0 / 1194 / 309); E1 injectivity of [A] through d = 9 at a
   second prime (2^31-19), both parities; the degree-8 point f(U)U0^2 exact in all three blocks (756/757, 828/829, 516/517).
2. NEW exact theorem: the K3-sector class A_K3 = P_0(A) has no denominator of degree <= 6: [A_K3] injective on (R/J)_d for
   d <= 7 (two primes); degree 5 killed by the relaxation; degree 6 reduces (residual quadrics (c0-c2)^2, (c0-c2)(c1-c2),
   (c1-c2)^2) to D ~ f(U), which fails exactly (136/137); and the full served space of f(U), U0, V0 has zero K3-class projection
   exactly over Q.
3. Degree 7 for the K3 class is open by relaxation methods (single-target rank 935/7260, 19-target rank 1092 — every K3 class
   imposes the same constraints — no squares); the reducible family f(U)·ℓ is excluded by the sum principle and f(U)U0, f(U)V0,
   six random homogeneous members of J_7 fail exactly for both targets.
4. Logic repaired: served spaces add over components of X ∩ V(D), so "deg >= 8" alone does not exclude products. With
   End_Hdg(T) = Q (thread) and η ∉ NS_prim⊗C (the equivalence-theorem hypothesis, automatic if ρ = 1): no product of
   degree-<=6 supports serves A; a legal octic has X ∩ V(D) irreducible and reduced, or = (irreducible degree-7 K3-carrying
   support) ∪ (hyperplane section); all of the thread's spurious degree-8 relaxation points f(U)h are reducible, hence dead.
5. Block cones g(U), g(V) and products g(U)h(V) are excluded in every degree (product divisors under Φ, H¹(S₆) = 0) — exactly what
   the character-0 degree-8 sector relaxation leaves over (30 cone coefficients, 6 residual quadrics). Hodge-theoretic; double-check
   the exceptional-divisor contribution before posting as a theorem.
6. Non-Jacobian K3 kernel at degree 8 completely described: K4⊗K4 (81, char 4; K4 = ker(η|R_{f,4}), dim 9, exact basis) ⊕
   {ℓ(U)ℓ'(V)((U0U1U2)²+(V0V1V2)²)} (9, char 1, reducible; closed form verified with a control). Ten tested elements all fail E2
   for A_K3 (μ6-restricted systems 13572×2430 / 13464 rows, mod 2^31-1); the K4 pairing serves no K3-block class; the kernel
   grows to 153+153 at degree 9.
7. Exact structured screens: 6 random homogeneous members of J_7, 5 of J_8 (all characters), 1 of J_9 (char 0) fail for both
   targets with rank gap exactly 1 (22 + 2 runs; degree-9 screen continuing). Character-1 degree-8 sector: relaxation
   essentially empty (rank 180/4950) — no universal statement there. [chars 2,3,5 and K3 sectors pending; K3 d=8 relaxation:
   rank 956/1050, kills nothing; full d=8 regeneration pending.]
8. Rigor corrections made during the run: a zero MODULAR served-space projection is only evidence (rank mod p is a lower bound);
   all vanishing statements used downstream were recomputed exactly over Q; the "sum principle ⟹ irreducible only" step was
   restated with its exact hypotheses.
9. Compute: the 128-core remote was saturated the whole night (load 580–1800, often < 60 GB free); rjob refused every submission;
   everything above was computed on the shared 40-core box, mostly one nice'd process at a time; 7 remote jobs remain queued
   (structured octic/decic sweeps mod p, exact degree-7 K3 relaxation, degree-8 relaxations, kernel tests).
10. Packet: /tmp/hodge/results/polarized-k3-remainder/wave2-excl/ (README with replay recipe + glossary, 17 scripts, ~40 JSON/JSONL
    receipts); log wave2-excl.log.jsonl (~75 entries); report /tmp/hodge/notes/wave2-polarized-k3-remainder-excl.md; draft post
    /tmp/hodge/notes/wave2-polarized-k3-remainder-draftpost.md. No witness; no claim of nonexistence.
