# Wave 2 — `polarized-k3-remainder` — lens EXCLUSION / RIGOR

STATUS: in progress — no witness (exclusion lens). Running summary (02:00 EST):
* Best rigorous bound for A (full target): deg D >= 8 (thread, degrees <= 6 now independently REPLICATED exactly; degree 7 replication running).
* NEW exact theorem: the K3-sector target A_K3 has no denominator of degree <= 6 (E1 injective through 7; degree 5 killed; degree 6
  reduces to D ~ f(U) which fails exactly); f(U), U0, V0 have zero K3-class projection of their full served space (exact over Q).
* Consequence (needs End_Hdg(T)=Q [thread] and eta not in NS_prim⊗C, automatic if rho=1): no product of denominators of degree <= 6 serves A;
  once degree 7 is done, a legal octic must have X ∩ V(D) irreducible and reduced (so the thread's spurious relaxation points f(U)h are all dead).
* Degree-8 fixed-D exact re-verification: f(U)U0^2 inconsistent in all 3 blocks (756/757 = thread, 828/829, 516/517), also for A_K3.
* Non-Jacobian K3 kernel at d=8: 51+39 dims; structured pairing over ker(eta|R_4) inconsistent (mod p, mu6-restricted 13572x2430).
* Remote saturated all night (load 700-1700 on 128 cores): everything above was computed locally; remote queue pending.

Started 2026-11-01 ~23:59 EST. Lane dir `/tmp/hodge/results/polarized-k3-remainder/`, packet dir `wave2-excl/`,
log `wave2-excl.log.jsonl` (one JSON line per experiment).

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

## DRAFT FINAL SUMMARY (to be finalized; numbers in [] pending)
1. Best rigorous bound for the instance target A: deg D >= 8 (thread); degrees <= 6 independently replicated exactly
   (p = 65521 and Q), degree 7 replicated [modular blocks + emptiness pending], degree-8 point f(U)U0^2 exact (3 blocks).
2. NEW: the K3-sector class A_K3 = P_0(A) has no denominator of degree <= 6 — exact theorem (E1 injective through 7;
   degree 5 killed; degree 6 ↦ D ~ f(U), which fails exactly 136/137). Degree 7: [relaxation survivors = f(U)·ℓ family
   (reducible) + ?; exact Q run queued].
3. The served space of f(U), U0, V0 has zero K3-class projection (exact over Q) — the join of the two copies of C
   carries no K3 class, as Ch. 3A predicts.
4. Logic repaired: additivity of served spaces over components means products can serve A even when no factor does; with
   the dichotomy (End_Hdg(T) = Q) and η ∉ NS_prim⊗C, items 2–3 exclude every product of degree-<=6 supports, and after
   degree 7: a legal octic must have X ∩ V(D) irreducible, reduced, singular, non-rational/non-normal.
5. Thread's degree-8 spurious family f(U)h: confirmed to lie in the relaxation (mod p) for both targets; all its members
   are reducible supports hence dead by item 4 (and f(U)U0^2 exact).
6. Non-Jacobian K3 kernel at d = 8: 51 + 39 dims (K4 = ker(η|R_4) has dim 9, large height); structured pairing
   Σ g_k(U)g_k(V) inconsistent (mod p, μ6-restricted 13572×2430); [more elements pending].
7. μ6-homogeneous octics (characters 0,1,2,3,5; <= 108 params): relaxation components = [pending] (goal: all reducible
   => no μ6-homogeneous legal octic, rigorous).
8. E1 tables re-verified at a second prime [pending chain3]; sum principle restated with its exact hypotheses.
9. Compute: remote saturated all night (load 700–1750 / 128 cores); all results above local; 4 remote jobs queued.
10. Packet: /tmp/hodge/results/polarized-k3-remainder/wave2-excl/ (README, 9 scripts, JSON receipts); draft post:
    /tmp/hodge/notes/wave2-polarized-k3-remainder-draftpost.md.
