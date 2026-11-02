# wave2 verification ledger: one-reflection-skew (adversarial referee)

STATUS: in progress — confirmed 0 / refuted 0 / weakened 0 / unverified 0 (updated continuously; see table)

Scope: every claim labelled rigorous/proved/exact/computed in `dossier-one-reflection-skew.md`, the lane logs under
`/tmp/hodge/results/one-reflection-skew/`, and any `attack1-`/`wave2-one-reflection-skew-*` notes and draft post as they appear.
Independent tooling: own polynomial code (`/tmp/hodge/verify/ors/vlib.py`, no engine import), python-flint `nmod_mat`
at primes 1000003 / 999983 (lane used 2^31-1, dossier used 32003), Singular 4.3.2 (`jac.sing`), msolve where needed.
JSON lines for every check: `/tmp/hodge/results/one-reflection-skew/wave2-verify.log.jsonl`.

Modular-logic reminders applied throughout:
* rank_p(M) <= rank_Q(M).  "Full rank mod p" => full rank over Q (OK).  "Kernel of dimension k mod p" => Q-kernel has dimension <= k; equality needs an exhibited Q-kernel of dimension k (e.g. J_d inside (J:A)_d).
* "Inconsistent mod p" alone does NOT prove Q-inconsistency (example: p x = 1).  It does when rank_p(M) = rank_Q(M), e.g. certified by an exhibited Q-kernel of the right size, or by an integer dual witness y with y^T M = 0 over Z and y^T b != 0.
* Good-prime condition for "independent modulo J over Q" from a mod-p computation: dim_p J_n must equal the char-0 value (Hilbert function of the complete intersection).

| # | claim | source | verdict | evidence | notes |
|---|---|---|---|---|---|
| 1 | F, A, q parsed from instance.rec equal the dossier's 10-monomial table | dossier 1.1 | CONFIRMED | own sympy parse == own table | |
| 2 | diagonal symmetry group: 144 elements of (Z/6)^6, 24 mod scalars, = <mu6(Y2), rho, sigma> ≅ Z6×Z2×Z2 | dossier 1.3 [computed] | CONFIRMED | own brute force; order statistics 1/7/2/14 | |
| 3 | only id and Y0<->Y1 preserve the support; no 6th-root diagonal twist makes Y0<->Y1 a symmetry | dossier 1.3 [computed] | CONFIRMED | own brute force over (Z/6)^6 with global 6th-root scalar | |
| 4 | no rescaled exchange with arbitrary complex scalars and F->κF | dossier 1.3 [proved] | CONFIRMED | hand: coupling/skew ratio gives λ0²=−λ3², sixth powers force κ=−κ | |
| 5 | residue characters of A: ζ (mu6), −1 (rho), +1 (sigma), +1 (tau), −1 (Y2-reflection) | dossier 1.3 [computed] | CONFIRMED | ζ^{<w,e>+Σw} recomputed | |
| 6 | Fermat-count Hodge numbers: target sector (0,44,77,3,0); conjugate (0,3,77,44,0); residual sector (0,126,492,126,0); tau-inv (0,208,866,208,0) | dossier 1.3 [computed] | CONFIRMED | own count over e∈[0,4]^6 | constancy along the diagonal-invariant smooth family is the standard semicontinuity+constancy-of-sum argument; fine |
| 7 | smoothness loads (23/7,30/7,1,18/7,18/7,0) < 6; V(F) and V(f) smooth | dossier 1.5 [proved] | CONFIRMED | own recomputation of loads; estimate checked by hand | |
| 8 | "the only two-element sign-change symmetries are {X1,Y1},{Y0,Y1}; they share Y1" | dossier 1.4 [computed] | REFUTED as stated / conclusion holds | own enumeration: pairs are {X1,Y0},{X1,Y1},{Y0,Y1} (three, pairwise intersecting) | no three disjoint pairs, so "no paired-sign cubic reduction" stands; the dossier itself lists (X1,Y0)=σρ in 1.3 |
| 9 | skew monomials have (5,1,5)-weights 4 and 2 mod 6 | dossier 2 / why.md | CONFIRMED | own | |
| 10 | Fermat cover: G(T,y)−G(y,T) = (αX0X1/7)(α²+1)yT(y²−T²) | dossier 2 / fqlx C5 | CONFIRMED up to sign convention | sympy gives the negative; magnitude/shape exact | harmless |
| 11 | rescaled exchange on the Fermat cover forces c⁶=1 and α²=−c², hence α⁶=−1 | dossier 2 / fqlx, DiamondJM | CONFIRMED | sympy on the skew term + hand on pure powers | |
| 12 | X-twist parity obstruction a+b≡3, 3a+b≡0 (mod 6) has no solution | DiamondJM C7 | CONFIRMED | enumeration | |
| 13 | Jacobian ring: vdim 15625, Hilbert function 1,6,21,…,1751,…,1; A ∉ J_6 | dossier 1.2 [computed] | CONFIRMED | Singular, char 1000003 and 32003 (different tool and primes) | char-0 std(J) still running |
| 14 | E1 kernel dims (J:A)_d = 6,36,135,480,1144,2196,3711,5762,8322,11502 (d=5..14); dim J_d = 6,36,126,336,756,1497,2682,4437,6882,10122 | dossier 4 table [computed mod 32003] | CONFIRMED | Singular quotient ideal (J:A), char 1000003 and 32003; own flint code p=1000003 agrees for d≤8 so far | also (J:A)_15=15448, (J:A)_16=20328 |
| 15 | 24-character decomposition of the E1 kernel at d=7 and d=8 (rows e5=0..5) | dossier 4 table | CONFIRMED | own flint code p=1000003: d=7 rows [30,20,19,20],[8,6,6,6],[2,1,1,1],0,0,[6,3,3,3]; d=8 rows [97,79,70,79],[30,20,19,20],[8,6,6,6],[2,1,1,1],0,[11,8,8,8] | exact match |
| 16 | 5-var (J_f:A)_d = 5,25,84,310,594,931,1330,1805 (d=5..12), J_f,d = 5,25,75,175,350,620,1000,1500; rho-parity split d=7: 47/39/8 and 37/36/1; d=8: 167/89/78, 143/86/57; d=9: 307/178/129, 287/172/115; d=10: 479/314/165, 452/306/146 | lane log 5var:e1-kernel-sweep; search.md e1-kernel5 | CONFIRMED | Singular (J_f:A) two primes; own flint p=1000003 by (rho,sigma) blocks | fqlx's (J_f:A)_7 = 84 also confirmed |
| 17 | five-variable reduction: certificate exists iff 6/5-equation A−u(f)+(6/5)f div u=0 has a rational solution u∈C(x)^5 | dossier 1.4 [proved up to last step] | CONFIRMED and STRENGTHENED | hand proof (see §A below) | the "last step" the dossier could not reproduce is fine: for the homogeneous problem div u has degree 0, so the Euler shift u'=u−(div u/5)E preserves the equation and gives the strict system u'(f)=A, div u'=0; for a Y2-free D one even keeps the same denominator (D | B_D(α) is forced) |
| 18 | lane 5-var fixed-D tests: for Y2-free D, "6/5 system inconsistent" ⟹ no 6-var certificate with that D | five.py docstring; all 310 engine-2 log rows | LOGIC CONFIRMED, RIGOR WEAKENED | hand proof of the direction 6-var(D) ⟹ 6/5(D) (§A); results are single-prime (2^31−1) mod-p inconsistencies | mod-p inconsistency alone is not a Q-proof (needs rank_p = rank_Q); treat as strong evidence, not exclusion; re-tests at p=999983 in progress |

## A. Hand proofs checked by the referee

**A.1 Even-D / parity reductions (dossier 1.2) are lossless.** If (D,N,Q) is a certificate then so is (DG, NG, QG²) for any G
(B_{DG}(NG) = G² B_D(N)); taking G = τ(D) gives a τ-even denominator.  For fixed τ-even D the equations are linear in (N,Q),
τ acts on solutions with N_i ↦ s_i N_i∘τ (s_i = −1 on Y-components) and A∘τ = −A, so (N − N')/2 is a solution with
N_0..N_2 Y-odd and N_3..N_5 Y-even.  Same argument for the ρ- and σ-sectors used in `five.py` (ρ, σ fix F and f; ρA = −A, σA = A).

**A.2 Euler normalization (q = 1) is lossless.**  Mod F: D(AD − N(F)) ≡ 0 and gcd(D,F) = 1, F irreducible ⟹ AD − N(F) = gF,
g ∈ S_d.  N' = N + (g/6)Z satisfies N'(F) = AD and B_D((g/6)Z) = (g/6)(D(d+6) − dD) = gD, so E(N') = E(N) − gDF + gDF = E(N).
(For q ≠ 1 the shift changes E by (1−q)gDF/q, so this step is special to q = 1 — correctly used by the lane, which only treats q = 1.)
Hence D ∈ (J : A)_d is necessary, and for such D the problem is: ∃ K Koszul (K(F)=0), F | B_D(N_0 + K).

**A.3 Five-variable reduction (dossier 1.4, `fringewidth`).**  F = f(x) + Y2⁶, A free of Y2.  μ6(Y2) fixes F and A; average
the field v = N/D (denominator becomes the norm D' = ∏ g·D, a polynomial in x and Y2⁶).  Write everything as polynomials in Y2⁶
over C(x) and reduce modulo F² using that D'(x, Y2⁶ = −f) ≠ 0 (coprimality): v ≡ (u(x), Y2·a(x)) + F·w mod F², and
K_A(F w) = F² div w, so K_A(v) ≡ K_A((u, Y2 a)) mod F² (denominators in C[x] only).  Then
K_A((u,Y2a)) = [A − u(f) + f(div u + a)] + Y2⁶[div u − 5a] has Y2⁶-degree ≤ 1 and is divisible by F² (degree 2 in Y2⁶) ⟹ both
brackets vanish ⟹ a = div u/5 and A − u(f) + (6/5) f div u = 0.  Conversely (u, Y2 div u/5) is a strict certificate; poles in
C[x] are coprime to F (F is monic in Y2).  CONFIRMED.

**The step the dossier could not reproduce (strict 5-variable system) is also correct**, for the homogeneous problem: u is
homogeneous of degree 1 in x, so g := div u/5 is homogeneous of degree 0 and E(g) = 0 (E = Σ x_i∂_i).  u' := u − gE has
u'(f) = u(f) − 6gf = A and div u' = div u − E(g) − 5g = 0.  So "6/5-equation solvable" ⟺ "N(f) = AD, D div N = N(D) solvable".
Moreover for a Y2-free denominator D the reduction above needs no norm, and the vanishing of the Y2⁶-bracket gives
a = B_D(α)/(5D²) with D·a polynomial, i.e. D | B_D(α); then u' = (α − (B_D(α)/5D)·x)/D keeps the denominator D.  Hence for
Y2-free D: 6-variable fixed-D consistency ⟺ 6/5 fixed-D consistency ⟺ strict 5-variable fixed-D consistency (same D).
The lane's `five.py` solves exactly the 6/5 system (checked column by column), and its 5→6 lift uses D² (correct; the docstring's
"(D, (N, Y2 div N/5))" is loose, the code uses D6 = D²).

**A.4 What "inconsistent mod p" proves.**  rank_p(M) ≤ rank_Q(M) and rank_p([M|b]) ≤ rank_Q([M|b]).  Inconsistent mod p
means rank_p([M|b]) = rank_p(M) + 1, which is compatible with Q-consistency when rank_Q(M) = rank_p(M) + 1 (bad prime for M).
A rigorous exclusion needs rank_Q(M) = rank_p(M), e.g. an exhibited Q-kernel of M of dimension (#cols − rank_p), or an integer
dual witness.  All 310 engine-2 fixed-D rows in `log.jsonl` are single-prime (2^31−1) statements; at a random large prime a
rank drop is astronomically unlikely, so they are strong evidence but not proofs.  (The thread's C1/C2/C4/C8 claims do give
integer witnesses or exact computations; see the ledger.)

**A.5 Mod-J necessary quadrics (used for C1 and my replays).**  For D = Σ g_i F_i and N = A g + K (K Koszul):
D div N ∈ J, and N(D) ≡ A Σ_{ij} g_i g_j F_ij mod J, so F | B_D(N) ⟹ A gᵀHess(F) g ∈ J_{2d}.  Rigor of a mod-p Gröbner
conclusion "V(I_p) = {0}": then I_p ⊇ (c)^k, and since (I_Q ∩ Z[c]) reduces to a space containing I_p degree-wise (saturated
lattice, equal ranks), (I_Q)_k = all degree-k forms, so V_C(I_Q) = {0}.  Needs only the good-prime check dim_p J_{2d} = dim_Q J_{2d}
(so that the F_p quadric space is the reduction of the Q quadric space).  Positive-dimensional mod-p results are only upper
bounds on the Q-variety's "size" in the other direction and are reported as such.
