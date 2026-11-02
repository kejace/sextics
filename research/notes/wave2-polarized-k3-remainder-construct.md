# Wave 2 — lane `polarized-k3-remainder`, lens CONSTRUCTION

STATUS: **no witness** (no certificate for the instance). Result of this lens: (A) a proof-level obstruction that no
quotient/cover mechanism available on X or its natural covers reaches the K3 sector of the target (so a construction for
the K3 part would be a new case of the generalized Hodge conjecture, as the dossier already argued — here made precise
with the *full* automorphism group, which is larger than previously recorded); (B) the complete explicit recipe,
with all identities verified and a finite-field incidence seed found, for the one part of the target that the book
supports only existentially — the six mixed (sixth-root-sector) monomials — i.e. the "book's recipe" item of the task.
(B) does **not** solve the instance (the K3 part A_K3 remains), but it is the first concrete step toward a
`HasCertificate.add` decomposition and is reusable for every polarized-split target with A(X,X)=0.

Log (one JSON line per experiment): `/tmp/hodge/results/polarized-k3-remainder/wave2-construct.log.jsonl`
Scripts: `/tmp/hodge/results/polarized-k3-remainder/{w23_seed_search.py,w23_seed_fastA.py,w23_seed_randomA.py,w23_job_B.sh,wave2_construct_util.py}`
Seed files: `construct-w23seedA-p5-kx4.jsonl` (16650 F_5 solutions), `construct-w23seedB-p5-kx4.jsonl` (profile B, all singular).
Remote job (queued in background, remote saturated at load ~300 by other agents when written):
`construct-polarized-k3-remainder-w23seed` = `w23_job_B.sh` (extra F_5 seeds for kx=1,2,3; random F_11 / F_13 seeds, outputs
`construct-w23seedA-p{5,11,13}-kx*.jsonl` in the lane dir; pull with `python3 -m hodge.rjob pull --lane polarized-k3-remainder`).

Conventions: U = X+Y, V = X−Y, F = ½(f(U)+f(V)), f = z0^6 + (5/7)(z0^5 z1 + 2 z1^5 z2 + 3 z2^5 z0 + z0^2 z1^3 z2) + z1^6 + z2^6,
Â(U,V) = 2^{-6} ∏(U_i² − V_i²); V-degree b ↦ residue character η^{b+3} under h(U,V) = (U, ηV).
A = A_K3 + A_mix with A_K3 = [(X0+Y0)²(X1+Y1)²(X2+Y2)² − (X0−Y0)²(X1−Y1)²(X2−Y2)²]/64 (b = 0,6) and A_mix the six
(4,2),(2,4) monomials (b = 2,4). Both are nonzero mod J(F)_6 (rank 18 → 19, p = 2^31−1) [log: target-split-AK3-Amix-jacobian-membership].

---------------------------------------------------------------------------------------------------------------------
## PART A — K3 sector: proven facts (all exact, local sympy / hand, logged)

A1. **Aut(X) ⊇ D_6 of order 12 — correction to the dossier.** h(U,V) = (U, ηV) preserves F exactly (f homogeneous of
degree 6); h^3 is the coordinate swap X ↔ Y; h^3τ is (X,Y) ↦ (−Y, X). The *diagonal-in-(X,Y)* group is ⟨τ⟩ (order 2,
as checks.json says) but the linear automorphism group is D_6 = ⟨h, τ⟩ (τ h τ = h^{-1}); with Aut(C_f) = 1 (probe) and the
uniqueness of the U/V splitting, Aut(X) = D_6 exactly. [log: symmetry-D6-planes]

A2. **X contains six planes** Π_λ = {V = λU}, λ^6 = −1 (in (X,Y): Y = μX, μ ∈ {±i, ±i(2±√3)}): f(U) + f(λU) = (1+λ^6)f(U) = 0.
They are the fixed planes of the six involutions h^kτ, k odd (two planes each). [log: symmetry-D6-planes, D6-fixed-loci-and-ages]

A3. **Character/Hodge table of H^4(X)_prim under D_6** (from Hilb(R_f) = (1,3,6,10,15,18,19,18,15,10,6,3,1), R_F = R_f ⊗ R_f,
character η^{b+3}; τ-invariant residues ↔ swap-antisymmetric polynomials) [log: D6-sector-hodge-table]:
| sector (D_6-rep) | h-characters | (3,1),(2,2),(1,3) τ-inv | τ-anti |
|---|---|---|---|
| K3 Λ²V_2 (χ: h ↦ −1, τ ↦ +1) | η^3 | 19, 172, 19 | 19, 191, 19 |
| sixth-root pair (ρ_1, rotation by 60°) | η^{±1} | 90, 261, 90 | 90, 261, 90 |
| third-root pair (ρ_2) | η^{±2} | 54, 333, 54 | 54, 333, 54 |
| joining Λ²H¹(C)(−1) (trivial) | 1 | 45, 100, 45 | 55, 100, 55 |
Matches the book's table (19,172,19), (90,261,90), (54,333,54), (45,100,45) exactly.

A4. **Fixed loci and ages (Kollár–Larsen / Reid–Tai) of every element of D_6 on X** [log: D6-fixed-loci-and-ages]:
h^k (k≠0): fixed curves C_U = {V=0, f(U)=0}, C_V = {U=0, f(V)=0}; tangent eigenvalues (1, η^k, η^k, η^k) at C_U ⇒ age k/2,
and (6−k)/2 at C_V. h^kτ, k odd: two planes, eigenvalues (1,1,−1,−1), age 1. h^kτ, k even (τ and conjugates): two curves,
eigenvalues (1,−1,−1,−1), age 3/2. **Junior elements (age < 1): only h and h^5 (age 1/2).** Consequently X/⟨h⟩ is uniruled
(indeed X/⟨h⟩ ≅_bir P²×P²: the map X ⇢ P²_U × P²_V has degree 6 with deck group ⟨h⟩), and H^4(X)^h = joining sector + Tate is
divisor-supported — a second proof of the joining support.

A5. **Averaging-ideal exhaustion (Theorem).** Let K ⊆ Aut(X) = D_6. X/K is uniruled iff K contains a junior element iff
K ∈ {⟨h⟩, D_6}. The K3 character χ (h ↦ −1, τ ↦ +1) has χ^K ≠ 0 iff K ⊆ ker χ = ⟨h², τ⟩ ≅ S_3, which contains no junior element.
Hence **no quotient of X by a subgroup of its full automorphism group (diagonal or not) serves the K3 sector**; the same
holds for the sixth-root (ρ_1^K ≠ 0 only for K ⊆ {1, h^kτ}, ages 1 or 3/2) and third-root sectors. Only the joining sector is
served (by ⟨h⟩), consistent with the book needing auxiliary varieties (W_f, W_{23,f}) for the companion sectors.
Intermediate quotients of the natural Galois closure X → P²×P² (Galois group ⟨h⟩ ≅ μ_6; the closure is X itself):
X/⟨h^3⟩ (X↔Y quotient), X/⟨h^2⟩ ≅_bir {w² = −f(U)f(V)} ≅_bir (S_2×S_2)/(ι,ι) (the K3 sector lives entirely here, as the
τ-invariant part of V_2⊗V_2), X/⟨h⟩ ≅_bir P²×P²; all but the last have Kodaira dimension 0 (ages ≥ 1). On S_2 × S_2 itself the
group ⟨(ι,1),(1,ι),flip⟩ ≅ D_4 has junior elements (ι,1),(1,ι) only, which act by −1 on V_2⊗V_2; g = (ι,1)∘flip (order 4)
acts by +1 on Λ²V_2 but has age 3/2. So (ii)/(iii) of the task are provably exhausted for the K3 sector.

A6. **Covers (task item (i)).** Every cover of X branched along an effective divisor has K > 0 (X is CY and simply connected),
so the age criterion is unavailable on covers; the §2.7 trick needs a Fano *quotient* of a cover by a group with a reflection.
On the natural covers S_6 × S_6 (general type, K = O(2,2)) the quasi-reflections (δ^a, 1) (fixed divisor C × S_6) act on
H_3 ⊗ H_3 by (−1)^a; those acting trivially ((δ², 1), (δ⁴,1)) give S_2 × S_6-type quotients which are not uniruled
(p_g ≥ 3 computed from ω_2 ⊗ H^0(K_{S_6})^{δ³-anti} = 3-dimensional). The double cover Z = {s² = f(U) = −f(V)} of X branched along
the joining divisor is singular along C_U ∪ C_V and of general type; its automorphism group μ_6 × μ_2 ⋊ Z/2 yields the same
quotients as above. No coordinate pairing or linear form produces a new exchange: Aut(C_f) = 1 kills every non-scalar M, and
(U,V) ↦ (aU, bV) preserves F iff a^6 = b^6, i.e. the elements of ⟨h⟩.

A7. **Shape of any K3-sector support (restating the obstruction).** By the dossier's equivalence theorem (End_Hdg(T) = Q,
MT(T) = SO(T), Λ²T irreducible, Hom_HS(Λ²T, T) = Hom_HS(Λ²T, Q(−2)) = 0), any legal D serving A_K3 cuts X in a divisor whose
resolution has H² ⊇ so(T)(1) (weight 2, h^{2,0} = 20−ρ ≥ 17, MT = SO(T)). Exceptional divisors of hypersurface sections are:
ruled over curves (ADE transversal types — algebraic H²), elliptic/genus-g fibrations over curves (H¹⊗H¹, Jacobian-type MT),
or the base surface of a cone (needs a surface S' with H²(S') ⊇ Λ²T(1) — a "Kuga–Satake surface"). The planes of A2 only
give (2,2) classes; hyperplane sections through a plane acquire 25 nodes (rational singularities) — nothing in (3,1).
Fibre products S ×_{P¹} S of pencils give H²(Z̃)^- ⊇ V_2 ⊗ (hyperplane) only. **No low-degree mechanism exists; the K3 part is
equivalent to GHC(coniveau 1) for Λ²T of a very general degree-2 K3 (ρ ≤ 3, End = Q) — new mathematics, not a construction.**

A8. (iv) q-raising / Gauss–Manin: §4.5 differentiates an existing generic certificate; nothing exists to differentiate in the
K3 sector. Along f_t = f + t·z0²z1²z2² the class η∧ω = (∇_t ω)∧ω is the derivative of the identically-zero section ω_t∧ω_t
in Λ²V_2 — no information (in Sym²V_2 it is ½∇_t(ω²), an unsupportable (4,0)-derivative). Special fibres with algebraic
Kuga–Satake (order-three family) do not propagate to the generic fibre. (v) There is no solved companion on this variety
(catalogue scan: the only other polarized-split instance is `contact46006`, a different F).

---------------------------------------------------------------------------------------------------------------------
## PART B — Mixed (sixth-root) sector: explicit construction data (all identities machine-verified)

B1. **Chart reduction for every V-degree b (formal identity, sympy, b = 0..6)** [log: chart-reduction-general-b].
Chart U = t(x,a,1), V = tλ(r,z,1): F = t^6 G/2, G = p + λ^6 q, p = f(x,a,1), q = f(r,z,1). For a target monomial of V-degree b,
A = t^6 λ^b s_b(x,a,r,z), the field w = (λ^b W_x, λ^b W_a, λ^b W_r, λ^b W_z, w_λ),
w_λ = (2λ^b s_b − λ^b W(p) − λ^{b+6} W(q))/(6λ^5 q), satisfies w(G) = 2λ^b s_b and has divergence (with chart term 2w_λ/λ)
divisible by G **iff**
    div W + ((b−3)/6)·W(p)/p − ((b+3)/6)·W(q)/q = ((b−3)/3)·s_b/p        (W ∈ C(x,a,r,z)^4).
Equivalently the cover field Z = u^{b−3} v^{−(b+3)} W on u^6 = p, v^6 = q has div Z = ((b−3)/3) s_b u^{b−9} v^{−(b+3)}.
Checks: b = 1 reproduces the accepted crossed-ternary equation (weights 4/3, 2/3 = u^4 v^2 in cube roots); b = 0,6 lands on
the double cover u^3v^3 = √(pq) (K3 sector = (S_2×S_2)/(ι,ι), consistent with A5); b = 3 has **zero** right-hand side.

B2. **Explicit strict certificate for the whole joining sector (b = 3)** [log: joining-sector-explicit-certificate, engine-exact].
For every (U,V)-bidegree-(3,3) target A (any parity): D = f(V) = f(X−Y) (degree 6), N = (A/3)·E_V, i.e.
N_{X_i} = A(X_i − Y_i)/6, N_{Y_i} = −A(X_i − Y_i)/6, Q = 0. Verified exactly with `core.verify_witness` on
A = U0²U1V0V1V2 − V0²V1U0U1U2 and on U0³V0³. (Explicit form of the book's "joining divisor" support; irrelevant for this target,
which has no b = 3 part, but it calibrates B1.)

B3. **Sources for the six mixed monomials.** s_2 = −(a²r² + x²z² + x²a²)/64 (b = 2), s_4 = (r²z² + a²r² + x²z²)/64 (b = 4).
By B1, b = 2 needs cover sources α·β with α = m_4(x,a)/u^7, β = n_2(r,z)/v^5 for (m_4,n_2) ∈ {(a²,r²),(x²,z²),(x²a²,1)};
b = 4 needs α = m_2/u^5, β = n_4/v^7 (mirror). Characters: u^{-7} ~ η^5 ↔ H_5^{1,1} (quartic × u^4), v^{-5} ~ η ↔ H_1^{2,0}
(quadric × u^0): exactly the book's H_5⊗H_1 ⊕ H_1⊗H_5 pair.

B4. **Second elliptic companion and its scalar conversion (formal identity, sympy)** [log: scalar-conversion-identity-(3,2)-companion].
W_{23,f}: H = 7f(x,a,y) + 7S² + 7T³ = 7x^6 + 5x^5a + 5x²a³y + 15xy^5 + 7a^6 + 10a^5y + 7y^6 + 7S² + 7T³ ⊂ P(1,1,1,3,2)
(smooth: cone gradient (∇f, 14S, 21T²) has no nonzero zero; −K = O(2)); affine companion p(x,a) + S² + T³ = 0 after radial descent
y ↦ 1. Cover S = u³α, T = u²β, u^6 = p, E_{23}: α² + β³ = −1 with invariant derivation δβ = α, δα = −(3/2)β²
(the differential dβ/α has character η; α ↦ −α, β ↦ ζ_3^{-1}β under ε). If a companion field V = (V_x,V_a,V_S,V_T) satisfies
V(H) = B and div V = 0, then with
    f = 6αV_T/u³ − 4βV_S/u^4,    C_x = −6V_x/u,    C_a = −6V_a/u :
    δf + D_x C_x + D_a C_a = B/u^7 − 3αB_S/u^4 − 2βB_T/u^5
(D_x = ∂_x + (p_x/6u^5)∂_u on the cover). Control: the same derivation reproduces the accepted crossed-ternary (2,2) identity
f = 3(b_eV_S − a_eV_T)/u², C = −3V/u, source B/u^4 − (a_eB_S + b_eB_T)/u³. Consequences: companion source B = m_4(x,a) gives
**exactly** δf + div C = m_4/u^7 (the b = 2 U-block source of B3); on the V-companion the source B' = n_2(r,z)·T gives
δg + div G = −(n_2/v^5)·β, i.e. the V-block source times θ = β (character η^4) — the θ-twisted primitive needed for the
contraction (crossed case: θ = b_e/a_e²). Weighted-degree reasons force θ ∝ β: polynomial B' with v^{-5} must be n_2·T.
Mirror for b = 4 (B = m_2·T on the U side, B' = n_4 on the V side).

B5. **Finite-field incidence seed on W_{23,f} (the task's "book's recipe" item)** [log: w23-companion-seed-profileA-p5].
Dimension count: conics (weighted degree e = 2) on W_{23} have −K·C = 4, so conics through a point form a **2-dimensional**
family; the crossed-ternary normalization (x = 4t, three endpoint parameters) is therefore not square here: profile B
(a, y linear; 13 unknowns) has det J = 0 for all 2320 F_5 solutions [construct-w23seedB-p5-kx4.jsonl]. Correct square system
("profile A"): x = 4t, a = c0+c1t+c2t², y = k0+c3t+c4t², S = k1 + c5t+…+c9t^5 + k3t^6, T = k2 + c10t+c11t²+c12t³ + k4t^4;
parameters (k0,k1,k2,k3,k4) = constant terms of y,S,T and the top coefficients S6,T4; 13 unknowns c0..c12; 13 equations
[t^0..t^12] of H(φ(t)). Over F_5: 16650 solutions, **4426 with invertible 13×13 incidence Jacobian and rank-12 freeness
matrices** (12×16 = H^0(O(1)^3 ⊕ O(5) ⊕ O(3)) → H^0(O(11)), and 12×14 without the x-columns, for the polynomial normal).
Chosen seed (independently re-verified with sympy over GF(5); H(φ) ≡ 0, det J = 3, ranks 12/12; dφ/d(params,t) has rank 4 =
dim of the affine cone, i.e. the family is dominant):
    x = 4t,  a = 1 + t + 2t²,  y = 2 + 3t²,  S = 3 + 4t + 4t^5 + t^6,  T = 1 + 2t² + t^4     (mod 5),
    parameters (k0,k1,k2,S6,T4) = (2,3,1,1,1); endpoint gradient (H_a,H_y,H_S,H_T) = (2,4,2,1) ≠ 0.
Newton/Hensel lifting of this seed in the *integral* system (true 5/7 coefficients cleared) converges modulo 5^6
[log: w23-seed-hensel-lift-5adic], as `PolynomialCurves.ArithmeticLifting.square_system_generic_lift` requires.
Weighted lines (e = 1) do not work: 14480 F_5 solutions, none with invertible 7×7 Jacobian [log: w23-companion-weighted-lines-e1-p5].
Random F_11/F_13 seeds (non-Fermat special fibre): remote job `construct-polarized-k3-remainder-w23seed` (pending capacity).

B6. **What a Lean proof of HasCertificate(F, 1, A_mix) would need beyond the shared library (precise gap list).**
(1) PolynomialCurves/HoldoutCover: the seed of B5 fits `square_system_generic_lift` (special-fibre solution, invertible
Jacobian, integral ring → generic field) and `anchored_polynomial_family_injective` (independent constants k0,k1,k2 (+k3,k4),
distinguished coordinate x = 4t); the freeness matrix supplies the polynomial normal (`full_polynomial_combination_one`).
Then `companion_residue_solution`-style integration gives residue primitives for every polynomial source on W_{23,f}
(this part is target-independent and mirrors MovingConic.lean).
(2) **New**: ScalarCalculus for weights (3,2) — the identity of B4 (verified here, not formalized).
(3) **New**: EllipticCompanions.Contraction is stated only for the Fermat cubic with its μ_3 structure (θ = b_e/a_e²,
kernels K0,K1 at r³ = 1728, infinity recurrence 5+2n). The sixth-root sector needs the same architecture on
E_{23}: α² + β³ = −1 (j = 0, isomorphic over C to the Fermat cubic but with the order-6 automorphism), θ = β (character η^4),
its θ-kernels / residue pairing R(fθ) = c·(m_4/u^7) − div ξ with its own infinity correction, and a translated-diagonal
contraction for diag μ_6 instead of μ_3 (book §3A.5). Pole separation/translation/trace descent are generic.
(4) Physical transport: the chart identity of B1 with b = 2 and b = 4 (verified here) replaces the crossed b = 1 identity;
legality via the same Eisenstein-type specialization.
Estimated size: comparable to PR #149 (≈ 13 Lean files); months, not a session. It would not solve this instance (A_K3 remains)
but is required by any eventual solve and is the first explicit sixth-root-sector construction for polarized-split sextics.

---------------------------------------------------------------------------------------------------------------------
## Ideas considered and rejected (with the reason)
* Using the six planes / the quadrics U_iV_j − U_jV_i (vanishing on all planes): hyperplane sections through a plane have
  25 nodes — rational singularities, no (3,1) contribution; products of the quadrics are excluded by the sum principle.
* Quotient by h² (order 3, acts trivially on the K3 sector): age exactly 1 at both fixed curves — Calabi–Yau-type quotient
  {w² = −f(U)f(V)}, not uniruled; its own natural quotients only reach the joining sector.
* Non-diagonal involutions h^kτ (k odd, fixed planes, age 1): act by −1 on the K3 sector; and their quotients are CY anyway.
* Fibre-product divisors S ×_{P¹} S of pencils on S_2: Gysin image in V_2 ⊗ NS only.
* Deformation to the order-three (Bolognesi–Laterveer) family: algebraic cycles do not propagate from a special fibre.

## Experiments (all in the log; wall times)
symmetry-D6-planes (11 s) · D6-sector-hodge-table / D6-fixed-loci-and-ages (<1 s) · joining-sector-explicit-certificate
(engine exact, d=6) · chart-reduction-general-b (sympy, b=0..6) · scalar-conversion-identity-(3,2)-companion (+ (2,2) control)
· target-split-AK3-Amix-jacobian-membership (p=2^31−1) · w23 seed enumerations: profile B (7 s, all singular),
profile A (16650 solutions, 4426 good), e=1 lines (32 s, none good), 5-adic Hensel lift (mod 5^6) · remote
construct-polarized-k3-remainder-w23seed (queued; extra F_5 kx and random F_11/F_13 seeds).
