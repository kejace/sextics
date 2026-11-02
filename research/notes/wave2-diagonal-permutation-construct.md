# Wave 2 — diagonal-permutation — CONSTRUCTION lens

STATUS: no witness. Mechanisms (i)–(v) analysed against the exact (F, q=1, A = X0^3Y0^3 − X1^3Y1^3); each fails for an identified structural reason (§2). New rigorous negative facts: every natural geometric degree ≤ 6 divisor tested (Shioda-type cones over the tau-fixed curves with any vertex, cones over the transposition-fixed sextic surface with equivariant or random vertex lines, the tau·s-cone, the sign-type quadric) supports NO primitive (3,1)-class (§3, X1/J5/J5b, 13 divisors, all kernels zero). Exact msolve tests for hidden reflections / §2.7 exchange partners on P^5 and on the (01)-pencil fibres were still running on the remote at hand-off (jobs `construct-dp-j1j2`, `construct-dp-fibre`; results will be appended to the REMOTE copy of wave2-construct.log.jsonl — retrieve with `python3 -m hodge.rjob pull --lane diagonal-permutation`); numerical screens found no candidate. A ready-to-run served-space job for D8 (w2c_d8_served.py) is prepared but not submitted (concurrency budget).

Log (one JSON line per experiment): /tmp/hodge/results/diagonal-permutation/wave2-construct.log.jsonl
Scripts: /tmp/hodge/results/diagonal-permutation/w2c_*.py, w2c_j1j2.sh, w2c_fibre.sh.

## 0. Summary (10 lines)
1. A spans the standard representation of the simultaneous S3; no quotient of X is uniruled (all ages ≥ 1, K ≡ 0), so support must come from an auxiliary variety or a correspondence, and the book's averaging-ideal criterion says a quotient-type correspondence serves A only through a subgroup fixing the (01)-invariant vector B = T0+T1−2T2.
2. (i) A §2.7-type exchange cover needs a point v ∉ X at which F = c y^6 + y G5 + G6 (y dual to v): equivalently rank Hess F(v) = 1 and the odd polars ∂_v^3F, ∂_v F vanish on the polar hyperplane ker ℓ_v. Pure-X or pure-Y partners are excluded by hand (the quadratic coupling matrices M(X)(v,v), N(Y)(v,v) are nonzero for all v ≠ 0); the mixed case is the exact job J2 (cubic polar condition, 6 charts). Lifted symmetries never create a reflection on a cyclic cover (fixed loci of lifts stay of codimension ≥ 2, and quotients by lifted groups containing the deck involution descend to X/s).
3. (ii) The unique Klein group containing a transposition and tau, G = <s,tau>, has no linear-section surface fibration with p_g = 0 fibres: p_g(S/G) = 1 for the five 1-parameter families of G-stable P^3's and 4 for the only 2-parameter family (exact character count, J3). Quotients by K ≤ S3 keep the 4-form (p_g = 1).
4. (iii) Natural cubic cover X/tau → P^2×P^2 with Galois closure Z and the commuting S3^sym: a subgroup K ≤ S3^mono × S3^sym serves A only if (V_std ⊗ V_std)^K ≠ 0, and for every such K the quotient Z/K is non-uniruled (K without monodromy reflections: Z → Z/K étale in codimension one, so Z/K uniruled ⇒ Q uniruled, false; K with one monodromy transposition: Z/K = X/(<tau> × K̄)). Extends book §4.7 to the symmetric case.
5. (iv) Deformation from the solved sub-families Fermat + (5/11)O2 / Fermat + (5/13)O3 (three reflections) is equivalent, degree by degree, to the E1 kernel bundle Ann(A)_d(λ) over the segment, whose generic rank is 0 for d ≤ 7 and 1 for d = 8 (the D8 line); so no deformation of a low-degree Fermat support reaches F, and the d ≥ 10 case is the search problem.
6. (v) No companion on the same variety.
7. Served-space facts (X1/J5/J5b): multiplication by [D] on R_6^- is injective for D = ΣX_i^6, ΣY_i^6, their sum, the cones over S_01, S_02 with the anti-invariant vertex line, the cones over Z_X with random vertex planes (both parities), the cone over S_01 with random vertex lines, the tau·s cones, the sign quadric u_X∧u_Y, X0^6+Y0^6 and the coupling part; since served(D) ⊆ Ann([D]) ∩ R_6^- these divisors support nothing. On the Fermat member ΣX_i^6 ∈ J and the join carries all bidegree-(3,3) classes (control: kernel 208).
8. Consequence: the (3,1)-class of A is transcendental for every "visible" divisor of degree ≤ 6; with the E1 bound deg D ≥ 8 (attack1), a support must be a divisor of degree ≥ 8 with non-rational singularities or a non-normal locus whose (2,0)-forms globalize (the join structure that does this at Fermat is destroyed by the couplings).
9. Correct test for ANY candidate divisor D (for the search lane): A is served by the S3-orbit product ∏_g D∘g iff
   A ∈ Q[S3]·served(D) + J_6^-, i.e. iff A lies in the span of the six S3-translates of served(D) = fixed_d.target_image(F,1,D)
   modulo J(F)_6^- (helper: w2c_served_orbit.py). Products of excluded denominators are NOT automatically excluded (served spaces
   add), and a (01)-stable divisor may serve the anti-invariant A without serving the invariant B (served spaces are (01)-stable
   subspaces, not necessarily (01)-invariant classes): the E1 facts A·D8 ∈ J, B·D8 ∉ J are consistent with V(D8) serving A.
10. No finite-field incidence seed can be proposed: there is no uniruled auxiliary variety on which the shared Lean modules (HoldoutCover, PolynomialCurves, …) could be instantiated; the first new object a construction would need is a (01)-stable divisor of degree ≥ 8 serving B, or a hidden fibre reflection on the (01)-pencil (job J2f).

## 1. Inherited rigorous facts (dossier + attack1; not re-derived)
- F = Fermat + (5/7) O1 + (5/11) O2 + (5/13) O3 with O1 = Σ_{i<j} P_ij(X) Y_i Y_j, O2 = Σ_k P_ij(X) Y_k^2, O3 = Σ_k X_k^2 P_ij(Y),
  P_ij(Z) = Z_i Z_j (Z_i^2 + Z_j^2), (i,j,k) a permutation. Equivalently F_{4,2} = Y^T M(X) Y with M_ii = (5/11) P_jk(X),
  M_ij = (5/14) P_ij(X), and F_{2,4} = (5/13) Σ_k X_k^2 P_ij(Y). [verified against instance.rec]
- Symmetry group S3 × <tau>; projective diagonal group {1,tau}; no reflection among monomial symmetries; all ages ≥ 1.
- Residue of A transforms in V_std (det = +1); 69 copies of std in (R_F)_6^- (dim 138 of 208).
- E1 kernel Ann(A)_d^±: 0 for d ≤ 7, one line D8 (Y-even) at d = 8, d = 9 unknown (search lane), 294/262 at d = 10.
- Prescribed-support theorem (docs/01-certificates.md, q = 1): a divisor X ∩ V(D) supports ω_A ⇔ the fixed-D system over D is
  consistent; "A·D ∈ J(F)" is necessary for support; served(D1 D2) = served(D1) + served(D2) (Gysin images add over components).
- benediamond (#177): A ⇔ B = T0 + T1 − 2T2 with (01)-equivariant data; X/(01) is a CY3-fibration over the t-line.

## 2. Mechanism-by-mechanism analysis

### 2.1 (i) Covers creating a reflection / exchange (book §2.7, §2.8)
Shape required by §2.7: coordinates (y, z_1..z_5) with F = c·y^6 + y·G_5(z) + G_6(z) and G_6 containing a binary sextic in two
of the z's that couple to y linearly. Let v be the point dual to y and H = {z} the complementary hyperplane. Then
∂_v^k F|_H = 0 for k = 2,3,4,5; ∂_v^5 F ∝ ℓ_v = Σ ∂_iF(v) z_i forces H = ker ℓ_v (so F(v) = ℓ_v(v)/6 ≠ 0), and
∂_v^4 F = z^T Hess F(v) z vanishing on ker ℓ_v together with Hess F(v)·v = 5∇F(v) forces Hess F(v) = λ∇F(v)∇F(v)^T (rank 1).
The decisive conditions are the odd ones: ℓ_v | ∂_v^3 F (cubic) and ℓ_v | ∂_v^2 F, ∂_v F. The same cubic condition with F(v) ≠ 0
characterises hidden reflections of F (F even w.r.t. (v, ker ℓ_v)), and with F(v) = 0 it is necessary for X ∩ T_vX to be a cone
with vertex v (the Fermat support type: cones over p_g = 10 sextic surfaces).
- Rank-deficient points do exist and are harmless: coordinate points (rank 1, off X) and (1, ±i, 0,0,0,0)-type points (rank 2, on X,
  where X0^6 + X1^6 and P_01(X) both vanish); both fail ℓ_v | ∂_v^3 F (exact check, log entry J1-numeric-interpretation).
- Pure-Y partner: ∂_v^2 F = 2 M(X)(v,v) + (Y-terms) with M(X)(v,v) = Σ_{i<j} P_ij(X)[(5/11) v_k^2 + (5/7) v_i v_j] ≠ 0 for v ≠ 0
  (the P_ij are independent and the three equations (5/11)v_k^2 + (5/7)v_iv_j = 0 force v = 0); a linear form through v dividing
  ∂_v^2F would have to divide this X-quartic after setting Y = 0, impossible (the ℓ^4-term argument). Pure-X partner: N(Y)(v,v) =
  (5/13)Σ v_k^2 P_ij(Y) ≠ 0 likewise. So only mixed v remain: exact job J2 (msolve, 6 charts, cofactor quadric, u·F(v) = 1) —
  RESULT: PENDING. Numerical screen (30 Powell starts on the sampled restricted-cubic norm): no candidate (log J2-numeric-reflection).
- Lifted symmetries on cyclic covers Z = {w^n = b(z)}: a lift of a transposition s has fixed locus of codimension ≥ 2 (b∘s = b:
  involution fixing π^{-1}(S); b∘s = −b: order 4 with s̃^2 = deck involution ι, fixed locus S̃ with eigenvalues (1,1,−1,±i), age 3/4
  but <s̃> ∋ ι so Z/<s̃> = X/s); lifts of tau, tau·s, 3-cycles have ages ≥ 1 (eigenvalue bookkeeping on T_xB ⊕ <w>); tau·c lifts act
  freely. Every automorphism of Z normalising the deck group descends to X, so a group L of lifts not containing ι acts with
  Z → Z/L étale in codimension one, and Z/L uniruled ⇒ Z uniruled ⇒ X uniruled, false. Ambient automorphisms of P(1^6,k) preserving
  the ideal (F, w^n − b) are lifts (w ↦ λw + h(z) forces h ≡ 0 mod F). Hence a cover can only help through a non-lifted BIRATIONAL
  automorphism, i.e. a fibration of X whose fibres acquire extra automorphisms after base change — which J2/J2f test in the two
  natural cases (P^5 and the (01)-pencil fibres).
- §2.8 (order-6 diagonal complementary quotient): the diagonal group is {1,tau}: excluded.

### 2.2 (ii) Fano / surface-fibration companions of quotients by subgroups of S3 × <tau>
- Uniruledness: ages ≥ 1 for all 11 nontrivial elements at all fixed strata (dossier table), K_X = 0 ⇒ every X/K has canonical
  singularities and K ≡ 0: none is uniruled (Reid–Tai, Kollár–Larsen). S3 preserves the 4-form (det = +1) ⇒ p_g(X/K) = 1 for K ≤ S3;
  X/tau, X/<s,tau>, X/<tau s>, X/<tau c> have p_g = 0 and Kodaira dimension 0.
- Hodge numbers of the quotients (from the dossier's equivariant Griffiths decomposition R_6^- = 39 triv + 31 sign + 69·2 std,
  R_6^+ = 46 triv + 28 sign + 72·2 std; residues of Y-odd numerators are tau-invariant, of Y-even ones tau-anti-invariant):
  | K | h^{4,0}(X/K) | h^{3,1}(X/K) = dim R_6^K | uniruled? |
  |---|---|---|---|
  | 1 | 1 | 426 | no (CY) |
  | <s> | 1 | 108 + 118 = 226 | no |
  | A3 | 1 | 70 + 74 = 144 | no |
  | S3 | 1 | 39 + 46 = 85 | no |
  | <tau> | 0 | 208 | no (κ = 0) |
  | <tau s> | 0 | 108 + 100 = 208 | no (κ = 0) |
  | <s, tau> | 0 | 108 | no (κ = 0) |
  | <tau c> (order 6) | 0 | 70 | no (free action, K trivial) |
  | <A3, tau> | 0 | 70 | no |
  | S3 × <tau> | 0 | 39 | no |
  The class of A lives in the 69·2-dimensional std-part of the 208; its (01)-fixed partner B in the 108 of X/<s,tau>.
- Surface fibrations with p_g = q = 0 fibres (book §3.2): impossible for K ≤ S3 (a nonzero 4-form restricts to a nonzero 2-form on a
  general fibre). For G = <s,tau>: the G-stable P^3's are sums of the eigenspaces E_{++} = <p,a>, E_{+-} = <q,b>, E_{-+} = <r>,
  E_{--} = <t> (p = X0+X1, r = X0−X1, q = Y0+Y1, t = Y0−Y1, a = X2, b = Y2); the 2-forms of the sextic surface S = X ∩ P^3 are
  Res(G·Ω_3/F|_{P^3}) with G a quadric, and the G-invariant ones are those with χ_G = det (J3, exact count):
  A = E_{++}⊕E_{+-} (Fix s): 6; B = E_{++}⊕E_{-+}⊕E_{--}: 1; C,D,E,F (1-parameter families): 1; G = <l'>⊕<l>⊕E_{-+}⊕E_{--}
  (the only 2-parameter family): 4. So no linear-section fibration over a surface has p_g = 0 fibres; non-linear fibres have
  K_S = O(d_1 + d_2), larger p_g. benediamond's fibration (fix t and λ = b/a) is the s-only version (p_g = 3), consistent.
- The (01)-pencil {r = μt} ⊃ P^3_0 = Fix(s): every fibre X_μ is a sextic threefold (K = O(1)) on which s is a reflection with fixed
  surface S (X_μ is smooth along S for every μ: T_vX ⊇ P^3_0 would need dF_v = 0). The quotient fibres W_t = X_μ/s are the CY3
  sextics in P(1,1,1,1,2) of #177. A Fano quotient of the fibres needs a SECOND reflection on X_μ (K = O(1) minus two hyperplane
  ramification divisors = O(−1), exactly the §2.7 fibre structure): exact job J2f (hidden reflections of F_μ at μ = 3/7; a random μ
  decides the geometric generic fibre) — RESULT: PENDING. If J2f is empty, no pencil of hyperplane sections of X has fibres with two
  reflections after base change (other pencils would need two hidden reflections).

### 2.3 (iii) Galois closures of the natural finite maps
Natural cubic cover Q = X/tau → B = P^2_X × P^2_Y, fibre {a + bt + ct^2 + dt^3 = 0}, Galois closure Z with group S3^mono. Since F
is S3^sym-invariant, S3^sym acts on B, Q, Z commuting with S3^mono, and the class of A lies in the V_std^mono ⊗ V_std^sym isotypic
part of H^4(Z) (book §4.7: the primitive part of H^4(Q) is the V_std^mono part). Averaging ideal: K ≤ S3 × S3 can serve A only if
(V_std ⊗ V_std)^K ≠ 0, which forces K ⊉ S3^mono × 1 and K ⊉ A3^mono × 1. Claim: Z/K is never uniruled for such K. Proof: the only
elements with a fixed divisor on Z are the monodromy transpositions (s_i, 1) (ramification over the discriminant); (s, g) and (1, g)
with g ≠ 1 have fixed loci over Fix(g) ⊂ B, of codimension 2 ({X0=X1}×{Y0=Y1}) or 0-dimensional (3-cycles: 9 points).
(a) K contains no (s_i,1): Z → Z/K is étale in codimension one; a covering family of rational curves on Z/K lifts to rational curves
on Z (general members avoid the codimension-2 branch locus), and Z → Q finite surjective would make Q uniruled — false (κ(Q) = 0).
(b) K ∩ (S3^mono × 1) = <(s,1)>: it is normal in K, Z/<(s,1)> = Q, and K/<(s,1)> embeds in 1 × S3^sym, so Z/K = X/(<tau> × K̄),
non-uniruled by 2.2. Hence the extra S3^sym does not enlarge the averaging ideal of any uniruled quotient of the Galois closure.
Other natural finite maps: F ≠ G(Φ) for a nonlinear polynomial map (no block/pullback structure, scalar Hessian symmetrizer), and
the non-Galois cubic cover X/C2 → X/S3 has Galois closure X itself (back to 2.2). (The canonical-class bookkeeping of the triple-cover
model needs the exceptional divisors over the tau-fixed curves; the étale-in-codimension-one argument above avoids it.)

### 2.4 (iv) q-raising / Gauss–Manin / deformation
§4.5 only raises q. On the segment F_λ = (Fermat + (5/11)O2) + λ((5/7)O1 + (5/13)O3) the start is served by three Y-reflections
(existence only, unbounded degree). A family of certificates (D_λ, N_λ, Q_λ) of fixed degree d would give D_λ ∈ Ann(A)_d(F_λ) for all
λ, and the generic rank of this kernel bundle equals its rank at λ = 1 (the instance): 0 for d ≤ 7, 1 for d = 8 (dossier/attack1,
rigorous mod p ⇒ over Q). So every low-degree Fermat-side certificate is obstructed at some finite order, and the d ≥ 10 case is
exactly the bilinear search problem (D in the 294-dimensional kernel). Not pursued further here.

### 2.5 (v) Companions: none (provenance has no sameVarietyInstances).

## 3. Experiments of this lens
| id | what | result | rigor |
|---|---|---|---|
| X1 | ker [D] on R_6^- (Y-odd targets), D ∈ {ΣX^6, ΣY^6, ΣX^6+ΣY^6, cone_01 = F|_{r=t=0}, cone_02, X0^6+Y0^6, F − Fermat} | all 0 (ring dim 208); Fermat control for ΣX^6: 208 | p = 2^31−1; zero kernel transfers to Q |
| J5 | ker [D_P] on R_6 (both parities), D_P = Σ(L_P Z)_i^6 cones over Z_X with random vertex planes (2 samples) | 0, 0 | mod p ⇒ Q |
| J5b | ker on R_6 (both parities) for the sign quadric u_X∧u_Y (sweep of the S3-orbit planes), cones over S_01 with random vertex lines (2), tau·s-cone over Fix(tau s)_+ with vertex E_- | all 0 (the second tau·s-cone was still running at hand-off; see log) | mod p ⇒ Q |
| J3 | p_g of <s,tau>-quotients of stable linear sections | min 1 (1-param), 4 (2-param) | exact |
| J1-num | 200 BFGS starts on the rank-≤2 Hessian objective | finds only the trivial rank-deficient points (coordinate points; (1,±i,0,..)) | numerical |
| J2-num | 30 Powell starts on the restricted-cubic objective | no candidate | numerical |
| J2 | ℓ_v ∣ ∂_v^3 F, F(v) ≠ 0 (hidden reflections / §2.7 partners), 6 charts, msolve F_p, 4 threads | RUNNING at hand-off (remote `construct-dp-j1j2`, chart 0 started 01:23 UTC) | mod p evidence |
| J2-onX | same with F(v) = 0 (cone tangent sections) | PENDING | mod p evidence |
| J1 | rank Hess F ≤ 2 locus (3×3 minors), 6 charts | PENDING (last in the job; expected to return the trivial points) | mod p |
| J2f/J1r1/J1f | hidden reflections of the fibre F_μ (μ = 3/7, 5 charts), then the rank-1 Hessian locus on P^5 (2×2 minors, 6 charts: order-3/6 quasi-reflections and §2.7 partners), then rank ≤ 2 on the fibre | RUNNING at hand-off (remote `construct-dp-fibre`) | mod p evidence |

Interpretation of X1/J5/J5b. For Y = X ∩ V(D) at q = 1 every served target A' satisfies A'·D ∈ J(F) (reduce the certificate
identity mod F), so served(D)/J_6^- ⊆ ker([D]: R_6^- → R_12^-); zero kernel ⇒ Y supports no nonzero tau-invariant primitive
(3,1)-class, whatever its singularities. Geometric reading: on the Fermat member the join {ΣX^6 = 0} ∩ X_0 = J(Z_X, Z_Y) (a cone over
Z_X × Z_Y) carries all 100 bidegree-(3,3) classes including A (Shioda); on X the divisor X ∩ {ΣX^6 = 0} still has the curve Z_Y of
transversal cone-over-Z_X singularities (exceptional surface Z_X × Z_Y, 100 candidate (2,0)-forms) but its tau-invariant Gysin image is
zero (rigorous). Interpretation (not proved here): either the candidate forms do not globalise to the resolution once the join
structure is broken (the d_2 differential H^1(Z_Y, R^1π_*O) → H^3(O_Y) of the Leray sequence need not vanish; along Z_X the
transversal type becomes the quadric cone Y^T M(x) Y = 0, which is rational), or they globalise and the Gysin map kills them. Likewise the cone over S_01 with vertex the anti-invariant line is the sweep of the s-equivariant plane
sextics C_p = X ∩ (p + E_−), p ∈ S_01, a family whose total space has H^{2,0} ⊇ H^{2,0}(S_01) (10-dimensional), and its cylinder map to
H^{3,1}(X)^tau vanishes. These are the first supported-sector computations on this variety and they show the obvious geometry is
"invisible" to the primitive (3,1)-cohomology.

## 4. What a construction would have to supply (precise)
- A divisor Y ⊂ X of degree ≥ 8 (E1) with non-rational singularities or a non-normal locus whose (2,0)-forms globalise, such that
  A ∈ Q[S3]·served(Y) (criterion in §0.9). All degree-≤ 6 natural candidates are excluded (§3).
- Or a second reflection on the (01)-pencil fibres after base change (J2f), which would make Z/<s̃,σ> fibrewise Fano and give
  H^4(X)^{(01)} ⊂ N^1 via the descent of §4.6 (with the exchange conjugating the deck action into a lift of a physical symmetry).
- Or a non-quotient correspondence from a uniruled 4-fold; none of the natural incidence varieties (lines, conics: 1-parameter
  families on a CY4; planes: parametrised by Σ ⊂ Gr(3,6) with p_g(Σ) > 0 — the cylinder map of the family {Π_σ ∩ X}) produced a
  nonzero served class in the equivariant cases tested.
- Lean: no instantiation of HoldoutCover/PolynomialCurves/EllipticCompanions/ResidueTransport/HigherResidues/TransitiveCubic is
  possible without a uniruled auxiliary variety; nothing to seed.
