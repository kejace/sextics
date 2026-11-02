# Wave 2 — diagonal-permutation — CONSTRUCTION lens

STATUS: in progress — no witness; mechanisms (i)-(v) analysed, two new rigorous negative facts (X1, see §3), remote msolve jobs J1/J2 queued (remote saturated at submission time), local J3/J5 running.

Log (one JSON line per experiment): /tmp/hodge/results/diagonal-permutation/wave2-construct.log.jsonl
Scripts: /tmp/hodge/results/diagonal-permutation/w2c_*.py, w2c_j1j2.sh (remote job `construct-dp-j1j2`).

## 0. One-paragraph summary
The target A = T0 - T1 (T_i = X_i^3 Y_i^3) spans the standard representation of the simultaneous S3; every quotient of X is
non-uniruled (all ages >= 1, K == 0), so the only known routes are auxiliary varieties (covers with hidden symmetry, Galois
closures, correspondences) or deformation.  This lens checks each route against the exact (F,q,A).  New facts: (a) the natural
"Shioda join" divisors — the cones over the tau-fixed plane sextics Z_X = X ∩ {Y=0}, Z_Y = X ∩ {X=0} with vertex the complementary
plane, which carry the bidegree-(3,3) classes (hence A) on the Fermat member — serve NO primitive tau-invariant (3,1)-class on this
X (multiplication by [ΣX_i^6], [ΣY_i^6] is injective on R_6^-; rigorous over Q by the modular-rank argument); (b) the cone over the
transposition-fixed sextic surface S_01 = X ∩ {X0=X1,Y0=Y1} with vertex the anti-invariant line (the s-equivariant cylinder
correspondence of plane sextics, p_g >= 10 on its resolution) also serves nothing; (c) the natural cubic cover's Galois closure, even
with the extra S3^sym action, has no uniruled intermediate quotient reaching V_std^mono ⊗ V_std^sym (étale-in-codimension-one
argument, §2.3); (d) no <s,tau>-stable linear-section fibration of X/<s,tau> has p_g = 0 fibres (minimum 1 for 1-parameter
families, 4 for the only 2-parameter family).  Pending: J1 (no point of P^5 with rank Hess F <= 2 — kills every §2.7-type exchange
cover and every cone hyperplane section), J2 (no hidden reflection).

## 1. Inherited rigorous facts (dossier + attack1; not re-derived)
- F = Fermat + (5/7) O1 + (5/11) O2 + (5/13) O3 with O1 = Σ_{i<j} P_ij(X) Y_i Y_j, O2 = Σ_k P_ij(X) Y_k^2, O3 = Σ_k X_k^2 P_ij(Y),
  P_ij(Z) = Z_i Z_j (Z_i^2 + Z_j^2)  ((i,j,k) a permutation).  [verified against instance.rec]
- Symmetry group S3 × <tau>; projective diagonal group {1,tau}; no reflection among monomial symmetries; all ages >= 1.
- Residue of A transforms in V_std (det = +1); 69 copies of std in (R_F)_6^- (dim 138 of 208).
- E1 kernel Ann(A)_d^±: 0 for d <= 7, one line D8 (Y-even) at d = 8, d = 9 unknown (search lane), 294/262 at d = 10.
- Prescribed-support theorem (docs/01, q = 1): a divisor X ∩ V(D) supports ω_A  <=>  the fixed-D system over D is consistent; hence
  "A·D ∈ J(F)" is a necessary condition for V(D) ∩ X to support A, and served(D1 D2) = served(D1) + served(D2).
- benediamond (#177): equivalent to serving B = T0 + T1 − 2T2 with (01)-equivariant data; X/(01) is a CY3-fibration over the
  t-line (fibres sextics in P(1,1,1,1,2)); its first surface fibration has p_g = 3.

## 2. Mechanism-by-mechanism analysis (what was tried, where it fails)

### 2.1 (i) Covers creating a reflection / exchange (book §2.7, §2.8)
The §2.7 construction needs, after a linear change of coordinates, F = c·y^6 + y·G_5(z) + G_6(z) (the partner coordinate y occurs
only in degrees 1 and 6; the base curve T^6 = h(Y_a,Y_b) then makes the fibre equation symmetric in y ↔ T).  Write v for the point
dual to y and H for the complementary hyperplane.  Then ∂_v^k F|_H = 0 for k = 2..5; ∂_v^5 F ∝ ℓ_v := Σ_i ∂_iF(v) z_i forces
H = ker ℓ_v, and ∂_v^4 F = z^T Hess F(v) z vanishing on ker ℓ_v together with Hess F(v)·v = 5∇F(v) forces
   **Hess F(v) = λ ∇F(v) ∇F(v)^T, i.e. rank Hess F(v) = 1, with v ∉ X.**
(Check on the holdout sextic: v = e_{Y1} gives Hess = diag(0,0,0,0,30,0).)  Similarly a hyperplane section X ∩ T_vX that is a cone
with vertex v ∈ X (the support type used by the Fermat solutions: cones over sextic surfaces, p_g = 10) forces rank Hess F(v) <= 2.
Job J1 (remote msolve, 6 affine charts, p = 1073741827) computes the ideal of 3×3 minors of Hess F.  Result: PENDING.
Hand check for v in the Y-span or X-span: ∂_v^2 F contains the X-quartic 2M(X)(v,v) with M(X)(v,v) = Σ_{i<j} P_ij(X)[(5/11)v_k^2 + (5/7)v_iv_j],
which is nonzero for every v ≠ 0 (the three P_ij are independent and the system (5/11)v_k^2 + (5/7)v_iv_j = 0 ∀(i,j,k) forces v = 0),
and no linear form through v divides it (the ℓ^4 term argument), so pure-Y and pure-X partner coordinates are excluded by hand; J1
settles the mixed case.  Consequence if J1 is empty: no coordinate pairing and no non-coordinate linear form admits a §2.7 cover.
The §2.8 complementary quotient needs an order-6 diagonal symmetry: the diagonal group is {1,tau}: excluded.
Lifting the symmetries of X to a cyclic cover never helps: a lift of a transposition s to a double cover w^2 = b(z) with b∘s = ±b has
fixed locus of codimension >= 2 (b∘s = b: involution, fixed locus π^{-1}(S); b∘s = −b: order 4 with s̃^2 = deck involution, so
Z/<s̃> = X/s), hence Z → Z/L is étale in codimension one for every lifted group L not containing the deck involution, and
Z/L uniruled ⇒ Z uniruled ⇒ X uniruled, impossible.  A cover can only help through an automorphism that is NOT a lift (the
§2.7 exchange), whose existence J1 tests.

### 2.2 (ii) Fano / surface-fibration companions of quotients by subgroups of S3 × <tau>
- Uniruledness: every nontrivial element has age >= 1 at every fixed point (dossier §1.6; the transposition has age exactly 1 along
  S_01, the 3-cycles age 2, tau·c acts freely).  With K_X = 0 the quotients X/K have canonical singularities and K ≡ 0: none is
  uniruled (Reid–Tai + Kollár–Larsen).  Elements of S3 preserve the 4-form, so X/K has p_g = 1 for K <= S3; only K ∋ tau-type
  elements kill it (X/tau, X/<s,tau>, X/<tau s>, X/<tau c>: p_g = 0, Kodaira dimension 0).
- Surface fibrations with p_g = q = 0 fibres (book §3.2): impossible for K <= S3 (a 4-form restricts to a nonzero 2-form on a
  general fibre).  For G = <s,tau> (the only Klein group containing a transposition and tau), the G-stable P^3's are sums of the
  eigenspaces E_{++} = <p,a>, E_{+-} = <q,b>, E_{-+} = <r>, E_{--} = <t> (p = X0+X1, r = X0−X1, q = Y0+Y1, t = Y0−Y1, a = X2, b = Y2),
  and p_g of the quotient of the sextic surface X ∩ P^3 is the number of quadrics G on P^3 with χ_G = det (J3, exact count):
  types A (fixed P^3 of s): 6; B: 1; C,D,E,F (1-parameter families): 1; G (the only 2-parameter family <l'> + <l> + E_{-+} + E_{--}): 4.
  So no linear-section fibration of X/<s,tau> over a surface has p_g = 0 fibres; non-linear fibres have larger p_g (K_S = O(d1+d2)).
  benediamond's fibration (t, λ = b/a fixed) is the type-C/D family with s only (p_g = 3), consistent.
- Non-diagonal elements add nothing: the ages are computed for all 12 elements (dossier table) and the fixed-point data are what the
  quotient singularities see.

### 2.3 (iii) Galois closures of the natural finite maps
Natural cubic cover Q = X/tau → B = P^2_X × P^2_Y, fibre {a + bt + ct^2 + dt^3 = 0}, Galois closure Z with group S3^mono.  Because
F is S3^sym-invariant, S3^sym acts on B, Q, Z commuting with S3^mono, and the class of A lies in the V_std^mono ⊗ V_std^sym
isotypic part of H^4(Z) (book §4.7: the primitive part of H^4(Q) is the V_std^mono part).  Averaging ideal: K <= S3 × S3 serves A
only if (V_std ⊗ V_std)^K ≠ 0, which forces K ⊉ S3^mono × 1.  Claim: Z/K is never uniruled for such K.  Proof: the only elements of
S3^mono × S3^sym with a fixed divisor on Z are the monodromy transpositions (s_i, 1) (fixed locus = a ramification divisor over the
discriminant); (s, g) and (1, g) with g ≠ 1 have fixed loci over Fix(g) ⊂ B, of codimension 2 (g a transposition: {X0=X1}×{Y0=Y1};
g a 3-cycle: 9 points).  (a) If K contains no (s_i,1): Z → Z/K is étale in codimension one, so a covering family of rational curves on
Z/K lifts to rational curves on Z, and Z → Q finite surjective would make Q = X/tau uniruled — false (κ = 0).  (b) If K ∩ (S3^mono × 1)
= <(s,1)>: it is normal in K and Z/K = Q/K̄ = X/(<tau> × K̄) with K̄ <= S3^sym, not uniruled by 2.2.  (c) K ⊇ A3^mono × 1 or S3^mono × 1
gives (V_std^mono)^K = 0.  Hence the extra symmetry does not enlarge the averaging ideal of any uniruled quotient of the Galois closure.
Other natural finite maps: F is not of the form G(Φ) for a nonlinear polynomial map (no block/pullback structure; Hessian symmetrizer
scalar), and the cubic cover X/C2 → X/S3 has Galois closure X itself (back to 2.2).

### 2.4 (iv) q-raising / Gauss–Manin / deformation
§4.5 only raises q.  The sub-families Fermat + (5/11)O2 (separately even in Y0,Y1,Y2) and Fermat + (5/13)O3 (separately X-even) serve
every odd target by three reflections (existence only, no explicit D; the book's reflection construction gives unbounded degree).  No
theorem transports support along the segment to F; the first-order deformation system (dossier E5) needs an explicit certificate at
the starting point, which does not exist in the literature (on Fermat + (5/11)O2 the primitive E1 kernel at d = 8 is 2-dimensional
(attack1 §2.6), so no low-degree witness there either).  Not pursued here.

### 2.5 (v) Companions on the same variety: none (provenance has no sameVarietyInstances).

## 3. New experiments (this lens)
| id | what | result | rigor |
|---|---|---|---|
| X1 | kernel of multiplication by [D] on R_6^- for D ∈ {ΣX^6, ΣY^6, ΣX^6+ΣY^6, cone_01 = F|_{r=t=0}, cone_02, X0^6+Y0^6, coupling part F−Fermat} | all kernels 0 (ring dim 208); Fermat control: kernel 208 | mod p = 2^31−1; zero kernel transfers to Q |
| J5 | same for cones over Z_X with a random (non-tau-stable) vertex plane, both parities | running | |
| J3 | p_g of <s,tau>-quotients of stable linear sections | min 1 (1-param), 4 (2-param) | exact character count |
| J1 | rank Hess F <= 2 locus (msolve, 6 charts) | PENDING (remote queue) | mod p evidence |
| J2 | hidden reflections (cubic polar condition) | PENDING | mod p evidence |

Interpretation of X1.  For a divisor Y = X ∩ V(D) at q = 1, every served target A' satisfies A'·D ∈ J(F) (reduce the certificate
identity mod F), so served(D)/J_6^- ⊆ ker([D] : R_6^- → R_12^-).  Zero kernel ⇒ Y supports no nonzero tau-invariant primitive
(3,1)-class.  Geometric reading: on the Fermat member the join {ΣX^6 = 0} ∩ X_0 = J(Z_X, Z_Y) carries all 100 bidegree-(3,3) classes
(Shioda), including A; along the deformation to F the divisor X ∩ {ΣX^6 = 0} keeps a curve Z_Y of transversal cone-over-Z_X
singularities (exceptional surface Z_X × Z_Y) but its tau-invariant Gysin image collapses to zero.  Hence the obvious "deform the
Fermat support" idea fails at the level of the support itself, not only at the level of the numerators.

## 4. Ideas not yet executed
- Degree-10 geometric members of Ann(A)_10 (294-dim): impose non-normality along S_01 ((r,t)^2) or a curve of non-rational
  transversal type inside the kernel — search-lane territory (attack1 item 4).
- The only remaining cover-type mechanism is a hidden (non-lifted) automorphism of a cover; J1/J2 test the two cheapest shapes.
