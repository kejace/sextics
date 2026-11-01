# Attack 1: `cyclic-triple-residual` (q = 1, A = X0^3 Y0^3) — resumed run

**No witness found.** The only "witness-candidate" file from the crashed first run
(`/tmp/hodge/results/cyclic-triple-residual/witness-candidate-1793572767.json`) is the *illegal control*
D = 6F (D must be coprime to F); it was produced on purpose to validate the solver and is not a certificate.

Work directory: `/tmp/hodge/results/cyclic-triple-residual/` (scripts, `log.jsonl` with one JSON line per tested
denominator, kernel elements `kernel_d*.json`). All linear algebra is exact over Q (python-flint `fmpz_mat`/`fmpq_mat`)
unless marked "mod p" (p = 2147483647, second prime 2147483629). The engine's sign bug in the F·N(D) term was fixed
before any result quoted here (all entries carry `"engine": 2`).

Headline results (details in sections 2-4):

1. **E1 is a sharp low-degree filter.** With R = Q[Z]/J(F) the Jacobian ring, the annihilator of [A] in R is
   zero in degrees d <= 7 and exactly one-dimensional in degree 8 (rigorous over Q: modular kernel dimensions are
   upper bounds for the kernel over Q, and dim J_d is the same over Q and F_p). Hence **every denominator of degree <= 7
   lies in J(F)** (in particular there is no denominator of degree <= 4 at all, since J_d = 0 there), and every
   denominator of degree 8 is congruent mod J(F)_8 to a multiple of one explicit polynomial `D8` (25 terms, sector
   (0,4) of the order-18 group). The degree-8 kernel is **not** forced by dimensions (sector dimensions 79 -> 97),
   so `D8` is a sporadic annihilator tied to the coupling structure; at degree 9 the kernel mod J has dimension 120 of
   which 116 are dimension-forced.
2. `D8` itself is excluded exactly over Q (sector-reduced system 2357 x 1188, rank 1155, augmented 1156), and random
   points of the coset D8 + J(F)_8 are excluded mod p (see table). The degree-9 sporadic kernel elements are excluded
   as listed below.
3. **Fermat calibration identifies the mechanism this class uses at t = 0**: on the Fermat sextic the target
   X0^3 Y0^3 is served by the *join divisor* D0 = X0^6 + X1^6 + X2^6 (equivalently Y0^6+Y1^6+Y2^6 = F - D0) with the
   strict block-Euler certificate N = A·(X0, X1, X2, 0, 0, 0)/6, Q = 0; it is **not** served by the (X0,Y0)-cone
   hyperplanes X0^6 + Y0^6 (the classical Shioda cones), nor by the reflection hyperplane X0^k, nor by any other
   Fermat-type sum tested. Geometrically: X ∩ V(D0) is the join of the plane sextic curves C_X, C_Y, whose resolution is
   a P^1-bundle over C_X x C_Y with H^{2,0} = H^{1,0}(C_X) ⊗ H^{1,0}(C_Y) (dimension 100), and the class of
   X0^3Y0^3 has Fermat character (4,1,1 | 4,1,1) of "join type" (both triples sum to 0 mod 6).
   The natural deformation of D0 to the coupled F, D = E_X(F) = Σ X_i ∂F/∂X_i (the ramification divisor of the
   degree-6 map X -> P^2 x P^2), fails (reduced Koszul system rank 370 of 388 columns, augmented 371, exact), and so do
   all twelve 0/±1 weighted Euler derivatives Σ m_i Z_i F_i (previous run) — the coupling terms of X-degree 2 reduce the
   multiplicity of T along the fixed curves from 6 (cone over a sextic, irrational) to 2 (transversal A_k, rational).

## 1. The obstruction and the geometric mechanism a certificate needs

**Data.** F = Z^2 + t Z C + H with t = 5/7, Z = X2^3, C = X0^3 + X0Y1Y2 + X1Y0Y1 (cubic, X2-free), H = the X2-free
sextic (Fermat part of the other five variables plus three couplings). A = X0^3 Y0^3, q = 1. The projective diagonal
symmetry group has order 18 (generators mu3: X2 -> z X2, g: (X0,X1,X2,Y0,Y1,Y2) -> (w, w^5, w^5, w^2, w^2, 1)·Z with
w = zeta6, containing tau = -1 on Y); the residue of A Ω/F^2 has mu3-character zeta3 (sector V_1, dim (107, 341, 69)),
and lies in the G-sector (0,3)+... of (3,1)-dimension 41 (dossier §1.6). The companion solve (A' = X0 X2^2 Y0^3) lives
in the mu3-invariant sector V_0, which is served by the Fano double cover W = X/mu3 ⊂ P(1,1,3,1,1,1); nothing in the
Galois closure of X -> P^4 (group Z/3 wr Z/2) gives a uniruled quotient seeing V_1 (dossier §1.7), and the stabiliser
of the residue character has no junior elements. So no quotient/Fano mechanism is available.

**What a certificate is, geometrically.** A fixed-D certificate exhibits the class of A as supported on T = X ∩ V(D):
it lies in the image of ⊕_i H^2(T̃_i)(-1) -> H^4(X) over the resolved components. Since X is a smooth hypersurface of
dimension 4, Pic(X) = Z·O(1) (Lefschetz), so *every* component T_i is itself X ∩ V(D_i). For any such Cartier divisor
H^2(O_{T_i}) = 0 (H^2(O_X) = 0 and H^3(O_X(-d)) = 0 by Kodaira), hence if T_i has rational singularities then
h^{2,0}(T̃_i) = 0 and T_i supports no (3,1)-class. **Every useful component must have non-rational singularities**
(non-normal counts, via H^2 of the conductor sheaf). Concretely the supply of (2,0)-forms comes from:

* a curve Γ ⊂ T along which T has multiplicity m >= 3 with transversal type a cone over a plane curve B of degree m
  (simple elliptic for m = 3): the exceptional surface E -> Γ is a B-fibration and contributes
  H^1(Γ, R^1) ∩ F^2 ≈ H^0(Γ, K_Γ ⊗ Hodge bundle), i.e. at least g(Γ) - 1 forms (g(Γ)·g(B) when the fibration is a
  product, as for the Fermat join where Γ = C_Y, B = C_X, m = 6, giving 100);
* an isolated point of T whose tangent cone is a cone over a surface of degree >= 4 (quartic: 1 form; sextic: 15);
* non-normal loci (double surfaces with irregular normalisation).

The dossier's heuristic bound (if the (67,238,67) Q-Hodge substructure is simple) requires Σ h^{2,0}(T̃_i) >= 67;
one simple-elliptic curve of genus 10 gives only ~9-10, so a successful D would need either several such curves /
high multiplicity, or the Hodge structure of the sector to split further (not known either way).

**The E1 filter (new, rigorous).** The identity splits (equivalently) into E1: Σ N_i F_i = A D and
E2: D div N - N(D) = F Q (shift N by R·Z/6 using Euler). E1 says [A]·[D] = 0 in the Jacobian ring R = S/J(F).
Computed by G-sector (mod p; kernel dimensions mod p bound those over Q from above and dim J_d is characteristic-free):

| d | dim S_d | dim J_d | dim Ann([A])_d mod J | dimension-forced part | sporadic |
|---|---|---|---|---|---|
| 0-4 | - | 0 | 0 | 0 | 0 |
| 5 | 252 | 6 | 0 | 0 | 0 |
| 6 | 462 | 36 | 0 | 0 | 0 |
| 7 | 792 | 126 | 0 | 0 | 0 |
| 8 | 1287 | 336 | **1** (sector (0,4)) | 0 | **1 = D8** |
| 9 | 2002 | 756 | 120 | 116 | 4 (sectors (0,0),(2,5),(4,1),(4,4)) |

(Degrees 10-12: see section 2.) The "dimension-forced" column is Σ_sectors max(0, dim R_d^χ - dim R_{d+6}^{χ+χ_A}),
computed from the Fermat monomial basis (the G-graded character of R is that of the Fermat ring). Consequences:
no denominator of degree <= 4 exists at all; every denominator of degree 5, 6, 7 lies in J(F) (so T is the tangency
locus of a polynomial vector field G of degree d - 5 with X: D = G(F)); every denominator of degree 8 is, mod J(F)_8,
a multiple of the explicit D8 (file `kernel_d8_0_4_0.json`). The degree-8 sporadic annihilator persists for every
coupling value t tested (t = 1, 2, 1/3, -5/7, 3/11), so it is a property of the coupling pattern, not of t = 5/7.
Reformulated by Gorenstein duality (socle degree 24): [A]·R_{18-d} = R_{24-d} for d <= 7 — in particular
[A]·H^{2,2}_prim -> H^{1,3} is onto — and [A]·R_10 has codimension exactly 1 in R_16.

## 2. Experiments (this resumed run; all entries also in `log.jsonl`, 83 lines incl. the first run)

Tools: `sect.py` (character-sector-reduced fixed-D solver, any q, exact fmpz/fmpq ranks or mod p; validated against fqlx's
exact exclusion of X2^6+Y0^6+Y1^6+Y2^6: 440 cols, rank 436, augmented 437, and against the illegal control D = F which
lifts and verifies), `kernel8.py`/`kernel9.py` (sector E1-kernel extraction with CRT + rational reconstruction),
`kdims.py`, `planes.py`, `fermat_cal.py`, `d8limit.py`, `sing.sing`. Machine load average was ~310 on 40 cores during
the run (other agents), so timings are 4-5x slower than nominal.

### 2.1 E1 kernel by degree (all G-sectors, mod p = 2147483647; rigorous upper bounds over Q)
| d | dim S_d | dim J_d | Ann([A])_d mod J | dimension-forced | sporadic | time |
|---|---|---|---|---|---|---|
| 0..7 | - | 0,0,0,0,0,6,36,126 | 0 | 0 | 0 | < 2 min total (first run, full space) |
| 8 | 1287 | 336 | 1 (sector (0,4)) | 0 | 1 | 44 s |
| 9 | 2002 | 756 | 120 | 116 | 4: (0,0) +1, (2,5) +1, (4,1) +1, (4,4) +1 | 74 s |
| 10 | 3003 | 1497 | 556 | 555 | 1 | 362 s |
| 11, 12 | - | - | (running at shutdown; forced parts 1020, 1325) | | | |

### 2.2 The degree-8 element D8 (sector (0,4), Y-even, 25 terms; `kernel_d8_0_4_0.json`)
```
D8 = X0^6 X1^2/2 - 54 X0^6 Y0 Y1/5 + 3 X0^6 Y1^2/2 - 25 X0^4 X1^2 Y0 Y2/588 + 6 X0^4 Y1^3 Y2 - 882 X0^4 Y2^4/25
   + 252 X0^3 X1^5/5 + 5 X0^3 X1^3 Y0 Y1/28 + 5 X0^3 X1^2 X2^3/28 - 6 X0^3 X1 Y0^2 Y1^2 + 5 X0^3 X1 Y0 Y1^3/7
   - 25 X0^3 X1 Y0 Y2^3/588 + 5 X0^3 X2^3 Y1^2/14 + 5 X0^2 X1^2 Y0 Y1 Y2^2/84 - 25 X0^2 X1^2 Y1^2 Y2^2/196
   - 25 X0^2 X1 X2^3 Y2^2/588 - 2 X0^2 Y0^2 Y1^2 Y2^2 + 5 X0^2 Y0 Y1^3 Y2^2/14 - 5 X0 X1^3 Y1^3 Y2/84
   - 25 X0 X1^2 X2^3 Y0 Y2/588 - 63 X0 X1 Y0 Y1 Y2^4/5 - 63 X0 X2^3 Y2^4/5 + X1^8/2 + 3 X1^6 Y1^2/2 + X1 X2^3 Y1^4
```
| test | system (sector-reduced) | rank / augmented | result |
|---|---|---|---|
| D8, exact over Q | 2357 x 1188 (N 939, Q 249) | 1155 / 1156, nullity 33 | inconsistent (62 s) |
| D8 + random J_8^(0,4) element, 6 trials (integer and rational coefficients), mod p, 2 primes | same | 1155 / 1156 each | inconsistent |
| membership of D8 in J + V for V = Z_i^k S_{8-k} (k=1..4, all i), all 20 cone subrings k[3 coords], 15 four-coordinate and 6 five-coordinate subrings, A S_2, C S_5, Z' S_5, H S_2, H' S_2, C^2 S_2, CZ' S_2, Z'^2 S_2 | exact (mod p ranks on 118-dim sector) | - | **no hit**: D8 is none of these mod J |
| sporadic degree-8 annihilator for F_t, t = 1, 2, 1/3, -5/7, 3/11 | mod p | - | present for every t in sector (0,4) (t = 2 has a second one in sector (4,3)); none in degrees 6, 7 at t = 1 |

### 2.3 Denominators from J(F) (first run, reduced Koszul system, exact; `j6.py`/`dstar.py`)
All thirteen weighted Euler derivatives Σ m_i Z_i F_i with m ∈ {0,±1}^6 (E_X F, E_Y F, single Z_i F_i, sign patterns): rank 370 of 388
columns, augmented 371 (X2 F_X2: 364/365; Y2 F_Y2: 369/370) — inconsistent. Control 6F: consistent, lifts (illegal). Dossier: all six polar
quintics, C, Z', H, H', X0Y0, A, CZ', HC, HZ', all 63 Fermat-type cones, random semi-invariant D of degree 1..5 in all 68 sectors — excluded exactly.

### 2.4 Fermat calibration (F_0 = Σ Z_i^6, A = X0^3Y0^3, sector solver with |G| = 7776, exact over Q)
| D | result |
|---|---|
| X0^6+X1^6+X2^6, X0^6+X1^6+Y1^6, X0^6+X1^6+Y2^6, X0^6+X2^6+Y1^6, X0^6+X2^6+Y2^6, X0^6+Y1^6+Y2^6 (the six join-type splits) | **consistent**, strict witness N = A·(coordinate vector of the triple)/6, Q = 0, verified |
| X0^6+X1^6+Y0^6, X0^6+X2^6+Y0^6, X0^6+Y0^6+Y1^6, X0^6+Y0^6+Y2^6 | inconsistent (7/8) |
| X0^6+Y0^6 (six cone hyperplanes), X0^6-Y0^6, X0^2+Y0^2, X1^6+X2^6, X1^6+X2^6+Y1^6+Y2^6, X0^6+X1^6+X2^6+Y0^6, X0^6+Y0^6+Y1^6, generic Σ c_i Z_i^6 (6 terms) | inconsistent |
| X0, X0^2, X0^3 (reflection hyperplane), X0Y0, X0^3Y0^3 | inconsistent |
| cubic/quartic cones over planes: X0^3+X1^3+X2^3, Y0^3+Y1^3+Y2^3, X1^3+Y0^3+Y2^3, X0^3+X2^3+Y1^3, X0^4+X2^4+Y1^4, X0^4+X1^4+X2^4+Y0^4 (K3 cone, vertex line), X0^2X2^2Y1^2, X0^3X2^3+Y1^6, X0^2X2^2+X2^2Y1^2+X0^2Y1^2 | inconsistent |

### 2.5 Deformation of the Fermat certificate (menu item 4)
D = E_X(F) (deformation of 6·D0) on the coupled F: exact rank 370/388, augmented 371 — inconsistent (see §4.1 for the geometric reason).

### 2.6 Plane / line multiplicity inside the E1 coset (degree 8, exact linear algebra over Q, `planes.py`)
For each coordinate plane Π (3 coordinates = 0) the maximal k with some D = D8 + Σλ_k J_k vanishing to order k along Π:
| Π | k_max | family | tests |
|---|---|---|---|
| {X0,X1,X2}, {X0,X1,Y0}, {X0,X1,Y1}, {X0,X1,Y2}, {X0,X2,Y0} | none (not even k = 3) | - | - |
| {X0,X2,Y1} (Γ = Fermat curve X1^6+Y0^6+Y2^6) | 3 | affine dim 18 | particular and random member: 1155/1156, inconsistent (exact) |
| {X0,X2,Y2} (Γ = Fermat curve X1^6+Y0^6+Y1^6) | 3 | affine dim 16 | particular member: 1155/1156, inconsistent (exact) |
| remaining 13 planes, 15 coordinate lines | running at shutdown | | |

### 2.7 Running at shutdown (results not available; scripts and partial outputs in the work directory)
`kernel9.py` (degree-9 sporadic kernel elements in sectors (2,5),(4,0),(4,1),(4,2),(4,4),(0,0), 4 primes, exact tests),
`kdims.py` (kernel dims d = 11, 12), `d8limit.py` (D8(t) at t = 5/7, 1/10, 1/100 for the t -> 0 limit),
`sing.sing` (Singular, char 32003: dimension/degree of Sing(X ∩ V(D8)), of Sing(X ∩ V(E_X F)), Fermat join control),
`batchB2.py` prepared (polar pencil, polar products, structural products of degree 8-12, mod p) — not launched.

## 3. Witness
None. (`witness-candidate-1793572767.json` is the labelled illegal control D = 6F.)

## 4. Mathematical ideas developed, partial constructions, what would finish them

**4.1 The Fermat mechanism for this class and why its deformation dies.** On F_0 = Σ Z_i^6 the class of X0^3 Y0^3
(character (4,1,1 | 4,1,1), both triples ≡ 0 mod 6) is a *join-type* class: it is served by D0 = X0^6+X1^6+X2^6 with
the strict certificate N = A·(X0,X1,X2,0,0,0)/6, Q = 0 (verified exactly, section 2). The identity behind it is
block-Euler: D0 = E_X(F_0)/6, so Σ N_i ∂_i F = A D0 (E1), and div(AG) = A, AG(D0) = A D0 (E2 = 0) because D0 is
X-homogeneous of degree 6 and A of X-degree 3. Geometrically X_0 ∩ V(D0) is the join of the plane sextics C_X, C_Y,
singular along both with transversal type the cone over a sextic curve; its resolution is a P^1-bundle over C_X × C_Y
with H^{2,0} = H^{1,0}(C_X) ⊗ H^{1,0}(C_Y) (dimension 100), and (4,1,1) is a character of H^{1,0} of each curve.
For the coupled F the same formula gives D = E_X(F) = Σ X_i ∂F/∂X_i: T = X ∩ V(D) is the **ramification divisor of
the degree-6 map X -> P^2_X × P^2_Y**, (X,Y) -> ([X],[Y]) (fibres = roots of the binary sextic F(sX, tY)). It still
contains the fixed curves Γ_X = X ∩ P^2_X = {X0^6 + X1^6 + X2^6 + tX0^3X2^3 = 0} and Γ_Y = {Y0^6+Y1^6+Y2^6 = 0}, but
only with multiplicity 2 (the couplings X0^2Y0Y1Y2^2 and X0X1Y1^3Y2 have X-degree 2; the Y-degree-2 couplings do
the same on the other side): along Γ_Y the quadratic part is t·X0·Y1Y2·(X0Y0Y2 + X1Y1^2), a rank-2 form, so the
transversal singularities are A_k (rational) and the exceptional surfaces carry no 2-forms. This is the concrete
reason the deformation of the Fermat certificate fails (exact rank 370/388, augmented 371), together with all
weighted Euler fields Σ m_i Z_i ∂_i F (= ramification divisors of the maps to weighted products): the E2 residual is
(A/6)·Σ_k k(6-k) F^{(k)} with F^{(k)} the X-degree-k part of F, i.e. (4tA/3)·(the five mixed couplings), and it is not
in F·S + B(Koszul corrections).

**4.2 Multiplicity along a plane as a linear condition inside the E1 coset.** For a plane Π ⊂ P^5 with
Γ = X ∩ Π a plane sextic (genus 10 when smooth) and a form D vanishing to order m along Π, T = X ∩ V(D) is locally
Γ × (cone over the plane curve B_P = {D_m(P;·) = 0} ⊂ P(N_{Π,P})) at smooth points of Γ (X is transverse to Π, so one
can solve one Π-coordinate). For m = 3 the transversal singularity is simple elliptic and the exceptional surface
E -> Γ is an elliptic fibration whose Leray piece H^1(Γ, R^1) ∩ F^2 has dimension >= g(Γ) - 1 = 9 (exactly
g(Γ)·g(B) = 10 if the fibration is a product, i.e. if the cubic leading form does not depend on P ∈ Γ). Since
"vanishing to order m along Π" is linear in D, it can be intersected exactly with the E1-admissible space
⟨D8⟩ + J(F)_8 (sector (0,4), 40 generators). Results (section 2.6): among the 20 coordinate planes only
Π = {X0 = X2 = Y1 = 0} admits order 3 (an 18-parameter affine family), with Γ = {X1^6 + Y0^6 + Y2^6 = 0} a Fermat
sextic; no coordinate plane admits order >= 4, and the tested members of the family are inconsistent (exact).
Non-coordinate planes (9 parameters) and the choice inside the 18-dimensional family remain open; the family is
small enough for a *parametric* treatment (the system is linear in (N, Q) with coefficients linear in the 18
parameters and right-hand side quadratic), which is the natural next computation.

**4.3 Why E1 is so restrictive here and what it says about T.** [A]·R_{18-d} = R_{24-d} for d <= 7 means every
(1,3)-class is A times a (2,2)-class (surjective Yukawa-type coupling) — the deformed Jacobian ring has no
"accidental" relations in low degree except the single sporadic one in degree 8 (and three more in degree 9).
Denominators in J(F) are tangency divisors of vector fields; denominators in the degree-8 coset D8 + J(F)_8 are the
only other option below degree 9. D8 is not congruent mod J to anything divisible by a coordinate power, to a cone
over any coordinate plane or line, or to a multiple of A, C, Z', H, H', C^2, CZ', Z'^2 (exact membership tests), so
it is not an obvious product of the structural pieces of F. Its t -> 0 limit in the family F_t (section 2.7) shows
which of the ~700 Fermat annihilators (X0^2, Y0^2)_8 survives the coupling to first order.

**4.4 What a complete argument would need.** Either (a) a D of degree >= 8 (necessarily in the E1-admissible
cosets computed here) whose divisor has non-rational singularities along curves or at many points, with
Σ h^{2,0}(T̃_i) large enough to contain the 41-dimensional (or 67-dimensional) sector — the plane/line families of
4.2 are the systematic way to generate such D compatibly with E1; or (b) a non-constructive argument as in the
companion solve (a family of Fano-fibred or rationally connected quotients whose transport hits the zeta3 sector),
for which the only structure specific to V_1 found so far is the isotrivial j = 0 elliptic family π^{-1}(L) over the
lines L ⊂ W (dossier §5.E5); or (c) a proof of non-existence of low-degree certificates: the E1 table plus a
parametric exclusion of the degree-8 coset would make "no certificate of degree <= 8" rigorous.

## 5. Ranked next steps (with estimated compute on an unloaded 40-core machine)

1. **Parametric exclusion/solve of the degree-8 coset D8 + J(F)_8^(0,4)** (39 parameters): treat λ as unknowns in the bilinear
   system (matrix linear, rhs quadratic in λ); a modular Gröbner/msolve run or a rank-drop scan along structured
   subfamilies (the 18- and 16-dimensional order-3 plane families of §2.6 first). Makes "no certificate of degree <= 8"
   rigorous or finds one. Hours of CPU; the sector systems are 2357 x 1188.
2. Finish §2.6 for all 20 planes and 15 lines and extend to non-coordinate planes spanned by eigenvectors of the
   order-18 group (keeps the sector reduction); test order-3 members with the leading cubic form forced to be
   P-independent (product fibration, 10 forms instead of 9). Minutes per plane.
3. Degree 9: the four sporadic kernel directions (sectors (0,0),(2,5),(4,1),(4,4)) and the coset geometry as in 1-2
   (`kernel9.py` already extracts them; exact test ~1-4 min each). The 116 dimension-forced directions are a
   120-dimensional space and need the parametric treatment, not basis-vector tests.
4. Compute Sing(X ∩ V(D8)) and the t -> 0 limit of D8(t) (jobs prepared) to identify the geometric meaning of the
   sporadic annihilator; if D8(0) is a recognisable Fermat join-type element, build the corresponding deformed divisor.
5. The dossier's E2 (fqlx's D32: 608k x 730k sparse modular system) and E3 (companion support divisor pulled back by
   Z = X2^3) remain the "informative but expensive" items; the sector solver reduces E3 to feasible size if the
   support divisor has degree <= 15 (table in dossier §5).
6. Theory: prove the h^{2,0} >= 67 bound (dossier E4) to settle that no divisor with a single simple-elliptic curve
   can work, and examine the j = 0 elliptic family over the lines of W (E5) for a descent argument in the zeta3 sector.
