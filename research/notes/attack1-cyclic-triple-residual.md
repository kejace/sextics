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

