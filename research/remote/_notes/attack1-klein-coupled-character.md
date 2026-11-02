# Attack 1: klein-coupled-character  (q = 1, A = X0^3 Y0^3)

**NO WITNESS FOUND.**  All tested denominators are inconsistent (exact mod-p linear algebra, p = 2147483647 unless stated).
Logs: `/tmp/hodge/results/klein-coupled-character/log.jsonl` (one JSON line per tested D / family);
code: `kcc.py` (weight-decomposed engine system, validated against 4 accepted witnesses), `kcc2.py` (Koszul
correction-map solver + E1-kernel machinery, validated by re-deriving 4 accepted witnesses exactly over Q),
`runA_kernel.py`, `runB_jacobian.py`, `runC_lift_served.py`, `runD_scans.py`, outputs `run*.out`.
This note continues (after a crash) a first session whose surviving results were the E1-kernel sweep (P5) and the
engine validation; everything below was recomputed or extended in the resumed session.

Notation: (x,y,z,t,r,s) = (X0,X1,X2,Y0,Y1,Y2).  F = K + M + G with K = x^5y + y^5z + z^5x,
M = (5/7)(x y^3 t r + 2 x^3 z r s + 3 y z^3 t s), G = t^6 + r^6 + s^6.  sigma = diag(zeta,zeta^2,zeta^4,1,1,1),
zeta^7 = 1, weights w = (1,2,4,0,0,0); F has weight 0, A has weight 3.  Projective symmetry group = <sigma, -1>
(order 14), no reflections, no block split (dossier §1.5).  "char k" = sigma-weight k mod 7.

## 1. The obstruction and the geometric mechanism a certificate needs

1. **Sector.**  The residue of A lives in the zeta^3-eigenspace of H^4_prim(X)^tau, a Q(zeta_7)-Hodge structure
   with Hodge numbers (h^{3,1}, h^{2,2}, h^{1,3}) = (26, 123, 37) (dossier §1.6; the sigma-trivial sector (19,128,19)
   is the only one the book's §4.3 cubic-surface quotient supports).  Averaging over any sigma^k, k != 0, kills A,
   so no quotient / Piola / cover mechanism of the book that factors through X/<sigma^k> can reach it.

2. **Support picture.**  A certificate with denominator D exists iff the class is supported on Z = X ∩ V(D)
   (ch. 3; at q = 1 existence over D <=> existence over any power D^m).  For *smooth* Z, Lefschetz gives
   H^2(Z) = Q·h, so the Gysin image is Q·h^2 and no primitive (3,1) class is supported: **Z must be singular, with
   singularities whose resolution carries (2,0)-forms of sigma-character 3 whose Gysin image contains [A]**.
   Isolated rational singularities (nodes, polar-type singularities) contribute nothing.  The natural sources are
   (a) non-isolated singularities with non-rational transversal type (cones: V(D) a cone with a plane vertex P;
   Z is then singular along the curve X ∩ P and fibred over the base curve of the cone in sextic *surfaces* S
   with p_g = 10, so the resolution has many (2,0)-forms), and (b) isolated non-rational singularities (cone over a
   surface with p_g > 0, vertex a line).  In case (a) with vertex P^2_X = {Y = 0}: Sing Z ⊇ B = {K = 0} (a smooth
   Z/7-symmetric plane sextic, H^{1,0}(B) = R(K)_3 carries character 3 twice: x^3, yz^2), fibres
   S_y = {K(x) + λ^2 M(x,y) + λ^6 G(y) = 0} ⊂ P^3_{[x:λ]}, and H^{2,0}(S_y) = R(S_y)_2 contains character 3 once (xy).
   So the geometry *does* produce character-3 (2,0)-forms; the question is whether their Gysin image hits [A].

3. **E1 kills everything of degree <= 8 outside J(F).**  A·D ∈ J(F) is necessary (F irreducible, F ∤ D:
   from A D^2 ≡ D N(F) mod F one gets F | (A D − N(F))).  The kernel of multiplication by [A] on the Y-even and the
   Y-odd Jacobian ring in every sigma-character is **zero for d <= 8** (table in §2.1; mod p, hence over Q by
   semicontinuity: the F_p-kernel bounds the Q-kernel from above, and both contain J_d).  Consequence
   (rigorous): **every legal denominator of degree <= 8 (either Y-parity) lies in J(F)**; in particular every cone,
   Fermat sum, monomial, companion (klein-hesse degree-15 and klein-rational degree-30 D's are excluded already by
   degree parity/weight of the served packets; their low-degree pieces p, r, K, UVW fail E1), block piece (K, M, G,
   F−M), and radial/Euler combination outside J is excluded through degree 8 **as a support** (not only as a polynomial).
   Geometric reading: the cones of item 2(a) with base curve of degree <= 8 do not support [A] even though their
   resolutions carry character-3 (2,0)-forms; the Gysin image misses x^3 t^3 (see §2.4 for the served-target test).

4. **What survives.**  (i) D ∈ J(F)_d, d <= 8 (dimensions 18/63/168 for Y-even d = 6/7/8): generic such Z is
   smooth (hence hopeless by item 2); special members would have to acquire non-rational singularities.  All
   structured members tested (polars F_i, Z_j F_i, products F_i F_j, Euler/radial combinations, random members
   per character, exhaustive GF(101)-scans of the 1- and 2-parameter polar families) are inconsistent (§2).
   (ii) d = 9: the E1-kernel mod J is 24-dimensional (Y-even; chars 0..6: 2,3,3,2,10,1,3) and 24-dimensional
   (Y-odd; chars 1,4,5,6: 5,9,6,4).  Duality explains the sizes (§4.1): the Y-even kernels in chars 1,2,4 are
   forced by dimension (93 → 90, 93 → 90, 93 → 83) and the "unforced" ones (chars 0,3,5,6: 2,2,1,3) are dual to
   the excess of the Y-odd kernels.  fqlx's Jhat is one point of the char-2 slice.  All 48 canonical kernel
   representatives and random points of kernel ⊕ J_9 in every character were tested: inconsistent (§2.2).
   (iii) d = 10: kernel dims 24,46,39,39,54,46,46 (dossier) — tested in §2.2 as far as time allowed.
   (iv) d >= 11: open; the structured candidates xyztrs·F_i (dossier A.5) and Hessian 3×3 block minors (§2.3) fail.
