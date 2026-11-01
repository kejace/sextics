# Attack 1: klein-coupled-character  (q = 1, A = X0^3 Y0^3)

STATUS LINE (updated at end of session): see section 3.

## 0. Attack plan (written before experiments)

Instance. F = K(x,y,z) + x^3 z (10/7) Y1 Y2 + x y^3 (5/7) Y0 Y1 + y z^3 (15/7) Y0 Y2 + Y0^6+Y1^6+Y2^6,
K = x^5 y + y^5 z + z^5 x (x,y,z = X0,X1,X2).  sigma = diag(zeta, zeta^2, zeta^4, 1,1,1), zeta^7 = 1,
fixes F; A has sigma-weight 3 (residue character zeta^3 since the volume form has weight 7 = 0).
Projective diagonal symmetry group = <sigma, -1>, order 14.  No coordinate reflection, no block split.

Known: the sigma-invariant sector (19 (3,1)-directions) is supported by the cubic-surface quotient
W = X/<sigma,tau> -> P^2_Y (book ch. 4.3).  Averaging over sigma kills A.  The discussion (fqlx) excludes
two specific "trace/Hamiltonian" families.  No D/N/Q known.

Plan.
 P0  Sanity + characters: A not in J(F) (done: D = 1 inconsistent).  Dimension of each sigma-weight
     sector of the Y-odd degree-6 Jacobian ring R_6^-; the weight-3 sector is the target's home.
 P1  Companion denominators (klein-hesse D_15, klein-rational D_30 pieces, p^3, r^3 K, K, r, UV/UW/VW
     monomials): E1 block test, then full weight-restricted fixed-D system.  Also what they serve.
 P2  Cone families: Y-only D (vertex P^2_X, base curve in P^2_Y, exceptional divisor B x C_D with
     H^{1,0}(B) carrying weight 3 twice: x^3, y z^2), four-variable D (vertex a line), degrees 2..6.
 P3  Structural pieces: partials, products, radial derivatives, F - coupling, coupling alone, Hessian
     minors, discriminants, perturbations D0 + t M.
 P4  Deformation from the split Klein member (coupling -> 0): klein-* D on this F, and rank/nullity
     distance.
 P5  E1-kernel sweep: dim{D even, deg d, weight k : A D in J(F)} / J(F)_d for d = 1..8, all k.
Mathematics in parallel: Leray/Higgs analysis of the sextic-surface fibration X -> P^2_Y in the
weight-3 sector; what singular divisors could carry a weight-3 (2,0)-form.

All tests exact mod p = 2147483647 (consistency); any consistent system is lifted to Q and verified
with core.verify_witness.  Log: /tmp/hodge/results/klein-coupled-character/log.jsonl.
