# Attack 1: quadric-contact-residual  (q = 1, A = X0^3 Y0^3)

Status: IN PROGRESS (resumed after a crash of the first run; sections are filled as experiments finish).
Work dir: /tmp/hodge/results/quadric-contact-residual/  (log.jsonl = every tested denominator)

## 0. Attack plan (written before executing)

F = S0^3 + 2 P0^2 S0 + 5 S0 Q1 + (5/7) P0 Q2 + X1^6+Y1^6+X2^6+Y2^6,
  S0 = X0^2+Y0^2, P0 = X0Y0, Q1 = P1^2 + P1 P2 + P2^2 (P_i = X_iY_i), Q2 = P1 (X2^2 + 2 Y2^2).
A = P0^3.  Visible symmetry group Gamma (order 96) = G_diag(24) x| <s0 = (X0<->Y0), r1 = (X1<->Y1)>;
the residue A*Omega/F^2 is invariant under G_aff (order 48, incl. tau_Y) and ODD under s0 and r1.

Known (fqlx, #125): deg D <= 6 impossible; (J:A)_d = J_d for d <= 7; at d = 8 a single new E1 direction
D8* = 4 G_S^2 - G_P^2 (invariant sector), which fails E2 alone; no certificate keeps (Y1/X1, Y2/X2) passive;
all quotients X/H (H in Gamma^+) have kappa = 0 and are not uniruled.

Reductions used here:
 (R1) Averaging: if (D,N,Q) certifies A then so does its image under any g in Gamma (up to the sign of A under
      tau_Y); HasCertificate.add then gives a Gamma-invariant denominator prod_g g(D).  So WLOG D is
      Gamma-invariant (at the cost of multiplying deg D by the orbit size).  All systematic closure runs below
      are over Gamma-invariant D; non-invariant D are only probed by the structural families.
 (R2) E1 kernel per V-sector: with a = X0-Y0, b = X0+Y0, u = a^2, v = b^2, F = G(u,v,X1,X2,Y1,Y2), the four
      V = <s0,sigma0>-sectors a^i b^j C[u,v,...] have Jacobian parts J_00 = (uG_u, vG_v, G_xi),
      J_10 = (G_u, vG_v, G_xi), J_01 = (uG_u, G_v, G_xi), J_11 = (G_u, G_v, G_xi); E1 for D = a^i b^j D' is
      D' in (J_ij : (v-u)^3).  Computed exactly with Singular (weighted degrees).
 (R3) E2 closure modulo J_{2d} (fqlx's lemma): for D = sum c_k D_k in the E1 kernel with lifts N_k,
      B_D(N) = sum c_k c_l B_{D_k}(N_l) must vanish in R_{2d}; Koszul and F*Q terms drop.  Quadratic system in c.

Experiments planned: (0) sanity/symmetry + D8* and its factors;  (1) E1 kernel profile, all degrees, 4 sectors;
(2) Gamma-invariant E1 kernel bases at d = 8..12 and full fixed-D solves of every basis element;
(3) quadratic E2 closure of the Gamma-invariant family at d = 8, 9, 10 (11, 12 if cheap);
(4) structural families (block pieces, derivatives, Hessian minors, cones over the (X0,Y0) P^3, products of
    invariant quadrics) filtered by E1 then solved;  (5) deformation from two-reflection-cubic-surface-d2.
