# Draft discussion comment (lane quadratic-cubic-pullback, #175) — EXCLUSION / RIGOR lens, wave 2

<!-- STATUS: DRAFT, to be completed when remote jobs finish; placeholders in [[...]] -->

### Degree bounds and certified exclusions for the residual target `X0^3 Y0^3` (q = 1)

All statements below are about the fixed instance `(F, 1, X0^3 Y0^3)`; none is a claim of nonexistence of a
certificate.  "Rigorous" means: exact rational/integer computation, or a modular rank computation whose
characteristic-zero validity is certified (modular rank <= rational rank; equality forced by a Hilbert value or by an
explicit rational kernel).  "Screened" means: a modular computation at one prime without such a certificate.

**1. E1 kernels (rigorous).**  Let `K_d = {D in S_d : A D in J(F)}` (the necessary condition obtained by reducing the
certificate identity modulo `F`).  Then `K_d = J(F)_d` for all `d <= [[9]]` in both Y-parity sectors (primes
2147483647 [wave 1, even d<=8 / odd d<=6] and 2147483629 [all d<=[[9]]]; the degree-(d+6) Jacobian blocks have their
Hilbert ranks, so the modular kernel dimension is an upper bound for the rational one and it equals `dim J_d`).
Consequence: every certificate has `deg D >= 10`, or `D = w . grad F` is of Jacobian ("polar") type.
Every D of degree <= 9 not in J(F) — in particular delta, all cubics in the Q_i except F, all joins, all quadric
sections, the coordinate-plane cones — is excluded at every pole order (prescribed-support theorem at q = 1).
At d = 10 the kernel modulo J has dimension [[294 (even) / 261 (odd)]], equal to the Gorenstein lower bound
`dim R_10 - dim R_16 = 555`, hence exact over Q.

**2. The polar family (new reduction).**  For `D = w . grad F` write `N = A w + H z/6 + c . grad F` (complete by the
Koszul syzygies).  All `H`-terms cancel identically and the certificate condition becomes
`Phi(w) + L(c,w) = F Q` with `Phi(w) = D (grad A . w + A div w) - A w . grad D` quadratic in `w` and
`L(c,w) = sum_{i,j} F_j [(d_i c_ij) D - c_ij d_i D]` bilinear.  Relaxing `W = w w^T`, `V = c (x) w` to independent
unknowns gives a linear system whose kernel projects onto the `W`-coordinates; admissible `w` satisfy `w w^T in pi_W(ker)`.
*Degree 5 (rigorous, exact integer nullspace):* `pi_W(ker) = 0` for the full 6-parameter family `w in C^6`
(no parity assumption).  **No certificate with `deg D = 5` exists; hence `deg D >= 6` for every certificate.**
*Degree 6:* [[result of the exact relaxation in the Y-even and Y-odd sectors; msolve candidate locus]].

**3. Structured Jacobian-type denominators of degree 6-10 (certified fixed-D tests).**  [[table: block Euler fields
E_X F (the deformation of the Fermat join E_X F_0/6 = X0^6+X1^6+X2^6, which serves 100 classes at the Fermat point),
weighted Euler, shape-symmetry fields rho/sigma, grad Q_k . grad F, grad(Q_kQ_l) . grad F, sum Q_i^2 F_i, |grad F|^2,
grad delta . grad F ... : served dimension modulo J, rank-certified]].  A fixed-D verdict is reported as rigorous only
when the modular nullity equals the explicit rational kernel `(G z, 0) + (g F, D div g - g.grad D)` (dimension
`129` at d=5, `244` at d=6, ...).

**4. Degree 10 (screened + structured).**  [[intersections of K_10 with: quintics in Q (252-dim), delta x quartics,
(Q_S)^2 cap (Q_S')^2, plane-order windows; dims; candidates tested; random kernel elements mod p]].

**5. Re-verification of earlier claims.**  h^{3,1}(X)^+ = 208, [A] and [delta] independent mod J_6, diagonal symmetry
group Z/2, F not a join (21 independent second partials), cubics in Q meet J_6 only in C.F, two-plane cancellation
nullity 0, degree-5 polars serve nothing — all reconfirmed at a second prime / exactly; the three degree-5 polar
results of wave 1 are upgraded to rigorous by the kernel certification.  One correction to the earlier notes:
`E_I F = sum_{i in I} z_i dF/dz_i` lies in `J(F)` for every index set `I`, so the block-Euler denominators are polar
candidates (tested above), not E1-failures.

**Reproduction.**  `/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl/` (scripts `w2x.py`, `wave2_excl_*.py`,
JSON receipts with primes, ranks, Hilbert values, structural nullities; README with the exact statements).
[[to be mirrored in the research branch]]
