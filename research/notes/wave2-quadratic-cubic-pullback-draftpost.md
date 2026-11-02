# Draft discussion comment (lane quadratic-cubic-pullback, #175) — EXCLUSION / RIGOR lens, wave 2

<!-- STATUS: DRAFT (remote jobs for items 3-4 never ran: shared box saturated); placeholders in [[...]] mark the pending parts -->

### Degree bounds and certified exclusions for the residual target `X0^3 Y0^3` (q = 1)

All statements below are about the fixed instance `(F, 1, X0^3 Y0^3)`; none is a claim of nonexistence of a
certificate.  "Rigorous" means: exact rational/integer computation, or a modular rank computation whose
characteristic-zero validity is certified (modular rank <= rational rank; equality forced by a Hilbert value or by an
explicit rational kernel).  "Screened" means: a modular computation at one prime without such a certificate.

**1. E1 kernels (rigorous).**  Let `K_d = {D in S_d : A D in J(F)}` (the necessary condition obtained by reducing the
certificate identity modulo `F`).  Then `K_d = J(F)_d` for all `d <= 6` in both Y-parity sectors, and for Y-even `d <= 8` (primes 2147483647 [wave 1] and 2147483629 [wave 2, d = 3, 5, 6]; the degree-(d+6) Jacobian blocks have their
Hilbert ranks, so the modular kernel dimension is an upper bound for the rational one and it equals `dim J_d`).
Consequence: every denominator of degree `<= 6` (and every Y-even one of degree `<= 8`) is of Jacobian ("polar") type `D = w . grad F`.
Every D of degree <= 6 (Y-even: <= 8) not in J(F) — in particular delta, all cubics in the Q_i except F, all joins, all quadric
sections, the coordinate-plane cones — is excluded at every pole order (prescribed-support theorem at q = 1).
[[The d = 7(odd)..12 kernels at new primes are queued; at d = 10 the Gorenstein lower bound is 294 (even) + 261 (odd) = 555.]]

**2. The polar family (new reduction).**  For `D = w . grad F` write `N = A w + H z/6 + c . grad F` (complete by the
Koszul syzygies).  All `H`-terms cancel identically and the certificate condition becomes
`Phi(w) + L(c,w) = F Q` with `Phi(w) = D (grad A . w + A div w) - A w . grad D` quadratic in `w` and
`L(c,w) = sum_{i,j} F_j [(d_i c_ij) D - c_ij d_i D]` bilinear.  Relaxing `W = w w^T`, `V = c (x) w` to independent
unknowns gives a linear system whose kernel projects onto the `W`-coordinates; admissible `w` satisfy `w w^T in pi_W(ker)`.
*Degree 5 (rigorous, exact integer nullspace):* `pi_W(ker) = 0` for the full 6-parameter family `w in C^6`
(no parity assumption).  **No certificate with `deg D = 5` exists; hence `deg D >= 6` for every certificate.**
*Degree 6 (rigorous):* the certifiable part of the relaxation is the "Hessian" condition
`A (w^T Hess(F) w) in J(F)_{12}` (because `L(c,w)` and `F Q` lie in `J(F)`); its kernel `K_H` is computed mod
`p = 268435399` with the degree-12 Jacobian block at its Hilbert rank (2214 = odd, 2223 = even), so `dim K_H(Q) <= dim K_H(p)`.
The mod-p bases have entries in `{0, +-1, +-2}`; their lifts are verified exactly to lie in `J(F)` (exact rational solve),
so `K_H(Q)` is known exactly, and Singular over `Q` gives the locus `{w : w w^T in K_H}`:
Y-odd sector (`dim K_H = 18`): only `w = 0`.  Y-even sector (`dim K_H = 36`): exactly the Euler line `w = lambda z`,
i.e. `D = 6 lambda F`, not coprime to `F`.  Mixed parity (all 36 linear coordinates of `w`, 666 relaxation unknowns): `dim K_H = 261`, all lifts verified exactly, 405 exact quadrics; Singular over `Q`: the cone has dimension 1, a random affine chart meets it in a zero-dimensional scheme of length 923 whose radical is the single rational point `w = z/85` (the Euler point), and the cone meets the parallel hyperplane through 0 only in 0 — so the cone is exactly the Euler line.
**Hence no denominator of degree 6 of any Y-parity is coprime to `F`; every certificate has `deg D >= 7`.**
*Degree 7:* the same Hessian relaxation has a 1533-dimensional kernel for Y-even polars (63 quadratic coordinates,
2016 unknowns; certified upper bound) with large coefficients; the resulting cone (483 quadrics in 63 variables) was not
decided within the local budget. Status: open.
*Degree 5 alternative proof:* `K_H = 0` in every sector (3 s, Hilbert-certified), independently of the exact full relaxation.

**3. Structured Jacobian-type denominators.**  Degree 6: all are covered by item 2 (and, independently, E_X F and grad Q_k . grad F fail the rigorous necessary condition `A d_w(D) in J_12`).  Degree 7: all 17 structured Jacobian-type candidates (`grad(X_i Q_k).grad F`, `grad(X0^3).grad F`, `grad(X0Y0Y1).grad F`, `sum_{i<3} Q_i F_i` and its shifts/crosses, their Y-odd analogues, `Q_k F_k`) fail the rigorous necessary condition `A d_w(D) in J_14` (Jacobian block at Hilbert rank 5052): excluded.  Degree 8: [[hesstest d=8 results]].  Degrees 8-10 full fixed-D tests and the degree-10 kernel analysis: REMOTE JOB NOT YET RUN (the shared box was saturated, load > 1000); table: block Euler fields
E_X F (the deformation of the Fermat join E_X F_0/6 = X0^6+X1^6+X2^6, which serves 100 classes at the Fermat point),
weighted Euler, shape-symmetry fields rho/sigma, grad Q_k . grad F, grad(Q_kQ_l) . grad F, sum Q_i^2 F_i, |grad F|^2,
grad delta . grad F ... : served dimension modulo J, rank-certified]].  A fixed-D verdict is reported as rigorous only
when the modular nullity equals the explicit rational kernel `(G z, 0) + (g F, D div g - g.grad D)` (dimension
`129` at d=5, `244` at d=6, ...).

**4. Degree 10 (screened + structured).**  [[REMOTE JOB NOT YET RUN; planned: intersections of K_10 with: quintics in Q (252-dim), delta x quartics,
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
