## 0. Setup and method (new in wave 2)

* Instance: `F = f(X) + g(Y) + (5/7) L(X)^2 M(Y)`, `A = X0^3 Y0^3`, q = 1.  Diagonal symmetry group `{1, tau}` only,
  so the only character sector is Y-parity: the denominator may be taken Y-even (lossless up to degree doubling);
  Y-odd families are searched too where cheap.  Prior work (dossier, attack 1): E1 forces every Y-even admissible
  denominator of degree <= 6 into J(F) (proved), kernels mod J of dimension 3/9/27 at d = 7/8/9, all natural
  low-degree denominators fail, a few members of the degree-7/8 families fail exactly (rank+1 pattern).
* Search machinery (`wave2/r1lib.py`, `wave2/e2search.py`): the **normalized E2 system** of book 5.2 instead of the
  cleared system.  For an E1-admissible family `D(a) = sum a_m D_m` with lifts `N0(a)` (`N0(F) = A D`), E2 holds iff
  `B_D(N0) + T_D(c) + F C = 0` for Koszul unknowns `c_ij` (degree d-4, parity-restricted) and `C` in `S_{2d-6}^-`
  (complete correction criterion, independent of the lift).  `C` is eliminated exactly by normal-form reduction
  modulo `F` (lex, LM = X0^6); the objective `f(a) = min_c |rho(B(a,a) + M(a) c)|^2` (a quartic in a) is minimized
  by L-BFGS with the exact variable-projection gradient (adjoint products), on the chart `a = Pc u/|u| + S s`
  where `S` spans the spurious F-multiples `F S_{d-6}` (lossless for denominators nonzero mod F).  Per-evaluation
  least squares: 2854x162 (d=6), 5160x420 (d=7), 8661x954 (d=8), 13717x1890 (d=9), 20736x3480 (d=10), versus
  16797x2602 ... 84903x18900 for the cleared system used by the engine's varpro.
* Residual measure `rel` = |r| / (|rho(D div N)| + sum_i |rho(N_i D_i)|) (cancellation-free scale), with
  |D| = 1 in the chart.  A true solution gives rel ~ 1e-12; candidates below 1e-6 are stored, below 1e-8 they are
  rationally reconstructed (continued fractions, several pivots/denominators) and re-solved exactly with the engine
  (`fixed_d.e1_solvable`, `solve_fixed_D_fast`, `verify_witness`, F-divisibility check).
* E1-admissible family in floats: SVD of the Jacobian-generator matrix in degree d+6, kernel of `Q_perp^T (A . m)`,
  lifts by pseudo-inverse (E1 residual ~1e-15).  Exact mod-p counterparts: `e1spaces.py` (kernel reps + lifts),
  `relax.py` (Jacobian-quotient relaxation rank: `B_D(N0) in J_{2d}`, for Jacobian-type D = dF(w) this is
  `A w^T Hess(F) w in J`), `structured.py` (intersection of structured subspaces with the admissible space),
  `e2exact.py` (single-D normalized E2 mod p).
* Validation (remote): accepted witness `coupled-pair-rational-d4-strict` (q=1, d=4): normalized residual 1e-16 at
  the witness point, Koszul-perturbed lifts recovered to 1e-16, random non-Koszul perturbations give O(1);
  gradient check analytic -0.0689627383 vs numeric -0.0689627386; a 12-restart search on that instance found the
  witness from 9 random starts and confirmed it exactly (strict).  On this lane `B_{eps_X F}(A X) = 8 t A L^2 M`
  is reproduced exactly and `D = eps_X F` is inconsistent (reduced rel 0.60, unreduced 0.091); d=6 family
  dimension 18 = dim J_6^+ and d=7 dimension 66 = 63 + 3 agree with the exact profile.
* Rigor of the exact exclusions.  All exact computations are modulo p = 2^31-1 (= `core.PRIMES[0]`); the Jacobian
  ring dimensions found there ((S/J)_10^- = 744, (S/J)_12^- = 866, (S/J)_14^- = 744, ...) equal the regular-sequence
  values, so every J-rank mod p equals its rank over Q.  Consequently (i) the Q-space of E1-admissible denominators
  reduces into the mod-p space, (ii) p-integral lifts reduce to lifts, (iii) normal forms modulo J commute with
  reduction.  If the relaxation map Phi: Sym^2(family) -> (S/J)_{2d} has full column rank mod p, its matrix over Q
  (in a saturated p-integral basis) has full column rank, so Phi(a (x) a) != 0 for every nonzero a over C: no
  denominator in that family passes E2.  For d <= 6 the family is the complete set of E1-admissible Y-even
  denominators (proved over C in the dossier), so a full-rank result there is an unconditional exclusion of that
  degree.  Numerical search results are evidence only; a failed multistart search excludes nothing.
