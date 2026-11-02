
## 8. Interpretation and what remains

* **Rigorous frontier moved from "deg D >= 7 for Y-even D (E1)" to "deg D >= 7 for any parity (E1 + E2
  relaxation)".**  Degrees <= 4: E1 forces D in J and J_d = 0 (dossier).  Degree 5 (both parities): the relaxation
  map on Sym^2(J_5) has full rank 6/6.  Degree 6, Y-odd: the relaxation ideal has a finite quotient (projectively
  empty).  Degree 6, Y-even: the relaxation variety is the single point [F] (the Euler/spurious multiple), length 184,
  at two primes; after slicing away the spurious direction (invariance lemma) the system is projectively empty.  All
  statements are mod p (2^31-1 and/or 1073741789, 1000000007) and lift to C by rank semicontinuity because every
  J-rank used equals its regular-sequence value.
* **Degree 7** is the first open degree.  The Y-even admissible family is 66-dimensional (63 polar + 3 kernel, the
  latter T-base-type per wave 1); its relaxation has 432 independent quadrics (rank Phi = 432 of 2211 pairs, 744
  available), which is far above the 63 unknowns of the sliced system, so emptiness is expected; msolve reached F4
  degree 3 in minutes, the degree-4 round decides (queued).  The numerical full-E2 search on the 66-dim family and on
  the 3-dim pure-kernel subfamily finds no residual below 1e-2 in the corrected metric (absolute residual
  0.011-0.021 at |u| = 1 over all restarts of both local methods); the degree-7 landscape has a floor around
  2e-2 .. 1e-1 with no near-solutions.  (Early LM runs reported rel 2e-6, an artifact of the first metric: the
  optimizer inflated the spurious F-multiple component to |s| = 89; see the metric note in section 0.)
* **Pure-kernel probes.**  The 3-dim degree-7 kernel (a projective plane) has no zero: 40 restarts, best absolute
  residual 0.02.  The 9-dim degree-8 Y-even kernel: 18 L-BFGS restarts all converge to one minimum, absolute residual
  0.121 -- no solution supported on the pure kernel directions at degree 8.  Degree-9 pure and mixed kernels: chain E.
* **Degree 8, 9, 10** are the "prime candidate" degrees of the task statement (non-Jacobian kernel elements 9 / 27 /
  >= 294).  The relaxation can no longer exclude by rank or Groebner (468 quadrics in 177 unknowns at d=8 has
  regularity far beyond reach; at d=9 the relaxation is underdetermined, 208 equations in 405 unknowns).  Only the
  full normalized-E2 search applies there; under the current machine load (load average ~1000 on 128 cores) the
  throughput is ~10x below nominal, so the restart counts in section 3 are what could be done, not what was planned.
* **Structured families** (block subrings, 3/4/5-variable cones, bidegree pieces, base-type P(X)M^a g^b and
  Q(Y)L^a f^b, (f,g,L,M)-monomials, Hessian minors; section 4).  Through degree 9 every structured family meets the
  admissible space only inside J_d, except the Y-base-type class P(X)M^a g^b (the packet's compact-Y class), which
  contains the whole degree-7 kernel (3), the whole degree-8 kernel (9) and 12 of the 27 degree-9 kernel dimensions.
  Hence the first E1-admissible denominators outside the packet's obstruction class are the 15 remaining degree-9
  kernel directions, spanned by K_7 . (S_2(Y) mod M) = {P_7(X) + L M beta_3(X)} x {Y-quadrics not proportional to M}
  modulo J + base type.  These are the target of the focused `kernel-mixed` searches (chain E).  At degree 10 new
  structured admissible directions appear: a 15-dimensional space of Y-odd bihomogeneous (5,5) denominators entirely
  outside J (pulled back from P^2 x P^2, i.e. still base type in the packet's sense), 3 new directions in each
  5-variable cone ring C[X_a X_b Y0 Y1 Y2], and 15 Y-base-type directions; members are tested exactly (e2exact).
* **No candidate** with residual below 1e-8 appeared on this instance; the only sub-1e-10 residuals of the whole
  run were the planted-witness validation on a solved instance.
