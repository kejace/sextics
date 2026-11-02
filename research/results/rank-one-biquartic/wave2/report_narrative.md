
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
  the 3-dim pure-kernel subfamily finds no residual below 1e-3 so far (see section 3); the degree-7 landscape has a
  floor around 5e-3 .. 1e-1, i.e. no near-solutions.
* **Degree 8, 9, 10** are the "prime candidate" degrees of the task statement (non-Jacobian kernel elements 9 / 27 /
  >= 294).  The relaxation can no longer exclude by rank or Groebner (468 quadrics in 177 unknowns at d=8 has
  regularity far beyond reach; at d=9 the relaxation is underdetermined, 208 equations in 405 unknowns).  Only the
  full normalized-E2 search applies there; under the current machine load (load average ~1000 on 128 cores) the
  throughput is ~10x below nominal, so the restart counts in section 3 are what could be done, not what was planned.
* **Structured families** (block subrings, 3/4/5-variable cones, bidegree pieces, base-type P(X)M^a g^b and
  Q(Y)L^a f^b, (f,g,L,M)-monomials, Hessian minors): see section 4 for the exact intersections with the admissible
  space; anything that intersects outside J is tested exactly (section 4, e2exact rows).
* **No candidate** with residual below 1e-8 appeared on this instance; the only sub-1e-10 residuals of the whole
  run were the planted-witness validation on a solved instance.
