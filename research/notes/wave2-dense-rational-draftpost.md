# Draft discussion comment — dense-rational (#59) — exclusion / rigor update (wave 2)

**No certificate found; no Lean claim.** (Draft at handoff: items marked PENDING await the remote msolve jobs; do not post before T5 is filled in.) This note reports (i) an independent re-verification of the published
degree bound, (ii) a new rigorous degree bound for tau-stable supports [PENDING: fill in once msolve finishes],
(iii) rigorous exclusions of structured polar-divisor families at degrees 7-8, and (iv) modular/numerical screens.
Scripts and receipts: `/tmp/hodge/results/dense-rational/wave2-excl/` (to be attached as a packet).

### 1. Set-up (same as fqlx's packet)
`E1: N.grad F = A D`, `E2: D div N - N.grad D = F C` after Euler normalization; `J = (F_0..F_5)`, `R = S/J`,
`(J:A)_d = {D : A D in J}`.

### 2. Re-verification of `deg D >= 7` (independent route, new primes)
* `(J:A)_d = J_d` for `d <= 9` in both Y-parities (p = 2147483587 here; p = 2147483647 wave 1; p = 101, 10007 dossier):
  every E1-admissible denominator of degree `<= 9` is a polar divisor `D = w(F)`.  Rigorous: the modular Koszul count
  is attained and the rational kernel reduces injectively.
* Degree 5: exact (rational) rank 21 of `Sym^2 C^6 -> S_4, S -> sum S_ij F_ij`: with `(J:A)_4 = 0` this excludes
  `D = v.grad F` for **every** direction `v` (the packet summary lists the six coordinate directions).
* Degree 6: the necessary condition `w^T Hess(F) w in J_6` (w linear) is a homogeneous quadratic system in each
  tau-sector (even 218 quadrics / 17 unknowns, odd 218 / 18, mixed 426 / 35).  msolve (`-g 1`, p = 1073741789) gives
  DRL leading ideals containing a pure power of every variable, so the only solution over `F_p`-bar is 0; a projective
  variety over Q with a C-point has `F_p`-bar points for every p at which the equations are p-integral, so there is no
  solution over C.  Hence `deg D >= 7` in all sectors (agrees with the packet's rank-one argument).

### 3. Degree 7 [PENDING]
Same reduction with w quadratic (`[A]` injective on `R_8`), 483 quadrics / 60 unknowns per pure sector
(p = 1073741789, 32 threads).  RESULT: [fill: leading ideal / wall time].  Consequence: no certificate whose denominator
has degree 7 and tau-stable support (a legal degree-7 D is squarefree, so tau-stable support <=> tau-pure D).
The mixed sector (951 quadrics / 120 unknowns) is not F4-decidable (generic emptiness only at Macaulay degree 5); it is
screened numerically ([fill: LM starts, best residual]) and rigorously excluded on structured sub-families
(weighted-Euler, block-diagonal, cross, gradient, two-component): [fill].

### 4. Degree 8 structured families [PENDING]
`D = w(F)`, w cubic, necessary condition `w^T Hess(F) w in (J:A)_10^+` (quotient dim 468).  The full sectors (156/159
unknowns) are not F4-decidable (Froeberg count positive through Macaulay degree 8), so we certify emptiness only on
families: [fill list with unknown counts and results].

### 5. Screens (modular at one large prime / floating point; not proofs)
[fill: served-target dimensions (all 18 or not), first-order Fermat deformation (12 seeds), E2 on degree-10 kernel
representatives, LM on degree-8 relaxation.]

### 6. What this does and does not say
Proved: `deg D >= 7` (re-proved); [PENDING: `deg D >= 8` for tau-stable supports]; the named families are excluded at
degrees 7-8.  Not excluded: degree-7 denominators with non-tau-stable support, degree >= 8 outside the named families,
everything at degree >= 10.  No claim of nonexistence.
