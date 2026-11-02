# dense-rational — wave-2 exclusion packet (lens: EXCLUSION / RIGOR), 2026-11-01

Instance: `dense-rational` (`F` = Fermat + (1..17)/9973 on all 230 non-pure even sextic monomials, `q = 1`,
`A = X0^3 Y0^3`), variables `X0,X1,X2,Y0,Y1,Y2`, `tau(X,Y) = (X,-Y)`.  Engine: `/home/kejace/devel/sextics-work/hodge`
(python-flint 0.9.0, msolve 0.6+).  Scripts: `/tmp/hodge/results/dense-rational/wave2_excl_*.py|.sh`.  Receipts, msolve
inputs/outputs and logs: this directory.  One JSON line per experiment: `../wave2-excl.log.jsonl` (field `exp`).
Nothing here is a certificate; everything is an exclusion or a screen.

## 0. Notation and the three reduction facts used

`S = C[X0..Y2]`, `J = (F_0..F_5)` (Jacobian ideal, `F_i = dF/dz_i`), `R = S/J` (Artinian Gorenstein, socle degree 24),
`S_d^{+/-}` Y-even / Y-odd forms, `(J:A)_d = {D in S_d : A D in J}`.  Certificate for `(F,1,A)`: nonzero homogeneous `D`
coprime to `F`, `N_i` of degree `d+1`, `A D^2 - D N(F) + F (D div N - N(D)) = F^2 Q`.

* **Euler normalization.** Adding `phi * Euler` (phi of degree d) to N leaves the cleared identity unchanged (direct
  expansion), and modulo F the identity reads `D (A D - N(F)) in (F)`; with `gcd(D,F) = 1` and `F` irreducible,
  `N(F) = A D + F Rho`, so after `N -> N - Rho Euler/6`: **E1** `N(F) = A D` and **E2** `B_D(N) := D div N - N(D) = F Q`.
* **E1 classification.** `E1` is solvable iff `[A][D] = 0` in `R` (Chapter 5).  If `D = w(F)` then `N = A w` solves E1,
  and every E1 solution is `N = A w + K grad F` with `K` skew (the syzygies of a regular sequence are Koszul).
* **Relaxation.** For `D = w(F)`: `B_D(A w) = D (w(A) + A div w) - A (w . grad w)(F) - A w^T Hess(F) w` and
  `B_D(K grad F) in J`, so **E2 implies `A w^T Hess(F) w in J_{2d}`**, i.e. `w^T Hess(F) w in (J:A)_{2d-6}`.
  (Checked exactly on a random quadratic field: `relaxation-sign-check` in the log.)  The Euler gauge
  `w -> w + l Euler` (`deg l = d-6`) changes `D` by `6 l F` and preserves the relaxation because `Hess(F) z = 5 grad F`;
  a legal `D` is never a gauge element (`D = 6 l F`), and `w(F) = 0` forces `w = 0` below degree 5.

Two rigor lemmas (both elementary):

* **(Sat)** If `V subset Q^n` is a linear subspace defined by integer equations, `Λ = V ∩ Z^n` is saturated, so
  `Λ/pΛ -> F_p^n` is injective and lands in the mod-p solution space.  Hence `dim_{F_p}(solutions mod p) >= dim_Q V`.
  Used for: `dim (J:A)_d(F_p) >= dim (J:A)_d(Q)`; and `dim J_d(F_p) <= dim J_d(Q) = 6 dim S_{d-5}` (Koszul, `d <= 9`),
  with equality when the modular count attains the Koszul value.  Thus `(J:A)_d(F_p) = J_d(F_p)` with the Koszul
  count attained proves `(J:A)_d(Q) = J_d(Q)`.
* **(Proj)** Let `I = (g_1..g_m) subset Z_(p)[c_1..c_n]` be homogeneous.  If `V(I)(C) != {0}` then `V(I)(F_p-bar) != {0}`:
  take a point over a number field, scale it to be integral at a prime above p with a unit coordinate, reduce.  So if
  the DRL leading ideal of `I mod p` contains a pure power of every variable (finite quotient => `V(I mod p) = {0}`),
  then `I` has no nonzero complex solution.  The msolve inputs are p-integral reductions of the rational systems: the
  normal form modulo `J_{2d-6}` (or `(J:A)_10`) is computed with a pivot minor that is a p-unit (the modular rank attains
  the exact dimension), so `NF_Q(v) mod p = NF_p(v mod p)`.

## 1. Exact statements proved here (over C)

* **T1 (E1 bookkeeping, both parities).**  `(J:A)_d = J_d` for `d <= 9`; at `d = 10` the kernel of
  `[A]: R_10^+ -> R_16^-` is exactly 294-dimensional and that of `[A]: R_10^- -> R_16^+` exactly 261-dimensional
  (lower bound `dim R_10 - dim R_16` by counting, upper bound by (Sat)).  Script `wave2_excl_e1kernel.py`, prime
  2147483587 (wave 1: 2147483647 even sector; dossier: 101 and 10007).  By socle duality `<Ax,y> = <x,Ay>`,
  `[A]: R_d -> R_{d+6}` has maximal rank in every degree.  Consequence: every E1-admissible denominator of degree
  `<= 9` is a polar divisor `D = w(F)`, `deg w_i = d-5`; non-polar E1-admissible denominators start at degree 10.
* **T2 (degree <= 4).**  `J_d = 0` for `d <= 4`, so no denominator of degree `<= 4` passes E1 (E1 splits by parity,
  so mixed-parity `D` is covered).
* **T3 (degree 5).**  `D = v . grad F` (`v in C^6`).  E2 forces `A v^T Hess(F) v in J_10`, and `(J:A)_4 = 0` gives
  `v^T Hess(F) v = 0` identically.  The map `Sym^2(C^6) -> S_4`, `S -> sum S_ij F_ij`, has exact rational rank 21
  (`wave2_excl_local_small.py`, fmpq_mat), so `v v^T = 0`, `v = 0`.  All directions `v` (the packet summary names the
  six coordinate ones).
* **T4 (degree 6, all sectors).**  `w` linear; relaxation `w^T Hess(F) w in J_6` (both parities) by injectivity of `[A]`
  on `R_6^{+/-}` (T1).  Homogeneous quadratic systems after the gauge: even 218 quadrics / 17 unknowns, odd 218 / 18,
  mixed 426 / 35.  msolve `-g 1` (p = 1073741789): leading ideals contain a pure power of every variable
  (`quadsys_d6_{even,odd,mixed}_p1073741789_ltcheck.json`).  By (Proj): **no denominator of degree 6**, so
  `deg D >= 7` for every certificate (independent re-proof of `fqlx`'s R3).
* **T5 (degree 7, tau-pure sectors).**  `w` quadratic, relaxation `w^T Hess(F) w in J_8^+` (`[A]` injective on
  `R_8^+`, T1); 483 quadrics / 60 unknowns per sector (`quadsys_d7_{even,odd}_p1073741789.ms`).
  RESULT: [PENDING msolve].  If both are projectively empty: no certificate with a Y-even or Y-odd denominator of
  degree 7.  Since a legal degree-7 `D` is squarefree (prescribed-support theorem at q = 1: `D_red` would be a legal
  denominator of degree `<= 6`), tau-stable support <=> tau-pure `D`, so **`deg D >= 8` for every certificate whose
  denominator support is tau-stable.**
* **T6 (structured families, rigorous when the leading ideal is m-primary).**  `wave2_excl_quadsys_fam.py` restricts the
  relaxation to a linear family `W_fam` of fields (d = 7: quotient `J_8`, both parities in the mixed sector; d = 8:
  quotient `(J:A)_10^+` from the T1 kernel basis at the same prime) and removes `W_fam ∩ gauge`.  A legal
  `D = w(F) + l F` with `w in W_fam` is excluded over C when the leading ideal contains a pure power of every variable.
  Results: `fam_d{d}_{sector}_{family}_p{p}_ltcheck.json` and the log (`msolve-leading-ideal-check`).  [table in the
  report `/tmp/hodge/notes/wave2-dense-rational-excl.md` §3]
* **T7 (serves nothing primitive, by E1 alone).**  For a denominator `D`, `served(D) ⊆ T_E1(D) := {A' in S_6^- : A' D in J}`,
  and `T_E1(D) ⊇ J_6^-` (dim 18).  If `dim T_E1(D) mod p = 18` then (Sat) gives `T_E1(D) = J_6^-` over Q, so `D` serves no
  primitive target over C.  With the Gysin additivity `served(D1 D2) = served(D1) + served(D2)` (Chapter 1/3), any
  product of such `D`'s serves nothing primitive either.  Script `wave2_excl_e1dual.py`; results in the log
  (`e1-dual-served-targets`).  [PENDING]

## 2. Modular / numerical screens (not proofs)

`target-image` (full E1+E2 served-target dimension at one prime; rigorous only through T7), `fermat-first-order-deformation`
(12 subset-Euler seeds, order-eps system mod p), `e2-normalized-fixed-D` (E2 for representatives of the degree-10 E1
kernel, mod p), `e1-kernel-structure`, `lm-hessian-relaxation-d7` (mixed sector, float LM), `lm-hessian-relaxation-d8`
(pure sectors, float LM).  A modular inconsistency is rigorous over Q only if the coefficient-matrix rank attains its
rational value, which is not certified for these fixed-D systems (the exact kernel of the (N,Q)-block is not known).

## 3. Replay

```
cd /tmp/hodge/results/dense-rational; export PYTHONPATH=/home/kejace/devel/sextics-work
python3 wave2_excl_local_small.py                            # T3 + A not in J_6 (seconds)
python3 wave2_excl_signcheck.py                              # relaxation identity check (minutes)
python3 wave2_excl_e1kernel.py d parD 2147483587 [--basis]   # T1 (d = 1..10; d = 10 with --basis)
python3 wave2_excl_quadsys.py d sector 1073741789            # T4/T5 systems (d = 6, 7)
python3 wave2_excl_quadsys_fam.py d sector family prime      # T6 systems (d = 7: any prime; d = 8: prime of the basis file)
msolve -v 1 -t 8 -g 1 -f IN.ms -o OUT.g1.out; python3 wave2_excl_ltcheck.py IN.ms OUT.g1.out
python3 wave2_excl_e1dual.py 2147483587 8                    # T7
```
