# Research dossier: `quadratic-cubic-pullback` (OPEN)

Prepared 2026-11-01 from the repository at `/home/kejace/devel/sextics` (branch `main`, HEAD `83fbf5c`), the saved
discussion thread `/tmp/hodge/discussions/quadratic-cubic-pullback.md` (GitHub discussion #175, node
`D_kwDOVLOOJc4AqVCW`; a single top-level comment, no replies), the instance folder
`data/quadratic-cubic-pullback/` (`instance.rec`, `why.md`, `provenance.json`, `checks.json`), the book chapters 01, 02,
03, 04 (`docs/`), the companion solves named below, the earlier local attack session (note
`/tmp/hodge/notes/attack1-quadratic-cubic-pullback.md`, work directory `/tmp/hodge/results/quadratic-cubic-pullback/`,
engine-2 log `log.jsonl`), and new computations made for this dossier (`/tmp/hodge/scratch/qcp_dossier_checks.py`,
`qcp_polar.py`, `qcp_polar6.py`; outputs quoted inline and appended to the same `log.jsonl`).

Revision note. The previous session crashed while a batch of degree-6 "served-space" computations (`listB_ours.txt`)
was running; only the first item (`join012`) finished. This dossier records exactly what finished, marks what did
not, and adds several cheap closures (parity of the E1 sweep, the Jacobian-ideal "polar" denominators, the symmetry
group, the monodromy correction). Two longer jobs were started at the end of the writing session and their status is
recorded in Appendix A (E1 kernel sweep at `d = 8, 9`; degree-6 polar denominators).

Conventions. Variables `Z = (X0,X1,X2,Y0,Y1,Y2)`; `tau : Y -> -Y`; "Y-odd/even" = parity of total `Y`-degree;
`J = J(F) = (F_0,...,F_5)` the Jacobian ideal, `R = S/J` the Jacobian ring, `R_k` its degree-`k` part;
`H^{3,1}(X)^+` the `tau`-invariant part. A certificate for `(F,q,A)` with denominator `D` of degree `d` means
`A D^2 - q D sum_i N_i F_i + F (D sum_i d_i N_i - sum_i N_i d_i D) = F^{q+1} Q`, `deg N_i = d + 6q - 5`,
`deg Q = 2d - 6` (at every `q`: `deg E = 2d + 6q`, `deg F^{q+1} = 6q + 6`), `gcd(D,F) = 1`.
"E1" denotes the necessary condition `A D in J(F)` (reduce the identity mod `F`; `F` is irreducible and coprime to
`D`). "Served space of `D`" denotes the space of Y-odd degree-6 targets `A'` for which the fixed-`D` system is
consistent (computed mod `p = 2^31 - 1`; by Chapter 1 at `q = 1` this is exactly the space of (3,1)-classes supported
on `V(D) ∩ X`, see 1.7). Its trivial part is `J_6^odd`, of dimension 18, so "imgdim 18" means "serves nothing".

---

## 0. Executive summary

* The target is `(F, q, A) = (sum_{i=0}^5 Q_i^3, 1, X0^3 Y0^3)` with six explicit Y-even quadrics `Q_i`
  (`Q_i = Z_i^2 + (5/p_i)(one X-pair product + one Y-pair product)`, `p_i = 7,11,13,17,19,23`). `X` is the inverse image
  of the Fermat cubic fourfold `C` under the finite flat degree-32 map `phi = (Q_0:...:Q_5)`. Its full linear symmetry
  group among monomial maps is `{1, tau}`; the projective diagonal symmetry group is exactly `Z/2` (new, Smith normal
  form). `[A]` is `tau`-invariant and lies in the 207-dimensional "residual" part of `H^{3,1}(X)^+` (dimension 208)
  complementary to the pulled-back class of `C`, which is the determinant density `delta = det(dQ_i/dZ_j)`.
* The only mechanism that applies to this sextic in the book (Ch. 3, §3.4 weighted-Fano / Piola) serves exactly the
  line of `delta`; `[A]` is independent of `[delta]` in `R_6` (checked mod 5 and mod `2^31-1`). No reflection, no
  paired signs, no block structure, no uniruled quotient carrying the residual sector (the natural quotients `X/tau`
  has `kappa = 0`, `C` carries the wrong sector).
* The discussion contains one comment (`benediamond`, 2026-10-23): an exact computation of the invariant eigensheaf
  splitting `(phi_* O)^+ = O + 6 O(-1) + 9 O(-2)` and a *conditional* obstruction: for any surface of lines on `C`,
  the pulled-back line cylinder `T = X ×_C U` has invariant (3,1) Gysin image inside the `delta` line, provided `T`
  is smooth or has rational singularities. No solve is claimed.
* Rigorous exclusions (mod-`p` ranks equal to the characteristic-zero Hilbert values, hence exact): every
  denominator of degree `<= 4` is impossible; every denominator of degree 5 or 6 (any parity) and every Y-even
  denominator of degree 7 or 8 must lie in the Jacobian ideal `J(F)` (`d = 8` closed during this session, 912 s;
  `d = 9` running). Consequently all "natural" degree-6 supports
  (`delta`, the ten join pullbacks `sum_{i in S} Q_i^3`, every cubic in the `Q_i` other than `F` itself, cones, quadric
  cubes, products of quadrics, Fermat-type blocks) are excluded *at every pole order* (Chapter 1, prescribed-support
  theorem at `q = 1`). The only degree-`<= 7` supports not yet fully closed are the Jacobian-ideal ("polar")
  denominators; the three tested at degree 5 serve nothing (new), the degree-6 test is pending.
* By Gorenstein duality in `R` the first degree at which non-Jacobian E1-admissible denominators can exist is
  `d = 10` (kernel of `A·: R_10 -> R_16` has dimension `>= 1506 - 951 = 555`, of which 294 Y-even if the pattern of
  `dense-rational` repeats), unless `A·: R_d -> R_{d+6}` fails to be injective at `d = 8` or `9` (sweep running).
  So any explicit witness needs `deg D >= 8`, almost certainly `>= 10` (`~19000` unknowns), unless the denominator is
  a polar.
* A concurrent local session (same machine, attack note `attack1-quadratic-cubic-pullback.md`, resumed after an OOM
  kill at 17:44) ran all ten joins `sum_{i∈S} Q_i^3` (every one: E1 fails for `A`; `join012` and `join345` are
  *consistent* for the target `delta`, so the single nontrivial class each join supports is exactly the `delta`
  line), proved that `F` is not a join in any coordinates (the 21 second partials are independent), and developed
  an adjoint-condition criterion `h^{2,0}(Ṽ) = h^1(X, I_{Z_adj}(d))` that reproduces every computed served space
  (Fermat join 100, deformed joins 1, cones 0) and argues that the Fermat two-curve mechanism is not realisable for
  this `F` in any degree (1.9; not independently verified here).
* Honest assessment: no construction in the book or in the accepted solves reaches this target; the residual sector
  of a non-Galois degree-16 cover `X/tau -> C` with monodromy inside `(Z/2) wr S_16` has no known uniruled carrier.
  The realistic path is computational: compute the E1-admissible kernel at `d = 10`, intersect it with structured
  subspaces, and run the full fixed-`D` system on structured candidates. A Lean proof would then be an explicit
  witness (`m3-coupled` already has an accepted degree-10 witness, 360 KB). A nonconstructive Lean proof has no
  available machinery.

---

## 1. Exact target

### 1.1 The record (`data/quadratic-cubic-pullback/instance.rec`)

`q = 1`, `A = X0**3*Y0**3`, and `F` is the 50-term Y-even sextic

```
F = X0^6 + (15/7) X0^4 X1 X2 + (15/7) X0^4 Y0 Y1 + (888750/10793861) X0^3 X1^3 + (1687250/16194277) X0^3 X2^3
  + (2599425/1002001) X0^2 X1^2 X2^2 + (75/289) X0^2 X1^2 Y0^2 + (375/2197) X0^2 X1^2 Y0 Y2 + (375/4913) X0^2 X1^2 Y1 Y2
  + (150/49) X0^2 X1 X2 Y0 Y1 + (375/12167) X0^2 X2^2 Y0 Y1 + (375/1331) X0^2 X2^2 Y1 Y2 + (75/529) X0^2 X2^2 Y2^2
  + (75/49) X0^2 Y0^2 Y1^2 + (15/11) X0 X1^4 X2 + (150/121) X0 X1^2 X2 Y1 Y2 + (15/13) X0 X1 X2^4 + (150/169) X0 X1 X2^2 Y0 Y2
  + (15/17) X0 X1 Y0^4 + (150/289) X0 X1 Y0^2 Y1 Y2 + (375/2197) X0 X1 Y0^2 Y2^2 + (375/4913) X0 X1 Y1^2 Y2^2
  + (375/12167) X0 X2 Y0^2 Y1^2 + (150/529) X0 X2 Y0 Y1 Y2^2 + (375/1331) X0 X2 Y1^2 Y2^2 + (15/23) X0 X2 Y2^4
  + X1^6 + (15/11) X1^4 Y1 Y2 + (900250/2352637) X1^3 X2^3 + (375/343) X1^2 X2^2 Y0 Y1 + (375/6859) X1^2 X2^2 Y0 Y2
  + (75/361) X1^2 X2^2 Y1^2 + (75/121) X1^2 Y1^2 Y2^2 + (375/343) X1 X2 Y0^2 Y1^2 + (375/6859) X1 X2 Y0^2 Y2^2
  + (150/361) X1 X2 Y0 Y1^2 Y2 + (15/19) X1 X2 Y1^4 + X2^6 + (15/13) X2^4 Y0 Y2 + (75/169) X2^2 Y0^2 Y2^2
  + Y0^6 + (15/17) Y0^4 Y1 Y2 + (1563750/4173281) Y0^3 Y1^3 + (1132000/15069223) Y0^3 Y2^3 + (33613425/55190041) Y0^2 Y1^2 Y2^2
  + (15/19) Y0 Y1^4 Y2 + (15/23) Y0 Y1 Y2^4 + Y1^6 + (780500/6539203) Y1^3 Y2^3 + Y2^6
```

**Sparse presentation (verified exactly, `F - sum Q_i^3 = 0`):**

```
Q0 = X0^2 + (5/7)  (X1 X2 + Y0 Y1)        Q3 = Y0^2 + (5/17) (X0 X1 + Y1 Y2)
Q1 = X1^2 + (5/11) (X0 X2 + Y1 Y2)        Q4 = Y1^2 + (5/19) (X1 X2 + Y0 Y2)
Q2 = X2^2 + (5/13) (X0 X1 + Y0 Y2)        Q5 = Y2^2 + (5/23) (X0 X2 + Y0 Y1)
F  = Q0^3 + Q1^3 + Q2^3 + Q3^3 + Q4^3 + Q5^3 .
```

Each quadric is `Z_i^2 + c_i (X-pair_i + Y-pair_i)` with `c_i = 5/p_i`, `p = (7,11,13,17,19,23)`; each of the three
X-pair products `X1X2, X0X2, X0X1` and each Y-pair product `Y0Y1, Y1Y2, Y0Y2` occurs in exactly two quadrics
(`X1X2` in `Q0,Q4`; `X0X2` in `Q1,Q5`; `X0X1` in `Q2,Q3`; `Y0Y1` in `Q0,Q5`; `Y1Y2` in `Q1,Q3`; `Y0Y2` in `Q2,Q4`). Modulo 5
every `Q_i` is `Z_i^2` and `F` is the Fermat sextic. `F` therefore sits in the six-parameter family
`F(c) = sum_i (Z_i^2 + c_i m_i)^3`, `m_i` the fixed cross forms, with `F(0) = Fermat` and the symmetric members
`c_i ≡ c` (see 1.3).

Monomial support: 50 monomials, all Y-even, including all six pure sixth powers; `F` is *not* dense (236 is the
Y-even count), and not a function of `X_i^2, Y_i^2` (there are odd-exponent monomials such as `X0^4 X1 X2`).

### 1.2 Degrees, parity, certificate dimensions

* `deg F = 6`, `q = 1`, `deg A = 6`, `A` is Y-odd (`Y0^3`), `F` is Y-even: the residue `A Ω/F^2` is `tau`-invariant
  (`A` and `Ω = iota_Euler dZ` are both `tau`-odd), so `[A] ∈ H^{3,1}(X)^+`.
* For a Y-even denominator of degree `d`: `N_0,N_1,N_2` are Y-odd and `N_3,N_4,N_5` Y-even of degree `d+1`; `Q` is Y-odd
  of degree `2d-6`; the identity lives in Y-odd degree `2d+6`. (For a Y-odd `D` the parities of the `N_i` flip; the
  identity is then Y-even.) Chapter 1 shows that parity can be imposed by averaging without loss.
* Monomial counts `S_k^{even}/S_k^{odd}`: `k=1: 3/3, 2: 12/9, 3: 28/28, 4: 66/60, 5: 126/126, 6: 236/226, 7: 396/396,
  8: 651/636, 9: 1001/1001, 10: 1512/1491, 11: 2184/2184, 12: 3108/3080`.
* Fixed-`D` linear system sizes at `q = 1` (Y-even `D`):

| `d` | unknowns in `N` | unknowns in `Q` | equations (Y-odd monomials of degree `2d+6`) | `dim S_d^{even}` (space of `D`) |
|---:|---:|---:|---:|---:|
| 5 | 1386 | 60 | 10152 | 126 |
| 6 | 2376 | 226 | 16797 | 236 |
| 7 | 3861 | 636 | 26532 | 396 |
| 8 | 6006 | 1491 | 40326 | 651 |
| 9 | 9009 | 3080 | 59332 | 1001 |
| 10 | 13104 | 5796 | 84903 | 1512 |
| 11 | 18564 | 10152 | 118608 | 2184 |
| 12 | 25704 | 16797 | 162248 | 3108 |
| 13 | 34884 | 26532 | 217872 | 4284 |
| 14 | 46512 | 40326 | 287793 | 5832 |

The engine solves these by random row compression to `(#cols + 40)` rows followed by `nmod_mat.rref` (single-threaded
FLINT); observed `~30-50 s` at `d = 6` (2602 columns), so a cubic scaling predicts `~3-5 h` per candidate at `d = 10`
(`~18900` columns) and `~10 h` at `d = 12` with the present code. Sparse elimination or a multithreaded rank would
cut this substantially (the matrices have `~50` nonzeros per column).

### 1.3 Symmetry (recomputed for this dossier)

* **Projective diagonal symmetry group:** exactly `Z/2 = {1, tau}`. The Smith normal form of the lattice of exponent
  differences of the 50 monomials, in the rank-5 degree-zero lattice, is `(1,1,1,1,2)` (the six pure sixth powers make
  the group finite; Chapter 1 "exponent description").
* **Monomial symmetries (permutations composed with sign changes, 46080 maps tested):** only the identity and `tau`.
  In particular there is no coordinate reflection, no paired-sign element, no exchange `X_i <-> Y_i`, no cyclic
  shift. `checks.json` agrees: `coordinate_reflections: []`, `coordinate_paired_signs: []`.
* **Continuous symmetries:** none; `checks.json` records `hessian_symmetrizer_dimension_mod101 = 1` (only scalars).
* **Residue character:** the only character available is the `tau`-parity; `[A]` is invariant. There is *no finer
  character decomposition* of the certificate problem: by Chapter 1, for this `F` the lattice `L(F)` has index 2 in
  the degree-zero lattice, so all Y-odd degree-6 targets lie in a single coset and no averaging beyond parity is
  possible.
* **Broken shape symmetries (useful for deformation ideas, not for the fixed target):** the cyclic shift
  `rho: (X0,X1,X2,Y0,Y1,Y2) -> (X1,X2,X0,Y1,Y2,Y0)` maps `Q_i` to the quadric of the same shape as `Q_{i+1}` and the
  exchange `sigma: X_i <-> Y_i` maps the shape of `Q_i` to that of `Q_{i+3}`; both would be symmetries if all `c_i`
  were equal. They are broken only by the six distinct constants `5/7, 5/11, 5/13, 5/17, 5/19, 5/23`. Under `sigma`
  (three transpositions, `det = -1`) the residue of `X0^3 Y0^3` would be *anti*-invariant; under `rho` (`det = +1`)
  `A` moves in the orbit `X0^3Y0^3, X1^3Y1^3, X2^3Y2^3`. So even at the symmetric point `c_i ≡ c` the target would lie
  in a non-trivial character, the same shape of obstacle as in `diagonal-permutation` and
  `symmetric-sextic-residual`.
* No non-diagonal symmetry was noted by anyone in the thread.

### 1.4 Smoothness

`checks.json`: good reduction to Fermat modulo 5. Modulo 5, `6 ≡ 1`, so `∂_i F ≡ Z_i^5`, whose only common zero is the
origin; the special fibre of the sextic over `Z_(5)` is smooth, hence so is the generic fibre (openness of
smoothness for the proper flat family). Equivalently (benediamond): mod 5 the quadrics are the coordinate squares,
so `phi` has no base points, hence is finite, and `X = phi^{-1}(C)` with `C` the smooth Fermat cubic is smooth
because `dF = dG(Q)·dQ` and `dQ` is invertible off `delta = 0`... (the clean statement is the mod-5 argument; the
instance rationale uses "smooth complete-intersection base change" for the independence statement in 1.6).

### 1.5 The finite cover and its invariants

* `phi = (Q_0 : ... : Q_5) : P^5_Z -> P^5_U` is finite flat of projective degree 32 (affine degree `2^6 = 64`, divided
  by `z ~ -z`); `phi^* O(1) = O(2)`. `X = phi^{-1}(C)`, `C = {sum U_i^3 = 0}`, `f = phi|_X : X -> C` finite flat of
  degree 32. Since the `Q_i` are `tau`-invariant, `f` factors as `X -> X/tau -> C` with `X/tau -> C` of degree 16.
* `delta = det(∂Q_i/∂Z_j)`: Y-odd sextic with 64 terms (`delta.txt`), `delta ≡ 4 X0X1X2Y0Y1Y2 = -X0X1X2Y0Y1Y2 (mod 5)`.
  `V(delta) ∩ X = R` is the ramification divisor, `R ∈ |O_X(6)|`. Its image `B_phi = phi(V(delta))` is a hypersurface of
  degree 48 in `P^5_U` (`phi_*(6H_Z) = 96 H_U` since a hyperplane maps to a degree-16 hypersurface, and `V(delta) -> B_phi`
  has degree 2 through `tau`); the branch divisor of `f` is `B_phi ∩ C`. (The attack note's "degree 48" is this; note
  `f^{-1}(B) = R ∪ R'` where `R'` collects the non-ramified points of ramified fibres, `R' ∈ |O_X(84)|`.)
* **Monodromy (correcting the first version of the attack note; its current version agrees).** The first draft
  asserted that the monodromy of `f` is `S_32` and that the residual local system is the standard representation
  of `S_32`. This is impossible: `tau` is a nontrivial deck
  transformation of `f` (`Q∘tau = Q`), and the deck group of a connected cover is the centralizer of the monodromy
  group in `Sym(fibre)`; the centralizer of `S_32` is trivial. Hence the monodromy group `G` lies in the centralizer
  of a fixed-point-free involution, `G ⊂ (Z/2) wr S_16`. Its image in `S_16` is the monodromy of `X/tau -> C`, which is
  generated by transpositions (meridians around the branch divisor, simple ramification generically, `C` simply
  connected) and transitive, hence equals `S_16` if the branch divisor of `X/tau -> C` has generically simple
  ramification. Then
  `f_* Q_X = Q ⊕ V_15 ⊕ W_16`, with `V_15` the standard representation of `S_16` (the residual part of `H^4(X/tau)`) and
  `W_16` the `tau`-anti-invariant rank-16 piece. The residual invariant `H^{3,1}` (207-dimensional, below) is the
  (3,1)-part of `H^4(C, V_15)`. Whether `G` is all of `(Z/2) wr S_16`, or whether `X/tau -> C` is imprimitive (which
  would produce an intermediate cover), has not been computed (experiment E7).
* **Invariant eigensheaves (benediamond, proved in the comment, valid for the actual quadrics):**
  `B = C[Z]/(Q_0..Q_5)` is a complete intersection with Hilbert series `(1+t)^6`, so
  `phi_* O = O ⊕ 15 O(-1) ⊕ 15 O(-2) ⊕ O(-3)`, and the equivariant Koszul trace `Tr(tau | B; t) = (1-t^2)^3` gives
  `(phi_* O)^+ = O ⊕ 6 O(-1) ⊕ 9 O(-2)` (ranks 1, 6, 9 = 16 ✓).
* **Hodge numbers (benediamond, reverified):** Y-odd degree-6 monomials: `3·C(7,2) + C(5,2)^2 + C(7,2)·3 = 226`;
  Y-odd degree-6 Jacobian relations `Z_i F_j` with exactly one index in the Y-block: 18, independent (rank 18 mod
  `2^31-1`); so `h^{3,1}(X)^+ = 208`, `h^{3,1}(X)^- = 236 - 18 = 218`, `h^{3,1}(X) = 426` (= number of moduli; `X` is
  Calabi–Yau, `K_X = O`). `H^{3,1}(C)` is one-dimensional and pulls back to the `delta` line; the residual invariant
  space has dimension 207, and `[A]` has nonzero projection to it (`A`, `delta` independent in `R_6`, 1.6).
* **Divisors of `X` that come for free:** the 405 planes of the Fermat cubic (`U_a = -ω^k U_b` for the 15 pairings and
  27 sign choices) pull back to 405 K3 surfaces `{L_1 = L_2 = L_3 = 0} ⊂ X`, complete intersections of three
  quadrics of degree 8, each rigid in `X` (`chi(N_{S/X}) = 4 - 96 < 0`). Lines `ℓ ⊂ C` pull back to curves
  `f^{-1}(ℓ)`, complete intersections of four quadrics (degree 16, genus 17). Hyperplane sections of `C` through an
  Eckardt point (`U_a + ω^k U_b = 0`) are cones over cubic surfaces (rational singularity); their pullbacks are the
  quadric sections `{Q_a + ω^k Q_b = 0} ∩ X`, singular only at the 32 preimages of the vertex (cone over a cubic
  surface pulled back through a local square root is still rational unless the ramification there has corank `>= 2`:
  weights `k/6 + (4-k)/3 <= 1` needs `k >= 2`, a finite set of points).

### 1.6 Jacobian-ring facts and the E1 ("annihilator") profile — the rigid part

Hilbert function of `R = S/J` (`J` generated by a regular sequence of six quintics):
`dim R_k = 1, 6, 21, 56, 126, 246, 426, 666, 951, 1246, 1506, 1686, 1751, 1686, 1506, 1246, 951, 666, 426, 246, 126, 56,
21, 6, 1` for `k = 0..24`; socle degree 24.

Necessary condition for any certificate (any `q`, any `D`): `A D ∈ J(F)` in degree `d + 6`. Equivalently `D` lies in the
annihilator `ann_R(A)` in degree `d`, lifted to `S`; `J_d ⊂ (J : A)_d` always. By Gorenstein duality,
`dim ann_R(A)_d >= dim R_d - dim R_{d+6}`, which is `0` for `d <= 9` and `555, 1020, 1325, 1440, ...` for `d = 10, 11, 12, 13`.

**Kernel sweeps (exact; the mod-`p` ranks of the Jacobian blocks equal the characteristic-zero Hilbert values, so the
mod-`p` kernel dimension is an exact upper bound that is attained by `J_d`):**

| `d` | Y-even `D`: `dim K_d`, `dim J_d^{even}` (file `kernel_sweep.jsonl`) | Y-odd `D`: `dim K_d`, `dim J_d^{odd}` (new, `qcp_dossier_checks.py`) |
|---:|---|---|
| 1 | 0, 0 | 0, 0 |
| 2 | 0, 0 | 0, 0 |
| 3 | 0, 0 | 0, 0 |
| 4 | 0, 0 | 0, 0 |
| 5 | 3, 3 | 3, 3 |
| 6 | 18, 18 | 18, 18 |
| 7 | 63, 63 (399 s) | not run |
| 8 | 168, 168 (912 s; `cols 5796, rowsJ 6006, rankJ 5052, rank_stack 5535`) | not run |
| 9 | running (Appendix A) | not run |

**Consequences (rigorous).** Write any `D = D_even + D_odd`. Since `J` is `tau`-stable, `A D ∈ J` iff
`A D_even ∈ J` and `A D_odd ∈ J`. Hence:

1. No certificate with `deg D <= 4` exists (both parities).
2. Any denominator of degree 5 or 6 lies in `J(F)`: `d = 5`: `D ∈ span(F_0,F_1,F_2) ⊕ span(F_3,F_4,F_5)`; `d = 6`:
   `D = sum ℓ_i F_i` with linear `ℓ_i`. Any Y-even denominator of degree 7 (resp. 8) is `sum q_i F_i` with quadrics
   (resp. cubics) `q_i`.
3. Because at `q = 1` support is independent of pole multiplicity (1.7), **every divisor `V(D_red) ∩ X` with
   `deg D_red <= 6` whose equation is not in `J` supports no part of `[A]` at any pole order**, and the same holds
   for Y-even `D_red` of degree 7 or 8. This disposes of `delta`, every cubic in the `Q_i` (the 56-dimensional space of
   cubics in `Q` meets `J_6^{even}` exactly in `C·F`, new), all ten joins `sum_{i ∈ S} Q_i^3`, all quadric
   sections and their cubes, products `Q_iQ_j`, cones in `<= 4` variables, binary blocks, Fermat-type blocks, and
   indeed every degree-`<= 6` polynomial not in `J`, whether or not it was ever run through the engine.
4. The target is not in `J_6` and is independent of `delta` modulo `J_6` (ranks 18/19/20 mod `2^31-1`; `checks.json`
   gives the mod-5 version).

### 1.7 Chapter-1 facts that control the search

* **Prescribed-support theorem (docs/01, "Prescribed support and pole multiplicity"):** for `1 <= q <= 4` the
  following are equivalent: the residue of `A` vanishes on `X \ V(D)`; a certificate exists with poles on `V(D)`; a
  certificate exists with denominator `D_red^q`. **At `q = 1`, existence over `D` is equivalent to existence over any
  `D^m`; an exact failure of the fixed-`D` system excludes every power of `D`.** So the engine's served space of
  `D` *is* the space of (3,1)-classes supported on `V(D) ∩ X`, and E1-failure of `D_red` kills the whole divisor.
* **Divisor support needs bad singularities (docs/03 §3.1):** on a smooth sextic fourfold a divisor supporting a
  nonzero (3,1) class must have a component that is non-normal or has non-rational singularities; exceptional
  components contracted to codimension `>= 2` contribute nothing to (3,1). Smooth or rationally-singular divisors
  (smooth quadric sections, generic polar divisors, generic line cylinders) serve nothing.
* **Combination:** certificates add (`HasCertificate.add`), with product denominators; served spaces add when
  supports are united. A bundle over several denominators can only help if the union of the served spaces contains
  `[A]`; so far every served space computed lies in `J_6 + C·delta`.
* **Jacobian reduction:** if `A ∈ J` the class is zero (for invariant `q = 1` targets); `A ∉ J_6` here.

### 1.8 The Fermat limit (`c = 0`) and what the deformation destroys

At `F_0 = sum Z_i^6`, the target `X0^3 Y0^3` has the strict subset-Euler certificate (docs/02 §2.2)
`v = A E_X / (6 (X0^6 + X1^6 + X2^6))`, `E_X = X0 ∂_{X0} + X1 ∂_{X1} + X2 ∂_{X2}` (verified exactly by
`core.verify_witness`; witness file `witness-candidate-fermat-92047.json`: `N = (X0^4Y0^3, X0^3X1Y0^3, X0^3X2Y0^3, 0,0,0)/6`,
`Q = 0`). Its support `V(X0^6+X1^6+X2^6) ∩ X_0` is the join of the two plane Fermat sextic curves
`Γ_X = {Y = 0} ∩ X_0`, `Γ_Y = {X = 0} ∩ X_0` (genus 10 each); a smooth model is a `P^1`-bundle over `Γ_X × Γ_Y`, with
`h^{2,0} = 100 = 10·10`, and the engine finds exactly `118 - 18 = 100` served classes, the products `m_1(X) m_2(Y)` of
cubic monomials (`fermatX/fermat`, engine 2, 408 s). On the cubic side this is the pullback of the join of two
Fermat plane cubics `{U_0^3+U_1^3+U_2^3 = 0 = U_3^3+U_4^3+U_5^3}`.

For the deformed `F` the analogous divisor is `f^{-1}(join) = V(Q_0^3+Q_1^3+Q_2^3) ∩ X` (singular along the degree-48
curve `C_1 = {Q_0=Q_1=Q_2=0} ∩ X` with transversal simple-elliptic type). The engine (`join012`, engine 2, 690 s)
finds served dimension `19 = 18 + 1`: the pullback supports a single nontrivial invariant class and **not** `[A]`.
The concurrent session identified that class: the fixed-`D` system with target `A := delta` is *consistent* for
`join012` and for `join345` (engine 2, rank 2358, nullity 244, 489 s / 329 s; `delta_served.py`), while E1 already
fails for `A`. So each deformed join supports exactly the `delta` line. (For `Fpart(6) = F_{60}` and for
`X0^6+X1^6+X2^6` on the actual `F`, the `delta` target is inconsistent: those divisors support nothing.) All ten
joins have the same compressed rank/nullity `2214/162` with E1 false for `A` (`out_listC_ours.txt`, seven of the
nine remaining joins finished at the time of writing). This is consistent with Chapter 3 §3.1: the extra holomorphic 2-forms on a
resolution live on exceptional divisors over the curves `C_1, C_2` and are contracted to codimension 3, so they
have zero (3,1) image; the join downstairs has `h^{2,0} = 1·1 = 1`. The attack note's hope that these forms land in
the residual sector is therefore not borne out, and by 1.6(3) the join is dead at every pole order anyway. The
subset-Euler structure itself is destroyed because `E_X(Q_i)` is not proportional to `Q_i` (each `Q_i` mixes
`X`-pairs and `Y`-pairs).

### 1.9 The adjoint-condition criterion of the concurrent session (reported, not independently verified)

The current attack note (§1.3–1.4) proposes an exact predictor of the served space from the singular locus of
`V = X ∩ V(D)`: since `V` is a complete intersection, `H^2(O_V) = 0`, so all of `H^{2,0}(Ṽ)` of a resolution comes from
the singularities, and by Serre duality plus Grauert–Riemenschneider

`h^{2,0}(Ṽ) = h^1(V, K_V ⊗ I_adj) = h^1(X, I_{Z_adj}(d))`, `K_V = O_V(d)`,

i.e. the failure of the adjoint subscheme `Z_adj` to impose independent conditions on `|O_X(d)|`:
`h^{2,0}(Ṽ) = h^0(O_{Z_adj}(d)) - (h^0(O_X(d)) - h^0(I_{Z_adj}(d)))`. Transversal cone over a plane curve of degree `m`
along a curve `Γ`: adjoint ideal `I_Γ^{m-2}` (`m = 2`: rational, no conditions; `m = 3`: simple elliptic, `I_Γ`;
`m = 6`: `I_Γ^4`); isolated singular points impose independent conditions and never help. The `tau`-invariant (3,1)
sector corresponds to the Y-odd part of the cokernel. Checks reported: Fermat join (`2·280` conditions, kernel
`span(G,H)/F` of dimension 1, `h^{2,0} = 560 - 460 = 100` in the odd sector) matches the solver's 118 − 18; deformed
join (singular exactly along `C_1 = {Q_3=Q_4=Q_5=0} ∩ X` and `C_2 = {Q_0=Q_1=Q_2=0} ∩ X`, curves of genus 145 in `|6H|`
on smooth degree-8 K3 surfaces, simple-elliptic transversal type) gives `h^{2,0} = 1`, the pullback of the resolved
join of two plane cubics in `C`, Gysin image `= delta` — matches 19 − 18 and the `delta`-consistency; cones
`G(Q_a,Q_b,Q_c)` with `G` not a sum of cubes: `h^{2,0} = 0`, so the planned Hesse/random cone lists were (correctly)
not run; cones with vertex a plane: `h^1(I_Γ^a(d)) = 0` for a single plane sextic.

Its §1.4 goes further: the Fermat mechanism needs `V` singular along *two* curves with transversal orders
`(m_1, m_2)`, `m_1 + m_2 >= d + 6`, which for the coordinate planes forces `D|_X` to kill a window of bidegrees that
`F_{42}, F_{24} ≠ 0` make impossible (12-unknown cancellation system `c_{20}F_{24} + c_{02}F_{42} = 0` has nullity 0 at
`d = 8`; the table in the note lists `d = 6, 8, 10` with odd-sector bound 0 in every realisable case), and any other
pair of complementary planes needs `F` to be a join in those coordinates, which it is not (21 second partials
independent, `jointest.py`; a join has rank `<= 12`). Conclusion claimed: **no denominator whose only non-rational
singularities are transversal cones along plane sections of `X` can serve `[A]`, in any degree.**

Skeptical remarks. (i) The criterion uses `h^{2,0}(Ṽ)` as an upper bound for the served dimension, which is right
(Gysin images of the strict transform); it also uses it as a *prediction* of the served dimension, which requires the
Gysin map to be injective — true in the computed cases but not in general. (ii) The "any degree" statement is about a
specific singularity configuration (two thickened plane curves) and does not exclude supports singular along
non-planar curves (e.g. `C_1`-type curves on the K3 surfaces, curves on quadric sections) or non-normal supports.
(iii) The inputs (`adjcount.py`, `struct1.py`) are exact mod-`p` rank computations; I did not rerun them.

---

## 2. Companion and related solved instances

`provenance.json` lists no `sameVarietyInstances`; no accepted instance shares this `F`. The related material is:

| Instance | Method / denominator | Why it does not transfer here |
|---|---|---|
| `cyclic-triple-fano` (`A = X0 X2^2 Y0^3`, same shape of theorem as why.md cites) | The exact application of docs/03 §3.4 (PR #153, `benediamond`, built on `fqlx`/`rbajaj5`): substitute `Z = X2^3`, weighted quotient `P(1,1,1,1,1,3)`, free weighted line `[0:s:0:is:t:it]` with Jacobian `-22500/2401`, polynomial-family integration, finite trace, legal clearing, cubic transport `D = φ(d)`, `Ñ_j = X2^2 φ(N_j)`, `Ñ_{X2} = φ(N_Z)/3`. **Nonconstructive**; no `witness.json`, no degree bound. | The §3.4 theorem serves only `A = delta·b(Φ)` with `wt(b) = 6(q+1) - M`. Here `M = 12`, so at `q = 1` `b` is a constant and the served sector is the single line `C·delta`; `X0^3Y0^3` is independent of `delta` in `R_6`. The open companion `cyclic-triple-residual` (`A = X0^3Y0^3` on the cyclic-triple sextic) is the same complementary-character situation. |
| `reflection-coupled-fano` (`A = X0^3 Y2^3`) | docs/03 §3.3 one reflection `u = Y2^2`, modular weighted quartic mod 5, 24×24 deformation matrix, étale lifting, integration + trace. Nonconstructive. | This `F` has no reflection (1.3). The open `reflection-coupled-residual` is the complementary (reflection-even) sector. |
| `holdout` (`A = X0^3Y0^3`, `F = Fermat + X0^3X1Y0Y1 + X1^3X2Y1Y2`), `holdout-*` family, `holdout-bundle`, `holdout-moving-branch`, `holdout-ternary-deformation` | `fqlx`: Fermat-curve cover creating a reflection, moving polynomial conic on the quotient, integration and trace; Lean library `HoldoutCover`. Nonconstructive. | Needs the specific holdout coefficient structure (monomials `X0^3X1Y0Y1`, …) and a cover on which a reflection appears; nothing analogous exists for `sum Q_i^3`. |
| `cycle-three`, `mixed-biweight`, `mixed-polar` (`A = X0^3Y0^3` or all-six) | Paired-sign cubic: with `R = X2^2`, `x = X0^2/R`, `y = X1^2/R`, the sextic is `R^3 U(x,y,1)` for a plane cubic `U` over passive ratios; Hessian covariants map to the Weierstrass Jacobian; integration + trace (`TransitiveCubic`, `ResidueTransport`). Nonconstructive. | This is the "three-square cubic-surface reduction used by the paired-sign constructions" that why.md explicitly says does not apply: `F` is not a polynomial in `X_i^2, Y_i^2` (e.g. `X0^4X1X2`, `X0X1Y0^4`), there are no paired signs (`checks.json`), and the `Q_i` mix the coordinate pairs. |
| `one-reflection`, `one-reflection-y2-sector` | One coordinate reflection, §3.3. | No reflection. |
| `family-block-cube` (`A = X0^3Y0^3`, explicit `d = 6`, `D = (2X0^2+9Y0^2)^3/243 + Y1^6 + (1681/1701) Y2^6`) and `family-hesse-web` (`A = X0^3Y0^3`, explicit `d = 6`, norm denominator `a^3 + t a^2 x + 4x^3`, `a = X0^2+Y0^2`, `x = X1^2`) | Binary-block contact / Hesse-pencil contact channels (docs/02 §2.3–2.5); strict degree-6 witnesses. | Require block structure `F ∈ C[X_i^2, Y_i^2]` with a binary or ternary cubic block; absent here. Moreover **every degree-6 denominator not in `J` is excluded** for this `F` (1.6), so no degree-6 contact denominator can exist. |
| Fermat (book, all targets) | Subset-Euler / Nambu fields, §2.2. | The subset-Euler condition `E_I^2(F) = 6 E_I(F)` fails for every block `I` once `c ≠ 0` (1.8). |
| `m3-coupled` (explicit `d = 10`, 360 KB `Witness.lean`), `klein-rational` (`d = 30`), `family-sylvester-norm` (`d = 12`) | Evidence that the Lean verifier accepts large explicit witnesses. | Relevant to feasibility (section 6), not to the mathematics. |

The other open "residual" lanes with the same target `X0^3Y0^3` (`cyclic-triple-residual`, `reflection-coupled-residual`,
`quadric-contact-residual`, `off-paired`, `rank-one-biquartic`, `four-plus-two-residual`, `short-power-presentation`,
`eight-plane-configuration`, `dense-rational`) share the pattern "target outside the character/density sector of the
only available construction"; their dossiers (`/tmp/hodge/notes/dossier-*.md`) contain no technique that transfers,
with one exception worth noting: for `dense-rational` the E1 classification was extended to `d = 9` (`K_d = J_d`) and the
`d = 10` kernel has dimensions 294 (even) / 261 (odd), exactly the Gorenstein lower bound `555` split by parity.
The same numbers are expected here.

---

## 3. Prior work, exhaustively

### 3.1 The discussion thread (#175)

**Comment `DC_kwDOVLOOJc4BJDaw`, `benediamond`, 2026-10-23T17:06:33Z (never edited; no replies).** Title: "The quadratic
cover's invariant eigensheaves, and a conditional obstruction for pulled-back line cylinders". Assisted by ChatGPT.
"No solve is claimed."

Proven in the comment (exact, valid for the actual quadrics):

1. `phi` is finite flat of degree 32, as is `f : X -> C`; `phi^*O(1) = O(2)`; properness + no base points mod 5.
2. `A = X0^3Y0^3` is independent of `delta` in `R_6` (mod 5 both are distinct Fermat basis monomials).
3. `h^{3,1}(X)^+ = 226 - 18 = 208`; residual complement of the `delta` line has dimension 207; `[A]` has nonzero
   projection to it. "This count alone gives no obstruction to divisor support."
4. `(phi_* O_{P^5})^+ ≅ O ⊕ 6O(-1) ⊕ 9O(-2)` via the Koszul/Hilbert-series argument of 1.5.

Conditional obstruction (proved under a hypothesis): let `S` be a smooth projective surface parametrizing lines on
`C` with generically finite parameter map to `Gr(2,6)`, `U = P_S(V)` the universal line, `h : U -> C`, and
`T = X ×_C U` the scheme-theoretic fibre product. Base change of item 4 gives `(π_* O_T)^+ = O_U ⊕ 6O_U(-1) ⊕ 9O_U(-2)`;
`H^2(U, O_U(-2)) = H^1(S, L^{-1}) = 0` by Kawamata vanishing (`L = det V^*` nef and big; Kodaira if finite), hence
`H^2(U,O_U) ≅ H^2(T,O_T)^+`. **If `T` is smooth or normal with rational singularities**, every invariant holomorphic
2-form on an equivariant resolution of `T` is pulled back from `U`, and push-pull `g_* π̃^* b = f^* h_* b` shows that
the invariant (3,1) Gysin image of the cylinder lies in `f^* H^{3,1}(C) = C·delta`. "Under the stated regularity
hypothesis, a pulled-back family of lines cannot serve this residual target." Explicitly left open: non-normal fibre
products, non-rational singularities (normalization/resolution changes `π_* O`), other families, correspondences
between sheets.

Links in the comment: the instance rationale
`https://github.com/hodge-computer/sextics/blob/main/data/quadratic-cubic-pullback/why.md`; Kawamata's vanishing
theorem `https://doi.org/10.1007/BF01456407`. No reproduction packet, no code.

Assessment. Items 1–4 are correct (recomputed: 226, 18, 208; `A`, `delta` independent mod `2^31-1`). The conditional
obstruction is sound as stated; its hypothesis is generically satisfied (for a generic surface of lines the cylinder
`T` is singular along a curve with transversal `A_1` type, which is rational), so line cylinders are effectively
closed for the residual sector. What the comment does not do: it does not address divisors of `X` that are not
pulled back from `C` at all, which — by 1.6 — is where any support must live.

### 3.2 The instance rationale and `checks.json` (project screening, `scripts/audit-structural-frontier.py`)

* why.md: the determinant density `delta` is the only direction served by the §3.4 map (`b` constant at `q = 1`);
  reparametrizing the quadratic system changes `delta` by a scalar; `A` is a different Jacobian basis element even
  mod 5; a rational inverse Jacobian does not help because its transformed source need not descend and trace keeps
  only the descending component; the map has more sheets than `tau` and is not a reflection-group quotient; suggested
  directions: ramification, the system of quadrics, correspondences between sheets.
* `checks.json`: `term_count 50`; smoothness by mod-5 reduction; target survives in `R` mod 5; `hessian_symmetrizer
  _dimension_mod101 = 1`; no reflections, no paired signs; the six quadrics; `known_density = det(dQ_i/dz_j) ≡ 64
  X0X1X2Y0Y1Y2 (mod 5)`; `independent_in_jacobian: true`; `q1_source_space_on_same_map: weight zero, constants`.
* Book (docs/04 §4.2, §4.7; docs/03 §3.4 "Six arbitrary globally Y-even quadrics composed with any cubic G illustrate
  the theorem: delta is an odd sextic and its class pulls back from a smooth cubic fourfold … No assertion that those
  densities exhaust all targets on F follows."). docs/04 §4.7 is the general statement that uniruled intermediate
  quotients serve only representations meeting their averaging ideal; `X/tau` has `kappa = 0` (fixed curves of age
  3/2, torsion `K`), so no quotient of `X` is uniruled except through `C`.

### 3.3 The local attack session (2026-11-01, this machine; engine `~/devel/sextics-work/hodge`, engine 2 only)

A sign bug in the `F·N(D)` term of the fixed-`D` system was fixed at ~17:25; all consistency/rank/served-space results
below are engine-2 recomputations (`log.jsonl`); engine-1 results are kept in `log_engine1_unreliable.jsonl` and are
not used, except that the E1 test (`e1_solvable`, independent code path) was correct in both engines.

Fixed-denominator tests on the actual `F` (all: E1 fails, served dimension 18 = nothing served):

| family | `D` | `d` | rank / nullity of the compressed system |
|---|---|---:|---|
| sanity | `X0^2 ± Y0^2` | 2 | 168 / 0 |
| cone3 | `X0^3+X1^3+X0Y0^2`, `X0^3+X1^3+X2^3`, `X0^3+X0Y0^2+X2^3`, `X0Y0^2+X1^3+X2^3`, `X0^3+X1Y0^2+X1^3`, `X0^3+2X1^3-X0Y0^2+3X1Y0^2+X0^2X1` | 3 | 378 / 0 |
| cone4 | `X0^4+X1^4+X2^4+Y0^4`, `X0^4+X1^4+Y0^4` | 4 | 747 / 9 |
| delta | `det(dQ)` (64 terms) | 6 | 2223 / 153 |
| quadric | `Q0`, `Q3`, `Q0±Q3`, `Q0+Q1`, `Q1+Q2`, `Q0+Q1+Q3`, `Q0+Q1+Q3+Q4` | 2 | 168 / 0 |
| qprod | `Q0Q3`, `Q0Q1`, `Q0(Q0+Q3)` | 4 | 747 / 9 |
| qcube (engine 1 only; E1 reliable) | `Q0^3+Q3^3` | 6 | E1 false |
| join012 (served space, engine 2, 690 s) | `Q0^3+Q1^3+Q2^3` | 6 | served dim **19**, `A` not served |
| join3, all ten splits (concurrent session, `listC`, 163–290 s each) | `Q_a^3+Q_b^3+Q_c^3` | 6 | E1 false for `A`, rank 2214 / nullity 162 (seven of nine remaining finished; served dims not computed) |
| target `:= delta` (concurrent session, `delta_served.py`) | `join012`, `join345` | 6 | E1 true, **consistent** (rank 2358, nullity 244): each join supports exactly the `delta` line |
| target `:= delta` | `Fpart(6) = F_{60}`, `X0^6+X1^6+X2^6` | 6 | E1 false, inconsistent (rank 2341 / 261): support nothing |
| polar (new, this dossier) | `F_0`, `F_0+2F_1+3F_2`, `F_3` | 5 | E1 true (automatic), served dim 18, `A` not served |
| structure (concurrent session) | `jointest.py`: rank of the 21 second partials = 21 (a join has `<= 12`), rank of first partials 6; `struct1.py`: Smith form `[1,1,1,1,2,0]`, bigraded parts `F_{60},F_{42},F_{24},F_{06}` with 10/15/15/10 terms, two-plane `(4,4)`-cancellation system nullity 0, smoothness of `F_{60}`, `F_{06}` as plane sextics, `F_{42}` Y-Hessian nonzero | — | — |

Comparison runs on Fermat `F_0` (for calibration; engine-1 served dimensions are unreliable and are *not* quoted):
`X0^6+X1^6+X2^6` at `d = 6`: consistent, exact witness verified, served dimension 118 (engine 2).

The first run of that session was killed by a machine-wide OOM at 17:44 (other agents' 13–42 GB jobs); the
concurrent session resumed it at ~18:14 with a reduced list `listC_ours.txt` (ten joins, `Fpart(6)`, `Fpart(0)`,
`Fpart(6)+Fpart(0)`, `Fpart(4)+Fpart(2)`, `X0^6+X1^6+X2^6`, `Y0^6+Y1^6+Y2^6`, block Euler images `EI([0,1,2])`,
`EI([3,4,5])`, `EI([0,3])`, `EI([0,1,2,3])`, `Q0^3+Q3^3`, `Q0^3+Q1^3`, `F - Q0^3`, `Q0^2`, `Q0^2+Q3^2`, `X0^6+Y0^6`,
`X0^4+Y0^4`, the `family-block-cube` denominator) of which the joins had finished at the time of writing, plus
`delta_served.py` (above) and the structural scripts. The cone lists `cones012.txt`, `cones345.txt` (Hesse-type and
random cubics in `Q_S`) were dropped after the adjoint criterion predicted `h^{2,0} = 0`; `bundle.py` (greedy
selection of denominators whose served spaces span `A`, joint modular solve + CRT) never had input. **All of these
are degree-6 polynomials not in `J` and are excluded as supports of `[A]` by 1.6(3) whether or not they finish**
(`EI(I)` of `F` is in `J` only when `I` is everything; `Fpart` sums are not in `J` unless equal to `F`).

E1 kernel sweeps: Y-even `d = 1..8` (`sweep.py`, `kernel_sweep.jsonl`; `d = 8` by this session, a duplicate `d = 8`
run by the concurrent session was also started), Y-odd `d = 1..6` (new). Results in 1.6.

Other sources of prior work: none. There is no DiamondJM-style sweep comment in this thread (the memory notes mention a
~1300-support sweep over the q=1 open lanes; it is not posted here and its coverage of this instance is unverified;
by 1.6 any degree-6 sweep would have found nothing).

### 3.4 Classification of what is settled

* Rigorous (exact linear algebra, Hilbert-function-certified mod-`p` ranks): 1.6 items 1–4; Fermat-limit witness;
  benediamond's items 1–4; the symmetry group (1.3).
* Screened-method failure, not a proof of nonexistence: served spaces (mod `p`, one prime) of the listed
  denominators; benediamond's cylinder obstruction (conditional on rational singularities).
* Conjectural / heuristic: `K_d = J_d` for `d = 8, 9` and `dim K_10 = 294 ⊕ 261` (running); monodromy of `X/tau -> C`
  equal to `S_16`; the identity of the single class served by `join012`.

---

## 4. What remains unexplored (concrete gaps)

1. **Polar (Jacobian-ideal) denominators of degree 6 and 7.** The only degree-`<= 7` supports not excluded. Degree 5
   is closed (three samples, all serve nothing; geometrically the ramification divisor of a projection has isolated
   singularities). Degree 6 (`sum ℓ_i F_i`, 18-dimensional Y-even family) is running (Appendix A); degree 7
   (`sum q_i F_i`, 63-dimensional) untested. Expected: nothing served (rational singularities), but it is the last
   low-degree door.
2. **E1 kernel at `d = 9, 10` (and the parity split at `d = 10`); Y-odd kernel at `d = 7, 8`.** `d = 8` (even) is closed
   (`K_8 = J_8`); `d = 9` running; `d = 10` not started. This fixes the minimal degree of any non-polar denominator
   and produces the candidate space itself.
3. **Structure of the `d = 10` kernel.** Nothing is known about whether `ann_R(A)_{10}` contains geometrically meaningful
   elements: quintics in the `Q_i` (pullbacks of degree-5 hypersurface sections of `C`), `delta · (Y-odd quartic)`,
   products of two non-Jacobian quintics, polynomials in `Q_S` for a triple `S` times something, elements of the form
   `H(Q) · F_i`-free combinations, etc. All of these are *linear* intersections once the kernel is known.
4. **No full fixed-`D` system has ever been run at `d >= 7` on this `F`** (the largest was `d = 6`).
5. **The joins are settled**: each is excluded as a support of `[A]` by 1.6, and the two checked (`join012`, `join345`)
   support exactly the `delta` line; the union/bundle over the ten splits is therefore dead. What is *not* known is
   which singular non-pullback divisors of degree `>= 10` could have `h^{2,0} > 0` in the odd sector; the adjoint
   criterion of 1.9 is the right tool to pre-screen singular configurations before any linear algebra.
6. **Monodromy / deck group of `phi` and of `X/tau -> C`.** Never computed. If `X/tau -> C` were imprimitive there would be
   an intermediate cover, possibly Fano, carrying part of `V_15`; if the deck group of `phi` were larger than `{±1, ±tau}`
   there would be an extra automorphism of `F` (impossible among monomial maps, not excluded among general linear
   maps, though the Hessian-symmetrizer screen makes it very unlikely).
7. **Deformation from Fermat.** Not attempted. The first-order obstruction to continuing the Fermat certificate
   `A E_X/(6 D_0)` along `F_t = sum (Z_i^2 + t m_i)^3` has not been computed; by 1.6 the degree-6 family cannot survive
   to `t = 5`, but the order at which it dies and the shape of the obstruction class are unknown.
8. **Deformation within the six-parameter family `c = (c_0..c_5)`.** No study of special members (symmetric point
   `c_i ≡ c`, half-deformed points `c_3=c_4=c_5=0`, points with a reflection) and of whether a certificate for `A`
   exists on a sub-family that could be continued.
9. **Non-pullback singular divisors.** No candidate divisor of `X` with non-rational singularities that is not a pullback
   from `C` has been proposed by anyone. Natural sources would be: unions of the 405 K3 surfaces (codimension 2, not
   divisors), cones over the K3 surfaces from points of `P^5` (not contained in `X`), degree-`>= 8` hypersurface sections
   of `X` singular along curves of positive genus.
10. **`tau`-anti-invariant / other sectors:** irrelevant for this target but a sanity structure never checked: the
    served spaces of `delta`-type denominators for *all* targets (what the `delta` line's certificate looks like
    explicitly; the §3.4 theorem is existential and no explicit `D` for `A = delta` is known).

---

## 5. Candidate attack plan

Ordered by (cost × information). Sizes refer to the engine's compressed modular solver.

**E1. Close the polar family (degrees 6, 7).** Idea: 1.6 says any degree-`<= 7` support is `V(sum ℓ_i F_i)`. Test served
spaces for random Y-even `sum ℓ_i F_i` (`d = 6`: 2602 columns, ~2 min each; 3–5 samples) and random `sum q_i F_i`
(`d = 7`: 4497 columns, ~10 min each; 3 samples). Expected outcome: served dimension 18 (polar divisors have isolated or
rational singularities: `V(F_0) ∩ X` is the ramification divisor of projection from `[1:0:0:0:0:0]`, singular only at
points where `e_0 - μ z ∈ ker Hess(F)`). Chance of a witness: `< 2%`. Value: turns "every support has degree `>= 8`" into a
screened statement; with the `d = 8, 9` sweeps, "every support has degree `>= 10`".

**E2. E1 kernel at `d = 10` with parity split, and its basis.** Matrix: Y-odd (resp. Y-even) monomials of degree 16 as
columns (`~10165` each), rows `3·S_11 = 13104` Jacobian rows plus 1512 (resp. 1491) rows `A·m`; rank mod `p`, then the
transpose nullspace restricted to the `A`-rows gives a basis of `K_10/J_10` (expected 294 / 261). Cost: ~1 h per
parity with the present pure-Python row builder; the rank itself is minutes in FLINT. Output: the complete list of
admissible denominators of the minimal possible degree. Chance that `K_10` alone yields a witness: nil (it is a
necessary condition only); value: prerequisite for E3–E5.

**E3. Structured intersections inside `K_10` (all linear, seconds each once E2 is done).** For each structured subspace
`W ⊂ S_10^{even}` compute `dim (W ∩ K_10)` and whether the intersection leaves `J`:
  * `W = C[Q]_5` (quintics in `Q_0..Q_5`, 252-dim; pullbacks of quintic sections of `C`);
  * `W = delta · S_4^{odd}` (60-dim; divisors containing the ramification divisor);
  * `W = {P(Q_S) : deg P = 5}` for the 20 triples `S` (21-dim each), and `{H(Q_S) H'(Q_{S'})}`;
  * `W = (Q_i) · K_8`-type products if `K_8 ≠ J_8` (it is expected to equal `J_8`, so skip);
  * `W = F_i · S_5` (polar-type, `⊂ J`, skip) and `W = {P_1 P_2 : P_1 ∈ fixed quintic}` for structured quintics `P_1`
    (e.g. `P_1 = F_0`, `P_1 = E_X(F)`, `P_1 = sum_{i<3} Q_i^2 ∂_{X0} Q_i`);
  * the `rho`/`sigma`-orbit-symmetrized subspaces (polynomials symmetric under the *shape* symmetries, in case the
    annihilator is "almost" symmetric).
Each nonzero intersection element not in `J` is a candidate for E4. Chance of finding a structured admissible `D`: maybe
30% (555 is a large kernel inside a 1512-dim space; random 252-dim subspaces meet it in dimension
`~ 252 + 555 - 1506 < 0`, so a nonzero intersection would itself be a structural signal). Chance that such a `D` then
serves `A`: low (`~5%`), since pullback-type supports should only serve the `delta` line by benediamond's argument.

**E4. Full fixed-`D` systems at `d = 10` for the E3 candidates and for 2–3 random elements of `K_10`.** Size 18900 unknowns
(13104 `N` + 5796 `Q`), 84903 equations; ~3–5 h per candidate with FLINT `rref` single-threaded, or far less with a sparse
eliminator (the matrix has `~50` nonzeros per column; structured Gaussian elimination / block Wiedemann mod `p` would
be minutes). The *served-space* variant (add the 226 target columns, compute the nullspace projection) costs the same and
answers "what does `V(D)` support" for every target at once; it is the right primitive. Random elements of `K_10` are
expected to serve nothing (their divisors are generic hypersurface sections, smooth by Bertini). Chance of a witness
from structured candidates: `5–10%` overall. This is the main computational experiment; an affirmative answer gives a
witness immediately (rational lift by CRT is already implemented, `solve_fixed_D_fast(want_exact=True)`).

**E5. Served space of the E1-admissible kernel "in bulk".** Instead of one `D` at a time, test whether `[A]` lies in the
*sum* of served spaces of a basis of `K_10/J_10` (a bundle certificate, `HasCertificate.add`, product denominator of
degree `10·k`). Cost: 294 served-space computations of E4 size — infeasible with the present code; feasible only if E4's
primitive is brought to minutes. Chance: low (`<5%`), but it is the only *exhaustive* statement available at degree 10:
if no element of `K_10` serves `[A]`, every support has degree `>= 11`.

**E6. First-order deformation from Fermat.** Along `F_t = sum (Z_i^2 + t m_i)^3`, write `D_t = D_0 + t D_1`, `N_t = N_0 + t N_1`,
`Q_t = t Q_1` with `D_0 = X0^6+X1^6+X2^6`, `N_0 = A E_X/6`. The order-`t` equation is linear in `(D_1, N_1, Q_1)`
(`2 A D_0 D_1 - D_1 N_0(F_0) - D_0 N_0(F_1) - D_0 N_1(F_0) + F_1 B_{D_0}(N_0) + F_0 B_{D_1}(N_0) + F_0 B_{D_0}(N_1) = F_0^2 Q_1`,
`F_1 = 3 sum Z_i^4 m_i`), size `2602 + 236` columns at `d = 6`: seconds. Continue order by order until obstructed; the
obstruction class lives in `R_6^{odd}` modulo the Fermat served space of `D_0` (100-dimensional). Expected: obstructed at
a low order (by 1.6 there is no degree-6 denominator at generic `t`). Value: tells whether the degree must grow already
at first order (then the deformation must be done in a bigger degree, `d = 10` or `12`, where the Fermat system itself
has huge nullity and the continuation becomes a Newton iteration in `~19000` unknowns — possible but a project). Chance
of leading to a witness: `~5%`; chance of structural insight: high.

**E7. Monodromy of `phi` and of `X/tau -> C`.** For random `u ∈ F_p^6` factor the zero-dimensional ideal `(Q_i(z) - u_i)` over
`F_p` (Singular/M2 in the `ddg-cas` Docker image, or FLINT univariate elimination via a random linear form): the residue
degrees give Frobenius cycle types on the 64 affine / 32 projective sheets. Fifty samples distinguish `(Z/2) wr S_16`
from imprimitive subgroups (block systems other than the `tau` pairs). Cost: hours of setup, minutes of compute.
Outcome if imprimitive: an intermediate cover `X/tau -> C' -> C` and a new quotient to test for uniruledness (a Fano `C'`
would open docs/03 §3.2 for the part of `V_15` it carries). Chance of imprimitivity: `~10%` (the mod-5 limit is the
Kummer cover with group `(Z/2)^5`, which *is* imprimitive; the deformation generically destroys it). Decisive either way.

**E8. Degree-9 polar-free candidates only if the running sweep shows `K_9 ≠ J_9`** (`K_8 = J_8` is now known). Otherwise
skip. (If `A·: R_9 -> R_15` failed to be injective, the kernel would be a very special structural fact — the first
thing to examine would be whether `delta·R_3` annihilates `A`.)

**E8b. Adjoint pre-screen of singular configurations at `d = 10` (uses 1.9; cheap).** Before any fixed-`D` system,
enumerate singular loci that a degree-10 hypersurface section of `X` can have with non-rational transversal type
(curves on the 405 K3 surfaces, intersections `S_Π ∩ S_Π'`, plane sections, the curves `C_S = {Q_i = 0, i ∈ S} ∩ X`),
compute the adjoint-condition bound in the Y-odd sector with `adjcount.py`-style counts, and keep only
configurations with positive bound. Each such configuration is a *linear* condition on `D` (`D ∈ I_Γ^{m}` restricted to
`X`), to be intersected with `K_10` (E3). This replaces blind sampling of the 294-dimensional kernel by a short list of
geometrically motivated subspaces. Chance that a positive-bound configuration exists at `d = 10`: unknown; the
two-plane analysis of 1.9 says the coordinate planes never work, but non-planar curves were not analysed.

**E9. Nonconstructive routes (assessment rather than experiment).** The accepted nonconstructive proofs all run through a
uniruled (Fano or elliptic/cubic-Jacobian) quotient with a polynomial curve family, integration and trace. Here the
candidate quotients are `C` (Fano, wrong sector: `f_*` kills the residual part) and `X/tau` (`kappa = 0`). docs/04 §4.7's
averaging-ideal criterion shows no intransitive quotient of the Galois closure of `X/tau -> C` can be uniruled (they all
dominate `X/tau` or carry no `V_15`-invariants). So no nonconstructive mechanism in the library applies unless E7 finds
new structure. A Hodge-theoretic support proof (e.g. via a Kuga–Satake-type or Shioda-type argument on the finite
cover) would not be formalizable in the present Lean library. Chance: `<2%` without new theory.

**E10. Degree 12 structured candidates (only after E4).** Products `H(Q_S)·H'(Q_{S'})`, `delta^2`, `delta·join`, products of
two polars — all E1-screenable linearly at `d = 12` (`25704 × 16797` Jacobian block, ~1 h) and then full systems of
42500 unknowns (needs the sparse eliminator). Chance similar to E4; cost 4× larger.

Not recommended: any further degree-`<= 7` sweep (closed by 1.6 except polars); line cylinders (benediamond); quadric
sections at any pole order (1.6 + 1.7); the nine remaining joins (excluded as supports of `[A]`).

---

## 6. Feasibility assessment

**Explicit witness.** Any witness has `deg D >= 8`, and `>= 10` unless the running sweeps surprise (or the denominator is a
polar, which E1 is closing). At `d = 10` a witness is `~19000` rational numbers; the accepted `m3-coupled` witness is
exactly this size (`d = 10`, `Witness.lean` 360 KB, accepted by the protected verifier), so the Lean side is known to
work at `d = 10`; `d = 12` (`~42500` coefficients, identity of degree 30 with 162248 Y-odd monomials) would be roughly
3–4× larger and has no precedent but no known hard limit either (`klein-rational` has `d = 30` with a sparse witness).
The bottleneck is the search, not the verification: there is no structured family of denominators with any
positive evidence, and the only exhaustive statement available at degree 10 (E5) needs a faster linear-algebra
primitive. With a sparse mod-`p` eliminator (or Magma/M4RI-class code on the 128-core box), E2–E5 are a few days of
compute; with the present FLINT-only code, E2–E4 on a handful of structured candidates is a weekend.

**Nonconstructive argument.** Nothing in `formal/lean/HodgeComputer/Constructions/*` applies: `HoldoutCover`,
`PolynomialCurves`, `EllipticCompanions`, `ResidueTransport`, `TransitiveCubic`, `HigherResidues` all need a quotient
with a polynomial curve family or an elliptic/cubic-Jacobian structure in the right character; `AveragingProjectors`
only transports support between quotients. The §3.4 weighted-Fano theorem (the one existing structure here) serves
exactly the `delta` line. A new mechanism would be needed, e.g. a proof that some explicit non-pullback divisor of `X`
with non-rational singularities carries `[A]`, converted into a polynomial-family integration on that divisor; no
such divisor is known.

**What a Lean proof would require, concretely.** For an explicit witness: `witness.json` with `D` (Y-even, degree `d`,
coprime to `F` — coprimality over `C` must be proved, e.g. by exhibiting a monomial order / specialization argument as
in the accepted `Witness.lean` files), `N_0..N_5` of degree `d+1` with the parities of 1.2, `Q` of degree `2d-6`, and the
exact identity `E = F^2 Q` checked by the generated `Witness.lean` (`hodge/leangen.py` + `leancheck.py` +
`sandboxcheck.py` replicate the verifier; a regenerated accepted witness passed the official sandbox in ~5 min). For a
nonconstructive proof: a `HasCertificate` existence theorem built from library primitives — currently unavailable for
this geometry.

**Bottom line.** Probability of a solve from the plan above, honestly: `10–15%` for an explicit witness (concentrated in
E3/E4 finding a structured degree-10 denominator), `<2%` for a nonconstructive proof. The highest-value next steps are
cheap: finish E1 (polars), E2 (`d = 10` kernel), E3 (structured intersections), and E7 (monodromy), in that order; only
then commit to the heavy E4/E5 computations.

---

## Appendix A. Reproduction and status of running jobs

* Engine: `~/devel/sextics-work/hodge/{core,fixed_d}.py` (engine 2; the `F·N(D)` sign fix of 2026-11-01 17:25 is in).
  Work module `/tmp/hodge/results/quadratic-cubic-pullback/qcp.py` (quadrics, `det`, logging), `sweep.py` (E1 kernels),
  `serves.py` (served spaces, union, `A` membership), `bundle.py` (greedy bundle + joint modular solve + CRT; never
  exercised), `run_list2.py` (list runner).
* Logs: `log.jsonl` (engine 2; each line `{D, d, family, e1, consistent, rank, nullity, seconds, note, engine, imgdim?,
  A_served?}`), `kernel_sweep.jsonl`, `out_*.txt`, `delta.txt`, `witness-candidate-fermat-92047.json` (Fermat only).
* New scripts for this dossier: `/tmp/hodge/scratch/qcp_dossier_checks.py` (F = ΣQ³ check; counts 226/18/208 and
  236/18/218; `A`, `delta` ∉ `J_6` and independent; SNF `(1,1,1,1,2)`; monomial symmetries; system sizes; Y-odd E1 sweep
  `d = 1..6`), `/tmp/hodge/scratch/qcp_polar.py` (cubics-in-`Q` ∩ `J_6^{even}` = 1-dimensional = `C·F`; polar `d = 5`
  served spaces: `F_0`, `F_0+2F_1+3F_2`, `F_3` all 18), `/tmp/hodge/scratch/qcp_polar6.py` (random Y-even element of
  `J_6` and `X0F_0 - Y0F_3`; output `out_polar6.txt`).
* Running at the time of writing (18:45): `python3 sweep.py 8 9 ours` (pid 3931550, started 18:20; `d = 8` finished
  at 18:35 with `dimK_mod_J = 0`, 912 s; `d = 9` in progress, output `sweep_8_9.out`, appends to `kernel_sweep.jsonl`);
  `qcp_polar6.py` (pid 3976314, started 18:31, output `out_polar6.txt`, appends to `log.jsonl`; no result yet — the
  machine is running four solver processes). Concurrent session's jobs: `run_list2.py ours listC_ours.txt` (pid
  3921344, eight of nine remaining joins done), a duplicate `sweep.py 8 8 ours` (pid 3933408). If `kernel_sweep.jsonl`
  shows `dimK_mod_J = 0` at `d = 9`, statement 1.6(2) extends to Y-even denominators of degree `<= 9` and, by the
  Gorenstein count, the first non-polar denominators have degree exactly 10 with a kernel of dimension
  `>= 555` (expected 294 Y-even). If `out_polar6.txt` reports served dimension 18 for both degree-6 polars, the
  degree-6 door is closed up to the usual sampling caveat (two samples of an 18-dimensional family).
* To rerun a served-space test: `python3 serves.py <listfile> <tag>` with lines `name | <expression in P, Q, core, delta,
  F, Fpart, EI>`; to run the E1 sweep at `d`: `python3 sweep.py d d ours` (writes `kernel_basis_d<d>_ours.txt` when the
  kernel exceeds `J_d`).

## Appendix B. Key formulas collected

* `F = sum_{i=0}^5 Q_i^3`, `Q_i = Z_i^2 + (5/p_i) m_i`, `p = (7,11,13,17,19,23)`,
  `m = (X1X2+Y0Y1, X0X2+Y1Y2, X0X1+Y0Y2, X0X1+Y1Y2, X1X2+Y0Y2, X0X2+Y0Y1)`.
* `F_i = ∂_i F = 3 sum_j Q_j^2 ∂_i Q_j`; `F = (1/6) sum Z_i F_i ∈ J_6`.
* `delta = det(∂Q_i/∂Z_j)` (degree 6, Y-odd, 64 terms); `phi^*(Ω_U/G^2) = delta Ω_Z/F^2`; the served sector of §3.4 at
  `q = 1` is `C·delta`.
* `h^{3,1}(X)^+ = 226 - 18 = 208 = 1 + 207`; `h^{3,1}(X) = 426`; `(phi_*O)^+ = O ⊕ 6O(-1) ⊕ 9O(-2)`;
  `f_*Q_X = Q ⊕ V_15 ⊕ W_16` (monodromy inside `(Z/2) wr S_16`).
* E1: `A D ∈ J_{d+6}`; `(J:A)_d = J_d` for `d <= 7` (even) and `d <= 6` (odd); Gorenstein bound
  `dim ann_R(A)_d >= dim R_d - dim R_{d+6}` = `0 (d<=9), 555 (d=10), 1020 (11), 1325 (12)`.
* Fermat witness: `v = X0^3Y0^3 (X0∂_0 + X1∂_1 + X2∂_2) / (6(X0^6+X1^6+X2^6))`, served classes `m_1(X)m_2(Y)`, 100 of them.
* Certificate degrees at `q = 1`: `deg N_i = d + 1`, `deg Q = 2d - 6`, identity in degree `2d + 6`; combination:
  `D_{A+B} = DE`, `N_{A+B} = EN + DM`, `Q_{A+B} = E^2 Q + D^2 R`.
