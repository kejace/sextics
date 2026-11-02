# Probe: `polarized-k3-remainder` — automorphisms and Picard number of the K3 double plane `w^2 = f`

Date 2026-11-01. Work directory `/tmp/hodge/work/probe-k3/` (scripts, raw outputs). Nothing in the repo was
modified; nothing was posted. Session was resumed after a crash of a previous attempt (which had only
written the scripts); all results below were recomputed in this session.

Conventions. `f(z) = z0^6 + (5/7) z0^5 z1 + (5/7) z0^2 z1^3 z2 + (15/7) z0 z2^5 + z1^6 + (10/7) z1^5 z2 + z2^6`,
`C = C_f = {f = 0} ⊂ P^2` (smooth plane sextic, genus 10), `S = {w^2 = f} ⊂ P(1,1,1,3)` (degree-2 K3),
`F = (f(X+Y) + f(X−Y))/2`, `s(q) = Σ_{z ∈ P^2(F_q)} χ_q(f(z))`, so `#S(F_q) = q^2 + q + 1 + s(q)` and
`Tr(Frob_q | H^2_prim(S)) = s(q)`, `Tr(Frob_q | H^2(S)) = q + s(q)`.

## 0. Summary (proven facts first)

**Proven in this session**

1. `F = (f(X+Y)+f(X−Y))/2` exactly (sympy, `verify_F.py`); `F` has 42 monomials.
2. The Lie algebra of infinitesimal linear automorphisms of `f` (`A ∈ gl_3` with `(Az)·∇f = c f`) is
   one-dimensional, spanned by the Euler field (`A = (c/6) Id`). So the group `{M : f(Mz) = f(z)}` is finite.
3. The 0-dimensional scheme `G̃ = {M ∈ Mat_3 : f(Mz) = f(z)}` (28 coefficient equations, 9 unknowns) has
   **length exactly 6 over F̄_65521**, and its points are exactly the scalar matrices `ζ·Id`, `ζ^6 = 1`
   (msolve F4 + rational parametrization, 9 min; eliminating polynomial `m33^6 − 1`, all off-diagonal
   entries `0`, `m11 = m22 = m33`). Since `G̃/μ_6 = Aut(C)` for a smooth plane curve of degree ≥ 4,
   **`Aut(C_f ⊗ F̄_65521) = 1`**. `C` has good reduction at 65521 (checked).
4. By injectivity of specialization of automorphisms for smooth curves of genus ≥ 2 (Deligne–Mumford,
   *Irreducibility of the space of curves*, Thm 1.11: the automorphism scheme of a stable curve over a DVR is
   finite and unramified; hence `Aut(C_{Q̄}) ↪ Aut(C_{F̄_p})` for every prime of good reduction),
   **`Aut(C_f ⊗ Q̄) = 1`: the plane sextic has no nontrivial projective automorphism at all** — in
   particular none of order 2 or 3. Consequently `f` is **not** projectively equivalent to
   `g_6(z0,z1) + z2^3 g_3(z0,z1) + a z2^6` (such a form has the order-3 homology `z2 ↦ ζ_3 z2`), and the only
   automorphism of `S` preserving the degree-2 polarization is the covering involution `ι`.
5. Point counts of `S` over `F_{3^n}`, `n = 1..11` (Zech-log C code with Frobenius-orbit reduction, validated
   against full enumeration and an independent Python count): `s(3^n) = 3, −3, 108, 153, 648, 1026, 8991,
   3321, 26244, 120042, 345546`. Newton identities + functional equation (the sign `ε = −1` is excluded
   because it produces roots with `|α| ≠ 3`, contradicting Deligne) give the characteristic polynomial of
   `Frob_3` on `H^2(S_{F̄_3}, Q_ℓ)`:

       χ(T) = (T − 3)^4 · Q_18(T),
       Q_18(T) = T^18 + 6T^17 + 33T^16 + 117T^15 + 405T^14 + 1134T^13 + 3402T^12 + 8748T^11 + 26244T^10
               + 72171T^9 + 236196T^8 + 708588T^7 + 2480058T^6 + 7440174T^5 + 23914845T^4 + 62178597T^3
               + 157837977T^2 + 258280326T + 387420489,

   `Q_18` irreducible over `Q`; normalized `Q_18(3z)/3^18 = (3z^18 + 6z^17 + 11z^16 + 13z^15 + 15z^14 + 14z^13
   + 14z^12 + 12z^11 + 12z^10 + 11z^9 + 12z^8 + 12z^7 + 14z^6 + 14z^5 + 15z^4 + 13z^3 + 11z^2 + 6z + 3)/3`,
   which is not a cyclotomic polynomial (non-integral coefficients), so `Q_18` has no root of the form
   `3·(root of unity)`. This **independently reproduces** the thread's `(z−1)^4 Q_18(z)/3`.
   Consequences:
   * Tate's inequality (unconditional): `ρ(S_{F̄_3}) ≤ 4`, hence **`ρ(S_C) ≤ 4`** (specialization is injective).
   * Tate conjecture for K3 surfaces in odd characteristic (Madapusi Pera): `ρ(S_{F̄_3}) = 4`, and since all
     four Tate eigenvalues are exactly `+3`, every geometric divisor class is already `F_3`-rational.
   * Artin–Tate (Milne; valid since Tate holds): `|Br(S_{F_3})| · |disc NS(S_{F_3})| = 3 · Q̃(1) = 211`.
     `211` is prime and `|Br|` is a square (Liu–Lorenzini–Raynaud), so **`|Br(S_{F_3})| = 1` and
     `NS(S_{F̄_3})` is a rank-4 lattice of discriminant `−211`** (signature `(1,3)`).
     It contains `⟨H, L_1⟩` (Gram `[[2,1],[1,−2]]`, disc `−5`) from the splitting line `a + b = 0`:
     verified `f|_{z0 = −z1} ≡ z2^6 (mod 3)`.
   * Predicted `s(3^12) = −1048302` and `s(3^13) = 8844228` from `χ`; the `3^12` count was still running when
     this report was written (see §5) — it is a genuine independent check, not an input.
6. Good primes `p ≤ 200` for the model `7f` (reduction of `C` smooth): all odd primes except 7 (the
   denominator prime; `7f mod 7` is a different form), plus `2` is bad (and char 2 is bad for `w^2 = f` anyway).
   `C` also has good reduction at 65521 and 1073741827.

**Relied on from the discussion thread (#115, `1337isnot1337`), not re-verified here**: the mod-9 lifting
obstruction for the splitting line (Elsenhans–Jahnel) giving `ρ(S_C) ≤ 3`; and `End_Hdg(T) = Q` (Frobenius
Galois-group argument at 3, 23, 61 + Charles/Zarhin). With those, `Aut(S)_tors = {1, ι}` (see §3).

**Heuristic / expected**: `ρ(S_C) = 1`; no elliptic fibration; no CM/RM; the Kuga–Satake route of the book
does not apply by any criterion the book states.

## 1. Automorphism group of `C_f` (part 1 of the task)

### 1.1 Method
* Lie algebra: solve the linear system `Σ_{ij} A_ij z_j ∂_i f = c f` (21 nonzero coefficient equations in 10
  unknowns): rank 9, nullity 1 — only scalars. (`verify_F.py`.)
* Finite part: the scheme `G̃ = V(coefficients of f(Mz) − f(z)) ⊂ A^9` over `F_p`. Any `M` with
  `f(Mz) = λ f(z)` can be rescaled to `λ = 1`, and `G̃ → PGL_3` has kernel `μ_6` and image `Aut(C)`
  (automorphisms of smooth plane curves of degree ≥ 4 are linear). `det M ≠ 0` is automatic (`f` smooth).
  So `|Aut(C_{F̄_p})| = (number of geometric points of G̃)/6 ≤ length(G̃)/6`.
* msolve (`aut_65521.ms`, generated by `gen_msolve.py` from the integral form `7f`): "Dimension of quotient: 6";
  rational parametrization = the six scalar matrices. Output `aut_65521.out`, log `aut_65521.log`.
  A second prime (`p = 1073741827`) and a char-0 Gröbner basis in Singular (`autQ.sing`) were launched as
  extra confirmation; they had not finished when the report was written (§5). The single-prime result is
  already a proof by the specialization argument of §0.4, because length 6 is the minimum possible (the six
  scalars are always there), so there is no room for a "lost" point.

### 1.2 Consequences for `S`
* An automorphism of `P^2` preserving `C` lifts to two automorphisms of `S` differing by `ι`; the
  Bolognesi–Laterveer hypothesis (non-symplectic order 3 with a fixed curve, in the book's normal form
  `g_6 + z2^3 g_3 + a z2^6`) is exactly an order-3 homology in `Aut(C)`. **It fails for this `f`.**
* The `why.md` remark "its displayed equation is not in the special coordinate order-three normal form" is
  now sharpened to: **no coordinate change brings it into that form**.
* Hidden symmetries of the plane model: none (not even an involution). The only symmetry of `S` of this
  type is the covering involution.

## 2. Picard number of `S` (part 2)

### 2.1 Frobenius at `p = 3` — see §0.5. Data table (`counts_p3.txt`):

| n | q=3^n | s(q) | #S(F_q) | Tr(Frob\|H^2) |
|---|---|---|---|---|
| 1 | 3 | 3 | 16 | 6 |
| 2 | 9 | −3 | 88 | 6 |
| 3 | 27 | 108 | 865 | 135 |
| 4 | 81 | 153 | 6796 | 234 |
| 5 | 243 | 648 | 59941 | 891 |
| 6 | 729 | 1026 | 533197 | 1755 |
| 7 | 2187 | 8991 | 4794148 | 11178 |
| 8 | 6561 | 3321 | 43056604 | 9882 |
| 9 | 19683 | 26244 | 387466417 | 45927 |
| 10 | 59049 | 120042 | 3486963493 | 179091 |
| 11 | 177147 | 345546 | 31381582303 | 522693 |

Reconstruction script `frob.py` (Newton identities for `e_1..e_11`, `e_{22−k} = ε 3^{22−2k} e_k`).
With `ε = +1` all 22 roots have `|α| = 3` (numerically, to 1e-6) and the polynomial factors as above;
with `ε = −1` the root moduli range over `0.74·3 … 1.36·3`, impossible.

### 2.2 Rigorous conclusions
* `ρ(S_C) ≤ ρ(S_{F̄_3}) ≤ 4` unconditionally (Tate inequality + injective specialization).
* `ρ(S_{F̄_3}) = 4`, `disc NS(S_{F̄_3}) = −211`, `Br(S_{F_3}) = 0` (Tate + Artin–Tate, both theorems for K3
  in odd characteristic).
* The thread's lift obstruction (not re-verified here) lowers this to `ρ(S_C) ≤ 3`. Note that if
  `ρ(S_C) = 4` held, `NS(S_{Q̄})` would be a finite-index sublattice of a lattice of prime discriminant 211,
  so `NS(S_{Q̄})` would have discriminant `−211·(square)`; the obstruction argument is what excludes this.
* The van Luijk two-prime trick cannot be run here by naive counting: the second prime `p = 5` would need
  `s(5^11)` (`≈ 10^14` field operations). Only partial traces at `5` and at `p ≤ 200` (`n = 1, 2`) were
  scheduled (§5); they give no further rigorous bound on `ρ`, only ordinariness tests (`p ∤ s(p)` ⟺ ordinary
  reduction; ordinary ⟹ not supersingular, which is already implied by the bound).

### 2.3 What would decide `ρ` exactly
Certified period computation (Lairez–Sertöz) of `w^2 = f` and LLL on the 22-dimensional period vector;
or a p-adic point-count (Kedlaya/controlled reduction) at a second prime to run van Luijk. Neither tool is
available in this environment. Expectation (heuristic): `ρ(S_C) = 1`.

## 3. Special structure (part 3)

* **Automorphisms**: `Aut(C) = 1` (proven). Given the thread's `ρ ≤ 3` and `End_Hdg(T) = Q`:
  any finite-order automorphism `g ≠ 1` of `S` acts on `T` by a Hodge endomorphism of finite order, so by
  `±1`; if `+1` it is symplectic, forcing `ρ ≥ 9` (Nikulin) — impossible; if `−1` then `gι` is symplectic
  hence trivial, so `g = ι`. Thus **`Aut(S)_tors = {1, ι}`** (conditional on the thread's two inputs; `Aut(C)
  = 1` alone already gives this for polarization-preserving automorphisms, unconditionally).
* **CM / RM**: excluded by `End_Hdg(T) = Q` (thread). Consistent with the data: `Q_18` is irreducible with
  Galois group `(C_2)^9 ⋊ S_9` per the thread; nothing in the present counts suggests otherwise.
* **Elliptic fibrations**: a degree-2 K3 has an elliptic fibration iff `NS` has an isotropic vector. With
  `ρ = 1` (expected) there is none; the thread excludes the low-degree ones (`E·H ≤ 3`, splitting lines/conics)
  over `C`. Nothing in the equation (`f` is not of the form `q^3 + c g^2`, not a sum of three cubic squares
  per the thread, Waring rank ≥ 10) suggests one.
* **Mod-3 geometry**: rank-4 Néron–Severi lattice of discriminant `−211` containing `⟨H, L_1⟩` from the
  splitting line `z0 + z1 = 0` (contact order 6 at a single point `(1:−1:0)`), `Br = 0`.
* **Flexes**: the flex scheme `(f, Hess f)` has length 72; the radical computation over `Q` (distinct
  flexes / hyperflexes) was still running when the report was written (§5) — informational only.

### Kuga–Satake route of the book (§3A.6)
The book's exterior-power theorem needs an *algebraic* Kuga–Satake correspondence for `S`. The criteria it
lists (and the standard ones): Kummer/Shioda–Inose (`ρ ≥ 17`), Paranjape six-line double planes (`ρ ≥ 16`),
non-symplectic order-3 automorphism (Bolognesi–Laterveer), RM/CM-type arguments. **None applies**:
`ρ ≤ 4` (proven here), `≤ 3` (thread), `Aut(C) = 1` (proven here), `End_Hdg(T) = Q` (thread). So the book's
route is closed for this `f` unless a criterion covering a *general* degree-2 K3 (no special automorphism,
`End_Hdg = Q`, `ρ ≤ 3`) is available. I am aware of recent work on algebraicity of Kuga–Satake
correspondences (Markman, O'Grady, Voisin's "Footnotes", Varesco) for hyperkähler/abelian-type situations;
whether any published result covers *all* degree-2 K3 surfaces I could not verify offline — **this is the one
literature question that would change the status of the instance** and should be checked before further
search. As things stand (and as the attack-1 note `attack1-polarized-k3-remainder.md` shows), the instance is
equivalent to a coniveau-1 statement for `Λ^2 T` of a K3 with `End_Hdg(T) = Q` and `ρ ≤ 3`, with no
applicable theorem.

**Obstruction, stated plainly**: the K3 `w^2 = f` has no symmetry, no (known) extra algebraic cycles and
a trivial Hodge endomorphism algebra; every known construction of the exterior-square support uses one of
those. A certificate for `polarized-k3-remainder` would be a new case of the generalized Hodge conjecture.

## 4. Reproduction
```
cd /tmp/hodge/work/probe-k3
python3 -I verify_F.py                 # F expansion + Lie algebra
Singular -q smooth.sing                # good primes <= 200
python3 -I gen_msolve.py 65521 && msolve -v 1 -t 4 -f aut_65521.ms -o aut_65521.out   # ~9 min
gcc -O3 -march=native -fopenmp -o count count.c && ./count 3 11 0                       # s(3^11)
python3 -I frob.py 3 counts_p3.txt     # chi(T), factorization, Artin–Tate
```

## 5. Still running when the report was written (results to be appended by whoever resumes)
* `count 3 12 0` → test `s(3^12) = −1048302` (`counts_p3.txt`).
* `counts_p5.txt` (`n ≤ 8`), `counts_small.txt` (`p ≤ 200`, `n = 1,2`; `n = 3` for `p ≤ 31`; `n = 4` for `p = 11, 13`).
* msolve at `p = 1073741827` (`aut_1073741827.log/.out`), char-0 Gröbner basis (`autQ.out`), flex radical
  (`flex.sing`, task output file).
None of these affects the proven statements above; they are confirmations and extra data.
