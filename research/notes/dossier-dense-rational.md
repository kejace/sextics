# Research dossier: `dense-rational` (OPEN)

Status date: 2026-11-01 (resumed after a crash of two background modular-rank jobs; their logs are referenced in 3.8 and Appendix C). Prepared from the repository at `/home/kejace/devel/sextics` (HEAD `83fbf5c`, read-only), the saved discussion thread `/tmp/hodge/discussions/dense-rational.md` (GitHub discussion #59, one comment, no replies), the full reproduction packet attached to that comment (downloaded to `/tmp/hodge/downloads/dense-packet/`, 167 files), the book chapters in `docs/`, and new independent computations in `/tmp/hodge/scratch/` (scripts and logs listed in Appendix C). Nothing in the repository was modified; nothing was posted.

Legend. **[proved]** = rigorous over C (exact kernel caps + nonzero modular minors, or symbolic identity). **[scoped]** = rigorous exclusion of an explicitly bounded family only. **[screen]** = modular or numerical evidence, not a proof. **[obs]** = my own remark. **[new]** = computed for this dossier, not in the packet or thread.

---

## 0. Executive summary

* Target: the 236-term "dense" even sextic, `q = 1`, `A = X0^3 Y0^3`. No certificate exists in the record; the only prior work is one AI-assisted report by `fqlx` (2026-10-20) with a 7.5 MB packet.
* The instance has **no structure that any accepted method uses**: projectively the only linear symmetry found is the global involution `tau`, there is no block splitting in any linear coordinates (Hessian centre is scalar), Waring rank is >= 56 (not a hidden Fermat or short-power form), and F is not even in any single coordinate (no reflection quotient). `why.md` says as much: "no singled-out contact configuration or favorable quotient".
* **Proved facts** (packet, independently re-verified in part here): F is smooth (derivative dominance, bound 2212/9973 < 6); every legal certificate has `deg D >= 7` (all parities, arbitrary quotient, complex coefficients); all degree-7 and degree-8 denominators that satisfy the first certificate equation are Jacobian-type `D = grad(G).w`; a nonlinear polynomial Piola atom with the exact density `A` cannot exist; the two-block-radial family `f E_X + g E_Y` fails at every degree; in the `tau`-quotient chart, base-only denominators and a simple `H_s` pole fail with arbitrary rational base coefficients.
* **New here**: the Jacobian-type classification extends to **degree 9** (multiplication by `A` is injective on `(R/J)_9` in both parity sectors, mod 101 and mod 10007); at **degree 10** multiplication by `A` is surjective, so the E1-admissible non-Jacobian denominators form spaces of exactly 294 (Y-even) and 261 (Y-odd) dimensions; and the two natural degree-6 "cone" denominators `F_{6,0}(X)`, `F_{0,6}(Y)` serve **no** primitive target whatsoever (not just not `A`). See Section 3.8.
* **Nothing has been excluded at degree >= 7 in general.** The degree-7 necessary condition (951 quadrics in 120 unknowns) was only probed numerically (31 starts) and by two tiny msolve runs (3.5 GB, <= 180 s). Degrees 8 and above were never searched.
* **Assessment.** For a Mumford-Tate-generic even sextic the invariant primitive Hodge structure (Hodge numbers 208/866/208) is irreducible, so any supporting divisor must carry >= 208 holomorphic two-forms on its resolution and supports *every* target at once (Chapter 4.4). No construction for a Zariski-open set of even sextics is known (Diamond's paper states that its methods reach about 17 projective dimensions of the 235-dimensional family of invariant sextics, and explicitly that a nonempty Zariski-open subset is out of reach). A constructive solve therefore needs either a new general mechanism or an unexpectedly structured divisor; a nonconstructive Lean proof has nothing to formalize yet. Realistic contributions are rigorous degree-bound extensions, structured-family exclusions, and cheap wide searches that have never been run (Section 5).

---

## 1. Exact target

### 1.1 The record

* `data/dense-rational/instance.rec`: id `dense-rational`, title "Dense rational · full coefficient support", `q = 1`, `A = X0**3*Y0**3`. `provenance.json`: `sameVarietyInstances: []`, discussion <https://github.com/hodge-computer/sextics/discussions/59>. No `checks.json`.
* Variables `Z = (X0, X1, X2, Y0, Y1, Y2)`; `tau(X, Y) = (X, -Y)`.
* `F` has all 236 degree-6 monomials that are even under `tau`. Equivalently (since total degree is 6) every monomial has even X-degree and even Y-degree. Split by Y-degree:

  ```
  F = F_{6,0}(X) + F_{4,2}(X,Y) + F_{2,4}(X,Y) + F_{0,6}(Y),   term counts 28 + 90 + 90 + 28 = 236.
  ```
* Sparse presentation **[proved, re-verified here]**: with `alpha` the exponent vector,

  ```
  F = sum_i Z_i^6 + (1/9973) * sum_{alpha not pure} c_alpha Z^alpha,
  c_alpha = 1 + ( sum_{i=0}^{5} (i+3)(alpha_i+1)^2  mod 17 )  in {1,...,17}.
  ```
  All 230 off-diagonal coefficients match this recipe (`/tmp/hodge/scratch/check_dense.py`). 9973 is prime. Histogram of numerators `c_alpha`: 1:15, 2:15, 3:10, 4:8, 5:11, 6:13, 7:13, 8:15, 9:15, 10:13, 11:12, 12:18, 13:10, 14:16, 15:17, 16:14, 17:15. The packet's "Fourier" rewriting: with `w = (3,4,5,6,7,8)`, `S(alpha) = sum w_i (alpha_i+1)^2 mod 17`, `h(r) = r+1`, and `zeta` a primitive 17th root of unity, `17 h(r) = sum_k B_k zeta^{k r}`, `B_k = sum_s (s+1) zeta^{-ks}`; the integral sextic is `G = 9973 F = C + 9965 X0^6 + 9968 X1^6 + 9971 X2^6 + 9957 Y0^6 + 9960 Y1^6 + 9963 Y2^6` where `C` is the generating-function expression and the uncorrected pure coefficients are `8,5,2,16,13,10`. None of the 16 nontrivial phase operators is an algebra homomorphism (it sends `z_i^2` to `zeta^{8kw_i} z_i^2`, not the square of `zeta^{3kw_i} z_i`), so the recipe yields no covering map or hidden order-17 symmetry **[proved, packet]**.
* Integer form used throughout the packet and here: `G = 9973 F`, `A0 = 9973 A`.

### 1.2 Degrees, parities, characters

* Certificate for `(F, 1, A)`: homogeneous `D != 0`, `gcd(D, F) = 1` over C, `N_i` of degree `d + 1` (`d = deg D`), polynomial `Q`, with the cleared identity

  ```
  E := A D^2 - D sum_i N_i dF/dZ_i + F ( D sum_i dN_i/dZ_i - sum_i N_i dD/dZ_i ) = F^2 Q.
  deg E = 2d + 6,  deg Q = 2d - 6  (Q = 0 forced for d <= 3: such certificates are strict).
  ```
* Euler-normalized form (`N -> N + (R/6) Z` keeps `D`): `E1: N . grad F = A D`, `E2: D div N - N . grad D = F C`, `deg C = 2d - 6`.
* Component parities: `v_X` Y-odd, `v_Y` Y-even. With `D` Y-even: `N_X` Y-odd, `N_Y` Y-even, `Q` Y-odd. A denominator of **mixed** parity is allowed at a fixed degree; `tau`-averaging replaces it by `D . tau(D)` (Y-even, degree `2d`). Hence "Y-even D of degree `<= 2d` excluded" would imply "any D of degree `<= d` excluded", but not conversely.
* Diagonal symmetry group of `F` **[proved, packet]**: among coordinate sign changes only `I, -I, tau, -tau` preserve `F`; among the 720 coordinate permutations only the identity; projectively the diagonal group is `{1, tau}`. (Direct argument in the packet: compare coefficients of `X0^6, X0^5 X1, X0^5 X2, X0^4 Y_j^2, X0^4 Y_j Y_k`; any diagonal map preserving `F` up to scalar has equal X-scalars and equal Y-scalars differing by a sign.) The infinitesimal projective automorphism system (462 x 37, rank 36, minor 68 mod 101) has only the Euler direction, so `Aut_lin(F)` is finite. **Not done**: a classification of the finite non-diagonal linear automorphisms (none are known; a generic sextic has only `{1, tau}`).
* Residue character. `det(tau) = -1`, so `tau* Omega = -Omega`; `A` is Y-odd; hence `omega_A = A Omega / F^2` is `tau`-invariant. In the diagonal-character formula of Chapter 3 (`sum_j (e_j + 1) w_j mod n`, weights `w = (0,0,0,1,1,1)`, `n = 2`): `(3+1) + 1 + 1 = 6 = 0 mod 2`. The class lies in `H^4(X)^tau_prim`, whose Hodge numbers are `(0, 208, 866, 208, 0)` (dimension 1282; `h^{3,1}(X) = dim R_6 = 426 = 218 even + 208 odd`; `b_4(X) = 2606`). The (3,1) class of `A` is **nonzero**: `A` is not in `J(F)_6` (462 x 37 matrix of rank 37, minor 50 mod 101 **[proved, packet]**); also nonzero at Fermat since all exponents are `<= 4`.
* Smoothness **[proved, re-verified here]**. Write `F = sum Z_i^6 + P`. The sums of absolute coefficients of `dP/dZ_i` are `(2113, 2103, 2212, 2126, 2167, 2191)/9973`, all `< 0.222`. At a nonzero point pick a coordinate `Z_i` of maximal modulus `M`; then `|6 Z_i^5| = 6 M^5 > (2212/9973) M^5 >= |dP/dZ_i|`, so `grad F != 0`. Smooth of dimension 4 implies irreducible. The perturbation's integer coefficients have gcd 1, so no prime reduces `F` to Fermat while keeping 9973 invertible (the PR #150-style modular Fermat reduction is unavailable) **[packet]**.
* Non-split **[proved, packet]**: Hessian centre `Z(F) = {C : Hess(G) C symmetric}` is `C . Id` (1890 x 36 system, 35 x 35 minor 76 mod 101), so no invertible linear change makes `F` a sum of forms in disjoint variable blocks (excludes split, crossed, binary and polarized block forms of Chapters 2a/3a). Middle catalecticant (56 x 56, `C[alpha,beta] = coeff(G, alpha+beta) prod (alpha_i+beta_i)!`) has determinant 93 mod 101, so the complex Waring rank is `>= 56`: `F` is not a six-term diagonal form in hidden coordinates, nor a short power presentation (compare `short-power-presentation`, catalecticant rank 8). The bidegree `(4,2)` and `(2,4)` coefficient matrices (15 x 6 and 6 x 15) have rank 5 with explicit quadric spaces `<Y0^2, Y0Y1+Y1^2, Y0Y2, Y1Y2, Y2^2>` and `<X1^2, X1X2, X2^2, X0X1, X0^2+X0X2>`; this excludes a hidden linear reflection acting on one block.
* `F` is not even in any coordinate (e.g. `X0^5 X1`, `Y0^5 Y1` occur), so Chapter 3.3's reflection theorem does not apply; the only finite-order automorphism available for Kollar-Larsen is `tau` with age `3/2` (fixed curves `C_X = {Y = 0, F_{6,0} = 0}`, `C_Y = {X = 0, F_{0,6} = 0}`), so `X/tau` is **not uniruled** **[proved, packet + Chapter 4.7]**.

### 1.3 Reductions valid for any even sextic (packet, Chapter 4.7)

* `tau`-quotient chart on `Y0 != 0`: `x_i = X_i`, `s = Y0^2`, `a = Y1/Y0`, `b = Y2/Y0`; quotient equation `H = P0(x) + s P1(x;a,b) + s^2 P2(x;a,b) + s^3 P3(a,b)` with X-degrees 6,4,2,0 (exact coefficients in the packet's `geometry/quotient-polynomial.json`); weights `(1,1,1,2,0,0)`; coordinate Jacobian `2/Y0`; transported source `B = x0^3 s^2 / 2` (so the target transports polynomially; the obstacle is the primitive, not a character). Piola pullback of a quotient field `w`: `v_Xi = (2/Y0) w_xi`, `v_Y0 = w_s/Y0^2`, `v_Y1 = a w_s/Y0^2 + 2 w_a`, `v_Y2 = b w_s/Y0^2 + 2 w_b`. The divergence-free quotient field `w_s = x0^3/(6 P3)` removes the `s^2` part of the source; the remainder is `-x0^3 (P1 + 2 s P2)/(6 P3)`; physically this is `v_Y = A Y/(6 F_{0,6})` with directional derivative `A [1 + (F_{4,2} + 2 F_{2,4})/(3 F_{0,6})]`.
* Natural cubic cover: `X/tau -> P^2_X x P^2_Y`, fibre equation `a + b t + c t^2 + d t^3 = 0` with `(a,b,c,d) = (F_{6,0}, F_{4,2}, F_{2,4}, F_{0,6})`, `t = (mu/lambda)^2`; discriminant `Delta~ = F_{4,2}^2 F_{2,4}^2 - 4 F_{0,6} F_{4,2}^3 - 4 F_{6,0} F_{2,4}^3 - 27 F_{6,0}^2 F_{0,6}^2 + 18 F_{6,0} F_{4,2} F_{2,4} F_{0,6}` (degree 24, a legal candidate denominator never tested). For this `F` the monodromy is the full `S_3`: an exact slice (`X1 = X2 = Y2 = 0, Y0 = 1, Y1 = t, s = X0^2`) has squarefree degree-12 discriminant (Bezout witness mod 101) **[proved, packet]**. Chapter 4.7: the primitive invariant classes live in the standard representation `V_std` of this cover, and no uniruled intermediate quotient of the Galois closure can serve it.

---

## 2. Companion and related solved instances

`sameVarietyInstances` is empty: no accepted certificate exists on this variety for any target. Related material:

| Instance | Status | Method / denominator | Why it fails here |
|---|---|---|---|
| Fermat, same target `A = X0^3 Y0^3` (Diamond, Chapter 2.2) | solved (not an atlas entry for this `A`) | Subset Euler `v = A E_I / E_I(F)`, needs `deg_I A + |I| = 6`. For this `A`: `I = {X0, Z_a, Z_b}` or `{Y0, Z_a, Z_b}` with `Z_a, Z_b` from `{X1,X2,Y1,Y2}` (12 choices); `D = X0^6 + Z_a^6 + Z_b^6` or `Y0^6 + Z_a^6 + Z_b^6`, degree 6, strict, `N = A (e_I)/6`. | `deg D >= 7` theorem; the deformed denominator `F_{6,0}(X)` fails even E1 (`F_{6,0}` is not in `J_6 = span(G_i)`); the whole radial family `f E_X + g E_Y` fails at every degree (Section 3.2). |
| `fermat-pair-rational-d4-strict` | solved (`benediamond`) | `D = X0^4 - X0^2 X1^2 + X1^4` (quartic cofactor of `X0^6 + X1^6`), `N_0 = B X1 (X0^2 - X1^2)/6`, `N_1 = B X0^3/6`, `B = X2 Y0`, strict. | Needs `A = x^3 y B` with `B` free of the pair and a split sum of two sixth powers; `A = X0^3 Y0^3` has no second linear pair coordinate, and `F` has no isolated binary block. |
| `fermat-q2` | solved (`benediamond`) | `D = E_I(F) = 6 (X1^6 + X2^6)`, `N = (A/2) E_I`, `q = 2`. | Same subset-Euler balance; needs separated powers. |
| `two-ternary-dense-two-odd-hessian-d48` | solved (`benediamond`, PR #41) | Dense-looking but `F = P(X0,X1,Y0) + Q(X2,Y1,Y2)`; binary-cubic covariants (Hessian `H`, discriminant `Delta`, `t^2 = 4H^3 - 27 delta f^2`), `D = 36 X0 Q_b G`, `G = H_P^3 Delta_Q Q_b^2 - H_Q^3 Delta_P P^2`, `deg D = 55`, non-strict (`C` of degree 104); complex coprimality by univariate specialization + Bezout. | Requires a ternary split (and coordinatewise squares); Hessian centre scalar excludes any split. Useful as the template for a large *structured* witness with a Lean proof (Section 6). |
| Diamond, arXiv:2603.29157 (Thm 3.16, Constr. 3.9) | paper | `F = f(X) - f(Y)` with `f` of Waring rank 3; Fermat via `dF(nu_I)`. | `F` has mixed `(4,2),(2,4)` terms and Waring rank >= 56. The paper states no Zariski-open set of even sextics is solved. |
| Holdout family (`holdout`, `-all-six`, `-p3`, `-p-y14`, `-bundle`, `-moving-branch`; `fqlx`, `benediamond`), target `X0^3 Y0^3` in `holdout` | solved | Fermat-plus-two-couplings sextic; adjoin `T^6 = Y0^6 + Y2^6`, reflection on the cover, polynomial curve family, trace; nonconstructive. | Needs few coupling monomials to write a cover with a reflection; dense `F` has all 236 monomials and no reflection on any cover is known. |
| `reflection-coupled-fano` (PR #150), `cyclic-triple-fano` (PR #153) | solved | Single-coordinate reflection, weighted Fano quotient, modular Fermat reduction to seed an incidence. | `F` not even in any coordinate; gcd-1 perturbation forbids the modular reduction. |
| `short-power-presentation`, `eight-plane-configuration` | OPEN, same target `X0^3 Y0^3` | Dense expansions with hidden structure (8 sixth powers; 8 planes `Y_i = +-i X_i`). | Same "no symmetry" class as this instance; neither has a certificate or a discussion result that transfers. |
| Conic theorem (Chapter 2.6): `F = H(Q1..Q4, r, s)` | theorem | `deg D = 10` over a degree-<= 27 extension (270 after trace). | Never tested for this `F` (Section 4); dimension count makes it implausible (family of such `F` has dimension <= 173 < 462). |
| Chapter 3A (split, crossed, polarized), Chapter 4.3 (cyclic-seven), Chapter 4.2 (paired signs) | theorems | all require block or diagonal symmetry | excluded by the Hessian-centre and symmetry computations. |

Chapter 4.4 is the relevant general statement: a `q = 1` certificate for *one* nonzero primitive target on a family dominating the even-sextic space would give all targets, and forces `sum_i h^{2,0}(E_i) >= 208` for the resolved support components. This instance is a single point; no dominating family is known.

---

## 3. Prior work (discussion #59, comment by comment)

### 3.1 The thread

Exactly one comment, no replies, no edits:

* `fqlx`, 2026-10-20T16:38:14Z (comment id `DC_kwDOVLOOJc4BIxgI`, discussion node `D_kwDOVLOOJc4AqE9a`): "Dense Rational search update — 20 October 2026". Explicitly: *no certificate found, no Lean `HasCertificate` claimed*, AI-assisted, with separate implementations of the main arithmetic checks. Five sections: (1) `deg D >= 7` for the full contract; (2) exclusion of a nonlinear polynomial Piola construction, with the restricted-gradient rank table; (3) degree-7 reduction to `w^T Hess(G) w in J_8` (951 quadrics / 120 unknowns), 31 numerical starts failed, two finite-field eliminations did not finish; (4) quotient-pole exclusions (simple `H_s` pole unrestricted; bounded `H_s^2`, `R`, `(sH_s)^2` spaces); (5) a K3 surface from a genus-4 trigonal slice falls under Bolognesi-Laterveer but no correspondence to the target. Attachment: `Dense-Rational-search-packet.zip`, <https://github.com/user-attachments/files/33728423/Dense-Rational-search-packet.zip> (start with `CONTINUATION.md`).

Everything below is from that packet (pinned input commit `d5881c13a7e785b6328cd823133f5245937b93e9`, instance blob `3b5f7734d70d400a1057bdb4a1666a2b3e5bd38c`, record SHA-256 `7337875384930427e148102a0f1c72de1f15f35d23308e1a7439af48a0d889f5`; refresh to main `55d7489d8578e3cb3c9cbf785ff9864fd4adc269`, 2026-10-20, no data change). Python 3.12.14, NumPy 2.3.5, python-flint 0.9.0, SymPy 1.14.0; core replay 26.6 s, full continuation replay 147 s, both logged as PASS. All "rank" statements use the pattern: explicit exact kernel vectors give the characteristic-zero upper bound, a nonzero modular minor gives the lower bound.

### 3.2 Rigorous results over C (all parities, arbitrary `Q`) **[proved]**

| # | Statement | Evidence (mod 101 unless noted) | Packet files |
|---|---|---|---|
| R1 | `deg D <= 4` impossible. Degree-4 E1 kernel is exactly the constant Koszul vectors (all with `D = 0`). | D even: 756 N-coeffs, 66 D-coeffs, Koszul kernel 9, rank 813. D odd: 756/60/6, rank 810. Smaller degrees by multiplying by `X0^(4-d)`. | `algebra/all_denominators_e1*.py/json`, `e1_coefficient_verification.json` |
| R2 | Degree 5: every E1 solution is `D = sum c_i G_i`, `N = A0 c + K grad G` (`c` constant, `K` skew linear). Each single `D = G_i` fails E2 (all six). | E1 rank 1464 each parity (1386 N, 126 D, 45 Koszul, 3 Jacobian-D). Per `G_i`: 1386 numerator + 60 quotient coefficients, 129 exact kernels (126 Euler + 3 fields `N = G W`), source rank cap 1317, augmented rank 1318. | `jacobian_poles.py/.json` |
| R3 | Degree 6: E1 forces `D = grad(G).w` with `w` linear; E2 reduced mod `J` gives `A0 w^T Hess(G) w in J`; the Hessian map on all symmetric quadratic `S_ij` has kernel exactly Koszul + the 36 Euler relations `S_ii = z_i u_i`, `S_ij = z_i u_j + z_j u_i` (from `Hess(G) z = 5 grad G`); rank-one `S` then forces `w_i = lambda z_i` with all `lambda_i` equal, so `D = 6 lambda G`, illegal. **Hence `deg D >= 7`.** | E1: even 3080 rows, 2376 N, 236 D, 162 Koszul, 18 Jac-D, rank 2432; odd 3108 rows, 226 D, 153 Koszul, rank 2431. Hessian: even-output 3108 x (2376 + 216), odd-output 3080 x (2376 + 225), rank 2421 both, minors 26 and 52; integer arithmetic bound `7.4e9 < 2^63`. Independence of kernels via gradient coefficient maps of rank 36 (linear), 126 (quadratic). | `verify_e1_coefficients.py`, `certificate_search/audit/audit_hessian_linear.py`, `hessian_obstruction.py` |
| R4 | Degree 7 and 8 E1 classification: `D = grad(G).w`, `N = A0 w + K grad G`, `w` quadratic (resp. cubic), `K` skew cubic (resp. quartic). Equivalently `(J : A)_7 = J_7`, `(J : A)_8 = J_8`: multiplication by `A` is injective `(R/J)_8 -> (R/J)_{14}`. | d=7: 4284 rows, 3861 N, 396 D, Koszul 420, Jac-D 63, rank 3774 (both parities). d=8: even 5796 rows, 6006 N, 651 D, Koszul 954, Jac-D 168, rank 5535; odd 5832/636/936, rank 5538. Gradient maps rank 336 (cubic), 756 (quartic). | `continuation/algebra/e1_degree7.py`, `e1_degree7.json`, `e1_degree8.json` |
| R5 | Degree 7 necessary condition: `w^T Hess(G) w in J_8`; after removing the 6 Euler gauge directions (`w -> w + ell z`), 951 homogeneous quadrics in 120 unknowns; a legal `D` needs a nonzero gauge class. | derived from R4 | `CONTINUATION_FINDINGS.md` |
| R6 | No two-block-radial field: `v = f E_X + g E_Y` with homogeneous degree-0 rational `f, g` (legal denominators of any degree) never certifies. Exact identity: with `delta = E_X`, `h = f - g`, `L = delta(F) - 3F`, mod `F^2` the defect is `8 F A (F_{4,2} + F_{2,4}) / L^2 != 0`. "At least one direction-ratio component is required." | symbolic, verified by `geometry/audit_geometry.py` | `geometry/geometry-findings.md` S5 |
| R7 | Quotient chart, base-only denominators: no certificate whose denominator depends only on `(a, b)`, with arbitrary rational coefficient functions in `C(a,b)`. | weight-7 reduction: 70 rows x 43 columns, 3 Euler identities cap rank at 40, augmented rank 41 at `(a,b) = (2,3)`, 41 x 41 minor 18. | `geometry/audit_quotient_modules.py --denominator base` |
| R8 | Quotient chart, simple `H_s` pole excluded with arbitrary `C(a,b)` coefficients and arbitrary `Q`. Necessary condition `B H_ss in (H_x0, H_x1, H_x2, H_s)_9 + H_a S_3 + H_b S_3` (obtained by differentiating the mod-`H` normal equation in `s`). | 125 rows x 126 cols (66 + 34 + 26), 3 Koszul relations cap 123, augmented 124; minors 6 (mod 101) and 260 (mod 503) at `(2,3)`; `(S/J)_9` has dim 28, obstruction quotient dim 2. For `H_s^k`, `k >= 2`, this obstruction carries no information. | `continuation/geometry/hs-simple-pole-proof.md`, `audit_hs_jacobian_obstruction.py` |
| R9 | No nonlinear dominant homogeneous polynomial `Phi: C^6 -> C^6` with `F = P(Phi)`, `A = det(dPhi) b(Phi)` (polynomial `P, b`). Proof: restricted gradients on `X0 = 0` and `Y0 = 0` have no syzygies of multiplier degree <= 4 and only the 15 constant Koszul ones in degree 5; a divisibility lemma forces `r = sum(d_i - 1) = 6`, all `d_i = 2`, then the all-quadratic case is killed by the chain rule. Finiteness and smoothness not needed. Independently reconstructed matrices. | `X0=0`: deg 4 715 x 420 rank 420 (minor 32), deg 5 1001 x 756 rank 741 (minor 48); `Y0=0`: 420 (94), 741 (1); 15 Koszul kernels, independence minor 92; integer bound `2.7e12`. | `continuation/geometry/polynomial_piola_obstruction.md`, `verify_piola_rank_premises.py`, `piola-rank-premises-review.json` |
| R10 | Structural: Hessian centre scalar; catalecticant rank 56; only `{+-I, +-tau}` sign symmetries, identity permutation, no infinitesimal automorphisms; `A` not in `J_6`; Fourier structure of the recipe; `(4,2)`/`(2,4)` blocks rank 5; restricted gradients independent even on `X0 = Y0 = 0` (6 x 6 integer minor `-162628800`); the six gradients independent (pure-fifth-power 6 x 6 determinant `45905248509902770925664702690`). | as listed | `research/audit_hessian_center.py`, `algebra/algebra_probe.py`, `symmetry_probe.py`, `target_nonjacobian.py`, `coefficient_structure_probe.py` |
| R11 | Smoothness and irreducibility of `F`; `X/tau` not uniruled (Kollar-Larsen Thm 2, age 3/2); cubic cover has `S_3` monodromy; the CM elliptic companion `z^3 = -P(u)` (Weierstrass `Y^2 = 4X^3 + Delta(P)`, `dX/Y = -du/(3 z^2)`) transfers nothing in the trace-free standard sector because the root 3-cycle acts on `H^1(E)` by a 3-torsion translation (trivially) while `det(C_3 - I) = 3` on `V_std`; the complementary block view gives the isomorphic curve (`v = 1/u`, `w = z/u`). | exact identities in `audit_cubic_resolvent.py`, `audit_elliptic_companion.py` | `continuation/research/` |

### 3.3 Scoped exclusions (bounded families) **[scoped]**

Quotient-chart higher poles with the full six-component field, `Q` arbitrary in `C(a,b)`, derivative-compatible coefficients (`H = alpha s^3 + beta s^2 + gamma s + delta`, `d = H_s`, `R = 3 alpha^2 s^2 + 2 alpha beta s + 4 alpha gamma - beta^2`, identity `Delta + H_s^2 R = H(27 alpha^2 H - 2(2 beta^3 - 9 alpha beta gamma + 27 alpha^2 delta))`; normalized equation `d grad H . N - H (d div N - k grad d . N) + H^2 Q1 = B d^{k+1}`, numerator weights `(4k+2, 4k+2, 4k+2, 4k+3, 4k+1, 4k+1)`):

| Pole | coefficient space (total degree `m` in `a,b`) | points | matrix | Euler kernel | rank cap | augmented | minor |
|---|---|---|---|---|---|---|---|
| `H_s^2` | polynomial, m = 0 | 3 | 2265 x 936 | 125 | 811 | 812 | 99 |
| `R^2` | m = 0 | 3 | 2265 x 936 | 125 | 811 | 812 | 67 |
| `H_s R` | m = 0 | 3 | 3525 x 936 | 125 | 811 | 812 | 87 |
| `H_s^2` | m = 1 | 6 | 4530 x 2808 | 375 | 2433 | 2434 | 38 |
| `R^2` | m = 1 | 6 | 4530 x 2808 | 375 | 2433 | 2434 | 7 |
| `H_s R` | m = 1 | 6 | 7050 x 2808 | 375 | 2433 | 2434 | 56 |
| `H_s^2` | m = 2 | 12 | 9060 x 5616 | 750 | 4866 | 4867 | 11 |
| `H_s^2` | `g/P3`, deg g <= 0 | 3 | 2265 x 936 | | 811 | 812 | 91 |
| `H_s^2` | `g/P3`, deg g <= 1 | 6 | 4530 x 2808 | | 2433 | 2434 | 63 |
| `(s H_s)^2` | `g/P3`, deg g <= 0 | 3 | 4263 x 2176 | | 1868 | 1869 | 31 |
| `(s H_s)^2` | `g/P3`, deg g <= 1 | 8 | 11368 x 6528 | | 5604 | 5605 | 13 |

Also: `P3^k` denominators (`k = 0, 1`) with polynomial coefficients of degree <= 8 (1800 coefficients, 135 Euler directions, augmented rank 1666 > 1665). Two "apparent consistencies" (8-point `H_s^2`, m = 2: ranks 4456/4456; 6-point `(sH_s)^2`: 5502/5502) disappeared with more evaluation points — a documented warning that finite-evaluation consistency is not a certificate. The physical polynomial behind `H_s` is `F_{4,2} + 2 F_{2,4} + 3 F_{0,6} = s H_s = Y0^2 H_s` (degree 6, already excluded by R3 as a physical denominator). Arbitrary rational coefficient functions with `H_s^k`, `k >= 2`, or other fibre poles remain **unresolved**.

### 3.4 Non-rigorous and incomplete searches **[screen]**

* Degree-7 Hessian system: 31 complex numerical starts (6 single Fermat-pair seeds `w_i = x_j^2, w_j = sqrt(-1) x_i^2`, 4 diagonal seeds `w_i = x_i ell_i`, all 15 perfect matchings with random diagonal perturbation, 6 random), scaled-residual threshold `1e-9`, 6 runs refined; none converged. Exact parity-even chart (`w0[X1^2] = 1`): 483 quadrics in 59 affine variables; msolve runs with 1 thread, 3.5 GB address space, 55 s (plain) and 180 s (`-l 42`) stopped at a degree-3 F4 matrix of 28,967 x 37,807. No modular emptiness established.
* Earlier `H_s` first-jet screen: 444 x 335, ranks 292/293 against cap 295 (inconclusive; superseded by R8). Pointwise `H_s^2` first-jet tests mod 503 at `(2,3)`: 825 x 1131 ranks 627/627 (consistent pointwise — but base-derivative compatibility unresolved); `(sH_s)` first jet 522/522 (825 x 690) vs compatible constant-base test 480/481 (minor 81).
* Monomial denominators of degree <= 3 (`algebra_results_p101_d3.json`) — superseded by R1.

### 3.5 Geometric leads (no construction) **[obs in packet]**

* Genus-4 trigonal slice `C: P(s; t0, t1) = 0` in `P(2,1,1)` (12 simple branch points); Galois closure `W` (genus 13), discriminant double cover `D` (genus 5), `J(W) ~ J(D) x J(C)^2`, `Prym(W/D) ~ J(C)^2` carrying the root-standard sector (Recillas-Rodriguez). Kummer twist `delta^6 = Delta` gives genus 25, compositum genus 73.
* Auxiliary K3: minimal resolution `S` of `z^3 + P(s; t0, t1) = 0` in `P(1,1,2,2)` (three `A_1` points), non-symplectic order-3 automorphism fixing `C`; Bolognesi-Laterveer's algebraic Kuga-Satake correspondence applies to `S`. **Missing**: any algebraic correspondence `Gamma in CH^3(T x Q)` from a surface `T` to the resolved quotient `Q` with nonzero action on the target class; the fixed-curve Gysin `H^1(C)(-1) -> H^3(S) = 0` gives nothing; blow-ups add Tate twists, not coniveau.
* Conditional observation: if the ramification divisor of `Q -> P^2 x P^2` (class `xi + 2H_X + 2H_Y`, ample) is smooth, weak Lefschetz gives `h^{2,0} = 0` for it and for the residual branch preimage, so discriminant-type supports alone could not serve a (3,1) class. Not proved for this `F`.
* Structural warning for curve-lifting routes: a bare six-coordinate degree-`d` polynomial curve seed has the reparametrization kernel `u_i = t phi_i'` in its square incidence Jacobian, so the square lifting lemma cannot be applied without fixing the reparametrization.
* Literature checked and found not applicable without new input: Voisin arXiv:1403.3904 S3 (poses the three-coordinate-sign sextic problem; propagation needs a cycle parametrization), Diamond arXiv:2603.29157 (Waring-rank-3 split family; Theorem 2.7/2.13 is the certificate criterion), Fu (coniveau of cup-product kernels; needs Lefschetz standard conjecture and a big class), Hu (complete intersections of quadrics with free involution), Laterveer (Enriques varieties), the Sept-2026 Kuga-Satake preprint (needs a lattice hypothesis and an actual correspondence).

### 3.6 Reproduction packet inventory

Top directory `dense-rational-20261020/`: `README.md` (lower-bound proof), `CONTINUATION.md` (overview), `REPLAY_LOG.md`, `verify_packet.py`, `manifest.sha256.json` (70 + 87 hashes), `algebra/` (E1 ranks d <= 6, symmetries, catalecticant), `certificate_search/` (Hessian audits, `.npz` matrices, kernel JSONs), `geometry/` (quotient polynomial, base/`H_s` modules, `inspect_hs_kernel.sing`), `continuation/algebra/` (d = 7, 8 E1, quadratic-`w` searches, msolve inputs/logs, Fourier/coupling probes, restricted gradients, target witness), `continuation/branch_poles/`, `continuation/geometry/` (simple-pole proof, Piola obstruction and its independent review), `continuation/research/` (K3/Prym addendum, Hessian centre, resolvent, elliptic companion, ramification Singular scripts), `continuation/verify_continuation.py [--full]`, `upstream/` (pinned copies of the record, two accepted solution.md files, docs). Replay commands are listed in each `.md`. Note: these are downloaded, untrusted scripts; they were **not executed** for this dossier — the checks below are independent re-implementations.

### 3.7 Reading of what is and is not excluded

* Rigorously excluded for the full contract: `deg D <= 6` (R1-R3).
* Rigorously excluded structured families (any degree): two-block radial fields (R6); polynomial Piola atoms (R9); quotient certificates with base-only or simple-`H_s` poles (R7, R8).
* Excluded only in bounded coefficient spaces: `H_s^2`, `R^2`, `H_s R`, `(sH_s)^2` quotient poles (3.3).
* **Not excluded**: every denominator of degree >= 7 outside the above families; in particular degree 7 (Jacobian-type, 951 quadrics / 120 unknowns) is open, degree 8 was never searched, and nothing at all is known for degrees >= 10 beyond E1 bookkeeping.

### 3.8 Independent re-verification and new computations for this dossier **[new]**

Scripts in `/tmp/hodge/scratch/` read `instance.rec` directly and use python-flint 0.9.0 `nmod_mat` ranks; none import packet code.

| Item | Result |
|---|---|
| Coefficient recipe, monomial counts (28/90/90/28), derivative sums `(2113,2103,2212,2126,2167,2191)/9973` | all reproduced (`check_dense.py`) |
| Degree-7 E1 ranks | 3774 in both parity sectors; N-block rank 3441 = exact `dim J_13` (Koszul count); multiplication by `A` has rank 333 = `dim (R/J)_7` in both sectors (`rank9/log_d7.txt`) — matches R4 |
| **Degree 9 E1 classification** | Rows 7752 (`S_15` of one parity), N columns 9009, D columns 1001. `rank[N | D] = 7752` (full), `rank[N] = 7129 = dim J_15` exactly; so multiplication by `A` has rank `623 = dim (R/J)_9` in **both** parity sectors mod 101, and in both sectors also mod 10007. Since a modular rank bounds the rational rank from below, `(J : A)_9 = J_9` over C: **every degree-9 E1 solution is `D = grad(G).w` with `w` quartic, `N = A0 w + K grad G`, `K` skew quintic.** (`rank9/mulA_rank.py`, `log_p101.txt`, `log_d9_p10007.txt`) |
| **Degree 10** ranks of multiplication by `A` (`(R/J)_10^+ -> (R/J)_16^-`, dims 762 -> 468; odd sector 744 -> 483) | Even: rows 10152, N 13104, D 1512, `rank[N | D] = 10152` (full), `rank[N] = 9684 = dim J_16^-`, so `x A` has rank 468 = `dim (R/J)_16^-`: **surjective**, kernel exactly `762 - 468 = 294`. Odd: rows 10197, N 13104, D 1491, ranks 10197 / 9714 (= `dim J_16^+`), `x A` rank 483 = `dim (R/J)_16^+`: surjective, kernel exactly `744 - 483 = 261`. Both rigorous over C (a modular rank attaining the target dimension is a lower bound). So at degree 10 the E1-admissible non-Jacobian denominators form spaces of the minimal possible dimensions 294 (Y-even) and 261 (Y-odd); multiplication by `A` has maximal rank in every degree 7-10 tested. (`rank9/log_p101.txt`) |
| E1 membership of structured degree-12 denominators (is `A D in J_18^-`?): the ten products of the Y-degree blocks `F_{6,0}, F_{4,2}, F_{2,4}, F_{0,6}` (`F_{6,0}F_{0,6}`, `F_{6,0}F_{2,4}`, `F_{4,2}F_{0,6}`, `F_{4,2}F_{2,4}`, `F_{6,0}F_{4,2}`, `F_{2,4}F_{0,6}`, and the four squares). Established so far: `rank N = 16589 = dim J_18^-` exactly (mod 101, 754 s), so the test is exact mod p. | **Still running at the time of writing** (`e1prod/log_joint_p101.txt`, joint matrix 16797 x 25714). Reading rule: `rank_joint - 16589 = 0` means all ten pass E1 mod 101 (then rigorous over C? No: membership mod p is only evidence; non-membership mod p is rigorous because rank N already attains its exact cap); `excess = k` means at least `k` independent candidates **fail E1 rigorously** (a modular rank increase lower-bounds the rational rank). Expect failures: at degree 12 the E1 condition is 208 linear conditions on an 885-dimensional space. |
| Served-target spaces `T(D)` of natural degree-6 denominators (which Y-odd sextic targets admit *any* certificate `N/D` with arbitrary `Q`): `dim T(D) = 226 + rank[N | Q] - rank[N | Q | A']` on a 16797 x 2828 system; `J_6^-` (dim 18) is always served. | `D = F_{6,0}` (X-block): ranks 2340 / 2548, `dim T(D) = 18`, i.e. **no primitive target at all** is served. `D = F_{0,6}` (Y-block): identical numbers, `dim T(D) = 18`. These are rigorous over C: the rational solution space reduces mod p into the modular one (saturated lattice), so `dim T_Q <= dim T_p = 18`, and `J_6^-` already gives 18. `D = F_{6,0} + F_{0,6}` (56 terms): ranks 2358 / 2566, `dim T(D) = 18`, again nothing primitive (372 s). Remaining candidates `F_{4,2}`, `F_{2,4}`, `F_{4,2}+F_{2,4}`, `F_{6,0} - F_{0,6}` were **still running at the time of writing** (`served/log_<name>.txt`); Fermat `X0^6+X1^6+X2^6`, each `G_i`, and a random control were not yet launched (run `python3 served/served_one.py 101 <name>`, ~3-6 min and ~0.8 GB each). Reading rule: `dim_T(D)_mod_J6odd = 0` is rigorous (no primitive target served over C, by lattice saturation); a positive value is modular evidence of a served primitive target and must be lifted to Q before any conclusion. Interpretation in 5.9: the two natural cone-type divisors, which plausibly carry ~100 two-forms each, have zero (3,1) Gysin image on the invariant primitive part — exactly what an irreducible invariant Hodge structure predicts. |

Caveat: these are modular computations; where a rank attains its exact characteristic-zero cap the conclusion is rigorous (lower-bound direction); where it does not, only the modular statement is established.

---

## 4. What remains unexplored

1. **Degree 7 is not excluded.** The 951 x 120 quadratic system was attacked with 31 local starts and two msolve runs limited to 3.5 GB and <= 180 s. No Groebner/Macaulay computation of realistic size was attempted; the mixed-parity sector was never set up as its own system.
2. **Degree 8 was never searched**, even though R4 gives the same Jacobian-type reduction: `w` cubic, Hessian condition `w^T Hess(G) w in (J : A)_10`. Per parity sector this is at most 468 quadrics in 156 (even) / 159 (odd) unknowns (Appendix A); mixed parity: at most 951 in 315.
3. **Degree 9** (new reduction above): `w` quartic, 350 unknowns per sector, but the mod-`J` Hessian condition has at most 208 equations (`(R/J)_18^-`), so it is *underdetermined* and no longer selective; the full E2 (linear in `K`, `C` for fixed `w`) has never been solved jointly. Nobody has run a full-E2 variable-projection search at any degree.
4. **Degrees >= 10**: E1 acquires a non-Jacobian kernel; at d = 10 it is exactly 294 (Y-even) / 261 (Y-odd) dimensional (3.8: multiplication by `A` is surjective there), and E1 is vacuous for `d >= 19`. The explicit kernel directions were never computed; no E2 test of any degree-10+ denominator exists.
5. **Structured denominators never tested**: block products `F_{6,0} F_{0,6}` etc. (degree 12), the cubic-cover discriminant `Delta~` (degree 24; the packet calls it "a legal candidate denominator" and stops), `det Hess(F)` (degree 24; Chapter 5.2 gives a universal E1 lift `N_0 = (5/6) adj(Hess F) grad A`), powers and products with hyperplanes (useless at q = 1: support unchanged), `D . tau(D)`-type even denominators of degree 14 built from degree-7 candidates.
6. **Deformation from Fermat**: `why.md` proposes expanding in the perturbation parameter; nobody did it. The Fermat degree-6 and degree-7 solution components for this exact `A` are explicit (Section 2 and 5.6), and the first-order system in direction `P = 9973(F - F_0)` was never written down.
7. **Symmetry sectors**: none exist (group `{1, tau}`), but the finite linear automorphism group was never classified exactly; a non-diagonal automorphism would be a game changer (nothing suggests one).
8. **Conic / Piola-type presentations** of `F` (Chapter 2.6: `F = H(Q_1..Q_4, r, s)`; Chapter 3.4 with rational `b` or sums of atom densities) were not tested.
9. **Cohomological side**: no attempt to decide whether `H^4(X, Q)^tau_prim` is an irreducible Hodge structure for this `F` (it governs whether any "small" divisor can help; see Section 5.9), nor to compute the Gysin images of the natural singular divisors `V(F_{6,0}) ∩ X`, `V(F_{0,6}) ∩ X` beyond the single target `A` (partly done in 3.8).
10. **Quotient-chart higher poles with unbounded rational coefficients** (`H_s^k`, `k >= 2`; other fibre divisors) remain open, but these are special physical denominators of degree <= 12 and the physical E2 test is the cleaner formulation.
11. **No Lean work of any kind** exists for this instance (no helper lemmas, no smoothness/coprimality infrastructure for this `F`).

---

## 5. Candidate attack plan

Sizes below use the exact dimension tables of Appendix A (`S_k^+/S_k^-` = Y-even/odd monomials of degree `k`; `R_k^+/R_k^-` = Jacobian ring of a smooth sextic).

### 5.1 Finish the E1 bookkeeping through degree 12 (rigorous, cheap)

Idea: for each `d`, the rank of `x A : (R/J)_d -> (R/J)_{d+6}` (both parity sectors) decides how many non-Jacobian denominators pass E1. Done for d <= 9 (injective) and d = 10 (surjective; kernels exactly 294 / 261). For d = 11, 12 the map cannot be injective (dims 843 -> 333, 885 -> 208); if it is again surjective the kernels are exactly 510 / 510 and 677 / 648. Matrix sizes: d = 11: 13167 x (18564 + 2184); d = 12: 16797 x (25704 + 3108). Cost: minutes each with `nmod_mat` (d = 10 took ~5 min per sector here). The next step, not yet done at any degree, is to extract the kernel vectors themselves (modular nullspace, then rational reconstruction or an exact recomputation over Q). Output: explicit bases of `(J : A)_d / J_d` (lift modular kernels to Q by rational reconstruction or recompute over Q with a Groebner basis of `J`). Chance of success: certain as a computation; it does not produce a certificate but is the prerequisite for 5.4.

### 5.2 Rigorous degree-7 exclusion (medium)

Idea: show the 951 quadrics (R5) have no common nonzero zero. Per parity sector the system is 483 quadrics in 60 variables; for generic quadrics the Hilbert series `[(1+t)^60 (1-t^2)^423]_+` becomes negative at `t^4` (`C(60,4) - 423 C(60,2) + C(423,2) = -171822`), so the degree-4 Macaulay matrix (483 x C(61,2) = 883,890 rows by C(63,4) = 595,665 quartic monomials) is expected to have full column rank. msolve F4 on this machine (40 cores, 125 GB) with no artificial memory cap is the natural tool; the packet's runs were capped at 3.5 GB / 180 s and only reached degree 3. The **mixed-parity** sector (both `w_e, w_o != 0`: equations `Q_+(w_e) + Q_+(w_o) = 0` in `R_8^+`, `B(w_e, w_o) = 0` in `R_8^-`, 120 unknowns) is harder: generic emptiness would show only at Macaulay degree 5 (288M x 225M), infeasible directly; one needs a reduction (e.g. eliminate `w_o` from the bilinear block, or exclude even denominators of degree 14 instead). Value: raises the proved bound to 8 (pure sectors first). Probability the sector computation finishes: good; probability of an actual solution appearing: low (overdetermined 8:1, Fermat seeds were tested).

### 5.3 Wide numerical search of the Hessian relaxation at d = 7 and d = 8 (cheap)

Idea: the packet used 31 starts. Run `10^4`-`10^5` Levenberg-Marquardt starts with analytic Jacobians on (a) d = 7: 951 x 120 (or 483 x 60 per sector); (b) d = 8, never tried: build `(R/J)_10` projection (dim 762 + 744; `S_10` = 1512 + 1491 monomials; `J_10` = 6 `S_5` - 15 `S_0` = 1497), unknowns `w` cubic: 6 x 56 = 336 minus 21 gauge = 315 (per sector 156/159), equations <= 468 per sector. Any residual below `1e-12` is then promoted exactly: solve E2 for `K` and `C` (linear) and test coprimality. Expected outcome: none, but it is the cheapest way to look at degree 8 and the first to use the mixed-parity space seriously.

### 5.4 Full-E2 variable-projection search at d = 7, 8, 9 (cheap to medium)

Idea (Chapter 5.3-5.4): since E1 forces `D = grad(G).w` through d = 9, parametrize the outer variables by `w` (gauge-fixed) and solve the inner *linear* system in `(K, C)`:

```
B_D(A0 w) + T_D(K) = G C,  T_D(K) = sum_{i<j} ( D[(d_i K_ij) G_j - (d_j K_ij) G_i] - K_ij (D_i G_j - D_j G_i) ),
```
unknowns per parity sector: d = 7: `K` cubic (15 x 28 = 420), `C` degree 8 odd (636), rows `S_14^-` = 5796; d = 8: 945 + 1491, rows 10152; d = 9: 1890 + 3080, rows 16797. Outer unknowns 60 / 156 / 350. Each inner solve is a dense least-squares of at most 16797 x 4970; the outer residual is a rational function of `w`; use LM with analytic derivatives of the bilinear map. This is the first genuine joint search at any degree and subsumes the Hessian relaxation. Expected outcome at d = 7: nothing (consistent with 5.3). At d = 9 the Hessian condition is underdetermined (<= 208 equations in 350 unknowns), so this is the first degree where a positive outcome is not already numerically disfavoured by the mod-`J` relaxation. Still, by Section 5.9 the geometric prior is strongly negative.

### 5.5 Bilinear Gauss-Newton at d = 10-12 including the E1 kernel (heavy)

Idea: `D = D_ker + grad(G).w` with `D_ker` in a complement of `J_d` inside `(J : A)_d` (from 5.1), `w` of degree `d - 5`; unknowns at d = 10: ~294 + (1512 - 126) outer, `K` sextic skew ~3465 + `C` degree 14 odd 5796 inner; rows `S_20^-` = 26532. Treat the whole cleared identity as one bilinear system (Appendix A: d = 10: 20412 unknowns, 84903 equations; d = 12: 45609 / 162248) and run sparse Gauss-Newton/LSQR from random and from structured seeds (5.6). Feasible on this machine; outcome unknown; the geometric prior (5.9) says a solution would require a divisor with a very singular locus, which random seeds will not find.

### 5.6 Deformation from Fermat along the pencil `F_t = F_0 + t P`, `t = 1/9973` (cheap, informative)

Idea (`why.md`): `P = 9973 (F - F_0)` is the integer perturbation (230 terms, coefficients 1..17). At `t = 0` the degree-6 strict certificates are the 12 subset-Euler fields of Section 2, and the degree-7 Hessian solutions include `w_i = x_i ell_i` (arbitrary linear `ell_i`, i.e. degree-6 solutions times a linear form) and the pair seeds `w_i = x_j^2, w_j = sqrt(-1) x_i^2` (for every pair `i != j`, since `30 x_i^4 w_i^2 + 30 x_j^4 w_j^2 = 0`). Write `D = D_0 + t D_1 + ...`, `N = N_0 + t N_1 + ...`, `Q = t Q_1 + ...` and impose the `t^1` coefficient:

```
2 A D_0 D_1 - D_0 N_1.grad F_0 - D_1 N_0.grad F_0 - D_0 N_0.grad P
 + F_0 (D_0 div N_1 + D_1 div N_0 - N_1.grad D_0 - N_0.grad D_1)
 + P (D_0 div N_0 - N_0.grad D_0) = F_0^2 Q_1 .
```
This is linear in `(D_1, N_1, Q_1)` for each Fermat solution `(D_0, N_0)`; inconsistency on every component proves that no analytic family of degree-`d` certificates passes through Fermat in direction `P` (then a degree-`d` certificate could exist only at finitely many special `t`, and `1/9973` would have to be one of them). Also use the generalized gauge `(D, N) -> (D + t ell F, N + t A ell E/6)` (`ell` of degree `d - 6`) to normalize limits: it preserves E1 exactly and E2 modulo `F`, by the identity `B_{D + t ell F}(N + t A ell E/6) - B_D(N) = t [F (ell div N - N(ell)) + ell (A D - N(F))] + t^2 ell^2 F A` (hand derivation checked symbolically on random data, `/tmp/hodge/scratch/gauge_check.py`). In particular a family of degree-`d` certificates for `F_t` can have a limit denominator divisible by `F_0` at `t = 0`, which is why degree-`d` denominators at Fermat should be studied modulo this gauge. Numerically: track solutions of the degree-7 (and 8) systems from `t = 0` to `t = 1/9973` by homotopy on the overdetermined system (Gauss-Newton path tracking). Expected: obstruction at first order for degree 6 (consistent with R3) and probably for 7; the computation explains the degree jump and tests whether the pencil's generic member has a low-degree certificate. Note the gcd-1 fact: no modular reduction maps `F` to Fermat, so this must be done over Q or numerically.

### 5.7 Structured denominators: cheap necessary filters, then one expensive E2 test

Candidates and the first test (E1 membership `A D in J`, free for `d >= 19`):

* `F_{6,0} F_{0,6}`, `F_{6,0} F_{2,4}`, `F_{4,2} F_{0,6}`, `F_{4,2} F_{2,4}` (degree 12): E1 test 16797 x 25714 mod p — running, see 3.8 for the reading rule. If one passes, run E2: fixed-`D` correction system with `c_ij` of degree 8 (15 x ~643) and `C` of degree 18 odd (16797), rows `S_24^-` = 59332 — a 59k x 26k dense modular rank, ~12 GB, feasible here.
* `Delta~` (cubic-cover discriminant, degree 24) and `det Hess(F)` (degree 24): E1 automatic; E2 needs `c_ij` of degree 20 (15 x 26.5k) + `C` of degree 42 (767k) unknowns and 1.43M rows — requires a sparse modular solver (or a Groebner basis of `J` mod p via msolve to work in normal forms). Geometric prior for `Delta~`: the packet's conditional weak-Lefschetz remark suggests the ramification/discriminant divisors carry no two-forms.
* Any candidate must satisfy Chapter 3.1's necessary condition: `V(D) ∩ X` must have a component that is non-normal or has non-rational singularities (a normal divisor with rational singularities on a smooth sextic fourfold has `H^{2,0} = 0` after resolution). Design `D` accordingly (e.g. cones whose vertex meets `X` in a curve, so that the transversal singularity is a cone over a curve of positive genus, as happens for `V(F_{6,0}) ∩ X` along `C_Y`), rather than random polynomials.

### 5.8 Groebner-basis infrastructure for the Jacobian ring (enabling)

Compute a degree-reverse-lexicographic Groebner basis of `J = (G_0, ..., G_5)` modulo a large prime with msolve (zero-dimensional, colength `5^6 = 15625`). Then every E1 question, every Hessian relaxation and every `(J : A)` kernel becomes a normal-form computation instead of a dense rank; this makes 5.1 and 5.3 routine up to degree 24 and allows the socle pairing `(R/J)_k x (R/J)_{24-k} -> (R/J)_{24}` to be used for compact obstruction tests.

### 5.9 Theory: decide whether small divisors can help at all

Observation **[obs]**: if `H^4(X, Q)^tau_prim` (dimension 1282) is an irreducible polarized Hodge structure for this `F` (expected for a Mumford-Tate-generic even sextic; nothing structural argues against it after R10), then the Gysin image of any divisor is either 0 or everything on the primitive invariant part; in particular a divisor whose resolution carries fewer than 208 holomorphic two-forms contributes **nothing** to `H^{3,1}(X)^tau`, and a successful support serves all 208 odd targets simultaneously. Consistency check: `E = V(F_{6,0}) ∩ X` is singular along `C_Y` with transversal type "cone over the plane sextic `C_X`" (genus 10), so its resolution plausibly carries about `g(C_X) g(C_Y) = 100` two-forms (as for the Fermat join), yet by R3 it serves no odd target — exactly what irreducibility predicts. The served-target computation of 3.8 tests this on several degree-6 denominators: any nonzero `T(D)/J` would disprove irreducibility and expose a Hodge substructure (a major structural discovery); zero is the expected outcome. Proving irreducibility is out of reach (Frobenius characteristic polynomials of degree 1282 are not computable here), but the search strategy should assume it: look only for divisors with huge singular loci, or for a general mechanism.

### 5.10 General mechanism (research-level, the only route to a real solve)

Chapter 4.4 says one generic certificate gives all; Diamond asks for a construction on a Zariski-open set. Candidates: (i) a divisor `E` fibred over a curve by surfaces with large `p_g` such that `H^{2,0}(E~)` has dimension >= 208 and nonzero Gysin image — e.g. `E = V(D) ∩ X` with `V(D)` a cone over a surface in `P^3` or `P^4` whose vertex meets `X` in a curve; compute `h^{2,0}` via the Leray spectral sequence for the resolution (`H^1(C, R^1 pi_* O)` terms) and test E2 only for candidates clearing 208; (ii) the Prym/K3 correspondences of 3.5, where the missing piece is an explicit `Gamma in CH^3(T x Q)` acting nontrivially on the root-standard sector. Neither has a concrete next computation beyond those already listed; the expected effort is large and the outcome uncertain.

---

## 6. Feasibility assessment

### 6.1 Explicit witness

* Any explicit witness must have `deg D >= 7` (R3), be non-Jacobian-type only from degree 10 on (R4 and 3.8), and its support must contain a component with non-normal or non-rational singularities (Chapter 3.1). The known constructive mechanisms all rely on structure this `F` lacks. Under the irreducibility expectation of 5.9 a witness of modest degree would have to carry >= 208 two-forms on a degree-`d` divisor — no example of such a divisor on a symmetry-free sextic is known at any degree.
* If a witness of degree `d` were found, the Lean obligations (`HodgeComputer.Verify`) are: `D != 0`, homogeneity of `D` and the six `N_i` (degree `d + 1`), `IsRelPrime F D` over C, and `E = F^2 Q`. Coprimality is convenient here: `F` is monic in every variable (coefficient 1 on each `Z_i^6`), so `relPrime_of_specialization` reduces it to a univariate Bezout identity after specializing five variables, exactly as in `two-ternary-dense` (`PDenseP0.lean`). The identity is the problem: `E` has degree `2d + 6` with up to `dim S_{2d+6}^-` terms (d = 7: 26,532; d = 9: 59,332; d = 12: 162,248), and `F` itself has 236 terms, so `F^2 Q` involves ~`6e3 x dim S_{2d-6}` term products (d = 7: ~4e6; d = 12: ~1e8). `ring`-style normalization inside the 60-minute / 14 GiB budget is plausible only for `d <= 8` or so unless the witness has a factored structure (the degree-55 `two-ternary-dense` proof worked because its covariant identities were proved symbolically, not by expansion). A dense unstructured witness of degree >= 10 would need new verified sparse-polynomial arithmetic in the shared library (a separate maintenance PR) or a reflective checker; `native_decide` is refused.

### 6.2 Nonconstructive argument

* The accepted nonconstructive proofs (holdout family, reflection/cyclic Fano, mixed-biweight) formalize a specific algebraic geometry: a cover with a reflection, a dominant polynomial curve family with polynomial normal, integration, trace, and descent (`HoldoutCover`, `PolynomialCurves`, `ResidueTransport`, `HigherResidues`, `TransitiveCubic`). None of these inputs exists for this `F`: no reflection on any known cover, no Fano quotient (`X/tau` is not uniruled), no polynomial Piola atom (R9), no split. A nonconstructive proof would first require a *mathematical* existence argument, which today would amount to a new case of the generalized Hodge conjecture for a symmetry-free even sextic; its formalization (Hodge theory, Gysin maps, mixed Hodge modules used in Chapter 1's prescribed-support theorem) is far beyond the current library.
* Chapter 1's descent lemma guarantees that any complex certificate yields a rational one (via the finite polynomial system and trace), so working over number fields or numerically is legitimate for discovery; the Lean statement is over C and needs only the final rational witness or an existence proof.

### 6.3 Bottom line

* Probability of a complete solve (Lean-accepted) by the methods in this dossier: very low. The instance is designed as the generic benchmark and sits exactly at the frontier the book describes (Chapters 4.4, 4.7).
* High-value contributions that are feasible now: (a) the degree-9 classification (new, this dossier) and degree-10-12 E1 kernels (5.1); (b) a pure-sector rigorous degree-7 exclusion via a real F4 run (5.2); (c) the first full-E2 searches at d = 7-9 (5.4) and the first degree-8 Hessian search (5.3); (d) the Fermat-deformation obstruction computation (5.6); (e) the structured-denominator E1/E2 tests (5.7) and the served-target computations (3.8), which also probe the Hodge-structure question (5.9). Each of these is a concrete, reproducible discussion post (with scripts and modular witnesses), not a certificate.

---

## Appendix A. Dimension tables (exact combinatorics, `/tmp/hodge/scratch/dims.py`)

Jacobian ring `R = S/J` of a smooth sextic in six variables, split by Y-parity (`+` even, `-` odd); `dim S_k^+ / S_k^-` = Y-even/odd monomials of degree `k`.

```
 deg | dimS+ dimS- | dimR+ dimR- | total
   0 |     1     0 |     1     0 |     1
   1 |     3     3 |     3     3 |     6
   2 |    12     9 |    12     9 |    21
   3 |    28    28 |    28    28 |    56
   4 |    66    60 |    66    60 |   126
   5 |   126   126 |   123   123 |   246
   6 |   236   226 |   218   208 |   426
   7 |   396   396 |   333   333 |   666
   8 |   651   636 |   483   468 |   951
   9 |  1001  1001 |   623   623 |  1246
  10 |  1512  1491 |   762   744 |  1506
  11 |  2184  2184 |   843   843 |  1686
  12 |  3108  3080 |   885   866 |  1751
  13 |  4284  4284 |   843   843 |  1686
  14 |  5832  5796 |   762   744 |  1506
  15 |  7752  7752 |   623   623 |  1246
  16 | 10197 10152 |   483   468 |   951
  17 | 13167 13167 |   333   333 |   666
  18 | 16852 16797 |   218   208 |   426
  19 | 21252 21252 |   123   123 |   246
  20 | 26598 26532 |    66    60 |   126
  21 | 32890 32890 |    28    28 |    56
  22 | 40404 40326 |    12     9 |    21
  23 | 49140 49140 |     3     3 |     6
  24 | 59423 59332 |     1     0 |     1
```
Total `dim R = 15625 = 5^6`; socle degree 24; `h^{3,1}(X) = 426`, `h^{2,2}_prim = 1751`, invariant primitive Hodge numbers `(208, 866, 208)`.

E1 for denominator degree `d` (target `A` Y-odd): `x A : R_d^{pd} -> R_{d+6}^{pd+1}`; forced kernel = `max(0, dim R_d - dim R_{d+6})`:

```
 d | pd | dim R_d | dim R_{d+6}' | forced kernel
 7 | 0,1|   333   |    843       |   0   (injective: proved, packet)
 8 | 0  |   483   |    744       |   0   (injective: proved, packet)
 8 | 1  |   468   |    762       |   0
 9 | 0,1|   623   |    623       |   0   (injective: proved here, 3.8)
10 | 0  |   762   |    468       | 294
10 | 1  |   744   |    483       | 261
11 | 0,1|   843   |    333       | 510
12 | 0  |   885   |    208       | 677
12 | 1  |   866   |    218       | 648
13 | 0,1|   843   |    123       | 720
14 | 0  |   762   |     60       | 702
...
18 | 0  |   218   |      0       | 218   (E1 vacuous for even D from d = 18; for all D from d = 19)
```

Hessian-relaxation sizes when `D = grad(G).w` (parity of `w_i` = parity(D) + parity(z_i); `w^T Hess(G) w` is always Y-even; Euler gauge `w -> w + ell z`, `ell` of degree `d - 6` with the parity of `D`):

```
 d | pd | #w coeffs | gauge | unknowns | #equations <= min(dim R_{2d-6}^+, dim R_{2d}^-)
 7 | 0,1|    63     |   3   |    60    | 483     (packet: 951 in 120 without parity split)
 8 | 0  |   168     |  12   |   156    | 468
 8 | 1  |   168     |   9   |   159    | 468
 9 | 0,1|   378     |  28   |   350    | 208     (underdetermined)
10 | 0  |   756     |  66   |   690    |  60
11 | 0,1|  1386     | 126   |  1260    |   9
12 |    |  2376     |       |          |   0     (relaxation vacuous)
```

Full cleared bilinear system (`D` Y-even; `N` standard parities; `Q` Y-odd of degree `2d - 6`; equations = `dim S_{2d+6}^-`; Euler redundancy = `dim S_d^-`):

```
 d | D unknowns | N unknowns | Q unknowns |  total  | equations | gauge
 6 |    236     |    2376    |     226    |   2838  |   16797   |   226
 7 |    396     |    3861    |     636    |   4893  |   26532   |   396
 8 |    651     |    6006    |    1491    |   8148  |   40326   |   636
 9 |   1001     |    9009    |    3080    |  13090  |   59332   |  1001
10 |   1512     |   13104    |    5796    |  20412  |   84903   |  1491
11 |   2184     |   18564    |   10152    |  30900  |  118608   |  2184
12 |   3108     |   25704    |   16797    |  45609  |  162248   |  3080
14 |   5832     |   46512    |   40326    |  92670  |  287793   |  5796
16 |  10197     |   79002    |   84903    | 174102  |  481194   | 10152
18 |  16852     |  127512    |  162248    | 306612  |  766843   | 16797
20 |  26598     |  197340    |  287793    | 511731  | 1174380   | 26532
24 |  59423     |  427518    |  766843    |1253784  | 2502990   | 59332
```

Fixed-`D` E2 correction system (Chapter 5.2): unknowns `c_ij` (15 of degree `d - 4`, roughly half by parity) + `C` (degree `2d - 6`, Y-odd); rows `S_{2d}^-`; second-syzygy redundancy `20 x S_{d-9}`:

```
 d | c_ij  |   C    |  rows  | redundancy
 7 |   420 |    636 |   5796 |     0
 8 |   945 |   1491 |  10152 |     0
 9 |  1890 |   3080 |  16797 |     0
10 |  3465 |   5796 |  26532 |    60
12 |  9645 |  16797 |  59332 |   560
18 | 87210 | 162248 | 374604 | 20020
24 |398475 | 766843 |1434680 |155040
```

## Appendix B. Formulas worth keeping at hand

* Cleared identity and normalized equations: Section 1.2. Euler shift `v -> v + phi E` (`phi` of degree 0, Y-odd) changes nothing; the degree-`d` gauge `(D, N) -> (D + t ell F, N + t A ell E/6)` (`deg ell = d - 6`) preserves E1 and E2 mod `F`.
* Fermat subset-Euler certificate for this `A`: `v = A E_I / E_I(F_0)`, `I = {X0, Z_a, Z_b}` or `{Y0, Z_a, Z_b}`, `Z_a, Z_b in {X1, X2, Y1, Y2}`, `D = X0^6 + Z_a^6 + Z_b^6` or `Y0^6 + Z_a^6 + Z_b^6`, strict.
* Fermat degree-7 Hessian seeds: `w_i = x_j^2`, `w_j = sqrt(-1) x_i^2` (`D = 6 x_i^2 x_j^2 (x_i^3 + sqrt(-1) x_j^3)`), and `w_i = x_i ell_i`.
* Universal E1 lift (Chapter 5.2): `D_0 = det Hess(F)`, `N_0 = (5/6) adj(Hess F) grad A`, from `Hess(F) Z = 5 grad F`.
* Chapter 5.2 correction map `T_D` and the E2 criterion `-B_D(N_0) in F S_{2d-6} + im T_D` (Section 5.4).
* Quotient-chart Piola transport and the cubic-cover discriminant: Section 1.3. Discriminant identities: `Delta + H_s^2 R = H(27 a^2 H - 2(2b^3 - 9abc + 27a^2 d))`; `H_s^{-1} = B/Delta mod H` with `B = (2ab^2 - 6a^2c) s^2 + (9a^2 d - 7abc + 2b^3) s + 3abd - 4ac^2 + b^2 c`.
* Two-block-radial obstruction: `3L - delta(L) = -9F + 8(F_{4,2} + F_{2,4})`, defect mod `F^2` equal to `8 F A (F_{4,2} + F_{2,4})/L^2`.
* Fermat-join two-form count: resolution of `V(X0^6+X1^6+X2^6) ∩ X_Fermat` is a `P^1`-bundle over `C_X x C_Y`, `h^{2,0} = 10 x 10 = 100`; the 208-form threshold of Chapter 4.4 needs at least three such components.

## Appendix C. Files produced for this dossier

* `/tmp/hodge/scratch/check_dense.py` — recipe, counts, derivative sums.
* `/tmp/hodge/scratch/dims.py` — Appendix A tables.
* `/tmp/hodge/scratch/rank9/mulA_rank.py` — E1 ranks and multiplication-by-`A` ranks; logs `log_d7.txt` (d = 7, p = 101), `log_p101.txt` (d = 9, 10), `log_d9_p10007.txt`.
* `/tmp/hodge/scratch/e1prod/e1_candidate.py` — E1 membership of block-product denominators (first run; died after the `rank_N` step, log `log_p101.txt`); resumed as `e1prod/e1_joint.py` (single matrix `[N | A D_1 .. A D_10]`, ten block products), log `log_joint_p101.txt`.
* `/tmp/hodge/scratch/served/served_targets.py` — served-target dimensions for degree-6 denominators (first run; completed `F_{6,0}`, `F_{0,6}`, `F_{6,0}+F_{0,6}`, log `log_p101.txt`); resumed per candidate as `served/served_one.py`, logs `served/log_<name>.txt`.
* `/tmp/hodge/scratch/gauge_check.py` — symbolic check of the degree-`d` gauge identity of 5.6.
* Packet extracted (untrusted, not executed) at `/tmp/hodge/downloads/dense-packet/dense-rational-20261020/`.
