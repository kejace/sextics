## 2. Experiments (all exact mod p unless stated; engine v2; machine shared with ~10 other agents, load ≈ 300/40 cores)

Conventions: `S_d^±` = Y-even/odd monomials of degree `d`; `R = S/J(F)`; `J_d^±`; `K_full(d) = {D ∈ S_d^± : A·D ∈ J}`;
"kernel mod J" = `dim K_full − dim J_d`. Prime `p = 2147483647` unless stated. A modular rank that attains its
characteristic-zero cap is a rigorous lower bound over `Q` (and over `C`); so **"kernel mod J = 0"** results are rigorous,
"inconsistent mod p" results are rigorous when the coefficient matrix rank mod p equals its rank over `Q`
(true whenever the modular rank equals the exact cap; otherwise they are modular evidence only).

### 2.1 Sanity (menu item 0)
| test | result |
|---|---|
| `A ∈ J(F)`? (`D = 1`, E1) | **no**: rank 18 / rank_aug 19 → `[A] ≠ 0` in `R_6^-`; `D = 1` inconsistent |
| diagonal symmetry / characters | group `{scalars} × ⟨τ⟩` only (full monomial support, exponent lattice rank 5): no sector decomposition; any `D` of pure τ-parity reaches the single Y-odd coset |

### 2.2 E1-kernel degree sweep (menu item 5): `dim ker([A] : R_d^+ → R_{d+6}^-)`
`direct`: Macaulay matrix `[A·S_d^+ | S_{d+1}·∇F]` with rows `S_{d+6}^-`, nullspace projected to the `D`-block.
`dual`: `ker_d = dim R_d^+ − rank([A] on R_{18−d}^+)`.  Script `e1kernel.py`, log `e1kernel_sweep.jsonl`.

| d | mode | rows × cols | K_full | dim J_d^+ | **kernel mod J** | time |
|---|---|---|---|---|---|---|
| 1 | direct | 396 × 66 | 0 | 0 | **0** | 0.1 s |
| 2 | direct | 636 × 180 | 0 | 0 | **0** | 0.1 s |
| 3 | direct | 1001 × 406 | 0 | 0 | **0** | 2 s |
| 4 | direct | 1491 × 822 | 0 | 0 | **0** | 5 s |
| 5 | direct | 2184 × 1512 | 3 | 3 | **0** | 11 s |
| 6 | direct | 3080 × 2612 | 18 | 18 | **0** | 32 s |
| 7 | direct | 4284 × 4257 | 63 | 63 | **0** | 89 s |
| 8 | direct | 5796 × 6657 | 168 | 168 | **0** | 437 s |
| 9 | direct | 7752 × 10010 | D9_KFULL | 378 | **D9_KER** | D9_TIME |
| 10 | dual (e=8) | 5832 × 6642 | – | 750 | 294 | 198 s |
| 11 | dual (e=7) | 4284 × 4257 | – | 1386 | 510 | 105 s |
| 12 | dual (e=6) | 3108 × 2602 | – | 2376 | 677 | 26 s |
| 13 | dual (e=5) | 2184 × 1512 | – | – | 720 | 10 s |
| 14 | dual (e=4) | 1512 × 816 | – | – | 702 | 3 s |
| 15 | dual (e=3) | 1001 × 406 | – | – | 595 | 1 s |
| 16 | dual (e=2) | 651 × 177 | – | – | 474 | 0.2 s |
| 17 | dual (e=1) | 396 × 66 | – | – | 330 | 0 s |
| 18 | dual (e=0) | 236 × 18 | – | – | 218 (= dim R_18^+, E1 vacuous from here) | 0 s |

Reading: **every τ-even denominator of degree ≤ 9 that passes E1 lies in `J(F)`** (`(J:A)_d^+ = J_d^+`, rigorous; d ≤ 8 new
here, d = 9 confirms the dossier's mod-101/10007 computation). For `d ≥ 10` the map `[A]: R_d^+ → R_{d+6}^-` is surjective
and the kernel has the minimal possible dimension `dim R_d^+ − dim R_{d+6}^-`. (The τ-odd sector was done by `fqlx`
for d = 7, 8 and by the dossier for d = 9, 10: same conclusion, kernel 261 at d = 10.)
Degree-10 kernel basis (both parities): see §2.7.

### 2.3 Degree ≤ 6 families (menu items 1–3): all excluded by §2.2
Because `ker = 0` for `d ≤ 9`, a candidate of degree ≤ 9 passes E1 iff `D ∈ J_d`. For `d ≤ 4`, `J_d = 0`: **every
denominator of degree ≤ 4 is excluded** (and `deg D ≥ 7` by fqlx's argument anyway). Membership `D ∈ J_d` is a small
rank test (`families.py: in_J`). All of the following are *not* in `J_6` and are therefore excluded (E1 fails, rigorous):
the 12 Fermat subset sums `Z_i^6+Z_j^6+Z_k^6` (one of `X0, Y0` in `{i,j,k}`) that are legal at Fermat; the full pure-power
sum `Σ Z_i^6`; the Y-degree blocks `a = F_{6,0}`, `b = F_{4,2}`, `c = F_{2,4}`, `d = F_{0,6}`; `a + d` (F minus
coupling), `a − d`; `X0^6 + Y0^6`; `(X0²+X1²+X2²)³`; `(X0²+Y0²)³`; `X0^5+X1^5+X2^5`, `X0^5+X0 Y0^4` (d = 5, not in `J_5`);
all degree-4/3/2 cones (`J_d = 0`). The three degree-6 elements of `J`:
`E_X(F) = Σ X_i F_{X_i} = 6a+4b+2c`, `E_Y(F) = 2b+4c+6d`, `E_X(F) − E_Y(F)` pass E1 trivially; full fixed-D systems
(16797 × 2602): DEG6B_RESULTS

### 2.4 Degree 7: the polar (Jacobian) sector, rigorous Gröbner computation (menu 3/5; dossier §5.2)
By §2.2 a τ-pure degree-7 denominator is `D = w(F) = Σ w_i F_i` with `w_i` quadratic of parity `pd + p_i`
(`pd = parity(D)`, `p_i = 0/1` for `X_i/Y_i`), modulo the gauge `w ↦ w + ℓ·ε` (`ℓ` linear of parity `pd`).
Necessary condition (§4.1, needs `[A]` injective on `R_8^±`, which is §2.2 / fqlx): `w^T Hess(F) w ∈ J_8^+`
(always Y-even), i.e. `dim R_8^+ = 483` homogeneous quadrics in `63 − 3 = 60` unknowns per sector
(`polar_quadsys2.py`, prime `1073741827`). Solved with msolve (F4, 8 threads):

| sector | system | mode | msolve result | time |
|---|---|---|---|---|
| τ-even (`pd = 0`) | 483 quadrics + random affine slice `L(c) = 1`, 60 vars | sliced | MS_EVEN_SLICE | MS_EVEN_SLICE_T |
| τ-odd (`pd = 1`) | 483 quadrics + slice, 60 vars | sliced | MS_ODD_SLICE | MS_ODD_SLICE_T |
| τ-even | 483 homogeneous quadrics, 60 vars | homogeneous (dimension of `V`) | MS_EVEN_HOM | MS_EVEN_HOM_T |
| τ-odd | 483 homogeneous quadrics, 60 vars | homogeneous | MS_ODD_HOM | MS_ODD_HOM_T |

MS_CONCLUSION

The **mixed sector** (`V(D)` not τ-stable, `w = w_e + w_o` both nonzero, 120 unknowns, 951 quadrics) is out of reach
of F4 (generic emptiness appears only at Macaulay degree 5: `C(124,5) ≈ 2.3·10^8` columns). Two substitutes:
(i) *support reduction*: if such a `D` is legal then `G = squarefree(D·τD)`, which is τ-pure of degree ≤ 14, is legal
(prescribed-support theorem); so excluding τ-pure denominators through degree 14 would exclude all degree-7 ones;
(ii) *numerical search* (menu 3/5.3; `lm_mixed7.py`): complex Levenberg–Marquardt on all 951 quadrics (both parities of
`R_8`), 126 − 6 = 120 gauge-fixed complex unknowns, random complex normalisation `L(w) = 1`, analytic Jacobian:
LM_MIXED_RESULTS

Cross-check of the engine on two explicit degree-7 polar candidates with the full fixed-D system
(rows `S_20^-` = 26532, cols 4497): DEG7X_RESULTS

### 2.5 Degree 8 (menu 3): structural `J`-elements and Hessian minors
Degree-8 τ-even candidates must again lie in `J_8` (§2.2). The Hessian minors `F_{X0X0}F_{Y0Y0} − F_{X0Y0}²`,
`F_{X0X0}F_{X1X1} − F_{X0X1}²`, the sums `X0^8+X1^8+X2^8`, `a(X)(X0²+X1²+X2²)`, `a(X)(Y0²+Y1²+Y2²)` are not in `J_8`
→ excluded. Polar candidates tested with the full fixed-D system (rows `S_22^-`, cols ≈ 7500): DEG8X_RESULTS
The degree-8 Hessian relaxation (`w` cubic, 156/159 unknowns per sector, 468 quadrics in `R_10^+/(J:A)_10`) cannot
be decided by F4 (generic emptiness only at Macaulay degree ≥ 6, `C(161,6) ≈ 2·10^10` columns); it needs the degree-10
kernel (§2.7) and a numerical search — not completed here.

### 2.6 Deformation from Fermat (menu 4; `deform.py`)
`F_ε = F0 + ε F1`, `F1 = 9973(F − F0)` (integer coefficients 1..17), `ε = 1/9973`. For each of the 12 Fermat
subset-Euler certificates `(D_I, N_0 = A E_I/6, Q_0 = 0)` the first-order equation of §4.2 is linear in
`(D_1 ∈ S_6^+, N_1 ∈ S_7 (standard parities), Q_1 ∈ S_6^-)`: DEFORM_RESULTS

### 2.7 Degree-10 E1 kernel basis and the first E2 tests of degree-10 denominators (`e1kernel2.py`, `e2norm.py`)
D10_RESULTS

### 2.8 Served-target dimensions of `J`-elements (simplicity probe, §1.3; `timage.py`)
Prediction under simplicity of `H^4(X,Q)^τ_prim`: every `D` passing E1 serves either only `J_6^-` (dim 18) or everything
(226). TIMAGE_RESULTS

## 3. Witness
**None found.** No candidate passed E2; `witness-candidate.json` was not written.
