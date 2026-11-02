# Referee report — klein-coupled-character, wave-2 exclusion packet (lane `klein-coupled-character`, discussion #116)

Referee: Claude (Anthropic) for `kejace`, 2026-11-01/02. No human review. Nothing Lean-verified. No certificate claimed.
Inputs: draft post, `wave2-klein-coupled-character-excl.md`, packet `wave2-excl/` (README, scripts, receipts), verify notes, aggregator status, dossier, saved thread.
Replays: own code in `wave2-excl/referee/` (ref_lib.py builds everything from the instance record; no dependence on `w2lib.py`).
Referee prime `p1 = 1000000009` (unused by the excl lens [1000003, 32003, 2^31-1] and the verify lens [10^9+7, 999999937]).
Exact = characteristic 0 (python-flint `fmpz_mat`/`fmpq_mat`, Singular over `Q`).

Verdicts: VERIFIED = replayed by the referee (exact, or modular at p1 with the transfer argument checked); PLAUSIBLE = modular-only
evidence or not replayed; UNSUPPORTED = no receipt / claim not established.

| # | claim (draft / excl report / README) | verdict | referee evidence |
|---|---|---|---|
| 1 | T1: `(J:A)_d = J_d` for `d <= 8`, both parities, all 7 characters (hence every legal denominator of degree <= 8 lies in `J(F)`) | VERIFIED | `ref_colon.py` at p1: kernel of `[A]` zero in all 14 sectors for d = 4..8 (`referee_colon_p1000000009.json`); `dim J_d` equals the CI value for all d; transfer over `Q` by semicontinuity needs flatness, see #3 |
| 2 | kernel tables d = 9: even `[2,3,3,2,10,1,3]`, odd `[0,5,0,0,9,6,4]`; d = 10: even `[24,46,39,39,54,46,46]`, odd `[36,39,39,36,39,36,36]` | VERIFIED (as mod-p upper bounds over `Q`; d = 9 even is exactly 24 by verify row 18) | same numbers at p1 |
| 3 | flatness at the lens's prime (`dim J_d mod p` = CI value), used for every mod p -> `Q` transfer | VERIFIED | Singular: Jacobian ideal of 7F is m-primary with `vdim 15625 = 5^6` and the CI Hilbert numerator at 1000003 and at p1 (`regseq_*.out`) |
| 4 | d <= 4: no denominator | VERIFIED | `(J:A)_d = J_d = 0` at p1 |
| 5 | d = 5: no denominator of any parity (relaxation empty, "Macaulay degree 2 full rank") | VERIFIED | exact rank 21 of the 21 quartics `d_i d_j (7F)` (`ref_hess21.py`), `(J:A)_4 = 0`; own relaxation at p1 m-primary for even/odd/mixed (`relax/d5_*`) |
| 6 | d = 6, Y-odd: relaxation empty (dim 0, vdim 57) | VERIFIED | own builder at p1: `gb_size 135 dim 0 vdim 57` |
| 7 | d = 6, Y-even: relaxation variety is the curve `C` inside the weighted-Euler family (`u_j^2` in `I_Q` for the 12 off-diagonal coordinates), `C` = CI of the 3 quadrics of `C_quadrics_exact.json` | VERIFIED (exact) | own builder over `Q`, Singular: 12 squares in `I`, restricted ideal equals the lens's 3 quadrics in both directions, dim 2 mult 8 (`relax/d6even_exactQ.sing`) |
| 8 | d = 6, mixed parity: all 30 odd/off-diagonal coordinates square into `I_Q` | VERIFIED (exact) | `relax/d6mixed_exactQ.sing`: 30/30 |
| 9 | second stage: 150 exact constant duals give 8 exact quadrics in `lambda`; `V` = gauge line `lambda = c(1,..,1)`; with `lambda_0 = 0` the ideal is m-primary (vdim 13) | VERIFIED | rebuilt from scratch at p1 with an unrestricted Koszul form and Y0-monic reduction (`ref_quadrics.py`): quadric space has dim 8 and contains the reductions of the 8 lens quadrics (rank 8), so the spaces coincide; Singular over `Q` on the lens quadrics: one component = gauge line, `(Q, l0)` has dim 0 vdim 13 |
| 10 | THEOREM `deg D >= 7` for every certificate denominator, any Y-parity and any character mixture | VERIFIED | chain 1-9; gauge argument and the Koszul-syzygy step re-derived by hand (see soundness notes) |
| 11 | T2 derivation (E2 mod J gives `w^T Hess(F) w` in `(J:A)_{2d-6}`; gauge `(D,N) -> (D + l'F, N + l'Az/6)`) | VERIFIED (logic) | re-derived; the gauge identity checked by direct expansion; `w -> D` injective for d <= 9 (Koszul syzygies start in degree 5) |
| 12 | NZD lemma: every nonzero `l` in span(X) or span(Y) satisfies `(J:l)_6 = J_6` | VERIFIED | own chart matrices at p1, Singular `quotient(N, freemodule(90)) = (1)` in all four charts (`nzd_*_p1000000009.sing`); the two remaining points `X2`, `Y2` and all six coordinate forms checked directly; chart covering of `P^2(F_p-bar)` and rank semicontinuity give the statement over `Q-bar` |
| 13 | use of the lemma: `served(l) = 0`, sum principle, prescribed-support theorem | VERIFIED (logic) | `served(l) = 0` is immediate from E1 for the denominator `l` (`A` in `(J:l)_6 = J_6`); the (3,1) sum principle follows from the Gysin description in docs/03 (only `c = 1` strata contribute to the pure part; components contracted to codim >= 2 give type (2,2)); docs/01 theorem applies to arbitrary homogeneous `D` ("without a field-parity requirement") |
| 14 | d = 7, sigma-eigen: relaxation EMPTY for Y-even k = 0,3,5,6 and Y-odd k = 1..6 | VERIFIED | own builder at p1: dim 0 in all ten sectors, same vdims (4, 16, 8) |
| 15 | d = 7, Y-even k = 1,2,4: relaxation variety = lifted curve `X_j . C` exactly over `Q` | VERIFIED (exact) | own builder over `Q`: all 6 vanishing coordinates and the 3 C-quadric images square into `I`; the lifted family lies in `V(I)` (`relax/d7_even_k*_exactQ.sing`) |
| 16 | hence no sigma-eigen Y-even denominator of degree 7 in characters 1, 2, 4 | VERIFIED | #13 + #15 + #10; independently, for k = 1, 2 the exact stage-2 duals (`ref_stage2_exact.py`) make the reduced system m-primary at p1 (dim 0 vdim 39), a transfer-safe proof not needing the sum principle; for k = 4 the stage-2 variety is still 2-dimensional, so k = 4 relies on #13 |
| 17 | Hilbert-function sandwich method; verdicts for even k = 1,2,4 | VERIFIED | inequality chain sound; `HS_L = (1+t)^3 (1-t)^9` (C is a CI of three quadrics, degree 8, genus 5) and `HF(S/in I_p)` = `8e - 4` from e = 5 on |
| 18 | d = 7, Y-odd k = 0: relaxation variety inside the lifted family `span(Y) . C` | PLAUSIBLE (modular only) | lens: radical membership mod 1000003 (GB 6358, dim 4). Referee: own leading ideal at p1 has the same dim 4, mult 72; since `deg L = 24` (Segre cone over `C x P^2`), `I_p` is non-reduced along `L`, so the Hilbert sandwich cannot close it; the exact msolve GB was killed (rc -15); exact stage-2 duals (68 quadrics) cut the family down to a curve of degree 64 but do not empty it. Not transferred to `Q-bar` |
| 19 | Y-quadric lemma `(J:q)_6 ∩ (Y-odd, char 3) = J_6^-[3]` for all Y-quadrics | PLAUSIBLE | not replayed (only needed for the pending degree-8 lifts); same chart method as #12 at 1000003 |
| 20a | polars `F_i` and products `Z_j F_i` excluded (hence no poles inside `V(Z_j) ∪ V(F_i)`) | VERIFIED (exact) | `ref_fixedD_exact.py`: 42/42 inconsistent over `Q` by exact `fmpz_mat` ranks (sector-restricted Koszul form); agrees with verify rows 7, 17 |
| 20b | README T4: "for each D listed in `duals/` a rational dual ..." | UNSUPPORTED (as a receipt) | no `duals/` directory exists in the packet; the exact duals live in the verify lens's files and in `kernel_d9_rational.json`-type files of wave 1 |
| 21 | degree-9 kernel representatives / random combinations inconsistent at 1000003 (`jobC_local`) | PLAUSIBLE | modular; only characters 0 and 1 (12 of 132 systems) finished in `jobC_local.out`; verify row 21 gives exact single-monomial duals for 5 of the 8 rational reps |
| 22 | draft: "`deg D >= 8` if the degree-7 analysis closes"; non-eigen d = 7 (60 unknowns) and d = 8 sectors | UNSUPPORTED (pending) | `local_d7full/d7_even_all.gb` is empty; `local_radical8/d8_even_k0.lead` empty; no receipt |
| 23 | Koszul-builder validation on 7 accepted witnesses | PLAUSIBLE | code validation of the lens's builder; not replayed (the referee used an independent builder instead) |
| 24 | re-verification of wave-1 kernel tables and the served-targets mislabel | VERIFIED | #2; verify row 12 |

Counts: 20 VERIFIED (1-17, 20a, 24 counting 16/17 separately), 4 PLAUSIBLE (18, 19, 21, 23), 2 UNSUPPORTED (20b, 22).

## Soundness of the transfer and support arguments (checked by the referee)

* **E1 / Koszul.** `N(F) = A D` with `D = sum w_i F_i` gives `N - A w` in `Syz(F_0..F_5)`; `F_i` is a regular sequence (vdim 15625), so the
  syzygies are Koszul: `N = A w + K(p)`, `p` antisymmetric of degree `d - 4`.  Koszul syzygies start in degree 5, hence `w -> D` is
  injective for `d <= 9` and the gauge slice ("`w_0` has no monomial divisible by `X0`") is a genuine complement of `{l z}`.
* **E2 mod J.** `D div N` lies in `J` because `D` does; `N(D) = sum N_i (sum_j (d_i w_j) F_j + sum_j w_j F_ij) = A w^T Hess(F) w` mod `J`
  (the Koszul part is killed by antisymmetry against the symmetric Hessian).  For mixed-parity `D` the condition splits by Y-parity
  because `J` is Y-graded; the lens's 35-unknown mixed system is exactly that.  So `w^T Hess(F) w` in `(J:A)_{2d-6}` is necessary.
* **Gauge.** `(D, N) -> (D + l'F, N + l'A z/6)`: E1 is preserved exactly (`z(F) = 6F`), and the change of `D div N - N(D)` modulo `F`
  is `6 c D - l' A D = 0` with `c = l'A/6` (uses `deg c = d`).  Coprimality with `F` is preserved.
* **Mod p -> characteristic 0.** (i) `(J:A)_Q` reduces into `(J:A)_p` and `dim_Q <= dim_p` because the degree-`(d+6)` Macaulay matrix of
  `J` has the same rank mod p as over `Q` (CI Hilbert function at 1000003 and 10^9+9); (ii) a `Q-bar`-point of the exact relaxation
  variety, scaled to a `p`-unit vector, reduces to an `F_p-bar`-point of the reduced system, so projective emptiness (ideal m-primary
  mod p) transfers; containment statements do not transfer, which is why the lens replays them over `Q` and why (7, odd, 0) remains
  open; (iii) the second stage uses exact integer duals, so the 8 quadrics are exact and their zero set was decomposed over `Q`.
* **Hilbert sandwich.** `HF(S/I(L)) <= HF(S/rad I_Q) <= HF(S/I_Q) <= HF(S/I_p)`; equality of the two ends in large degree forces
  `rad I_Q = I(L)` in large degree, hence `V(I_Q) = L`.  Correct; it needs `I_p` reduced along `L` (true for the even k = 1,2,4 sectors,
  false for odd k = 0 where the multiplicity is 3).
* **Support.** `served(l) = 0`: a certificate with denominator `l` has E1 `A' l` in `J`, so `A'` in `(J:l)_6 = J_6`, i.e. zero class.
  Sum principle for (3,1)-classes: docs/03's Gysin description gives `N^1 H^4` as Gysin images from the strict transforms of the
  components of an embedded resolution (c = 1 strata; deeper strata have higher weight; components contracted into codimension >= 2
  contribute type (2,2) only), and each strict transform maps into one component of `V(D_1)` or `V(D_2)`.  Prescribed-support theorem
  (docs/01): at `q = 1` a class supported on `V(D)` has a certificate with denominator `D_red` (so with `D` itself), for arbitrary
  homogeneous `D` ("without a field-parity requirement"); powers of `D` never help at `q = 1`.
* **Farkas duals.** `y^T M = 0`, `y^T b != 0` is a valid exact inconsistency certificate; the referee replaced it by exact rank
  comparisons (`rank [M|b] = rank M + 1`) for the 42 polar-type denominators.

## Reconciliation with the verify lens (parallel session)
* Its 21 rows concern the dossier / wave-1 claims.  Its exact confirmations (Jhat, polars, `Z_j F_i`, the degree-9 even kernel being
  exactly 24, the single-monomial functionals for 5 of the 8 rational reps) are cited in the post; its REFUTED/WEAKENED/UNVERIFIED items
  (served-targets run in character 4; only one GF(101) scan; Hessian 3x3 minors never tested) are wave-1 internal claims that do not
  appear in the post.  No contradiction between the lenses: "`deg D >= 9` unless `D` in `J(F)`" (verify) and "`deg D >= 7` for arbitrary
  `D`" (excl) are both stated.  The aggregator's independent finding that the mixed degree-6 relaxation is not m-primary (pure powers
  missing for the `Y0 d/dY0`, `Y1 d/dY1` coefficients) matches the curve `C` inside the weighted-Euler family.
* Changes forced on the draft post: the conditional "`>= 8` if the degree-7 analysis closes" is dropped; (7, odd, 0), non-eigen degree 7,
  degree 8 and degrees 9-10 are reported as screened/open; the missing `duals/` receipt is replaced by the referee's exact replay;
  the k = 1, 2 sectors get an additional transfer-safe proof (exact stage-2 duals, m-primary mod p).
