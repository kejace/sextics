# Referee report — `dense-rational` (#59), wave-2 exclusion lens — 2026-11-01

Referee: Claude (Anthropic) for `kejace`; no human review; nothing Lean-verified; no certificate claimed by anyone.
Inputs: draft post `wave2-dense-rational-draftpost.md`, report `wave2-dense-rational-excl.md`, packet
`results/dense-rational/wave2-excl/`, verify-lens notes `wave2-dense-rational-verify.md`, dossier, saved thread (#59, one
comment by `fqlx`). Referee replays: `results/dense-rational/wave2-excl/referee/` (log `referee.log.jsonl`); gaps for a
stranger are listed in `results/dense-rational/wave2-excl/README-referee.md`.

Conventions. VERIFIED = receipt located and the statement re-derived by the referee (exact over Q, or mod p with the
rank caps that make the mod-p -> C transfer valid, at a prime not used before). PLAUSIBLE = modular evidence only, or
consistent evidence on record that the referee did not replay. UNSUPPORTED = no finished receipt.

## A. Claim table

| # | claim (draft post / excl report) | receipt | referee check | verdict |
|---|---|---|---|---|
| 1 | `A = X0^3 Y0^3 ∉ J_6`, `dim J_6^- = 18` | log `A-not-in-J6` (p = 2147483587, 18 -> 19) | exact `fmpq_mat` rank 18 -> 19 | VERIFIED |
| 2 | `(J:A)_d = J_d` for `d <= 6`, both parities (T1, T2) | `e1kernel_d{1..6}_par{0,1}_p2147483587.json` | logged numbers imply `rank_p(N) = dim J_{d+6}` at every d; own sweep d = 4, 5, 6 at p = 1000000009 with both caps: kernel 0 | VERIFIED |
| 3 | `(J:A)_d = J_d` for `d = 7, 8`, both parities (T1) | **none at 2147483587** (chain never ran); wave 1 even at 2147483647; aggregator both parities at 1000000007 | own sweep at 1000000009 with both caps: kernel mod J = 0 for d = 7, 8 and both parities; rank_p(N) = 3441, 3441, 5052, 5070 = exact dim J_13^∓, J_14^∓; J_d caps 63, 63, 168, 168 attained (69-221 s each) | VERIFIED |
| 4 | `(J:A)_9 = J_9`, both parities (T1) | none in packet; dossier 3.8 (101, 10007), aggregator (1000000007), wave 1 even (2147483647) | not replayed; three independent primes with `rank N = dim J_15` recorded in the dossier | VERIFIED (by prior independent runs; not re-run here) |
| 5 | degree-10 kernels exactly 294 (even) / 261 (odd) (T1) | none in packet; wave 1 dual sweep (even, 2147483647), dossier (both parities, 101, with `rank N = dim J_16` attained) | not replayed | PLAUSIBLE (two primes agree; cap argument sound) |
| 6 | draft §2 wording "`(J:A)_d = J_d` for `d <= 9` ... (p = 2147483587 here)" | packet has that prime only for d <= 6 | — | UNSUPPORTED as worded; the statement itself is covered by #2-#4 |
| 7 | degree 5 excluded for every direction v (T3): rank `Sym^2 C^6 -> S_4` = 21 | log `degree5-hessian-kernel` (exact) | exact `fmpq_mat` rank 21 of 126 x 21; logic needs #2 at d = 4, 5 (holds) | VERIFIED |
| 8 | relaxation identity `B_D(Aw) ≡ -A w^T Hess(F) w (mod J)` and `B_D(K grad F) ∈ J` | log `relaxation-sign-check` (quadratic w, p = 2147483587: + in J, - not) | hand derivation (section B) + own mod-p check on a random linear field at 1000000009 (rank J_12 = 4437; + stays 4437, - gives 4438) | VERIFIED |
| 9 | degree 6 excluded in even/odd/mixed sectors (T4) => `deg D >= 7` for every certificate, all parities | `quadsys_d6_{even,odd,mixed}_p1073741789.g1.out` + `_ltcheck.json` (135/153/495 leading monomials) | own parser on the agent's outputs (m-primary); own generator + msolve `-g 1` at 1000000009: 135/153/530 leading monomials, pure powers of all 17/18/35 variables; also verify lens (1000003) and aggregator (1000000007) | VERIFIED |
| 10 | degree-7 mixed structured families `weuler` (30 unknowns), `blockdiag` (36), `cross` (36) excluded over C | `fam_d7_mixed_{weuler,blockdiag,cross}_p1073741789.g1.out` (315/562/562) | own parser (m-primary); own generator at 1000000009 incl. exact gauge-in-family (6/0/0): 315/562/562 leading monomials, all variables | VERIFIED |
| 11 | degree-7 mixed family `gradient` (56 unknowns) excluded over C | `fam_d7_mixed_gradient_p1073741789.g1.out` (5360) | own parser (m-primary); own generator at 1000000009: 5267 leading monomials, pure powers of all 56 variables (276 s, 4 threads) | VERIFIED |
| 12 | degree 7, τ-pure sectors excluded (T5) => `deg D >= 8` for τ-stable supports | `quadsys_d7_{even,odd}_p1073741789.ms` only; both `.g1.out` empty (even job killed in F4 degree 4; odd still in that step after 2.5 h, 39 GB) | nothing to replay | UNSUPPORTED (pending) |
| 13 | degree-7 mixed families `xonly`/`yonly`, 15 pairs, 20 triples (T6) | `xonly` msolve unfinished (0-byte output); others never ran | nothing to replay | UNSUPPORTED |
| 14 | degree-8 structured families, even/odd (T6) | never ran (needs the d = 10 kernel basis at 2147483587, also never computed) | nothing to replay | UNSUPPORTED |
| 15 | T7 "serves nothing primitive by E1 alone" + sum principle | never ran | logic sound if the `J_{d+6}` rank cap is checked (section B) | UNSUPPORTED (not run) |
| 16 | §5 screens (served-target dims, Fermat deformation, E2 on degree-10 representatives, LM) | never ran in wave 2 | — | UNSUPPORTED (drop from the post) |
| 17 | "mixed sector / full d = 8 sectors are not F4-decidable" | Fröberg-count heuristic | heuristic, not a theorem | wording: "not decided here" |
| 18 | "a legal degree-7 D is squarefree" | docs/01 prescribed-support theorem | the theorem gives: legal D => `D_red` legal (q = 1, arbitrary homogeneous D by its last sentence); so "WLOG squarefree", then τ-stable support <=> τ-pure | wording corrected; logic VERIFIED |

Counts: VERIFIED 9 (#1, 2, 3, 4, 7, 8, 9, 10, 11); PLAUSIBLE 1 (#5);
UNSUPPORTED 7 (#6 as worded, #12, 13, 14, 15, 16 and the wave-1 "degree-7 τ-pure is rigorous" statement that the excl
report already retracts); wording fixes 2 (#17, 18).

## B. Soundness of the transfer arguments (as used here)

* **mod p -> characteristic 0 for E1 (saturated lattice).** `V = {D : A D ∈ J_{d+6}}` over Q; `V ∩ Z^n` is saturated so it
  reduces injectively mod p. Its image lies in the mod-p solution space **only if** the generator matrix of `J_{d+6}`
  keeps its Q-rank mod p (otherwise `A D = Σ N_i F_i` may need `p` in the denominators of `N`; counterexample `(p)` vs
  `(1)`). With that cap and `rank_p(J_d) = dim J_d`, `dim_p K_full = dim J_d` proves `(J:A)_d = J_d` over Q, hence over C.
  The agent's script checks only the second cap; the first holds numerically in every logged line (README-referee §1) and
  is explicit in the referee's sweep. Same remark for T7 (served-target test by E1): sound with the cap.
* **m-primary leading ideal => no projective zero over C.** Homogeneous `I ⊂ Z_(p)[c]`; a nonzero complex zero gives a
  nonzero `Q̄`-zero, scale it `𝔭`-integral with a unit coordinate for a prime `𝔭 | p`, reduce: a nonzero `F̄_p`-zero of
  `I mod p`. A DRL leading ideal containing a pure power of every variable makes `F_p[c]/I` finite, so no such zero.
  Correct. The mod-p quadrics must be reductions of p-integral rational quadrics equivalent to the rational relaxation:
  true because the normal form uses a pivot minor of the `J_{2d-6}` generator matrix that is a p-unit (rank attained:
  18/36 at degree 6, 336 at degree 8, checked in both the agent's and the referee's generators).
* **Hessian relaxation.** For `D = w(F)` and `N = A w + K grad F` (`K` skew), `B_D(N) = D(w(A) + A div w) - A (w·∇w)(F)
  - A w^T Hess(F) w + (terms in J)`, so E2 (`B_D(N) = F C ∈ J`) forces `A w^T Hess(F) w ∈ J_{2d}`; injectivity of `[A]` on
  `R_{2d-6}` (d = 6, 7) gives `w^T Hess w ∈ J_{2d-6}`; at d = 8 the condition is `w^T Hess w ∈ (J:A)_10^+` (both pure
  sectors give a Y-even form). Gauge `w -> w + l·Euler` adds `10 l D + 30 l^2 F ∈ J`; a legal D is not a gauge element
  (`D = 6 l F`), and `w(F) = 0` forces `w = 0` for `deg w <= 4`. Sector projection: for a Y-even (odd) legal D the
  parity projection of w still gives D and differs from w by a syzygy, so pure-sector systems cover pure-parity D and the
  mixed system covers the rest. Correct.
* **Family restriction.** Exclusion of `w ∈ W_fam ⊗ C` needs `W_fam ∩ gauge` of the same dimension over Q and mod p
  (else the gauge fixing drops a direction that is not gauge over Q). Verified exactly for the four finished families.
* **Prescribed support (q = 1).** docs/01: for arbitrary homogeneous D, a certificate with denominator D exists iff one
  with `D_red` exists. So degree bounds may assume D squarefree; "legal => squarefree" is not literally true.
* **Sum principle `served(D1 D2) = served(D1) + served(D2)`** (T7, not run): sound for coprime `D1, D2` — the supported
  classes are images of Borel–Moore homology `H_4(Z1 ∪ Z2)`, and `H_4(Z1) ⊕ H_4(Z2) -> H_4(Z1 ∪ Z2)` is onto because
  `H_3^{BM}(Z1 ∩ Z2) = 0` for the curve `Z1 ∩ Z2`.
* **Farkas duals** are not used in this lane's wave-2 material; nothing to check.

## C. Agreement with the verify lens (`wave2-dense-rational-verify.md`)

Agree on every item that overlaps: claim 17 (served-target "rigorous over C" is only modular as argued in the dossier),
claim 20 (wave-1 "degree-7 τ-pure exclusion is rigorous" is unproven), claims 21-22 (degree 5 and degree 6 closures),
the six wording fixes for the draft. Additions here: the implicit `J_{d+6}` rank cap in the E1 sweep (holds), the exact
gauge-in-family check (holds), the `xonly` output is empty (not a certificate), the rig chain never started, and the
degree-7 even msolve job was killed (not merely slow). No disagreement.

## D. What the post may say

Proved (all parities, arbitrary Q, over C): `deg D >= 7`; `(J:A)_d = J_d` for `d <= 9` (so every E1-admissible
denominator of degree `<= 9` is a polar divisor `w(F)`); no degree-7 denominator `D = w(F) + lF` with `w` in the four
named mixed-parity families. Not proved: anything at degree 7 outside those families (including both τ-pure sectors),
anything at degree `>= 8`, any served-target statement, any screen. No claim of nonexistence.
