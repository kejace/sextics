# Wave 2 verification — lane `reflection-coupled-residual` (adversarial referee)

STATUS: IN PROGRESS — confirmed 9 / refuted 1 (typo-level) / weakened 0 / unverified (pending) — see table; last update see bottom.

Scope: every claim labelled rigorous/proven/exact/computed in `dossier-reflection-coupled-residual.md`,
`attack1-reflection-coupled-residual.md`, the lane logs under `/tmp/hodge/results/reflection-coupled-residual/`
(no `wave2-reflection-coupled-residual-*.md` and no draft post existed at 00:30 UTC; re-checked periodically).
Independent code (no engine dependency) lives in `/tmp/hodge/results/reflection-coupled-residual/verify/`
(`vlib.py`, `v_symmetry.py`, `v_sectors.py`, `v_e1.py`, Singular scripts `sing_*.sing`); every check is appended
to `/tmp/hodge/results/reflection-coupled-residual/wave2-verify.log.jsonl`.
Primes used by the lane: 2147483647, 2147483629 (sibling, python-flint), 32003 (dossier, NumPy). Primes used here:
1000000007 (own flint code), 1000003 (Singular), plus exact char-0 Singular where stated.

## Logic conventions checked

* **Modular-rank direction.** For an integer matrix, `rank_p <= rank_Q`. The E1-admissible space is
  `K = {D : A·D ∈ J}` with `dim K = dim S_d − rank[J;A·S_d] + rank J`, so `dim K_Q <= dim K_p` holds **only if**
  `rank_p J_{d+6} = rank_Q J_{d+6}`. That equality is certified in every run below by comparing `dim (S/J)` mod p with the
  Fermat count (all exponents <= 4) in both degrees (`hilbert_ok`). With it, `ker_p = 0 ⇒ ker_Q = 0` (zeros are
  rigorous) and `ker_p = 1 ⇒ ker_Q <= 1` (an upper bound; the exact rational element must be verified separately).
* **"E1 false mod p" is rigorous** under the same Hilbert-function check. **"E1 true mod p" is not** a proof of E1
  over Q; K8's exact membership is checked separately (Singular, char 0).
* **"Inconsistent mod p" for a fixed-D (N,Q) system is NOT automatically rigorous** (needs `rank_p = rank_Q` of the
  homogeneous system, whose kernel depends on D). The dossier labels those results [reported]/[computed, mod p] and
  never calls them rigorous — correct. The E2-mod-J obstruction system, by contrast, has a kernel of structural
  (Koszul) dimension, so a modular inconsistency there does transfer (see claim table).
* **Character restriction is lossless for E1** (J and A are Klein-homogeneous, so `A·D ∈ J` iff each character
  component of D satisfies it); restricting to Y-even D is lossless for *existence* (replace D by D·τ(D)) but not for
  degree bounds (a Y-even exclusion up to degree k only excludes arbitrary D up to degree k/2 via the product; the
  all-character computation is what gives the unrestricted bounds).

## Claim table

| # | claim | source | verdict | evidence | notes |
|---|---|---|---|---|---|
| 1 | Projective diagonal symmetry group of F is the Klein group {1, r, ρ, τ}; no permutation symmetry; F(Y,X) ≠ F; no μ3/μ6 action | dossier §1.3, attack1 §1.1 | CONFIRMED | `v_symmetry.py`: exhaustive μ6 weight search (Z_i^6 ∈ F forces z_i ∈ μ6 after normalising z_0 = 1) gives exactly 4 weight vectors; 720 permutations: identity only | — |
| 2 | Residue characters: A = X0³Y0³ is τ-invariant, r-anti, ρ-anti; companion X0³Y2³ is (+,+,+) | dossier §1.3 | CONFIRMED | `v_symmetry.py` (numerator sign × det) | formula direction checked: Ω ↦ det(g)Ω |
| 3 | Sector Hodge numbers: target (126,492,126); companion (82,374,82); τ-odd (144,529,144), (74,356,74); totals 426/1751/426; τ-inv total 208, τ-anti 218 | dossier §1.3, App. A | CONFIRMED | `v_sectors.py` equivariant Hilbert series Π(1−χ_i t⁵)/(1−χ_i t); equals Fermat monomial counts in every character | valid for every smooth Klein-symmetric F |
| 4 | Character counts: R_8^(+,+)=301 → R_14^(+,−)=426; R_8^(−,−)=182 → R_14^(−,+)=318; first forced kernel d=9: 381 → 360 (≥21); dim R_16^odd = 468, 274 in (+,−); 112 invariant-sector unknowns (= dim J_8^(+,+)) | attack1 §1.3–1.4 | CONFIRMED | `v_sectors.py` | — |
| 5 | "R_6^odd has dim 866 = 208 (Y2-even) + 658" | attack1 §1.3 | **REFUTED (as stated; typo-level)** | `v_sectors.py`: R_6^odd = 208 = 126 + 82; R_12^odd = 866 = 492 (Y2-even) + 374 (Y2-odd); 658 matches nothing | harmless: no downstream use |
| 6 | Smoothness: 5/7 ≡ 0 mod 5 so F ≡ Fermat mod 5; good reduction ⇒ smooth over Q, C; A ∉ J | dossier §1.5, checks.json | CONFIRMED | 5·7⁻¹ ≡ 0 (mod 5); A ∉ J over Q follows from A ∉ J mod 5 because rank_Q J_6 = rank_5 J_6 = 36 (no syzygies below degree 10) | — |
| 7 | 138 Y-odd/Y2-even sextic monomials; 158 Y-even/Y2-even (= one-reflection coefficient dimension) | dossier App. A | CONFIRMED | `v_sectors.py` | — |
| 8 | E1-admissible space: zero in every character for d ≤ 4; equal to J_d for d = 5, 6 in every character (3,2,1,0 / 14,12,6,4) | dossier §3.1, App. A (p = 32003) | CONFIRMED (new prime 1000000007, own code) | `v_e1.py`, `e1_p1000000007.jsonl`; `hilbert_ok = true` in every cell so zeros are rigorous over Q and C | — |
| 9 | Y-even E1 kernel mod J = 0 for d = 1..7 (dims 3,12,28,66,123,218,333) and 1 at d = 8 (dim 483) | attack1 §1.3, sweep_e1.out, k8_d8.out (primes 2³¹−1, 2³¹−19) | CONFIRMED (3rd prime) and **EXTENDED**: at d = 7 the two Y-odd characters also have kernel 0 (dims 135, 120), so "deg D ≤ 7 ⇒ D ∈ J(F)" now holds for every character, not just Y-even | `v_e1.py` p = 1000000007: d=7 kernels (0,0,0,0), d=8 (0,0)-sector kernel 1 (dimK 113 = 112 + 1) | the d = 8 statement is an upper bound ker_Q ≤ 1; exactness of K8 is row 11 |
| 10 | K8 is irreducible over Q; r-even and ρ-even; gcd(K8, F) = 1 | attack1 §1.3 | CONFIRMED | Singular `factorize` in char 0: one factor, multiplicity 1; `gcd(K8,F) = 1`; all 37 monomials have character (0,0) | absolute irreducibility not claimed and not needed for coprimality (F is absolutely irreducible, companion Lean) |
| 11 | K8 (rationally reconstructed from two primes) is an exact E1 class over Q: A·K8 ∈ J(F) | attack1 §1.3 (implicit; only mod-p evidence existed) | PENDING → see update below | A·K8 ∈ J mod 1000000007 confirmed (`v_k8props.py`); exact char-0 Singular reduction running (`sing_exactK8.sing`) | until this is done, "kernel = 1 at d = 8" is only "≤ 1" over Q |
| 12 | K8 is not proportional mod J_8 to H_u², Y2⁸, X0⁴Y0⁴, X0⁸ | attack1 §1.3 | CONFIRMED | `v_k8props.py` p = 1000000007 | — |
| 13 | K8 is not proportional mod J_8 to Y2⁴H_u | attack1 §1.3 | WEAKENED (vacuous) | Y2⁴H_u = Y2³·(∂F/∂Y2)/2 ∈ J_8, its class is zero | the sentence only restates K8 ∉ J |
| 14 | E1-admissible degree-8 Y-even denominators form the 168-dim affine family K8 + Σ g_i F_i, g_i ∈ S_3 | attack1 §1.3 | CONFIRMED | dim J_8^{Y-even} = 168 = 112 + 56 (`v_sectors.py`, `v_shapes.py`), no syzygies below degree 10 | — |
| 15 | No member of that family has any of the 55 listed shapes (divisible by a variable, cone with vertex a line, few-variable, bounded X/Y-degree, …) | attack1 §2.3 (repr.py, one prime) | CONFIRMED mod p at a 2nd prime; **rigorous over Q only for 20 of the 55** (those with rank T = 168 = #columns: all 14 "cone" shapes except (X0,Y2), "only even exponents", (X0,Y0)⁴, (X0,Y0)⁶); the other 35 are mod-p statements (rank_Q T = rank_p T not certified) | `v_shapes.py`, `shapes_p1000000007.json` | the note states them flatly; should say "mod p" for the 35 |
| 16 | ob(D) := [D·div N − N(D)] ∈ R_{2d}^odd is independent of the E1 solution N and vanishes for every certificate; ob(D+T) = ob(D) + [A(D div g − 2g(D))] − [A g(g(F))] | attack1 §1.4 | CONFIRMED (derivation checked) | N − N' = aF (Koszul, F_i regular); div(aF) = Σ_j(Σ_i ∂_i a_ij)F_j by antisymmetry + Hessian symmetry; (aF)(D) ∈ J; (F) ⊂ J by Euler; the affine/quadratic formula re-derived term by term | E1&E2 ⇔ certificate uses F irreducible (companion Lean) and gcd(D,F)=1 |
| 17 | "The companion target X0³Y2³ has character (+1 under r2, −1 under s)" | attack1 §1.1 | **REFUTED (s-part; typo-level)** | `v_symmetry.py`: X0³Y2³ is s-invariant and det(s) = +1 ⇒ residue character +1 under s; dossier §1.3 has (+,+,+) correctly | irrelevant for the descent to X/r, which only needs the r-character |
| 18 | Fermat control (c = 0): Φ12(X0,Y0), X0⁴+Y0⁴, X0²Y0², X0²±Y0² pass E1; X0Y0 and X0⁴+X1⁴+X2⁴+Y2⁴ fail E1 | attack1 §2.2 | CONFIRMED (inspection: J_Fermat = (Z_i⁵) is monomial) | — | the "inconsistent" verdicts for the passing ones are mod-p engine results (not recomputed); geometrically expected since those supports pair X0 with Y0 (a_0 + a_3 = 8 ≢ 0 mod 6 in Shioda's inductive structure) |
| 19 | 42 degree-4 denominators fail E1; 61 degree-6 DiamondJM supports fail | attack1 §2.1, dossier §3.1 | CONFIRMED by implication (deg ≤ 4) / UNVERIFIED [reported] (DiamondJM) | row 8 | the dossier labels these [reported], correctly |
| 20 | dim J_d = 6·dim S_{d−5} for d ≤ 9 (no Koszul syzygies below degree 10) | dossier App. A | CONFIRMED | `v_e1.py`: rank J_7 = 45+40+23+18 = 126 = 6·21; rank J_8 = 112+104+64+56 = 336 = 6·56 | — |
| 21 | Linear-system size table (§1.2), deg Q = 2d − 6, N_i ∈ S_{d+1}; §4.1 parameter counts (d=5: c 6, K 90, m 6, Q 126) | dossier §1.2, §4.1 | CONFIRMED | binomial counts rechecked by hand; deg N_i = d + 6q − 5 per contract | — |
| 22 | Ages: r junior (1/2, Fano quotient W); ρ age 1 along the smooth Fermat sextic surface X ∩ {Y0=Y1=0}; τ age 3/2 along two plane sextic curves ⇒ X/ρ, X/τ not uniruled (Kollár–Larsen) | dossier §1.4 | CONFIRMED (arithmetic of ages; criterion cited from book §4.7) | X is CY (K_X = O), so Kollár–Larsen applies | consistent with checks.json `juniorElementsFixingTarget = []` |
| 23 | W = {H_c = 0} ⊂ P(2,1,1,1,1,1) Fano with −K_W = O(1); B_W = {H_c(x,0)=0} ⊂ P⁴ smooth | dossier §1.4 | CONFIRMED | weights sum 7, degree 6; B_W ≡ Fermat sextic threefold mod 5 (c ≡ 0) ⇒ smooth by good reduction | — |
