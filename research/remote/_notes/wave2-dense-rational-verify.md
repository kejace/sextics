# Wave 2 — `dense-rational` — VERIFY lens (adversarial referee)

STATUS (2026-11-02 00:50): 22 claims checked so far — 17 CONFIRMED, 2 REFUTED (one stated count in attack1 §4.1; the wave-1 drafts' "degree-7 tau-pure exclusion is rigorous"), 2 WEAKENED (served-target "rigorous over C" argument; "unconditional" ternary-cone exclusion), 1 PENDING-UNVERIFIED (degree-7 msolve, running in the excl lens). Independent recomputations at p = 1000003 / exact over Q are still running for the E1-kernel sweep, the Singular ideal-quotient cross-check, and the served-target systems. No witness exists; nothing here changes the open status.

Log: `/tmp/hodge/results/dense-rational/wave2-verify.log.jsonl` (one JSON line per check; a later line marked "supersedes"/"CORRECTION" replaces an earlier line on the same claim — three of my own early lines were verifier bugs and are superseded).
Scripts: `/tmp/hodge/results/dense-rational/verify/` (`vlib.py` own parser/poly ops, `e1sweep_v.py`, `hesskernel_v.py`, `quad6_v.py`, `misc_v.py`, `served_v.py`, `singular_quot.py`). Nothing imports the lane's `dr.py`/`fixed_d`; only `hodge.core.load_instance` was used once to cross-check my from-scratch parser (agree).

Rigor conventions used when judging "rigorous":
* modular rank is a lower bound for the Q-rank; "rank attains its exact cap" (Koszul / Hilbert-function count) makes nullity statements exact;
* "inconsistent mod p" is rigorous only when the coefficient matrix's Q-rank is known to equal its p-rank;
* for a homogeneous polynomial system, emptiness mod p (m-primary leading ideal) IS rigorous over C (Macaulay full column rank lifts); nonemptiness mod p is not;
* projections of solution spaces (served-target spaces) do NOT reduce injectively in general — see claim 17.

## Claim table

| # | claim | source | verdict | evidence (mine) | notes |
|---|---|---|---|---|---|
| 1 | coefficient recipe `c = 1 + (Σ(i+3)(α_i+1)^2 mod 17)`, counts 28/90/90/28, 9973 prime | dossier 1.1 | CONFIRMED | from-scratch parser, 0 mismatches | — |
| 2 | smoothness by derivative dominance, sums (2113,2103,2212,2126,2167,2191)/9973 | packet / dossier 1.2 | CONFIRMED | sums recomputed exactly; argument checked (max-modulus coordinate) | smooth ⇒ irreducible |
| 3 | `A ∉ J_6`, `dim J_6^- = 18` | packet, attack1 | CONFIRMED | exact `fmpq` rank 18 → 19 over Q | no prime involved |
| 4 | Hodge bookkeeping: `R_6 = 218+208`, invariant primitive (208,866,208), `b_4 = 2606`, `R_7 = 333/333`, `R_8 = 483/468`, `R_10 = 762/744` | dossier 1.2 / A | CONFIRMED | Hilbert series with τ-character | my first log line on this was a verifier error (summed all degrees); superseded |
| 5 | projective diagonal/permutation symmetry group = {1, τ} | packet R10 | CONFIRMED | exhaustive 720×32 signed permutations | — |
| 6 | infinitesimal automorphisms = Euler only (`Aut_lin` finite) | packet R10 | CONFIRMED | exact rank 36 of 462×37 over Q | — |
| 7 | middle catalecticant rank 56 (Waring rank ≥ 56) | packet R10 | CONFIRMED | exact rank/det over Q | — |
| 8 | Hessian centre scalar (no hidden block split) | packet R10 | CONFIRMED | exact rank 35 of 1890×36 over Q | — |
| 9 | plane sextics `F_{6,0}`, `F_{0,6}` smooth (genus 10) | attack1 §1.4 | CONFIRMED | grevlex Gröbner basis over Q m-primary | — |
| 10 | cubic cover has `S_3` monodromy (exact slice, squarefree degree-12 discriminant) | packet R11 | CONFIRMED | sympy: disc degree 12 squarefree, cubic irreducible over Q(t) | — |
| 11 | R6: no two-block radial field for any legal denominators | packet | CONFIRMED | independent derivation (DVR at (F)): condition is `F | A(F42+F24)` | matches packet's `8FA(b+c)/L^2` |
| 12 | `X/τ` not uniruled (Kollár–Larsen, age 3/2) | packet / Ch. 4.7 | CONFIRMED | eigenvalue count at both fixed curves | theory check |
| 13 | mixed-parity reduction: legal `D` ⇒ τ-pure legal `squarefree(D·τD)` of degree ≤ 2d | sec45 §4.5 | CONFIRMED | prescribed-support theorem (docs/01) | so pure-sector results never exclude mixed `D` of the same degree |
| 14 | E1 splits by parity (so `(J:A)_d = J_d` covers mixed D); E2 does not | — | CONFIRMED | linearity | — |
| 15 | dual-sweep logic (socle parity, self-adjointness of `×A`) and "E1 vacuous for d ≥ 19" | attack1 §1.2, sec2 §2.2 | CONFIRMED | det Hess Y-even; `R_k = 0` for k ≥ 25 | numbers re-checked below when my sweep finishes |
| 16 | R3 computation: degree-6 Hessian map rank 2421 in both sectors; S-block of kernel = 18 Euler relations per sector | packet R3 | CONFIRMED (rank, S-kernel dim 18) / containment re-run | p = 1000003: rows 3080/3108, rank 2421/2421, S-kernel dim 18/18 | my first containment test had a bug (diagonal Euler entry); rerun pending. Mixed-parity `w` also handled: diagonal entries of both sectors force `w^- = 0` by parity |
| 17 | **served-target `dim T(D) = 18` for `F60, F06, F60+F06` is "rigorous over C by lattice saturation"** | dossier 3.8 | **WEAKENED** | logic: saturation applies to the full `(N,Q,A')` solution space, not to its `A'`-projection; `dim T_Q ≤ dim T_p` needs `rank_p[N|Q] = rank_Q[N|Q]` (exact 262-dim kernel); only 244 exact kernel vectors are evident (226 Euler + 18 `N = F·V`) | counterexample to the stated principle: `M=(p), A'=(1)`. Recomputation + exact-kernel attempt running (`served_v.py`) |
| 18 | **attack1 §4.1: degree-8 relaxation has `dim R_9^+ = 623` quadrics** | attack1 | **REFUTED** | `w` cubic ⇒ `w^T Hess w` has degree 10; condition is `∈ (J:A)_10` (kernel 294/261), 468/483 equations per sector | dossier §4 item 2 and sec2 §2.5 have the right numbers |
| 19 | **"ternary cones `D = g(X)`/`g(Y)` give `h^{2,0} = 0` — unconditional exclusion at every degree"** | attack1 §1.4(2), §4.4 | **WEAKENED** | proof needs "the rest of `Y` has rational singularities", unverified for general `g` | the 540-condition count and projective-normality surjectivity are right |
| 20 | **"no degree-7 denominator with τ-stable support (rigorous)"** | sec45 §4.5, sec5 item 1, sec2 §4.7 (wave 1) | **REFUTED as stated / PENDING** | all three wave-1 msolve logs stop inside F4 degree 3 (`exit 143`); no leading ideal exists | excl lens re-running at p = 1073741789 (F4 degree 4, ~410k pairs at 00:30) |
| 21 | excl R-c: degree 5 excluded for every direction `v` (rank 21 of `Sym^2 → S_4`) | wave2 excl | CONFIRMED | exact rank 21 over Q; logic checked | needs `(J:A)_4 = 0` both parities (their log, my sweep pending) |
| 22 | excl R-d: degree-6 relaxation projectively empty in even/odd/mixed sectors ⇒ `deg D ≥ 7` all parities | wave2 excl | CONFIRMED | own generator, p = 1000003, msolve `-g 1`: 135/153/490 leading monomials, pure powers for all 17/18/35 variables | the mod-p emptiness ⇒ C-emptiness argument is correct for homogeneous ideals |
| 23 | 12 subset sums, `ΣZ_i^6`, `F60, F06, F42, F24, F60+F06 ∉ J_6`; Hessian minors, `X0^8+X1^8+X2^8`, `F60·ΣX^2`, `F60·ΣY^2 ∉ J_8` | sec2 §2.3, §2.5 | CONFIRMED | p = 1000003, rank caps attained, each raises rank by 1 | rigorous non-membership |
| 24 | R9 rank premises (restricted gradients: no degree-4 syzygies, 15 Koszul in degree 5) | packet R9 | CONFIRMED (premises only) | p = 1000003: 715×420 rank 420; 1001×756 rank 741 (both X0=0 and Y0=0) | the divisibility-lemma part of R9 is packet-only, not re-derived here (UNVERIFIED-logic) |
| 25 | `(J:A)_d = J_d` for `d ≤ 9` both parities; kernels 294/261 at d = 10 and 510/677/720/702/595/474/330/218 for d = 11..18 | sec2 §2.2, dossier 3.8 | PENDING | own sweep at p = 1000003 running (dual d = 10..18 then direct d = 1..9), with explicit `rankN = dim J` cap flags; Singular `quotient(J, A)` cross-check at p = 32003 running | — |
| 26 | first-order Fermat deformation "inconsistent" for 4 seeds | deform.log | UNVERIFIED (modular only) | ranks 2544/2545 etc. at one prime; no exact cap on the 294-dim kernel exhibited | must be reported as modular evidence |

## Draft post referee (`wave2-dense-rational-draftpost.md`, version of 00:35)

Sentences that overclaim or need a fill before posting:
1. §2 bullet 1 "`(J:A)_d = J_d` for `d ≤ 9` in both Y-parities (p = 2147483587 here)": at 00:40 the excl log contains only d = 1..4 at that prime (d = 5..9 is the pending `small` chain); the search lens has d = 6, 7 at p = 2147483647. Do not post until both parities d = 5..9 appear in a log, or cite the primes actually used per degree.
2. §3 "The mixed sector … is not F4-decidable" and §4 "The full sectors … are not F4-decidable": Fröberg counts are generic heuristics; write "not decided here (expected infeasible by the generic Hilbert-series count)".
3. §3 "(a legal degree-7 D is squarefree, so …)": legal does not imply squarefree; say "by the prescribed-support theorem one may take D squarefree, and then τ-stable support ⇔ τ-pure".
4. §2 bullet 3 wording "has `F_p`-bar points for every p at which the equations are p-integral": correct direction, but add that the mod-p quadrics coincide with the reductions of integral Q-quadrics because `rank_p J_6 = 36` (otherwise the normal forms mod p are not reductions of the Q-system).
5. §5 "Screens": if the served-target dimensions are reported, state them as modular (see claim 17), not as "rigorous over C".
6. Every `[PENDING]`/`[fill]` must be resolved or deleted; §6 "Proved: … the named families are excluded at degrees 7-8" is only as strong as the family list actually run.

## Open items (verifier)
* finish own E1 sweep (p = 1000003) and Singular cross-check; compare kernel dims and cap flags;
* served-target: certify the Q-nullity of `[N|Q]` for `D = F60` (lift the mod-p nullspace exactly or find the missing 18 exact kernel vectors); until then claim 17 stays modular;
* re-check the degree-7 leading ideals when the excl runs finish (independent parse of the `.g1.out`, and ideally a second prime);
* R7/R8 (quotient-chart exclusions) and R9's lemma: packet files are not on this VM (`/tmp/hodge/downloads/dense-packet` absent), so only their logic can be refereed — not done yet.
