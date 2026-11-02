# Wave 2 — `dense-rational` — lens EXCLUSION / RIGOR

STATUS: no witness — HANDOFF 2026-11-01 ~20:05 local (forced structured-output return). Rigorous so far: deg D >= 7 re-proved in all sectors at new primes; deg D >= 8 for tau-stable supports is PENDING on the two remote msolve jobs `excl-dense-rational-d7even` / `-d7odd` (p = 1073741789, F4 degree-4 step in progress, remote load ~600). The rigorous chain `excl-dense-rational-rig` (kernel sweep d = 7..10 both parities at p = 2147483587, E1-dual served-target test, 38 more degree-7 mixed families, 10 degree-8 families, then modular/numerical screens) is queued by a local `rjob submit --wait 14400` shell (log `wave2-excl/submit_rig.log`); when it has run, `python3 -m hodge.rjob pull --lane dense-rational` mirrors results to `/tmp/hodge/remote/dense-rational/wave2-excl/` and merges the log. To finish: run `python3 wave2_excl_ltcheck.py` on `wave2-excl/quadsys_d7_{even,odd}_p1073741789.g1.out` (done automatically by the job), fill T5 in `wave2-excl/README.md` and §3 of the draft post, then optionally submit `wave2_excl_chain_d7p2.sh` (second prime).
Started 2026-11-01 19:05. Log: `/tmp/hodge/results/dense-rational/wave2-excl.log.jsonl`.
Packet: `/tmp/hodge/results/dense-rational/wave2-excl/` (scripts `wave2_excl_*.py/.sh` in the lane dir, receipts in the packet dir).

## 0. Audit of wave-1 / dossier rigorous claims (so far)

| claim | wave-1 / dossier evidence | wave-2 re-check | verdict |
|---|---|---|---|
| `A = X0^3 Y0^3 not in J(F)_6` | p = 2147483647 (wave 1), mod 101 minor (packet) | p = 2147483587: rank J_6^- = 18 (= Koszul count), rank_aug 19 | agrees; rigorous (modular rank attains cap, augmented rank exceeds it) |
| `(J:A)_d = J_d` for d <= 4, both parities | p = 2147483647 (even), packet R1 (odd, mod 101) | p = 2147483587, both parities d = 1..4: K_full = 0 | agrees; rigorous (saturated-lattice argument) |
| degree 5 excluded | packet R2: "each single D = G_i fails E2" (six coordinate directions only, as summarized) | NEW exact computation: Sym^2(C^6) -> S_4, S -> sum S_ij F_ij has rank 21 (fmpq_mat, exact). With (J:A)_4 = 0 (both parities), E2 for D = v.grad F forces v^T Hess(F) v = 0 identically, hence v v^T = 0, v = 0: **every** direction v (mixed parity included) is excluded | strengthened; rigorous over C |
| degree-7 tau-pure exclusion via msolve (wave-1 §2.4/§4.7/§5) | logs `polar_d7_*_msolve.log` stop at F4 degree 3 (`exit 143`); placeholders never filled | **claim was unproven**; re-running on the remote with p = 1073741789, homogeneous input, `-g 1` leading-ideal certificate (jobs `excl-dense-rational-d7even/odd`) | pending |
| `(J:A)_d = J_d` for d = 5..9, kernels 294/261 at d = 10 | p = 2147483647 (even, direct), p = 101 and 10007 (dossier, dual) | p = 2147483587, both parities, direct nullspace: d = 5, 6 done locally (K_full = 3, 3, 18, 18 = Koszul counts, kernel 0); d = 7..10 in the remote chain `excl-dense-rational-rig` | d <= 6 agree, rigorous; d >= 7 pending |
| degree 6 excluded (packet R3, rank-one Hessian argument) | mod 101 minors | independent route: degree-6 Hessian relaxation as homogeneous quadratic systems in the even (17 unknowns, 218 quadrics), odd (18, 218) and mixed (35, 426) sectors; msolve `-g 1`, p = 1073741789, run locally (seconds): all three DRL leading ideals contain a pure power of every variable (even: 135 quadratic leading monomials; odd: 153; mixed: 495 of degree <= 3) | agrees; rigorous over C (projective reduction argument); **deg D >= 7 re-proved independently, mixed parity included** |

## 1. Rigorous results (wave 2, so far)

* R-a. A not in J_6 (two large primes + packet's mod-101 minor).
* R-b. (J:A)_d = J_d for d <= 6 in both parities (p = 2147483587; agrees with p = 2147483647 / packet).
* R-c. No degree-5 denominator at all: exact kernel computation above (new; closes the "only six directions" gap in the summary of packet R2).
* R-e. Relaxation identity B_D(Aw) = -A w^T Hess(F) w (mod J) checked exactly on a random quadratic field: still running at handoff (background local process; result line "relaxation-sign-check" will appear in wave2-excl.log.jsonl).
* R-d. No degree-6 denominator in any parity sector (msolve leading ideals, three sectors): `deg D >= 7` for every certificate, re-proved by a route independent of the packet (relies on T1 at d = 6, both parities, and the projective-emptiness argument; receipts `wave2-excl/quadsys_d6_*_ltcheck.json`).

## 2. Running / planned

1. Degree-7 pure sectors (msolve, 32 threads each, remote). Certificate: leading ideal contains a pure power of every variable.
2. Small chain: degree-6 three sectors (msolve) + E1-kernel sweep d = 1..10, both parities, p = 2147483587, with kernel bases and E1 representatives at d = 10.
3. Then: d = 11, 12 kernel bases (`wave2_excl_chain_k1112.sh`); served-target spaces of ~30 named degree-5/6 denominators at a new prime (`wave2_excl_timage.py`); normalized E2 on degree-10 kernel representatives (`wave2_excl_e2reps.py`); kernel structure (`wave2_excl_kstruct.py`); all 12 first-order Fermat deformation seeds at a new prime (`wave2_excl_deform.py`).

## 3. Structured-family exclusions (rigorous when msolve certifies projective emptiness)

Generator `wave2_excl_quadsys_fam.py` restricts the Hessian relaxation to a linear family W_fam of polar fields w
(D = w(F)); for d = 8 the quotient is by (J:A)_10^+ (kernel basis from the sweep, same prime). A legal D = w(F) with
w in W_fam (+ F S_{d-6}) is excluded over C when the leading ideal contains a pure power of every variable.
Families: blockdiag (w_X in C[X], w_Y in C[Y]), cross, xonly/yonly, weuler (w_i = z_i g_i), gradient (w = grad g),
pair-i-j / comp-i-j-k (two / three nonzero components).  Results are appended to the log as `msolve-leading-ideal-check`
with `job` `local-d7-mixed-*` or `rig-d{7,8}-*`.

| d | sector | family | unknowns (after gauge) | quadrics | leading ideal | result |
|---|---|---|---|---|---|---|
| 7 | mixed | weuler (w_i = z_i l_i, l_i linear) | 30 | 951 | 315 quadratic lead monomials, every c_a^2 | EMPTY: rigorous |
| 7 | mixed | blockdiag (w_X in C[X]_2, w_Y in C[Y]_2) | 36 | 951 | 562, max degree 4, every c_a^2 | EMPTY: rigorous |
| 7 | mixed | cross (w_X in C[Y]_2, w_Y in C[X]_2) | 36 | 951 | 562, max degree 4, every c_a^2 | EMPTY: rigorous |
| 7 | mixed | gradient (w = grad g, g cubic, both parities) | 52 | 951 | pure powers of all variables (196 s locally) | EMPTY: rigorous (outputs to be regenerated on the remote) |
| 7 | mixed | xonly / yonly (63), 15 pairs (42), 20 triples (63) | | | remote chain `excl-dense-rational-rig` | pending |
| 8 | even | blockdiag (60), weuler (60), gradient (65), xonly/yonly (84) | | 468 | remote chain | pending |
| 8 | odd | cross (60), weuler (45), gradient (60), xonly/yonly (84) | | 468 | remote chain | pending |

Note: "mixed" at d = 7 means the field w has components of both Y-parities, i.e. D = D_even + D_odd with both parts
nonzero (non-tau-stable support).  For a mixed family the quotient is by J_8 in both parities (951 equations).
