# Referee report — `short-power-presentation`, wave 2, exclusion/rigor packet

Referee: Claude (Anthropic) for `kejace`, 2026-11-01/02. No human review. Nothing Lean-verified. No certificate claimed.
Scope: every claim in `wave2-short-power-presentation-{excl,draftpost,verify,search,construct}.md` and this packet's
`receipts.jsonl`, reconciled with the parallel verify lens and the aggregator (`wave2-status.md` §2.15).
Own receipts: `referee/referee.jsonl`, `referee/*.out` (jobs `referee-spp-*`, also appended to the lane log
`/tmp/hodge/results/short-power-presentation/wave2-excl.log.jsonl`). Fresh prime 1073741741 (`PRIMES30[2]`, unused before).

Target (fixed): `F = e(X)+e(Y)+c((L+M)^6+(L-M)^6)`, `c = 5/7^7`, `L = X0+2X1+3X2`, `M = Y0+3Y1+5Y2`, `A = X0^3 Y0^3`, `q = 1`,
variables `X0,X1,X2,Y0,Y1,Y2`, `tau` = global Y-sign, `K_d^± = {D in S_d^± : A D in J(F)}`.

## 1. Replays performed by the referee (all under `ulimit -v 6 GB`, one process at a time)

| # | replay | result | receipt |
|---|---|---|---|
| R1 | E1 kernels `d = 5,6,7`, both parities, `p = 1073741741`, own-code + engine `mult_kernel`, Hilbert check (R1) | `K_5/J = K_6/J = 0`, `K_7^±/J = 1`; every `rank J_{d+6}` and `dim J_d` at the Hilbert value | `referee/e1k_p1073741741.out` |
| R2 | exact `D7^±` certificates from `excl_kernel_D7.json` (remote job `exclB-spp-kernel2`): `A*D = sum N_i dF/dZ_i` over Q with `fmpq_mpoly` | both identities exact; `D7^+ in Q[X0,X1,X2,M]` (killed by `d/dY1-3d/dY0`, `d/dY2-5d/dY0`); `D7^± not in J_7` mod `p` with `rank_p J_7 = 63` = Hilbert | `referee/referee.jsonl` (`referee-D7-exact`) |
| R3 | exact lifted certificates `L*D7^+ in J_8`, `M*D7^- in J_8` (`exact_membership`, verified with `Fractions`) | both exact over Q | `referee-exact-membership-J8` |
| R4 | search-lens odd cone representative `D7odd_cone.json` | rational, in `Q[Y0,Y1,Y2,L]`, and `D_cone + D7^- in J_7` exactly over Q (lifted certificate) | `referee-D7odd-cone` |
| R5 | mixed `d = 6` annihilator-functional systems (`excl_mixed.py 6 --pidx 2`), both sides, plus the exact dependency step | y-side: 36 quadrics / 18 unknowns, Singular dim 0 vdim 4096; x-side: 18 / 12, vdim 729; **exact dependency certificates 114/114 on both sides** (`dim_Q L^+ mod J = 39`, `dim_Q L^- mod J = 57`, lifted with 2 primes, verified with `fmpq_mpoly`) | `referee/mixed6_p1073741741.out`; lane log `mixed-annihilator-exact-certificate` ×2 (job `referee-spp-mixed6`) |
| R6 | Euler-sliced relaxations `d = 5` (e/o/m) at `p = 1073741741`, msolve + Singular | all proj-empty: Hilbert functions `[1,3,2,0]`, `[1,3,2,0]`, `[1,6,13,8,0]`; `rank_p J_10^-` = 747 = Hilbert | `referee/relax_p1073741741.out` |
| R7 | Euler-sliced relaxations `d = 6` (e, o) at `p = 1073741741`, single thread | even: proj-empty, Hilbert function `[1,17,114,358,501,225,0]`, Singular vdim 1216 (178 s); odd: proj-empty, `[1,18,114,272,129,6,0]`, vdim 540 (68 s); `rank_p J_12^-` = 2214 = Hilbert; identical to the `1073741789` receipts | same file |
| R8 | `A not in J_6` at `p = 1073741741`; Hilbert reference values `R_d^±` from `(1+t+..+t^4)^6`, `Theta = (1+t^2+..+t^8)^3` | augmented rank 19 vs 18 = Hilbert; library values agree | `referee-A-notin-J`, `referee-hilbert-values` |

## 2. Transfer arguments (mod p -> characteristic 0), checked line by line

* **E1 kernels.** `dim K_d = dim S_d - rank[J_{d+6} | A S_d] + rank J_{d+6}`. `rank_p <= rank_Q` on both matrices gives no
  inequality by itself; it does once `rank_p J_{d+6} = rank_Q J_{d+6}` (Hilbert value, `F_i` regular sequence) and
  `dim_p J_d = dim_Q J_d`: then `dim_Q K_d/J_d <= dim_p K_d/J_d`. Both equalities are recorded (`hilbert_ok`). Lower bounds
  need exact elements: provided at `d = 7` by R2. (Agrees with verify-lens row 4.)
* **Relaxation (R3).** `A*D in J` and `D = w(F)` force `A*w^T Hess(F) w in J_{2d}` (re-derived: with `N(F) = A D` exactly,
  `N - A w` is a Koszul syzygy, and `D div K`, `K(D)`, `D div(Aw)`, `A D div w` all lie in `J`; `E2` mod `F` implies mod `J`).
  The quadrics are `phi(Xi(w))`, `phi` in `Ann(J_{2d})`; when `rank_p J_{2d}` is the Hilbert value the mod-`p` annihilator is
  the reduction of the saturated integer one, so the mod-`p` ideal is the reduction of the rational ideal and
  `h_{S/I_p}(t) >= h_{S/I_Q}(t)`; a zero of the mod-`p` Hilbert function therefore proves `V(I_Q) = {0}` in the Euler
  slice, i.e. `w = G*E`, `D = 6 G F`, not coprime. The slice `{w_{0,m} = 0 : X0 | m}` is a transversal of the gauge.
  All `d <= 6` receipts record `rank_p J_{2d}` at the Hilbert value (`parts`).
* **Mixed parity (R4).** Odd part of `Xi` has no cross terms (checked: `w^{+T} H w^-` is Y-odd, so `A*` of it is Y-even).
  Functionals killing `J_12^- + L^+` are reductions of rational ones iff `dim_Q (L^+ + J)/J = dim_p`: this needs the exact
  dependency certificates, which were **absent from every log** (local and remote) when the packet was closed — the remote
  job `exclB-spp-relax-mixed` skipped the exact step. R5 now supplies them (114/114 on each side). With them, `y = 0` over
  `Q-bar` for every mixed solution, then `x = 0` by the pure even case. Independently, the aggregator's full 35-unknown
  mixed system has an m-primary leading ideal at `10^9+7` with the Hilbert check (`hessrelax6_mixed_mprimary`).
* **Block elimination (R5 of the excl note).** Logic is sound, but each step's functionals need the same exact-rank
  certificate; the certification job (`exclB-spp-certify7e`) ended `INCOMPLETE` (rational lift failed with 12 primes) for
  both steps it attempted. All degree-7 block statements are therefore single-prime modular evidence.
* **Parity coverage.** `J` is `tau`-stable, so `A D in J` splits into `A D^± in J`; for `d <= 6` this gives `D = w(F)` with
  `w = w^+ + w^-`, and the relaxations cover pure and mixed `w`. The only symmetry is `tau` (verify rows 15–16), so no further
  character splitting exists and no symmetry is imposed on `D`. Coverage through degree 6 is complete.

## 3. Claim ledger

**VERIFIED** (replayed here, or complete receipt chain):
V1 the sparse form of `F`, `A`, `q` (two independent exact receipts; not re-parsed here) · V2 `A not in J` (R8) ·
V3 `K_d/J_d = 0`, `d <= 6`, both parities (R1 for `d = 5,6`; receipts at `2^31-1`, `10^9+7`, `1073741789`, `2147483629`
for `d <= 6`; `d <= 4` is `J_d = 0`, `K_d = 0`) · V4 `dim_Q K_7^±/J_7^± = 1` exactly (R1 + R2) · V5 `D7^+ in Q[X0,X1,X2,M]`;
the odd class has an exact representative in `Q[Y0,Y1,Y2,L]` (R2, R4) · V6 `L*D7^+`, `M*D7^- in J_8` exactly (R3) ·
V7 no coprime denominator of degree 5 (any parity mix) (R6; also elementary: `Phi(w) = 30 sum l_k^4 (w.grad l_k)^2`,
`J_4 = 0`, `K_4 = 0`, the eight `l_k^4` independent) · V8 no Y-pure degree-6 denominator (receipts at `1073741789`,
msolve + Singular; R7) · V9 no mixed degree-6 denominator (R5 with exact certificates; aggregator's independent route) ·
**V10 theorem: every certificate denominator has degree >= 7** (V3 + V7 + V8 + V9) ·
V11 `K_8^-/J_8^- = S_1 * K_7` over Q, dimension 6 (receipt `kernel-generation d=8 parD=1`, Hilbert-checked, + R2 + V6 logic).

**PLAUSIBLE** (modular evidence only, or chain incomplete):
P1 `d = 7` Y-even: `lambda = 0` and the 18 coefficients `w_{k,X_iY_j}` (`j != k`) vanish — one prime, exact certificates
INCOMPLETE; the 42-unknown residual Gröbner was OOM-killed (msolve exit 137 at F4 degree 6, matrix 2533891 x 2261343) ·
P2 `d = 7` Y-odd mirror (`lambda = 0`, X-blocks) — one prime; residual run ended without a verdict ·
P3 `d = 7` mixed: `lambda^± = 0` mod `p`, residual 120 unknowns / 636 quadrics ·
P4 fixed-`D` systems of `D7^±` inconsistent mod `p` (nullity exceeds the trivial count) ·
P5 `lambda`-only `d = 8` sub-families (mod-`p` representatives; family not canonical over Q) ·
P6 the twelve degree-8 kernel basis elements individually inconsistent mod `10^9+7` (aggregator) ·
P7 `dim K_8^+/J_8^+ = 6` (upper bound rigorous, `4 <= dim_Q` rigorous via V6; the two extra classes have no exact E1
certificate in the receipts) and `dim K_9^±/J = 18` (upper bounds) · P8 verify-lens "characteristic 0" Singular run of the
pure `d = 6` relaxation (artifact present on the remote; characteristic not confirmed here; not needed for V8).

**UNSUPPORTED** (as stated in the inputs; corrected above):
U1 excl report / packet README: "114/114 exact dependencies verified, remote job `exclB-spp-relax-mixed`" — no such record
existed anywhere; the numbers were right but the attribution was not (now supplied by R5) ·
U2 draft §3 "partly proven" and "Gröbner basis ... still running" — nothing at `d = 7` is proven; the run died ·
U3 draft §4 "covers `D7^+ + J_7` entirely once it finishes" — it did not finish · U4 draft "[TODO]" placeholder ·
U5 construct note: "dims mod p are lower bounds for Q-dims" — wrong direction (mod-`p` dimensions are upper bounds).

## 4. Open items
Degree 7 (all parities): rigorous closure needs either exact-rank certificates for the eliminating functionals (lift failed
at 12 primes; try more primes or a different functional basis) or a Gröbner basis of the 42-unknown residual (memory-bound).
Degree 8: ~160 unknowns; not semi-regular. No claim of nonexistence of a certificate is made; all results are degree bounds.
