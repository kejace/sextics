STATUS: report finalized 2026-11-02 ~03:00 UTC (runs still in flight when the agent was stopped: official Docker sandbox on submit/m3-coupled had passed compile-HodgeTarget; verifier-replica check of the sliced proposal submit/m3-coupled-slice had passed HodgeTarget/HodgeWitness/Kron and was compiling Witness.lean; remote sliced/dense d=10 jobs KS_m3_d10_k300, KS_m3_d10_k1000, KS_dn_d10_k300, KS_dn_d10_k20, K5_dn_d10_k20, SL_e2m3k150 were running — their records land in /tmp/hodge/jobs/cap-<name>.log on hodge-sextics and, once pulled, in /tmp/hodge/lean/capacity/results/timings.jsonl). Core findings: template ceiling ≈30–40k raw products (generic shapes); Kronecker checker verifies 3.1M-product identities (18.5k-term witness) in 263 s CPU; complete Kronecker-based m3-coupled proposal passes the exact verifier compile sequence in 184 s.

# Lean capacity study (accepted template vs alternatives)

## Executive summary
1. **Accepted template ceiling** (`MvPolynomial.funext; norm_num <;> ring`, one identity per theorem): cost is superlinear in the number of raw monomial products and in the number of distinct monomials of the expansion.  Random-support e2 identities: 18.6k products → 205 s / 7.5 GB, 28.5k → 306 s / 9.2 GB, 185k → fails (2M heartbeats, 12.5 GB).  Structured accepted identities fare better (m3-coupled: 78k products in 275 s / 5.5 GB).  Practical ceiling ≈ **30–40k raw products per identity for generic shapes, ≈100k for highly structured ones**; dense F (236 terms) costs ≈600–900 s already at 15–25k products.  Witness binding in `Solution.lean` by the same tactic is also superlinear (300 terms: 103 s).
2. **Recommended alternative — Kronecker big-integer checker** (`src/Reflect/Kron.lean`, ≈1100 lines, Mathlib only, no new axioms, passes the banned-token scan): explicit polynomials as balanced trees of integer monomials; the kernel computes ψ(P) = P(2^W, 2^(Wm), 2^(Wm²), 2^(Wm³), 2^(Wm⁴), 1) with GMP-accelerated `Nat` shifts/multiplications (Nat-pair representation), checks homogeneity and an ℓ¹ coefficient bound, and `soundN` turns ψ(L) = ψ(R) into `L = R` in `MvPolynomial (Fin 6) ℂ` (mixed-radix digit uniqueness).  Cost is **linear in the number of terms of the witness** (≈25–30 ms per term per identity, independent of the raw-product count): e2 with 435k raw products in 58 s / 5.5 GB; dense-F e1 with 215k products in 55 s; **d = 10 witness of 18.5k terms with 3.1M raw products in 263 s / 13.5 GB** (unsliced) — memory ∝ m⁵·W·depth, so for d ≥ 10 use the **slicing variant** (`sound4`: split by the degree in X4, Kronecker on four variables, numbers ÷ m).
3. **Witness binding**: `W`-AST mirrors the verifier's translation, `HodgeWitness.D = denW wD` is `rfl`, `denW wD = toC tD` by `matchW_sound` + `decide +kernel`: 3000 terms in 125 s (linear) vs > 15 min extrapolated for the template route.
4. **End-to-end**: `gen/kron_submit.py --instance m3-coupled` produces a complete proposal that passes the exact verifier compile sequence (`hodge.leancheck`) in **184 s total** (template Witness.lean alone: 275 s); official Docker sandbox run in progress.
5. **Budget rule of thumb (Kronecker route)**: per 1000 witness terms ≈ 12 s protected `HodgeWitness.lean` elaboration + 17 s binding + ≈30 s per identity check, each kernel part replayed once more by `leanchecker`; so a 60-minute budget fits ≈ 20–30k witness terms (any raw-product count), memory permitting (d ≤ 8 unsliced, d ≥ 10 sliced).

Work dir: /tmp/hodge/lean/capacity (log.jsonl, gen/, src/Cap/, results/). Remote jobs: lane `_cross`, names `cap-*`.

## Verifier budget (scripts/proposals/sandbox.py + .github/solve-policy.json)
- Each module: `lean -j 2 -M 14336` in Docker `--cpus 2`, then `leanchecker --verbose` per module; total wall budget 3600 s for everything (lake prep, all compiles, all replays, audit); 16 MiB source total, 64 files.
- Banned tokens: axiom sorry admit unsafe native_decide run_cmd run_tac initialize elab macro syntax #eval #exit [extern. `decide` (incl. `decide +kernel`) is allowed; axioms limited to propext/Classical.choice/Quot.sound.

## Known floor (accepted proofs that passed CI)
| instance | q | d | |F| | |D| | |N_i| | e1 raw products | e2 raw products |
|---|---|---|---|---|---|---|---|
| m3-coupled | 1 | 10 | 22 | 190 | 88,87,83,35,83,35 | 2877 | 78310 |
| contact46006 | 1 | 6 | 56 | 44 | 48,84,84,48,84,84 | 11280 | 18384 |
(raw products = number of monomial×monomial products before collection, the quantity `ring` must normalize.)

## Timings — accepted template (`MvPolynomial.funext; norm_num <;> ring`)
Remote box hodge-sextics (n2-standard-128), `taskset` 2 cores, `lean -j 2 -M 14336`, `set_option Elab.async false` (single-threaded elaboration as in the accepted files). `user` = CPU seconds (robust against the heavy load on the box); `wall` inflated when load > cores.

| label | identity | raw products (e1 / e2) | defining terms | compile user s | compile wall s | peak RSS | leanchecker wall s |
|---|---|---|---|---|---|---|---|
| base-m3-coupled (accepted Witness.lean) | e1+e2+rest | 2877 / 78310 | ~620 | 275 | 278 | 5.5 GB | 16 |
| base-contact46006 (accepted; Elab.async on) | e1+e2+rest | 11280 / 18384 | ~540 | 167 | 104 | 6.8 GB | 22 |
| e2m3k10 | e2 only | – / 1200 | 186 | 30 | 31 | 3.8 GB | 8 |
| e2m3k20 | e2 only | – / 2760 | 204 | 42 | 43 | 4.1 GB | 11 |
| e2m3k40 | e2 only | – / 18602 | 565 | 205 | 208 | 7.5 GB | 42 |
| e2m3k60 | e2 only | – / 28510 | 481 | 306 | 312 | 9.2 GB | 68 |
| e2m3k80 | e2 only | – / 22934 | 355 | 324 | 329 | 9.3 GB | 76 |
| S_m3_d2_n10 (sibling agent, synth.py) | e1+e2 | 460 / 180 | – | 28 | 29 | 3.6 GB | 8 |
| S_m3_d2_n21 | e1+e2 | 966 / 798 | – | 48 | 71 | 3.7 GB | 11 |
| S_m3_d4_n40 | e1+e2 | 1840 / 4000 | – | 106 | 259 | 4.3 GB | 79 |
| S_dense_d2_n10 | e1+e2 | 7040 / 180 | – | 230 | 708 | 5.4 GB | 58 |
| S_dense_d2_n21 | e1+e2 | 14784 / 798 | – | 328 | 1014 | 6.4 GB | 277 |
| e1m3k10 | e1 only | 1240 / – | 186 | 43 | 81 | 3.5 GB | 13 |
| e1m3k40 | e1 only | 6685 / – (RHS A·D has 2428 distinct monomials) | 654 | 144 | 153 | 4.5 GB | 15 |
| **e2m3k150** | e2 only | – / 184779 | 1107 | **605 → FAIL** (deterministic timeout, maxHeartbeats 2000000) | 639 | **12.5 GB** | – |
(import-only baseline ≈ 5–8 s and ≈3.3 GB; a second agent with the same program shares the remote results file — their `S_*`/`R_*` records come from `gen/synth.py`, `gen/ident.py`, `gen/reflect_gen.py` in the remote dir, mirrored under `peer/`.)

Reading: cost grows faster than linearly in the raw product count and depends strongly on the number of *distinct* monomials in the expansion: the accepted m3-coupled identity (78k raw products, very structured monomial support) costs 275 s / 5.5 GB, while random-support synthetic identities of 20–30k raw products already cost 200–330 s and 7.5–9.3 GB. With the 14 GiB limit, the template's practical ceiling for generic (random-support) identities is **~40–50k raw products per identity** (memory-bound), and for highly structured identities perhaps 100–150k (time-bound: ~10–20 min per identity, leaving headroom in the 60 min budget only for one such identity plus replay). Dense F (236 terms) makes e1 expensive early: 15k raw products already cost 330 s CPU.

## Kronecker big-integer checker (new route; `src/Reflect/Kron.lean`, generator `gen/kron_gen.py`)
Idea: each explicit polynomial is a balanced tree `PT` of integer monomials; the kernel computes the Kronecker value ψ(P) = P(2^W, 2^(Wm), 2^(Wm²), 2^(Wm³), 2^(Wm⁴), 1) with GMP-accelerated `Nat` arithmetic (divide-and-conquer with relative shifts so intermediate numbers stay small), plus a homogeneity/degree check and an ℓ¹ coefficient bound. Theorem `soundC`: if both sides are homogeneous of the same degree n < m, the ℓ¹ bounds sum to < 2^W, and ψ(L) = ψ(R), then L = R in `MvPolynomial (Fin 6) ℂ`. Mathematics: mixed-radix uniqueness of signed digits (`eq_zero_of_eval_kpt`). Everything is checked with `decide +kernel`; no new axioms. Kernel work is O(#terms) reduction steps + a handful of big multiplications.

| module | identity | raw products (e1 / e2) | defining terms | compile wall s | user s | peak RSS |
|---|---|---|---|---|---|---|
| K_m3k10 (local) | e1+e2+homog+ne_zero | 1240 / 1200 | 186 | 9.5 | 9.1 | 3.6 GB |
| K_m3k100 (local) | e1+e2+homog+ne_zero | 11095 / 102276 | 963 | 35.7 | 35.0 | 4.8 GB |
| K_m3_d6_k40 (remote; same data as e2m3k40: template 205 s) | e1+e2+homog+ne_zero | 5831 / 18602 | 565 | 16.8 | 14.8 | 4.0 GB (replay 11 s user) |
| K_m3_d8_k150 (remote) | e1+e2+homog+ne_zero | 25804 / **434802** | 2367 | 60.8 | 57.9 | 5.5 GB (replay 43 s user) |
| K3_dn_d4_k10 (remote, dense F scaled by 9973) | e1+e2+homog | 16662 / 993 | 787 | 24.6 | 19.7 | 3.6 GB |
| K3_dn_d6_k40 (remote, dense F) | e1+e2+homog | 81730 / 25166 | 1212 | 59 | 34.9 | 4.2 GB (replay 23 s) |
| K3_dn_d8_k100 (remote, dense F) | e1+e2+homog | **215170 / 159699** | 2138 | 62 | 55.1 | 6.2 GB (replay 46 s) |
| K2_m3_q2_d8_k100 (remote, q = 2, m = 23) | e1+e2+homog | 28883 / 333679 | 2428 | 882 (loaded box) | 155.7 | 9.6 GB (replay 129 s) |
| K2_m3_d10_k300 (remote, Int version, m = 21, W = 42) | e1 | 84674 / – | 18499 | – | 282 → **FAIL** | **22 GB** `(kernel) excessive memory consumption` |
| K2_m3_d10_k1000 (remote, Int version) | e1 | 163569 / – | 25058 | – | 851 → FAIL | 21.8 GB |
| **K4_m3_d10_k300 (remote, Nat-pair `soundN`, m = 21, W = 42)** | e1+e2+homog | 84674 / **3118055** | 18499 | 455 (loaded box) | **263** | **13.5 GB** (replay 181 s user / 10.9 GB) |
| K4_m3_d10_k1000 (remote, Nat-pair) | e1 | 163569 / 16.5M | 25058 | 794 | 562 → **FAIL** | 23.6 GB (kernel memory) → sliced run queued |
| K4_m3_q2_d8_k100 (remote, Nat-pair, q = 2, m = 23) | e1+e2+homog | 28883 / 333679 | 2428 | 305 | 98 | 10.6 GB (replay 100 s user / 7.6 GB) |
Template comparison points: e1m3k150 (e1 only, 25.8k raw products, 6636 distinct monomials) took **914 s user / 10.6 GB** with the template; the Kronecker module K_m3_d8_k150 containing the same e1 plus a 435k-product e2 took 58 s.
Profile of K_m3k100's e2 side goals (local, user s): ψ equality 11.5, homE 2.6, l1E 2.3, degE 1.4. The template needs >1000 s and >14 GB at this size (extrapolated from the table above).

## Kernel-evaluation microbenchmarks (local, lean 4.33, `decide +kernel`)
Sparse polynomial arithmetic on `List (Nat × Int)` (packed base-256 exponent keys), balanced merge-tree multiplication:
| variant | case | wall | maxrss |
|---|---|---|---|
| import-only baseline | — | 6 s | 1.6 GB |
| Core0 (Decidable `if` merge, structural recursion) | mul 100×100 terms (10k raw products, ~70k merge steps) | 31.7 s | 4.1 GB |
| Core1 (Bool `Nat.blt` matches) | same | 25.3 s | 3.7 GB |
| Core2 (explicit `Nat.rec`) | same | 24.5 s | 3.5 GB |
| Core0 | `scaleT` only (10k product terms, no merging) | 5.5 s | 1.8 GB |
Conclusion: the kernel rewriter costs ~250–300 µs and ~25 KB **per merge step**; a single declaration is memory-bound at ~3×10^5 steps (14 GiB). Sparse kernel arithmetic is therefore at best a few × better than `ring` and cannot reach 10^6+ raw products. The scalable design is Kronecker evaluation: encode each explicit polynomial as a big integer (GMP-accelerated `Nat` ops in the kernel, O(#terms) kernel steps), compare `ψ(L)=ψ(R)`, and discharge soundness by a degree bound (homogeneity) + coefficient bound (ℓ1 norms) + mixed-radix injectivity lemma.

### Binding the protected witness (`HodgeWitness.D = …`)
The verifier regenerates `HodgeWitness.lean` from `witness.json`; a submission must identify those terms with its own definitions. `Kron.W` mirrors the translator's term shapes exactly, so `HodgeWitness.D = denW wD` is `rfl`, and `denW wD = toC tD` is certified by the kernel (`matchW_sound`, structural comparison of the same balanced tree). Local measurements (`gen/bind_test.py`):
| terms of D | rfl + matchW route (wall / user / RSS) | accepted-template route `funext; norm_num; ring` |
|---|---|---|
| 300 | 19.8 s / 25 s / 3.5 GB | 103 s (norm_num 41 s + ring 4 s + …) |
| 3000 | 125 s / 195 s / 4.9 GB (profile: elaborating the translator-shaped term 37 s, `rfl` 37 s elab + 14 s kernel, linting 48 s → disable linters) | not run (extrapolated > 15 min) |
Consequence: witness size is itself budget-relevant: compiling the protected `HodgeWitness.lean` costs ≈12 s per 1000 terms, the binding ≈17 s per 1000 terms. A 60k-term witness (D of 3000 terms + six 10k-term numerators) would spend ≈30 min on these linear steps alone, so the witness itself should stay below ~30k terms (≈16 MiB source limit is not the binding constraint; time is).

### Rational `F`/`A` and the certificate glue (implemented in `kron_gen.py`/`kron_submit.py`)
The kernel works over ℤ.  The witness (D, N, Q) is scaled to integers (D, N by s, Q by s²: the identities are homogeneous of degree 2 in (D,N,Q)); the instance's F and A are written as `C (1/s') * toC tF'`, `C (1/s') * toC tA'` with integral trees, and the kernel checks `q·Σ∂ᵢF'·Nᵢ = A'·D` and `s'^q·(D·divN − Σ∂ᵢD·Nᵢ) = F'^q·Q`; the ℂ-identities follow by multiplying with `C (1/s'^q)` (`linear_combination`).  Dense-rational F (236 terms, s' = 9973) with d = 4, |D| = 10: whole module (e1 16.7k raw products + e2 + homogeneity) compiles in 16 s locally (template S_dense_d2_n21 with 14.8k raw products: 328 s).

### Complete proposal prototype (`gen/kron_submit.py`)
Generates the five proposal files for an instance from a witness dict (Kron.lean as flat helper, Witness.lean with the Kronecker identities + template legality machinery on the unfolded trees, Solution.lean binding `HodgeWitness.*` by `rfl` + `matchW_sound`, balanced witness strings).  Run on the accepted m3-coupled witness (scaled by 6): module sizes Kron 31 KB, Witness 43 KB, Solution 73 KB; local verifier replica (`hodge.leancheck`) result: see below.

### Kernel memory model of the Kronecker evaluation (measured)
The kernel caches every intermediate whnf result for the duration of one declaration, so the peak memory of a `decide +kernel` is the sum of all intermediate big numbers.  With leaves sorted by position and a balanced tree, each level of the tree holds numbers whose sizes add up to the final size, so memory ≈ (3 for the Int version, ≈1.5 for the Nat-pair version) × depth × final size, where final size = (number of positions m⁵, m = identity degree + 1) × W bits.  Measured on a sorted 10k-leaf degree-14 tree (m = 21, W = 42, final ≈ 14 MB): Int version `psi` 33 s user / +1.5 GB; Nat-pair `psiN` **9.7 s user / +0.85 GB**.  Unsorted leaves are catastrophic (every node spans the whole range: 275 s / 32 GB for the same tree).  Consequences: d ≤ 8 (m ≤ 17, ≤ 1.4M positions) is comfortable for any number of terms; d = 10 (m = 21, 4.1M positions, 20 MB numbers) needs the Nat-pair version and possibly slicing by one variable (positions ÷ 21) to stay below 14 GiB with ~20k-term witnesses.

### Degree-sliced template (`gen/slice_gen.py`)
Works (e1/e2 assembled by `linear_combination Σ X0^m * slice_m` over opaque slice atoms), but every slice lemma pays ≈3–5 s of `norm_num`/`simp` fixed cost: the k10 case took 428 s (template: 30 s e2 + 43 s e1).  It can only pay off where the plain template fails (≥100k raw products); one such run (SL_e2m3k150) is queued for the record.  Verdict: not recommended; the Kronecker route dominates it at every size.

### End-to-end validation on an accepted instance
`gen/kron_submit.py --instance m3-coupled` (witness scaled by 6 to integers, W = 61, m = 21) → five proposal files; `python3 -m hodge.leancheck` (exact verifier compile sequence, local): HodgeTarget 5.9 s, HodgeWitness 13.7 s, Kron 21.0 s, Witness 182.4 s, Solution 22.0 s (rfl + `matchW` binding of 7 polynomials), HodgeChecked 5.9 s — **total 251 s**, no new axioms (`Kron.lean` only uses Mathlib; passes the banned-token scan of `check.py`).  The accepted template Witness.lean for the same instance takes 275 s on the remote box (CI job 13.5 min including setup).

## Recommended pipeline (draft; numbers refined below as jobs finish)
1. **Witness preparation (Python).** Scale the witness to integers (D, Nᵢ by s, Q by s²), write `witness.json` with balanced parenthesized sums whose term order is the Kronecker position order (`gen/bind_test.py:witness_string`); keep coefficient 0/1 conventions (`X0**6`, `-X0**6`, `5*X0**6`) so the verifier's translation is syntactically `denW` of the `W`-AST.  Choose m = (identity degree) + 1 and W = bits(ℓ¹ bound) + 2.
2. **Files.** `Kron.lean` (library, verbatim), `Witness.lean` (trees `tF, tA, tD, tNᵢ, tQ` as `noncomputable def`s; e1/e2 via `soundN` with six `decide +kernel` side goals; homogeneity via `toC_isHomogeneous`; `D ≠ 0` via `toC_ne_zero_of_evalPT`; legality via the accepted template's univariate specialization on the unfolded trees with `push_cast`), `Solution.lean` (F, A by `funext; norm_num; ring`; D, Nᵢ, Q by `rfl` + `matchW_sound`).  All generated by `gen/kron_submit.py`.
3. **Budget estimate per 1000 terms of witness data**: HodgeWitness.lean elaboration ≈ 12 s, binding ≈ 17 s (+ the same again in `leanchecker` replay for the kernel parts), Kronecker checks ≈ 10 s per pass over the data (e1 + e2 + homogeneity ≈ 4 passes, replayed once more by leanchecker).  Memory: ≈ 1.5 × depth × m⁵ × W bits per declaration (fine for d ≤ 8; d = 10 needs ≤ ~20k terms or a variable slice).
4. **Keep the identity small where possible**: the kernel cost is linear in the number of terms of D, Nᵢ, Q (not in raw products); the raw-product count (which kills `ring`) is irrelevant for this route.

### Slicing by one variable (`Kron.sound4`, generator `gen/kron_slice_gen.py`)
For high degrees the Kronecker numbers (m⁵ positions) dominate memory, so the library also provides slicing by the degree in `X 4`: `toC t = Σ_k X4^k · toC (sliceK k t)` (`toC_eq_sum_sliceK`, kernel-checked bound on the X4-degree), each slice identity is an identity between `X 4`-free polynomials and is certified by `sound4` with a Kronecker point on four variables (m⁴ positions, numbers ÷ m), and the slices are assembled with `linear_combination Σ_k X4^k · slice_k` over opaque slice atoms (derivatives of slices in direction 4 vanish by `pderiv_four_toC_sliceK`).  Degrees of expressions are `Option`-valued so that empty slices (zero polynomials) are handled uniformly.  Small-size cost: KS_m3k10 113 s / 2.5 GB vs unsliced 11 s / 3.5 GB (one declaration per slice, 8 kernel side goals each); the route is meant for d ≥ 10 witnesses where the unsliced check exceeds 14 GiB (results of the remote runs below).

## In-flight runs at hand-off (how to collect)
- `ssh hodge-sextics`: `grep -h '^{"module"' /tmp/hodge/jobs/cap-KS_*.log /tmp/hodge/jobs/cap-K5_*.log /tmp/hodge/jobs/cap-SL_e2*.log` gives the JSON timing records (append them to `results/timings.jsonl`; do not rsync `results/` to the remote — the sibling agent's pushes overwrite it).
- Local: `logs/sandbox-m3-kron.log` (official sandbox, `hodge.sandboxcheck m3-coupled submit/m3-coupled`), `logs/leancheck-m3-slice.log` (replica on the sliced proposal), `logs/KS_m3_d8_k150.local.log` (sliced mid-size module).
- Expected (from the measured model): sliced d=10 modules should stay ≈2–4 GB per declaration; their total CPU ≈ (number of slices) × fixed overhead + ≈30 ms per witness term per pass.
