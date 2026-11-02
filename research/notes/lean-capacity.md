STATUS: in progress — template ceiling located (random-support e2: 28k raw products OK at 306 s/9.2 GB; 185k raw products times out at 2M heartbeats after 605 s/12.5 GB); Kronecker checker verifies 435k-product e2 in 58 s/5.5 GB (remote) and binds verifier-shaped witness terms by rfl+kernel (3000 terms: 125 s); more Kronecker scaling jobs running on the (overloaded) remote; degree-slicing generator in progress.

# Lean capacity study (accepted template vs alternatives)

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
