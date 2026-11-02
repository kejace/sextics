STATUS: in progress — 19 template timing jobs running on hodge-sextics (baseline m3-coupled/contact46006 + synthetic e1/e2/full at 1k..500k raw products); kernel sparse-arithmetic microbenchmarks done (too slow/memory-hungry, see below); now building the Kronecker big-integer checker (Lean scaffolding under /tmp/hodge/lean/capacity/src/Reflect/).

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

## Timings
(to be filled)

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
