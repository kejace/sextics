---
name: hodge-lean-capacity
description: "Lean verifier capacity facts and the Kronecker big-integer checker for large certificate identities (hodge.computer), 2026-11-02"
metadata:
  node_type: memory
  type: project
  originSessionId: 2f213184-0ffd-40a3-8b2f-58f1332191ff
---

Report: /tmp/hodge/notes/lean-capacity.md; scaffolding: /tmp/hodge/lean/capacity (README.md, src/Reflect/Kron.lean, gen/*.py).

- Verifier budget: each module `lean -j 2 -M 14336` in Docker (2 cpus), 60 min total incl. leanchecker replay of every module; 16 MiB source. `decide +kernel` allowed (axioms stay propext/choice/Quot.sound); `native_decide`, `macro`, `syntax`, `elab` banned. `matches` is a Lean keyword (don't name a def that).
- Accepted template (`funext; norm_num <;> ring`) ceiling: random-support e2 identities ~30k raw monomial products (300 s, 9 GB); 185k raw products fails (2M heartbeats, 12.5 GB); structured accepted identities up to ~80k ok (m3-coupled 275 s).
- Kernel facts: sparse-list polynomial arithmetic in the kernel costs ~280 us + ~25 KB per merge step (useless at scale). `Nat.pow` kernel acceleration stops above exponent 2^24 (use `Nat.shiftLeft`); Nat mul/shift on 10^7-bit numbers are ~ms. The kernel caches every intermediate for a declaration, so memory = sum of intermediates. Big data `def`s must be `noncomputable` (the code generator blows up on 10k-leaf trees: 49 s + GBs, "maximum recursion depth"). Disable linters (`set_option linter.all false`) in big files.
- Kronecker checker (src/Reflect/Kron.lean): balanced `PT` trees of integer monomials sorted by position, `soundN`: homogeneity + ℓ1 bound + Nat-pair Kronecker value equality ⇒ identity in MvPolynomial (Fin 6) ℂ. Verified 435k-product e2 in 58 s/5.5 GB; dense F (236 terms) e1 215k products in 55 s. Memory ≈ 1.5 × depth × (m^5 · W bits): d ≤ 8 easy, d = 10 (m = 21) tight at 14 GiB with ~20k-term witnesses.
- Binding verifier-generated `HodgeWitness.*`: `W` AST mirrors check.py's `polynomial()` shapes, `HodgeWitness.D = denW wD := rfl`, then `matchW_sound` by `decide +kernel`: 3000 terms ≈ 125 s (template funext/norm_num/ring: 103 s for only 300 terms). Protected HodgeWitness.lean elaboration itself ≈ 12 s / 1000 terms.
- `gen/kron_submit.py --instance <id>` builds a complete proposal; m3-coupled passes `hodge.leancheck` in 251 s total (first failure: Int.cast coefficients need `push_cast` before `polynomial`).
- A sibling agent ran the same program with sparse-list reflection (`Cap/Reflect.lean`, `gen/synth.py`) in the shared remote dir; results merge in results/timings.jsonl.
