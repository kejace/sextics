# Lean capacity study — scaffolding (2026-11-01/02)

Report: /tmp/hodge/notes/lean-capacity.md.  Log: log.jsonl (append-only).  Timings: results/timings.jsonl
(merged from the remote box; records carry `label`, compile/kernel user+wall+RSS, and generator meta).

## Lean sources (src/)
- `Reflect/Kron.lean` — **Kronecker big-integer identity checker** (the recommended route).  Trusted interface is only
  `MvPolynomial`; proves `soundC : homE L → homE R → degE L = degE R → degE L < m → l1E L + l1E R < 2^W → psiE m W L = psiE m W R → denC L = denC R`
  plus `toC_isHomogeneous`, `toC_ne_zero_of_evalPT`, and the witness-binding lemmas (`W`, `denW`, `matchW_sound`).
  Copy it verbatim as the flat helper `Kron.lean` of a proposal (`import Proposal.Kron`); no banned tokens.
- `Reflect/Core0.lean`, `Core1.lean`, `Core2.lean` + `Bench*.lean` — sparse-list kernel arithmetic microbenchmarks (rejected route).
- `Reflect/KronRec.lean` — recursor-style variant experiment (no gain).
- `Cap/M3.lean`, `Cap/Contact.lean` — accepted Witness.lean files used as template baselines.
- `Cap/*.lean` — generated test modules (template `e1*/e2*/full*`, Kronecker `K_*`, binding `Bind*`, sliced `SL_*`).

## Generators (gen/)
- `gen_capacity.py` — exactly-true synthetic certificate identities (Koszul construction), accepted-template modules.
- `kron_gen.py` — same data, Kronecker-checked modules (handles rational F/A by a common scale).
- `slice_gen.py` — degree-in-X0 slicing of the accepted template with `linear_combination` assembly.
- `bind_test.py` — binding test: verifier-shaped term (`check.polynomial`) = `denW` by `rfl`, then `matchW` by kernel.
- `kron_submit.py` — **complete proposal generator** (Kron.lean, Witness.lean, Solution.lean, proof.json, witness.json)
  from an instance id + witness dict; check with `python3 -m hodge.leancheck --instance <id> --dir <out>`.
- `batch1.sh`, `batch2_remote.sh`, `batch3_remote.sh`, `kjob.sh` — remote job launchers (run_module.sh = verifier-like timing harness).

## Harness
- `run_module.sh <Module> <label> <cpuset>` — `lean -j 2 -M 14336` pinned to 2 cores + `leanchecker`, appends to results/timings.jsonl.
- `push.sh` / `pull.sh` / `rsubmit.sh` — sync with the remote box (hodge-sextics).
