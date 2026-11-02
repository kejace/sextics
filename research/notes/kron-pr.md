# Kronecker checker → shared library module (maintenance PR draft)

Status: IN PROGRESS (started 2026-11-01). Log: /tmp/hodge/lean/kron-pr/log.jsonl. PR draft: /tmp/hodge/lean/kron-pr/PR.md.

## Layout
- Worktree: /home/kejace/devel/sextics-kron, branch `maint/kronecker-checker` from origin/main 83fbf5c9 (never pushed).
- New files: formal/lean/HodgeComputer/Reflection/{Kronecker.lean, KroneckerTest.lean, README.md}; paragraph in formal/README.md.
- Source: /tmp/hodge/lean/capacity/src/Reflect/Kron.lean (namespace renamed HodgeComputer.Kron → HodgeComputer.Kronecker;
  `set_option linter.unusedSimpArgs false` removed, unused simp args to be fixed from build warnings).
- Generators in /tmp/hodge/lean/capacity/gen/ still say `import Proposal.Kron` / `open HodgeComputer.Kron`; a solve PR
  would instead `import HodgeComputer.Reflection.Kronecker` and `open HodgeComputer.Kronecker` (one-line change in kron_submit.py).

## Decisions
- Import root formal/lean/HodgeComputer.lean imports only Verify + CertificateOperations; Constructions are not registered
  there either (built on demand by name; scripts/check-library.py builds changed modules by name). So the new modules are
  not added to the root; README documents `lake build HodgeComputer.Reflection.Kronecker`.
- Test module keeps linters on; `#print axioms` lines live in the test module only.

## Progress
- [x] worktree created, .lake copied (in background)
- [x] Kronecker.lean written
- [ ] README.md, KroneckerTest.lean, formal/README.md paragraph
- [ ] lake build (-j 4) + fix warnings; lake env lean; leanchecker replay; axioms
- [ ] timings table, PR.md, commit
