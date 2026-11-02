# Kronecker checker → shared library module (maintenance PR draft)

Status: DONE 2026-11-02 (commit 575e666b on maint/kronecker-checker, not pushed, no PR opened). Log: /tmp/hodge/lean/kron-pr/log.jsonl. PR draft: /tmp/hodge/lean/kron-pr/PR.md.

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
- [x] worktree created, .lake copied
- [x] Kronecker.lean written; 5 lint warnings fixed (2x eval_C, map_C, Finsupp.single_apply unused simp args; unused hyp hQ' dropped from abs_coeff_mul_le + 2 call sites); six 268-char simp lines wrapped
- [x] README.md (Reflection/), KroneckerTest.lean (small identity + m3-coupled e1/e2 via soundN + homogeneity + ne_zero + #print axioms), formal/README.md paragraph
- [x] lake build --wfail passes for both modules (Lake 5 has no -j; used taskset -c 0-3). Kronecker 21 s; KroneckerTest 93-97 s
- [x] axioms (all six #print axioms): [propext, Classical.choice, Quot.sound]
- [x] lean -j 2 -M 14336 KroneckerTest: 95 s wall / 86 s user / 10.9 GB with async elaboration; added `set_option Elab.async false` -> lake build 97 s / 8.6 GB
- [x] leanchecker replay Kronecker: 18.5 s wall / 10.6 s user / 2.9 GB
- [x] leanchecker replay KroneckerTest: 97 s wall / 83 s user / 8.2 GB; lean -j 2: Kronecker 32 s / 3.6 GB, KroneckerTest 97 s / 8.5 GB
- [x] content checks: not runnable locally (mistune missing); formal/README.md is not website content
- [x] PR.md written; commit 575e666b
- [x] check-library.py --base origin/main: README under formal/lean => full library rebuild + Docker recheck of all applications (CI will do this); stopped locally after 88 modules

## Findings worth keeping
- Memory for d=10 identities is dominated by full-span Kronecker numbers per tree occurrence (m^5*W bits, ~31 MB at m=21,W=61), ~independent of term count: e2 of m3-coupled (13 tree occurrences) = 8.6 GB alone. d<=8 is ~4x smaller.
- Lake 5.0 (Lean 4.33 toolchain) has no `-j`; limit CPU with taskset; `lake build --wfail` is the warnings-as-errors switch.
