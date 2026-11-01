---
name: hodge-sextics-project-setup
description: "Workspace layout, tooling and compute options for the hodge.computer sextics challenge (open-lane attack started 2026-11-01)"
metadata:
  node_type: memory
  type: project
  originSessionId: 2f213184-0ffd-40a3-8b2f-58f1332191ff
  modified: 2026-11-01T22:03:00.684Z
---

Working on the 16 open lanes of hodge.computer (repo ~/devel/sextics, GitHub hodge-computer/sextics; user account `kejace`, gh CLI authenticated). Started 2026-11-01.

- Scratch engine lives OUTSIDE the repo at ~/devel/sextics-work/hodge/ (core.py = exact polys/parsing/verify; fixed_d.py = fixed-denominator solver mod p + rational lift, target_image, e1_solvable). python-flint 0.9 installed via `pip install --user --break-system-packages`.
- Research notes/dossiers/logs under /tmp/hodge/ (notes/, discussions/<id>.md, results/<id>/log.jsonl, ext/hodgeHoldout = DiamondJM sweep repo).
- Lean toolchain v4.33.0 + mathlib cache already built in ~/devel/sextics/formal/lean (.lake present); check files with `lake env lean` from that dir.
- Machine: 40 cores, 125 GB RAM, no GPU. User offers a gcloud box (128 cores, 512 GB) on request. Docker image `ddg-cas:latest` has Singular, Julia, python-flint.
- Witness-to-Lean pipeline verified 2026-11-01: hodge/leangen.py (generate Witness.lean/Solution.lean in the accepted template), hodge/leancheck.py (local replica of the verifier compile sequence), hodge/sandboxcheck.py (official Docker verifier); a regenerated accepted witness passed the official sandbox (kernel replay + axiom audit) in ~5 min.
- Shell aliases to avoid: `du` is `dust`, `gp` is `git push`. `rm -rf` is blocked by a hook (ask user).
- Open instances (2026-11-01): cycle-q2, cyclic-triple-residual, dense-rational, diagonal-permutation, eight-plane-configuration, four-plus-two-residual, klein-coupled-character, off-paired, one-reflection-skew, polarized-k3-remainder, quadratic-cubic-pullback, quadric-contact-residual, rank-one-biquartic, reflection-coupled-residual, short-power-presentation, symmetric-sextic-residual.

**Why:** setup took an hour; next session should resume directly from the engine and notes.
**How to apply:** reuse ~/devel/sextics-work and /tmp/hodge (if /tmp was cleared, dossiers must be regenerated from discussions). See [[hodge-sextics-math-insights]] for mathematical findings.
