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
- Remote compute (created 2026-11-01): gcloud instance `hodge-sextics`, project gluck-proof, zone us-east1-b, n2-standard-128 (128 vCPU/503 GB), 35.229.16.59; reach with `gcloud compute ssh hodge-sextics --zone us-east1-b`. Fully provisioned mirror: fork clone at ~/devel/sextics, research worktree, engine, /tmp/hodge, Lean built, gh authenticated. Work is versioned on branch `research` of fork kejace/sextics via research/bin/sync.sh (run every 10 min by cron) and restored with restore.sh; HANDOFF.md there has the restart recipe. Don't `pkill -f` a script name inside a gcloud ssh --command (kills the session).
- Remote job runner: `python3 -m hodge.rjob {up|submit|status|log|wait|pull|sh}` (engine file hodge/rjob.py) pushes the engine + lane dirs to hodge-sextics over ssh (key ~/.ssh/google_compute_engine) and runs detached jobs there; pull merges remote log.jsonl lines into the local lane log. Verified 2026-11-01: a degree-6 fixed-D solve takes 7 s there.
- Wave 2 launched 2026-11-01 ~19:00 EST from this session (single runner on brunos; heavy compute on hodge-sextics via rjob): workflows hodge-wave2-panels-A (wf_56479152-f22), panels-B (wf_0a9f2713-e92), search (wf_c55cd42a-6ac), cross (wf_e1e9e78c-568); scripts under ~/.claude/projects/-home-kejace-devel-sextics/2f213184-0ffd-40a3-8b2f-58f1332191ff/workflows/scripts/ (mirrored in research/workflows). Reports land as /tmp/hodge/notes/wave2-<lane>-{excl,construct,search,draftpost}.md plus wave2-status.md, wave2-trio-theory.md, lean-capacity.md.
- Shell aliases to avoid: `du` is `dust`, `gp` is `git push`. `rm -rf` is blocked by a hook (ask user).
- Open instances (2026-11-01): cycle-q2, cyclic-triple-residual, dense-rational, diagonal-permutation, eight-plane-configuration, four-plus-two-residual, klein-coupled-character, off-paired, one-reflection-skew, polarized-k3-remainder, quadratic-cubic-pullback, quadric-contact-residual, rank-one-biquartic, reflection-coupled-residual, short-power-presentation, symmetric-sextic-residual.

**Why:** setup took an hour; next session should resume directly from the engine and notes.
**How to apply:** reuse ~/devel/sextics-work and /tmp/hodge (if /tmp was cleared, dossiers must be regenerated from discussions). See [[hodge-sextics-math-insights]] for mathematical findings.
