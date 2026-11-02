# Handoff: hodge.computer open-lane attack

Last updated 2026-11-01 (local machine `brunos`; remote `hodge-sextics`).

## Where things are

- Upstream repo: `hodge-computer/sextics` (clean checkout at `~/devel/sextics`, branch `main`). Nothing is ever edited there except in a solve-PR branch.
- Fork: `kejace/sextics` (remote name `fork`). Branch `research` = this directory, checked out as a worktree at `~/devel/sextics-research`.
  - `research/engine/` — the certificate engine (`hodge/core.py`, `fixed_d.py`, `leangen.py`, `leancheck.py`, `sandboxcheck.py`, tests). Live copy at `~/devel/sextics-work`.
  - `research/notes/` — book technique maps, accepted-proof catalogue, per-lane dossiers (`dossier-<id>.md`), attack reports (`attack1-<id>.md`), probes. Live copy at `/tmp/hodge/notes`.
  - `research/results/<id>/` — per-lane experiment logs (`log.jsonl`, scripts). Live copy at `/tmp/hodge/results`.
  - `research/discussions/` — saved GitHub discussion threads for the 16 open lanes.
  - `research/workflows/` — the Workflow scripts used by Claude Code (resumable with `resumeFromRunId`).
  - `research/claude-memory/` — Claude Code memory files for this project.
  - `research/bin/sync.sh` — snapshot everything above into this branch and push. Run it from cron every 10 minutes:
    `*/10 * * * * $HOME/devel/sextics-research/research/bin/sync.sh >> /tmp/hodge/sync.log 2>&1`
  - `research/bin/restore.sh` — unpack the branch onto a fresh machine.

## Fresh machine bootstrap (what was run on `hodge-sextics`)

```sh
# packages: build-essential git rsync python3-pip docker.io msolve singular nodejs npm gh tmux
python3 -m pip install --user --break-system-packages python-flint sympy numpy scipy pyyaml markdown-it-py
curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
sudo npm install -g @anthropic-ai/claude-code
gh auth login            # or copy ~/.config/gh/hosts.yml
gh ssh-key add ~/.ssh/id_ed25519.pub --title hodge-sextics   # after ssh-keygen on the machine
git clone git@github.com:kejace/sextics.git ~/devel/sextics && cd ~/devel/sextics
git remote rename origin fork && git remote add origin https://github.com/hodge-computer/sextics.git && git fetch origin && git checkout -B main origin/main
git worktree add ~/devel/sextics-research research
~/devel/sextics-research/research/bin/restore.sh
cd ~/devel/sextics/formal/lean && lake exe cache get && lake build HodgeComputer     # ~10 min
```

## Resuming the Claude Code session itself

Session transcripts live in `~/.claude/projects/-home-kejace-devel-sextics/` (this session: `2f213184-0ffd-40a3-8b2f-58f1332191ff.jsonl`, with `subagents/` and `workflows/`). Copy that directory to the same path on the new machine (same user `kejace`, same checkout path `~/devel/sextics`), set the Anthropic token, then from `~/devel/sextics`:

```sh
claude --resume 2f213184-0ffd-40a3-8b2f-58f1332191ff
```

Background workflows that were running are not carried over; resume them with `Workflow({scriptPath, resumeFromRunId})` using the ids in `research/workflows/` (completed agents return cached results from the journals).

## State of the work (2026-11-01 evening)

- Engine verified against accepted witnesses; witness-to-Lean pipeline verified by the official Docker verifier.
- No lane has a certificate. Rigorous partial results so far: polarized-k3-remainder is conditionally equivalent to an open case of the generalized Hodge conjecture; rank-one-biquartic has no degree-6 denominator family through the split limit; diagonal-permutation has one exceptional degree-8 E1-admissible denominator under test; at q=1 products of excluded denominators are excluded (Gysin additivity over components).
- Next: finish wave 1 (16 lane reports), run wave 2 (rigorous degree bounds with msolve, joint bilinear searches, theory probes), draft discussion posts per lane with reproduction packets, and sync this branch after each wave.

## 2026-11-02 04:30 UTC — session state (brunos, session 2f213184-0ffd-40a3-8b2f-58f1332191ff)
- Twelve refereed comments gated and ready at /tmp/hodge/posts/<lane>/comment.md (links pinned to research commit ae801ece…); the user posts them with `bash /tmp/hodge/posts/post_all.sh`.
- Workflows in flight: wave 4 routes `wf_11d84a5e-70d` (11 routes + verifiers; notes wave4-<lane>-<route>.md, logs results/<lane>/wave4/log.jsonl), referee D `wf_2fec8a28-1af` (polarized-k3-remainder, short-power-presentation), wave 3 closure+synthesis `wf_2ad30444-807` (closure table notes/wave3-closure.md; synthesis notes/wave3-construction-ranking.md COMPLETE), referee C `wf_8c1c2178-f69`, search `wf_c55cd42a-6ac`, construct relaunch `wf_17d5e026-da0`, trio theory `wf_85383db5-341`, Kron PR prep `wf_cd5e71da-4bb` (worktree ~/devel/sextics-kron, branch maint/kronecker-checker, never pushed).
- Resume: `claude --resume 2f213184-0ffd-40a3-8b2f-58f1332191ff` from ~/devel/sextics; re-launch a workflow with Workflow({scriptPath: ~/.claude/projects/-home-kejace-devel-sextics/<session>/workflows/scripts/<name>-<runid>.js, resumeFromRunId}).
- rjob gate now fails closed and admits jobs only at remote load <= 300, >= 110 GB free, <= 200 allocated cores (engine CHANGELOG 2026-11-02).
