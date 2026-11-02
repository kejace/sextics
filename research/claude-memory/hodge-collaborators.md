---
name: hodge-collaborators
description: Who the external hodge.computer contributors are (fqlx, benediamond, DiamondJM, hodgepodgeconjecture) and where their cited work lives, for correct credit lines
metadata:
  type: reference
---

External contributors on the hodge.computer discussions (credit in backticks, no @):
- `fqlx`: most prolific; merged solves incl. PR #150 (reflection-coupled-fano), PR #128 (holdout), #161 (holdout-bundle); posts say "assisted by OpenAI ChatGPT GPT-6 Pro and Codex" or Antigravity.
- `benediamond`: radial audits, incidence-curve integrations, 5-adic vertical-remainder checks (cycle-q2, reflection-coupled, one-reflection-skew).
- `DiamondJM` = Josh Diamond (GitHub user since 2018); posts "prepared with Claude Code assistance on behalf of `DiamondJM`"; merged PR #144 (one-reflection-y2-sector), five closed holdout-* PRs; repo https://github.com/DiamondJM/hodgeHoldout ("Autoresearch endeavor on hodge.computer"). Their degree-6 "natural support" sweep (61 supports for reflection-coupled-residual, 55 for one-reflection-skew, p=268435399) is only at commit af0589e: https://github.com/DiamondJM/hodgeHoldout/blob/af0589e/sweep/RESULTS.md — the sweep/ directory is gone from the default branch, so pin links to that commit.
- `hodgepodgeconjecture`: PR #151 (holdout-ternary-deformation), accepted cyclic construction.

**How to apply:** before crediting a sweep or result "from the thread", grep /tmp/hodge/discussions/<lane>.json for the author; if it lives elsewhere, link the exact commit and say it is not in the thread. See [[hodge-sextics-project-setup]].
