---
name: feedback-subagent-incremental-logging
description: "User requires every subagent to log/write results incrementally to disk, not only at the end"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 2f213184-0ffd-40a3-8b2f-58f1332191ff
  modified: 2026-11-01T23:10:51.696Z
---

Every subagent/workflow agent must log incrementally: append results to files as they are produced (JSONL logs, running report files updated section by section), never hold findings until a final report.

**Why:** a Claude Code session crash on 2026-11-01 lost ~45 min of 16 attack agents' work that had not yet been written out; resumed runs could only build on what was on disk.
**How to apply:** in every agent prompt, require an append-as-you-go log path and a report file updated after each experiment; prefer resumable workflows (resumeFromRunId) and tell agents to read their own prior partial files on restart.
