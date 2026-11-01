---
name: feedback-single-session
description: "Never run two copies of the same Claude session; check `pgrep -fa \"claude --resume\"` before resuming workflows"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 2f213184-0ffd-40a3-8b2f-58f1332191ff
  modified: 2026-11-01T23:51:39.694Z
---

On 2026-11-01 a `claude --resume <id>` was started while the earlier resume of the same session (pts/3) was still alive and running wave 1. The new copy received "workflow didn't finish before the previous session ended" notifications and almost re-launched everything, which would have duplicated 50+ agents on an already saturated box.

**Why:** background workflows and agents belong to the process that launched them; a second process cannot see them and misreads them as crashed.
**How to apply:** before resuming any workflow after a "stopped" notification, run `pgrep -fa "claude --resume"` and `ps` for recent jobs in /tmp/hodge/results; if another copy is alive, coordinate with the user to keep exactly one session, and never launch duplicate workflows. Heavy compute goes to the remote box via `python3 -m hodge.rjob` (see [[hodge-sextics-project-setup]]).
