# Engine maintainer report (wave 2)

STATUS: active on the VM since 2026-11-02 00:13 UTC.  Baseline tests ALL OK.  rjob local mode shipped.
Now benchmarking the dense/compressed path on real d=8 systems to design the sparse elimination.

Work dir: /tmp/hodge/results/_engine/ (log.jsonl, bench_*.log).  Changelog: ~/devel/sextics-work/CHANGELOG.md.
Inbox: /tmp/hodge/notes/engine-maintainer-inbox.md, bugs: /tmp/hodge/notes/engine-bugs.md (polled every 15 min).

## Done

1. rjob runs locally on the VM (`python3 -m hodge.rjob submit --lane L --name N -- cmd` works here;
   before, the ssh-to-self failed with "Permission denied (publickey)").

## In progress

2. Sparse modular elimination for large fixed-D systems (task 2): measuring structure first.

## Queue

3. served-space tool (target_image for a denominator list, one pass per lane).
4. varpro improvements (on request).
5. proof-packet tool (verify -> leangen -> leancheck -> sandboxcheck).
