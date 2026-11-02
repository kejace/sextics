# Wave 2 — `dense-rational` — lens EXCLUSION / RIGOR

STATUS: in progress (no witness; this lens does not search for witnesses). Started 2026-11-01.

Log: `/tmp/hodge/results/dense-rational/wave2-excl.log.jsonl` (one JSON line per experiment).
Packet: `/tmp/hodge/results/dense-rational/wave2-excl/` (scripts, receipts, README).

## Findings so far

1. (audit) The wave-1 degree-7 msolve runs (`polar_d7_*_msolve.log`) never finished: every log stops at the F4
   degree-3 step (one with `exit 143`). The wave-1 section drafts (§2.4 placeholders MS_*, §4.7, §5 item 1) treat the
   degree-7 tau-pure exclusion as established; it is NOT. Re-running on the remote with a new prime.
