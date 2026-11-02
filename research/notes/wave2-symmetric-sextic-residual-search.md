# Wave 2 search report: `symmetric-sextic-residual`

STATUS: NO WITNESS FOUND (wave 2 started 2026-11-01 ~19:10 EST; report rewritten after every batch).

Log: `/tmp/hodge/results/symmetric-sextic-residual/wave2-search.log.jsonl` (one JSON line per batch).
Scripts: `/tmp/hodge/results/symmetric-sextic-residual/w2lib.py` and `w2_*.py` (remote jobs via `hodge.rjob`).

## 0. Corrections to the wave-1 record

* `scan_served.py` (wave 1, last 260 log lines, family `quot-curve6X`): the lines "A ENTERS THE SPAN" and
  "A-system consistent=True" are artefacts. The candidate `1EX(F)+1EY(F)` is `6F` (Euler), so its served
  space is everything and the cumulative span is polluted from that point on; and `S.solve()` was called on
  the system that still contained the target columns, so "consistent" is vacuous. Real content of that
  scan: every tested D of degree 3-6 serves nothing beyond `J(F)_6^G`.

## 1. Plan (wave 2)

1. E1-admissible denominator spaces per degree d = 4..14 in the character sectors triv/sgnX/sgnY/sgnsgn
   of S3xS3 (both Y-parities), plus S2xS2 / S3xS2 / S2xS3 types at d = 10 and the unrestricted kernel at d = 9, 10.
2. Single-D exact (mod p, CRT-lifted) solves of every kernel class + random combinations (never done for
   the non-invariant classes).
3. Character-reduced variable-projection (varpro) joint search over the full E1-admissible family
   (kernel + J^chi_d) per (chi, d, pD), thousands of restarts on the remote.
4. Structured families (cones over plane curves with the three G-stable vertex planes, line-vertex cones,
   discriminant/Hessian-block products, products of kernel elements) via restricted E1 kernels + solves.
5. Rational reconstruction + exact re-solve of every candidate with residual < 1e-10.

(Results sections are appended below as batches finish.)
