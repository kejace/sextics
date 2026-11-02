# Coordination note (wave 2, lane symmetric-sextic-residual)

Two agents appear to be running searches on this lane concurrently (remote jobs `search-symmetric-sextic-residual-{struct,vp-d8triv0,vp-d8sgnX1}`
and `excl-ssr-relax1` were submitted by another process at 00:24 UTC; this agent's jobs use the prefix
`search-symmetric-sextic-residual-w2b-*` and the scripts `w2lib.py`, `w2_varpro.py`, `w2_kernels.py`, `w2_single*.py`,
`w2_struct.py`, outputs in `w2/*.log.jsonl`, merged into `wave2-search.log.jsonl`).

Scope of this agent (to avoid duplication): invariant (`triv`) varpro sectors d = 6..14, the unrestricted
Jacobian w-space (`jac`) at d = 5..8, subgroup-type (S2xS2, S3xS2, S2xS3) kernels/solves/varpro at d = 6..10,
single-D solves of all 1-dim-character kernel classes d = 8..12, structured vertex families (7 types) d = 7..12,
kernel dimensions d = 4..14 and the Fermat-deformation scan.  Please do not reuse these job names or overwrite
`w2/` files; please append your own batch records to `wave2-search.log.jsonl` with a distinct `job` field.
