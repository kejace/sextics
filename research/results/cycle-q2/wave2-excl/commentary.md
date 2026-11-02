## Summary of established results (hand-written; updated after each experiment)

- E0 (validation): independent solver `xlib.py` (no engine code except parsing) reproduces the accepted `block-q2` and `fermat-q2` witnesses exactly (identity verified in own arithmetic; sector solver returns SOLVED-EXACT), and its matrix agrees with direct polynomial evaluation on a random (N,Q) for `cycle-q2`, `D = F_{X2}`.
- Rigor note: a mod-p inconsistency (rank/aug = r/r+1) is NOT by itself a characteristic-zero proof (p could divide all r x r minors of M). The packet therefore attaches an exact Farkas certificate (integer vector Y with Y^T M = 0, Y^T b != 0, checked in exact integer arithmetic) to every exclusion, which proves inconsistency over Q and hence over C.
