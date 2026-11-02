# Wave 2 — lane `quadratic-cubic-pullback` — lens EXCLUSION / RIGOR

STATUS: in progress (no witness; re-verification and degree-bound jobs being launched)

Started 2026-11-01 (wave 2).  Log: `/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl.log.jsonl`.
Packet: `/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl/`.

## Results so far
- Local check: F = sum Q_i^3 exactly; the fixed-D kernel of a polar denominator D = dF/dX0 (d=5) has the explicit
  rational basis (G z, 0), G in S_5^odd (126) and (g F, D div g - g.grad D), g constant on the Y-slots (3): 129 elements,
  equal to the modular nullity (engine reports 128 = rank of augmented matrix convention).  Consequence: the three
  wave-1 degree-5 polar inconsistencies are rigorous over Q (modular rank of M = rational rank, sandwich argument).
- Engine note written to /tmp/hodge/notes/engine-bugs.md (reported nullity off by one for inconsistent systems).
- Dossier correction: E_I F = sum_{i in I} z_i dF/dz_i lies in J(F) for EVERY index set I (not only I = all), so the
  block-Euler denominators EI([0,1,2]) etc. listed in wave-1 `listC` were never tested; they are natural
  degree-6 polar candidates (the Fermat join X0^6+X1^6+X2^6 is exactly E_X F_0 / 6) and are queued now.
