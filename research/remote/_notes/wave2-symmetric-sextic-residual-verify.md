# Wave 2 — `symmetric-sextic-residual` — lens VERIFY (adversarial referee)

STATUS: 0 confirmed / 1 refuted / 0 weakened / 0 unverified (started 2026-11-02 00:20 UTC; table rewritten after every check).
Log: `/tmp/hodge/results/symmetric-sextic-residual/wave2-verify.log.jsonl` (one JSON line per check).
Scratch: `/tmp/hodge/verify-ssr/` (independent code; no engine imports for the load-bearing checks).

## HEADLINE (read first)

* **REFUTED — dossier §1.4 sector table is mislabelled (polynomial character used as residue character).**
  `mult.py` counts polynomial characters of R_d (fixed monomials of the permutation action). The residue
  character is χ_poly·det, and det(π_X,π_Y) = sgn·sgn. So the residue sector of [A] (sgn,sgn,τ=+) is the
  **polynomial (triv,triv) Y-odd** column: Hodge numbers **(17, 52, 17)**, not (1,10,1). The (1,10,1) row is the
  polynomial (sgn,sgn) sector = residue (triv,triv) sector, spanned at (3,1) by Δ_XΔ_Y, i.e. exactly what the
  symmetric Fano quotient serves. Consequences: "[A] is the unique (3,1)-class of its character", "K3-type rank-12
  Hodge structure", "τ-invariant part of H^4/U is (1,10,1)", "H^4(V)_− is (1,10,1)⊕(0,6,0)" are all wrong; the
  attack-1 note ("17-dimensional sgn⊗sgn residue sector") is right. The construct lens §0 has copied the wrong
  numbers and should be corrected. (Independent recomputation below.)

## Claim table

| # | claim | source | verdict | evidence | notes |
|---|---|---|---|---|---|
