# Wave 2 verification — lane `reflection-coupled-residual` (adversarial referee)

STATUS: IN PROGRESS — confirmed 9 / refuted 1 (typo-level) / weakened 0 / unverified (pending) — see table; last update see bottom.

Scope: every claim labelled rigorous/proven/exact/computed in `dossier-reflection-coupled-residual.md`,
`attack1-reflection-coupled-residual.md`, the lane logs under `/tmp/hodge/results/reflection-coupled-residual/`
(no `wave2-reflection-coupled-residual-*.md` and no draft post existed at 00:30 UTC; re-checked periodically).
Independent code (no engine dependency) lives in `/tmp/hodge/results/reflection-coupled-residual/verify/`
(`vlib.py`, `v_symmetry.py`, `v_sectors.py`, `v_e1.py`, Singular scripts `sing_*.sing`); every check is appended
to `/tmp/hodge/results/reflection-coupled-residual/wave2-verify.log.jsonl`.
Primes used by the lane: 2147483647, 2147483629 (sibling, python-flint), 32003 (dossier, NumPy). Primes used here:
1000000007 (own flint code), 1000003 (Singular), plus exact char-0 Singular where stated.

## Logic conventions checked

* **Modular-rank direction.** For an integer matrix, `rank_p <= rank_Q`. The E1-admissible space is
  `K = {D : A·D ∈ J}` with `dim K = dim S_d − rank[J;A·S_d] + rank J`, so `dim K_Q <= dim K_p` holds **only if**
  `rank_p J_{d+6} = rank_Q J_{d+6}`. That equality is certified in every run below by comparing `dim (S/J)` mod p with the
  Fermat count (all exponents <= 4) in both degrees (`hilbert_ok`). With it, `ker_p = 0 ⇒ ker_Q = 0` (zeros are
  rigorous) and `ker_p = 1 ⇒ ker_Q <= 1` (an upper bound; the exact rational element must be verified separately).
* **"E1 false mod p" is rigorous** under the same Hilbert-function check. **"E1 true mod p" is not** a proof of E1
  over Q; K8's exact membership is checked separately (Singular, char 0).
* **"Inconsistent mod p" for a fixed-D (N,Q) system is NOT automatically rigorous** (needs `rank_p = rank_Q` of the
  homogeneous system, whose kernel depends on D). The dossier labels those results [reported]/[computed, mod p] and
  never calls them rigorous — correct. The E2-mod-J obstruction system, by contrast, has a kernel of structural
  (Koszul) dimension, so a modular inconsistency there does transfer (see claim table).
* **Character restriction is lossless for E1** (J and A are Klein-homogeneous, so `A·D ∈ J` iff each character
  component of D satisfies it); restricting to Y-even D is lossless for *existence* (replace D by D·τ(D)) but not for
  degree bounds (a Y-even exclusion up to degree k only excludes arbitrary D up to degree k/2 via the product; the
  all-character computation is what gives the unrestricted bounds).

## Claim table

| # | claim | source | verdict | evidence | notes |
|---|---|---|---|---|---|
