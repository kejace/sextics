# Packet `wave2-excl` — lane `cycle-q2`, lens EXCLUSION / RIGOR (agent session sextics-76, 2026-11-01)

Instance: `F = Σ(X_i⁶+Y_i⁶) + ½(X0³X1Y0Y1 + X1³X2Y1Y2 + X0X2³Y0Y2)`, `q = 2`, `A = X0²X1²X2³Y0²Y1²Y2` (`data/cycle-q2/instance.rec`).
Certificate identity (formal/README.md): `E = A D² − 2D Σ N_i ∂_iF + F(D Σ ∂_iN_i − Σ N_i ∂_iD) = F³Q`, `deg N_i = deg D + 7`, `deg Q = 2 deg D − 6`.
Engine: `hodge` v2 (`/home/kejace/devel/sextics-work`, sign fix of 2026-11-01), python-flint 0.9.0, Singular 4.3.2, msolve. Primes: P1 = 1000000007,
P2 = 998244353, 2^31−1 (Singular). All systems are built with `2F` (integer coefficients); consistency is unchanged.

A second EXCLUSION agent (`xlib-v1`: `xlib.py`, `a_*.py`, `b_*.py`, `c_*.py`, `commentary.md`, `farkas/`, `receipts/`) shares this directory; its files are
independent re-implementations and are not described here.

## Exact statements established (see `/tmp/hodge/notes/wave2-cycle-q2-excl.md` for the full report)

1. **No certificate with `deg D ≤ 5`.** `Ann_d([A]) = 0` for d ≤ 5 (`ann.sing` → `ann.out`, Singular at 2^31−1 and P1; `vdim J = 15625`) and the
   21 classes `[H_ij] ≡ −[A F_ij]` are independent in `R_16` (`hess5.py` → `hess5.receipt.json`, P1 and P2; `dim J_16 = 19398`).
   Transfer to characteristic 0: Jacobian ranks equal their regular-sequence values at these primes, and `rank_p ≤ rank_Q`.
2. **Exact Farkas certificates** (`w2lib.farkas`, replay with `w2lib.cert_from_file`): for every excluded fixed denominator, `certs/<tag>.json` holds the
   row monomials, the rational vector `y` (entries as strings), `yᵀb`, SHA-256, and the pivot data; `yᵀM = 0`, `yᵀb ≠ 0` is a prime-free proof that the
   character-sector system (lossless by averaging) is inconsistent over any field of characteristic 0. Candidates and ranks: `reverify.results.jsonl`.
3. **Degree 6, single character.** `hess6.py --deg 6` → `hess6_p1000000007.receipt.json` (dim K_6 = 62, 36 in J, 46 sectors) and the 426 quadrics
   (`.quadrics.json`, `.ms`, `.sing`). `sectors.py --deg 6` → log records `sector-d6`, `sector-d6-element`, Singular scripts `sectors/d6_*.sing`:
   every nontrivial sector with `dim K ≥ 2` has empty projective relaxation variety; every 1-dimensional sector element has a Farkas certificate.
   Exact (reduction-free) confirmation of the ten nontrivial sectors with dim K ≥ 2: `sector_exact.py --sector k1,k2,k3` (E1 kernel, Koszul-unique
   E1 solutions, normal forms modulo J_18, quadrics — all over ℚ; Singular dim 0). Trivial character: relaxation surface `V₀` (three explicit
   quadrics), reduced Koszul system (`redsys.py`, mod p; `redsys_lift.py`, exact): 60 exact rational degree-1 polynomial Farkas vectors whose cubic
   forms have the Euler point `D = 6F` as their only common zero in P⁵ (`polyfarkas/red_trivial6_exact.receipt.json`, `..._Q.sing`); `D = F` is
   inconsistent (`certs/D_equals_F.json`) and illegal. Mixed characters: `subsets6.py --size s` restricts the 426 quadrics to unions of s sectors
   (s = 2, 3 done: all unions of nontrivial sectors empty, unions with the trivial sector contained in V₀; receipts `sectors/subsets_d6_size*.json`).
4. Degrees 7+: `sectors.py --deg 7`, `cones.py --deg 7..9` (sector kernels, relaxation dimensions, element certificates) — see the log.

## Replay
```
cd /tmp/hodge/results/cycle-q2
python3 wave2-excl/hess5.py                       # 20 s
Singular -q wave2-excl/ann.sing                   # minutes
python3 wave2-excl/reverify.py --workers 8        # hours (resumable; writes reverify.results.jsonl, certs/)
python3 wave2-excl/hess6.py --deg 6 && python3 wave2-excl/gb_check.py wave2-excl/hess6_p1000000007.ms --threads 16   # Groebner (heavy)
python3 wave2-excl/sectors.py --deg 6 --workers 8 # sector analysis (resumable)
python3 wave2-excl/redsys.py --family trivial6 --degree 1 --check-degree1 ; python3 wave2-excl/redsys_lift.py --nvec 60
python3 -c "import sys; sys.path.insert(0,'wave2-excl'); import w2lib; from hodge import core; F,A,q,G=w2lib.load(); print(w2lib.cert_from_file('wave2-excl/certs/B6.json',F,A,q))"
```
Log of every experiment: `/tmp/hodge/results/cycle-q2/wave2-excl.log.jsonl` (lines with `"engine": 2`); remote jobs append to
`wave2-excl.remote.log.jsonl` on the compute box, merged by `merge.py`.

Note on certificate files: `certs/<tag>.json` stores the full rational vector `y` unless its text exceeds the storage cap (raised to 60 MB on
2026-11-01 22:20; a few early degree-6/7 certificates were stored with `y = null` and `y_omitted_size`); such a certificate is regenerated
deterministically by `w2lib.farkas(F, A, q, D, p, G, save=path)` with the stored `D` and prime, and all stored vectors replay with `cert_from_file`.
