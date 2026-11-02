# Referee replay guide — lane `quadratic-cubic-pullback`, wave-2 EXCLUSION packet

Instance: `F = sum_{i=0}^{5} Q_i^3` (six Y-even quadrics, `data/quadratic-cubic-pullback/instance.rec`), `q = 1`,
`A = X0^3 Y0^3`.  Variables `X0 X1 X2 Y0 Y1 Y2`; "Y-parity" = parity of the total Y-degree (A is Y-odd, F is Y-even).
Certificate (formal contract): nonzero homogeneous `D` coprime to `F`, `N_i` homogeneous of degree `deg D + 1`, and
`E := A D^2 - D N(F) + F (D div N - N(D)) = F^2 Q`.  `J = J(F) = (F_0..F_5)`, `R = S/J`, `K_d := {D in S_d : A D in J}`.

This file tells a stranger how to replay every rigorous claim of the packet (`README.md`, `../wave2-excl.log.jsonl`,
`receipts/`) and the referee's own independent replays (`referee/`).  The referee's verdicts are in
`/tmp/hodge/posts/quadratic-cubic-pullback/referee-report.md` (mirrored in the research branch under `posts/`).

## Pinned engine

* python 3 with **python-flint 0.9.0** (`nmod_mat`, `fmpq_mat`, `fmpz_mat` ranks / rref / nullspace),
* **msolve 0.6.5** (`/usr/bin/msolve`) — used only for characteristic-p screening; note the two pitfalls below,
* **Singular 4.3.2** (`/usr/bin/Singular`) — all candidate-locus statements over Q use `std`, `dim`, `vdim`, `radical`
  (`primdec.lib`) — the packet says 4.3; the referee's box has 4.3.2,
* the lane helper `qcp.py` and the engine `hodge.core` (`~/devel/sextics-work`) are used only for parsing the instance
  and for dict-polynomial arithmetic; every rank/kernel computation in `referee/` builds its own matrices.

msolve pitfalls (both confirmed by the lens; documented in `/tmp/hodge/notes/engine-bugs.md`): write squares as `w^2`
(msolve parses `w*w` as `w`), and keep the characteristic below `2^30` (primes near `2^31` are silently mishandled).

## Transfer principles used (checked by the referee)

1. **Modular rank is a lower bound.**  For an integer (p-integral) matrix, `rank_p <= rank_Q`.  The Jacobian block
   `J_k^{par}` has rational rank equal to the Hilbert value `dim S_k^{par} - dim R_k^{par}` with
   `HS(R) = (1-t^5)^3 (1-u t^5)^3 / ((1-t)^3 (1-u t)^3)`, `u^2 = 1` (regular sequence; attaining the value mod p proves it).
   When `rank_p(J_k) = Hilbert`, then for any p-integral `T`: `T notin J mod p  =>  T notin J over Q`; and any kernel
   computed as `#columns - (rank_p[J ; T-block] - rank_p[J])` is an **upper bound** for the rational kernel.
2. **E1 kernel.**  `dim K_d(Q) <= dim K_d(p)` (principle 1) and `J_d subset K_d(Q)` with `dim J_d(Q) = Hilbert`; so
   `dim K_d(p) = dim J_d` proves `K_d(Q) = J_d`.  The lower bound `dim K_d/J_d >= dim R_d - dim R_{d+6}` is rank–nullity.
   Parity: `J` is tau-stable and `A` is Y-odd, so a mixed-parity `D` is in `K_d` iff both parity parts are; E1 facts
   therefore cover mixed-parity denominators.
3. **E1 => polar.**  Every certificate has `A D in J` (reduce `E = F^2 Q` modulo `F`; `D` coprime to `F`).  If `K_d = J_d`
   then `D = w . grad F` for a polynomial vector field `w` of degree `d - 5`.
4. **E2 => Hessian condition.**  With `N = A w + H z/6 + c . grad F` (`c` antisymmetric, degree `d - 4`; completeness
   uses that `dF` is a regular sequence so syzygies are Koszul), all `H` terms cancel and `E = F (Phi(w) + L(c,w))`
   with `Phi(w) = D (grad A . w + A div w) - A w(D)`, `L(c,w) = sum_{i,j} F_j [(d_i c_ij) D - c_ij d_i D]`.  Since `L` and
   `F Q` lie in `J`, a certificate forces `A . (w^T Hess(F) w) in J_{2d}`.  Trivial solutions: `w = lambda z` gives
   `D = 6 lambda F`, which is **not coprime to F** and so is not a certificate; the condition is invariant under
   `w -> w + h z` because `z^T Hess(F) w = 5 w(F) in J` and `z^T Hess(F) z = 30 F in J`.
5. **Candidate cone.**  Relaxing `W = w w^T` gives a linear condition whose kernel `K_H(p)` is an upper bound for
   `K_H(Q)`; exact lifts of a mod-p basis that are verified exactly to satisfy the condition pin `K_H(Q)`.  The cone
   `{w : w w^T in K_H(Q)}` is then decided by Singular over Q (`dim std = 0` means `{0}`; `dim 1` means a finite union of
   lines through 0; a random affine chart `L = 1` with radical a single rational point plus `dim std(I + (L)) = 0`
   means exactly one line).  Equivalent cheaper form (used by the aggregator and by the referee): when `K_{d'} = J_{d'}`
   for `d' = 2d - 6`, multiplication by `[A]` is injective on `R_{d'}`, so the condition is `w^T Hess(F) w in J_{d'}`
   (both parity parts) and the quadrics have exact rational coefficients after an exact rref of the `J_{d'}` blocks.
6. **m-primary leading ideal mod p => no projective zero over C.**  The Macaulay matrices of the exact integer quadrics
   reduce mod p to those of the mod-p ideal (the reduced quadrics span the reduction of the exact span once the `J_{d'}`
   block is at its Hilbert rank), so `HF(Q[w]/I_Q)(e) <= HF(F_p[w]/I_p)(e)`; a leading ideal containing a pure power of
   every variable makes the right side `0` for large `e`, hence `I_Q` is m-primary.  Together with principle 4's
   translation invariance, `dim std(I + (c_E)) = 0` on a slice transverse to the Euler direction proves the cone is
   exactly the Euler line.

## Replaying the packet's own computations (lens scripts)

    cd /tmp/hodge/results/quadratic-cubic-pullback
    python3 wave2-excl/wave2_excl_reverify.py                        # wave-1 structural facts at 2147483629 / exact
    python3 wave2-excl/wave2_excl_kernel.py k 2 "5:0:1,5:1:1,6:0:1,6:1:1"   # K_d = J_d (writes wave2-excl.k.jsonl)
    python3 wave2-excl/wave2_excl_relax.py r5 5 all                 # exact integer relaxation, degree 5, all parities (~15 s)
    python3 wave2-excl/wave2_excl_hess.py h5 5 all 268435399         # K_H = 0 at degree 5 (3 s)
    python3 wave2-excl/wave2_excl_hess.py h6 6 1 268435399           # degree 6, Y-odd: dim K_H = 18, writes the KH json
    python3 wave2-excl/wave2_excl_hess.py h6 6 0 268435399           # degree 6, Y-even: dim K_H = 36
    python3 wave2-excl/wave2_excl_hess_exact.py x6 6 1 268435399     # exact lift + Singular over Q (needs file hessL6 -> rename h6)
    python3 wave2-excl/wave2_excl_hesstest.py 7 1                    # 17 structured degree-7 candidates (J_14 reducer ~4 min)
    Singular -q wave2-excl/receipts/wave2-excl-hess-exact-hessX6-d6-parall-charts2.sing   # mixed-parity cone, charts + radical

Receipts: `wave2-excl.log.jsonl` (master, 90 lines), `receipts/wave2-excl.hessL6.jsonl`, `receipts/wave2-excl.hessX6.jsonl`,
`receipts/wave2-excl.hesstest.jsonl`, `receipts/hessX6_all_charts.out`, KH basis files `receipts/wave2-excl-hess-hessL6-d6-par*-KH.json`.

## Referee's independent replays (`referee/`, fresh prime 1000000009, exact where stated)

    cd /tmp/hodge/results/quadratic-cubic-pullback/wave2-excl/referee
    python3 referee_replay.py exact      # F = sum Q_i^3; rank of 21 second partials = 21 (exact); A, delta vs J_6^odd exact ranks 18/19/19/20
    python3 referee_replay.py hilbert    # own Koszul series; J_k ranks k = 6,10,11,12 both parities = Hilbert
    python3 referee_replay.py kernels    # K_d for d = 0..7, both parities (dim K_d = dim J_d, J blocks certified)
    python3 referee_replay.py deg5       # Hessian relaxation at degree 5: K_H = 0 (joint, both parity blocks) + two-line hand argument
    python3 referee_replay.py khodd      # the lens's Y-odd degree-6 K_H basis has no diagonal entries (=> w = 0 by hand)
    python3 referee_replay.py deg6       # exact 426 quadrics in 36 unknowns (own build) -> referee_deg6_exact.sing
    Singular -q referee_deg6_exact.sing  # FULL36 dim 1, SLICE_c35 dim 0, EVEN_SECTOR dim 1 (+slice 0), ODD_SECTOR dim 0
    python3 referee_hess7.py hess7 kernel8   # 3 structured degree-7 candidates (J_14^odd rank 5052 certified); K_8 both parities
    python3 referee_hess8.py             # remote: 2 structured degree-8 candidates (J_16^odd rank 9684 certified)

Output: `referee/referee-replay.jsonl` (one JSON line per check), `referee/referee_replay.out`, `referee/referee_hess7.out`,
`referee/referee_deg6_singular.out`, `referee/referee-replay-remote.jsonl` (degree 8, if the remote job ran).
