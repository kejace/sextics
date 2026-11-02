# Wave 2 — `one-reflection-skew` joint search (report rewritten after every batch)

STATUS: NO WITNESS FOUND (searches running on the remote; the shared box is saturated by ~65 other jobs,
load ~520 on 128 cores, so one 4-CPU chain was forced through the capacity gate and three 8-CPU chains are
queued behind it). Pipeline validated end-to-end on an accepted witness (family-hesse-web): exact 6→5
reduction, operator (rel 2.6e-15 at the known D), LM variable projection finds the known D from random
starts (4/8 restarts, full E1-kernel family), snap-to-rational + exact confirm recovers exactly the accepted
denominator and verifies the 6-variable lift (`core.verify_witness`).

## 0. Setting (exact reduction used throughout)

`F = f(X0,X1,X2,Y0,Y1) + Y2^6`, `A = X0^3 Y0^3` (free of Y2), `q = 1`. Averaging a certificate over
`μ6(Y2)` and expanding in powers of `F` (details in §6) gives the **exact equivalence**

> a certificate exists  ⇔  a rational 5-variable field `u = N/D` (deg D = d, deg N_i = d+1) satisfies
> `A D² − D N(f) + (6/5) f (D div N − N(D)) = 0`     (*)

and then `(D², (N_i D)_{i<5}, Y2 (D div N − N(D))/5, Q = 0)` is a strict 6-variable certificate of degree
`2d`. Necessary condition (E1): `A D ∈ J_f`. `D` can be taken in any `(ρ, σ)` character
(`ρ`: sign of (Y0,Y1), `σ`: sign of (X1,Y1)); `N` and the rows then live in determined sectors.
Verified on `family-hesse-web` (accepted degree-6 witness, same target): the reduced field satisfies (*)
exactly, the float operator gives relative residual `3.7e-15` at it, `exact_confirm` recovers and
verifies the 6-variable certificate.

Code: `/tmp/hodge/results/one-reflection-skew/w2lib.py` (operator linear in D over monomial coordinates,
variable projection with Tikhonov normal equations, full Golub–Pereyra Jacobian LM, rationalisation +
LLL, exact modular confirm with CRT), `w2_search.py` (multiprocess restarts), `w2_dims.py`,
`w2_struct.py`, `chains/chain{1,2,3,4}.sh`. Log: `wave2-search.log.jsonl` (one line per batch).

## 1. E1-admissible denominator spaces (task 1)

(pending: `w2_dims.py` on the remote; wave-1 mod-p numbers for the 5-variable ρ-even sector:
kernel = 3, 13, 47, 167, 307, 479, 677 at d = 5..11 with J_d = 3, 13, 39, 89, 178, 314, 506, i.e.
non-Jacobian directions 0, 0, 8, 78, 129, 165, 171.)

## 2. Optimizer calibration (needle = accepted degree-6 witness of family-hesse-web, same target A)

* First attempts (L-BFGS, LM-Kaufman, LM with the full Golub–Pereyra Jacobian — the latter checked
  against finite differences to 3e-9) all converged to recurring local minima (rel 3.056e-2, 4.256e-2,
  2.3e-3, …), 0/26 hits. Basin probe: LM started 5°–45° from the needle still ended at rel 5e-2.
* Cause: the family had the f-multiple subspace `f·S_{d−6}` projected out (as the engine's varpro does
  with `exclude_F_multiples`). The true `D*` lies in the E1 kernel (projection residual 8e-16) but has a
  component along `f`, so the restricted span did not contain it. **Removing the F-multiple subspace is
  wrong**: it removes solutions `D0 + F·G` that are not divisible by F. Recorded in
  `/tmp/hodge/notes/engine-bugs.md`; the lane driver now searches the full kernel (in the strict
  5-variable problem there are no spurious solutions: `(fG, fW)` solves iff `(G, W)` does).
* With the full kernel (K = 32): **4/8 LM restarts converge to the needle** (rel 2e-13 … 6e-8,
  |cos(u,u*)| = 1.0000); the other restarts end at the usual local minima. So random-restart LM has a
  ~50% per-restart hit rate on a 32-dimensional family with an 8-term solution. Hit rates at larger K
  are measured by `calib-hesse-d6` (400 restarts, LM vs L-BFGS) and `calib-xtype-d6` on the remote.

## 3. Search batches (filled from the log)

(pending)

## 4. Structured families (task 4)

(pending: `w2_struct.py`, cone families)

## 5. Rational reconstruction attempts (task 5)

Finding from the calibration needle: a numerical zero (rel 2e-13) does **not** rationalise directly — the
zero set is positive-dimensional (LM Jacobian at D* has two near-null directions beyond scaling, singular
values 2.7e-10 and 7.9e-6; e.g. direction `X0⁴X1² − X1²Y0⁴`), so a generic numerical zero is an irrational
point of a family of denominators. Remedy implemented (`w2lib.snap_confirm`): drop monomials with
|coefficient| < τ (τ ladder 1e-6 … 1e-1), refit by LM inside the sub-family, rationalise
(`limit_denominator` ladder + LLL), confirm exactly (E1 mod p, fixed-D rref mod p, CRT lift, exact 5-variable
identity, `core.verify_witness` on the 6-variable lift). On the needle: τ = 1e-4 keeps 8 monomials, refit
rel 4e-15, exact solve → the accepted D, verified. Skew-lane attempts: none yet (no restart below 1e-7).

## 6. Derivation of the reduction (for the record)

Let `v = N/D` with `K_A(v) := A − v(F) + F div v = F² Q/D²`. Average over `g: Y2 ↦ ζ^k Y2` (F, A invariant;
`K_A` is natural): `v̄` is invariant, `K_A(v̄) = F² · (avg of g·(Q/D²))`, denominator `∏ g(D)²` coprime to F.
Invariance: `v̄_i = u_i(x, w)` for i < 5 and `v̄_5 = Y2 a(x, w)`, `w = Y2^6`. With `F = f + w`, regard
everything as rational in `(x, F)`: `K_A = A − ũ(f) − 6(F − f)ã + F(div_x ũ + ∂_F ũ(f) + ã + 6(F − f)ã_F)`.
`F² | K_A` ⇒ coefficients of `F⁰, F¹` vanish: `A − u₀(f) + 6 f a₀ = 0` and `div u₀ − 5 a₀ = 0`
(`u₀ = ũ(x,0)`), hence (*) for `u₀`. Conversely `v = (u, Y2 div u/5)` gives `K_A(v) = A − u(f) + (6/5) f div u = 0`
with denominator `D²` coprime to `F` (F is monic in Y2, D is Y2-free).
