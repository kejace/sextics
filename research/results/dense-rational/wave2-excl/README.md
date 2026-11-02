# dense-rational — wave-2 exclusion packet (lens: EXCLUSION / RIGOR)

Instance: `dense-rational` (`F` = Fermat + (1..17)/9973 on all 230 non-pure even sextic monomials, `q = 1`,
`A = X0^3 Y0^3`), variables `X0,X1,X2,Y0,Y1,Y2`, `tau(X,Y) = (X,-Y)`.  Engine: `/home/kejace/devel/sextics-work/hodge`
(python-flint 0.9.0).  All scripts live in `/tmp/hodge/results/dense-rational/wave2_excl_*.py|.sh`; receipts, msolve
inputs/outputs and logs live in this directory.  One JSON line per experiment: `../wave2-excl.log.jsonl`.

## Exact statements proved here (over C)

Notation: `S = C[X0..Y2]`, `J = (F_0..F_5)` the Jacobian ideal, `R = S/J`, `S_d^{+/-}` the Y-even / Y-odd forms,
`(J:A)_d = {D in S_d : A D in J}`.  A certificate for `(F,1,A)` with denominator `D` of degree `d` is
`E1: N.grad F = A D`, `E2: D div N - N.grad D = F C` (Euler-normalized form), `gcd(D,F) = 1`.

* **T1 (E1 bookkeeping).** `(J:A)_d = J_d` for `d <= 9` in both Y-parities; for `d = 10` the kernel of
  `[A]: R_10^{+} -> R_16^{-}` is exactly 294-dimensional and of `[A]: R_10^{-} -> R_16^{+}` exactly 261-dimensional.
  Proof: direct Macaulay nullspace mod `p` (`wave2_excl_e1kernel.py`, p = 2147483587; wave 1 used p = 2147483647,
  the dossier p = 101, 10007).  Rigor: `(J:A)_d(Q) ∩ Z^n` is a saturated lattice, so it reduces injectively into
  `(J:A)_d(F_p)`; `dim J_d(F_p) = dim J_d(Q)` whenever the modular count attains the Koszul value `6 dim S_{d-5}`
  (`d <= 9`); and for `d = 10` the rational kernel is at least `dim R_10 - dim R_16` by dimension count.  Hence
  every E1-admissible denominator of degree `<= 9` is a polar divisor `D = w(F) = sum_i w_i F_i`, `deg w_i = d - 5`.
* **T2 (degree <= 4).** `J_d = 0` for `d <= 4`, so no denominator of degree `<= 4` passes E1 (both parities, hence also
  mixed parity, since E1 splits by parity).
* **T3 (degree 5).** `D = v.grad F` (`v in C^6` constant).  E1 numerators are `N = A v + K grad F` (Koszul), so E2 forces
  `A v^T Hess(F) v in J_10`, and `(J:A)_4 = 0` (T1) gives `v^T Hess(F) v = 0` as a polynomial.  The linear map
  `Sym^2(C^6) -> S_4`, `S -> sum_ij S_ij F_ij`, has rank 21 (exact rational rank, `wave2_excl_local_small.py`), so
  `v v^T = 0`, `v = 0`.  (The packet of `fqlx` tested the six coordinate directions; this covers every direction.)
* **T4 (degree 6, all three sectors).** `D = w(F)` with `w` a linear vector field (36 coefficients; sectors even /
  odd / mixed by Y-parity of the components).  E2 forces `w^T Hess(F) w in J_6` (both parities), because `[A]` is
  injective on `R_6^{+/-}` (T1).  After the Euler gauge (`w -> w + c z`), these are homogeneous quadratic systems
  (even: 218 quadrics / 17 unknowns; odd: 218 / 18; mixed: 426 / 35).  msolve (`-g 1`, p = 1073741789) returns DRL
  leading ideals containing a pure power of every variable (`quadsys_d6_*_ltcheck.json`), so each system has only the
  trivial solution over `F_p`-bar, hence over `C` (a projective variety defined over `Q` with a `C`-point has an
  `F_p`-bar point for every `p` at which the equations are `p`-integral; the mod-`p` equations are the reductions of the
  rational ones because the pivot minor of `J_6` is a `p`-unit).  So **no denominator of degree 6 exists**, in any
  parity sector, and `deg D >= 7` for every certificate — an independent re-proof of `fqlx`'s R3.
* **T5 (degree 7, tau-pure sectors) — PENDING msolve** (`quadsys_d7_{even,odd}_p1073741789.ms`: 483 quadrics / 60
  unknowns each).  Same reduction with `w` quadratic and `[A]` injective on `R_8^{+}` (T1).  If both leading ideals
  contain pure powers of all variables: no certificate with a Y-even or Y-odd denominator of degree 7, i.e. **no
  denominator of degree 7 whose support is tau-stable** (a legal degree-7 `D` is squarefree by the prescribed-support
  theorem and `deg D >= 7`, so tau-stable support <=> tau-pure `D`).  The mixed sector (951 quadrics / 120 unknowns)
  is out of F4 reach (generic emptiness only at Macaulay degree 5, 2.3e8 columns) and is screened numerically only.

## Modular / numerical screens (not proofs)

See `../wave2-excl.log.jsonl` entries with `exp` in {`target-image`, `fermat-first-order-deformation`,
`e2-normalized-fixed-D`, `e1-kernel-structure`, `lm-hessian-relaxation-d7`, `lm-hessian-relaxation-d8`}.

## Replay

```
cd /tmp/hodge/results/dense-rational; export PYTHONPATH=/home/kejace/devel/sextics-work
python3 wave2_excl_local_small.py                       # T3 + A not in J_6
python3 wave2_excl_e1kernel.py d parD 2147483587 [--basis]   # T1 (d = 1..10)
python3 wave2_excl_quadsys.py d sector 1073741789       # T4/T5 systems; then
msolve -v 1 -t 8 -g 1 -f wave2-excl/quadsys_d{d}_{sector}_p1073741789.ms -o OUT; python3 wave2_excl_ltcheck.py MS OUT
```
