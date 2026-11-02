# README-referee.md — diagonal-permutation, wave-2 exclusion packet (referee notes)

Written by the referee pass (kejace's team, Claude/Anthropic assistance, no human review) on 2026-11-02 UTC.
The agent's files in this directory were not modified; everything the referee produced lives in
`../wave2-referee/` (scripts `ref_*.py`, receipts `ref_*.json[l]`, copies of remote-only receipts).
In the research branch of the fork these paths are
`research/results/diagonal-permutation/wave2-excl/` and `research/results/diagonal-permutation/wave2-referee/`.

## 0. Verdict in one paragraph

The packet's **Theorem (no certificate with `deg D ≤ 6`, any Y-parity, over C)** is sound and was replayed
independently at a fresh prime (999999937, never used in this lane) with code that does not import the `hodge`
engine or `w2lib`: E1 kernels `ker([A]·: R_d → R_{d+6}) = 0` for `d ≤ 7` (both parities, good reduction in degrees
`d` and `d+6` checked), the degree-5 step by an **exact rational** rank (the 21 second partials `F_ij` are linearly
independent over Q), and the degree-6 Hessian relaxation rebuilt from scratch in all three sectors (17/134, 18/153,
35/404 unknowns/independent quadrics — the same counts as the packet), m-primary by msolve 0.6.5 **and** Singular 4.3.2
(`dim 0`, `vdim` 39/39/924).  Degree 7 is **open** (the remote msolve run was killed after 7079 s, `relax_d7e.json`:
`msolve_status -9`); degree 8 is open except for the canonical `D8` line (see §4); degree 9 kernels are exactly 21+21
and (01)-invariant.  Modular-only items are listed in §5.

## 1. Where the receipts are (gap in the agent's README)

The agent's README cites `e1k_*.json`, `d8_*.json`, `reptheory.json`, `fam_*.json`, `d8_u8_exact.json`.  Only the
local receipts (`relax_d6*`, `relax_d5_mixed.ms`, `ann_ann-local.json`, `singular_d6.out`, `d8check.out`) are in this
directory.  The remote-produced receipts were pulled to `/tmp/hodge/remote/diagonal-permutation/wave2-excl/`, which
`sync.sh` does not mirror; the referee copied them to `../wave2-referee/remote-receipts/`:

| file | content |
|---|---|
| `e1k_e1kF.json` | E1 kernels d = 1..9 at P1 = 1000000007 and d = 1..8 at P2 = 998244353, both parities, good-reduction flags, mod-p kernel bases at d = 8, 9 |
| `e1k_e1kS.json` | (01)-sector kernels d = 8..12 at P1 (good-reduction flag is `null` there; see §3) |
| `d8_d8G.json` | D8: exact lift of `u8` **failed** (4 primes); relaxation `Ψ(D8;u8) ∉ J` at P1 and P2; reduced E2' inconsistent at 3 primes; rank certification certified **0** of 159 null vectors |
| `reptheory.json` | symmetry/representation re-checks |
| `relax_d7e.json` | degree-7 Y-even relaxation: 60 unknowns, 480 independent quadrics, msolve killed (`-9`) after 7079 s, no basis |

`fam_*.json`, `d8_u8_exact.json`, `d9_*.json` **do not exist** (not computed).  The wave-1 rational `D8` (414 terms)
is `../D8_e1kernel.json` (key `D`, a sympy-parsable string).

## 2. Pinned environment

python-flint 0.9.0, msolve 0.6.5 (Ubuntu package 0.6.5-1build2), Singular 4.3.2, numpy 2.5.2, sympy 1.12, Python 3.
The `hodge` engine is the directory `research/engine/` of the research branch (it is not a git repository of its own;
the research-branch commit pins it).  Caveat: `Singular --version` does not exit (it waits on stdin); use
`echo 'quit;' | Singular -q`.  `w2lib.py` hard-codes `sys.path.insert(0, '/home/kejace/devel/sextics-work')`,
`LANE = '/tmp/hodge/results/diagonal-permutation'` and the instance record under `~/devel/sextics`; recreate those
paths (symlinks) or edit the three constants before running the agent's scripts.  The agent's scripts append to
`../wave2-excl.log.jsonl`.

## 3. How to replay (stranger's checklist)

Agent's route (needs the engine):
```
python3 w2_e1kernels.py all          # e1k_<job>.json  (remote-sized at d >= 8)
python3 w2_relax.py 6 e; python3 w2_relax.py 6 o; python3 w2_relax.py 6 m   # relax_d6*.{ms,gb,json}; msolve -g 2
python3 w2_singular_degbound.py      # Singular cross-check -> singular_d6.out
python3 w2_ann.py                    # ann_<job>.json
python3 w2_d8.py                     # d8_<job>.json (remote-sized; exact lift fails)
```
Referee's route (no engine; `cd ../wave2-referee`):
```
python3 ref_hess5.py                    # exact rank 21 of the second partials (degree-5 step)         -> ref_hess5.json
python3 ref_sectors.py                  # (tau x (01))-sector dimensions, forced kernels                 -> ref_sectors.json
python3 ref_e1.py 999999937 1 9         # E1 kernels d=1..9, good reduction, A*D8 in J mod p at d=8     -> ref_e1_p999999937.jsonl
python3 ref_relax6.py 999999937 e,o,m   # degree-6 relaxation from scratch + msolve -g 2                  -> ref_relax6_p999999937.json, .ms/.gb
python3 ref_singular_check.py ref_relax6m_p999999937.ms   # Singular dim/vdim of the same system          -> ref_singular_check.jsonl
python3 ref_ann.py 999999937            # 11 annihilators of the sum-principle catalogue                 -> ref_ann_p999999937.jsonl
python3 ref_quick.py                    # A not in J_6 (degree-0 denominators)                            -> ref_quick.json
python3 ref_d8.py 999999937             # D8: u mod p, Psi(D8;u) in J_16 mod p?  (~1 GB, 10-20 min)      -> ref_d8_p999999937.json
```
Reading the receipts: in `relax_d6*.json` / `ref_relax6_*.json`, `zero_dimensional` means every variable has a pure
power among the grevlex leading monomials of the reduced basis (for a homogeneous ideal this is `V(I) = {0}`);
`independent_equations` is the rank of the quadric system; `good_reduction_2d` / `good_reduction_12` is
`dim R_12 mod p = Koszul value` (866 even, 885 odd), the hypothesis of the transfer argument.  Singular's `dim(std(I)) = 0`
is the same statement.  In the E1 receipts, `kernel_dim` is `dim{D : A D ∈ J} − dim J_d`; the transfer
`kernel mod p ≥ kernel over Q` needs `good_reduction` in degree `d+6` (recorded).  In `e1k_e1kS.json` the
good-reduction flag is `null`; it can be confirmed from the data because the two sector values of `jac_dim_d6` add
up to the Koszul `dim J_{d+6}` (e.g. d = 10, odd target: 4788 + 4926 = 9714 = 10197 − 483).

## 4. Soundness of the transfer arguments (checked by hand)

* **Lemma "E2 modulo J".**  Rederived from the contract identity `E = A D² − D·N(F) + F(D div N − N(D)) = F² Q`
  (formal/README.md, q = 1): `N(F) = A D + F H`, `N = u + K + (H/6)E` by the regular sequence, the H-terms cancel,
  `D div K − K(D) ∈ J` for Koszul `K`; `Ψ(D;u) mod J` is independent of `u` and of `D mod F·S`; for `D = w(F)`,
  `Ψ ≡ −A·wᵀHess(F)w`.  Correct.
* **Zero kernel mod p ⇒ zero over Q (and C).**  Needs `dim J_{d+6} mod p = Koszul`; recorded.  Correct.
* **m-primary mod p ⇒ m-primary over Q.**  The quadric space mod p is the reduction of the saturated lattice of the
  rational quadric space when `dim J_{2d} mod p` is the Koszul value (the equations are the image of `(J_{2d})^⊥`);
  then `dim (S/I_p)_k ≥ dim (S/I_Q)_k` for all k.  Correct; over C by base change.
* **Euler gauge.**  `A·H(w+gE, w+gE) ≡ A·H(w,w) + 10 g A·w(F) (mod J)` and `w(F) ∈ J`, so the relaxation is
  invariant under `w ↦ w + gE`; for linear `w`, `w(F) ∈ (F)` iff `w ∈ Q·E` (no Koszul fields of degree < 5).  Correct.
* **D8 without an exact lift (README §1.3).**  Correct with one precision: in case (b) the mod-p test refers to the
  canonical rational representative *with respect to the reducer's standard monomials* (which is p-integral and reduces
  to `D8 mod p`), together with its F-multiples; case (a) should read "the rational kernel is zero", not
  "A·D8 ∉ J for the lifted polynomial".  Either way the conclusion "no certificate with `D = c·D8 + gF`" holds for that
  representative, and nothing is said about `D8 + w(F)` with `w ≠ gE`.
* **Farkas / rank certification.**  Not available for the D8 E2' system (0 of 159 null vectors certified), so the
  three-prime inconsistency is a screen.
* **Sum principle.**  `served(D)` depends only on `X ∩ V(D)` by the prescribed-support theorem (docs/01, which also
  covers arbitrary Y-parity); additivity over components follows from the weight argument (the image of `H_4(Y)` in
  `H^4(X)` equals the image of `H_4` of a resolution, which is a disjoint union over components); `served(D)/J_6 ⊆ Ann_6(D)`
  is E1; `Ann_6` computed mod p is an upper bound.  Rigorous for every product of catalogued factors (coprime to F).
* **Prescribed-support powers.**  Excluding `D^m` from an *exact* failure for `D` is covered by the same theorem; the
  wave-1 power statements rested on modular inconsistencies and remain screens.

## 5. What remains modular-only or open

* A·D8 ∈ J over Q **exactly**: not established (lift failed); verified modulo six primes (2147483647, 2147483629,
  1000000007, 998244353 by the agents; 999999937 by the referee).  Screened.
* Reduced E2' inconsistency of D8 (three primes, no rank certificate).  Screened.
* Degree 7 (all sectors) and the degree-8 family `λ D8 + w(F)`: **open** (no Gröbner basis obtained).
* Structured families d = 8..12 (`w2_families.py`) and individual degree-9 kernel elements (`w2_d9.py`): not run.

## 6. Referee D8 replay (added after the run)

`ref_d8_p999999937.json`: for the wave-1 rational `D8` reduced mod 999999937, `A·D8 ∈ J_14` (rank of `J_14^-` = 5052 =
Koszul, E1 nullity 954 = Koszul), a solution `u` of `u(F) = A·D8` was recomputed from scratch, and
`Ψ(D8;u) = D8·div u − u(D8)` (degree 16, Y-odd, 10067 terms) is **not** in `J_16` (rank of `J_16^-` = 9684 = Koszul, so good
reduction holds).  Third prime for this exclusion, independent code; scope as in §4.
