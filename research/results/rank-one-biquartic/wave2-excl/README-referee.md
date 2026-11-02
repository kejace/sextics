# README-referee.md — independent replay guide for the `wave2-excl` packet (lane `rank-one-biquartic`)

Referee pass of 2026-11-01/02 (Claude, Anthropic, for `kejace`; no human review; nothing Lean-verified; no
certificate claimed).  This file tells a stranger how to re-derive every rigorous statement of the packet from the
receipts in this directory, and which statements are modular evidence only.  The lane's own description of the
scripts is `README.md`; the claim-by-claim verdicts are in `referee-report.md` (copy of
`/tmp/hodge/posts/rank-one-biquartic/referee-report.md`).

## 0. What is proved (characteristic zero) and what is not

Target: `F = f(X) + g(Y) + (5/7) L(X)^2 M(Y)` (54 terms, see `data/rank-one-biquartic/instance.rec`),
`A = X0^3 Y0^3`, `q = 1`.  Variables `X0,X1,X2,Y0,Y1,Y2`; `J = (F_0..F_5)` the Jacobian ideal (a regular
sequence: `F` is smooth, `F == sum Z_i^6 (mod 5)`).

* **PROVED over Q (replayed by the referee):** there is no certificate `(D, N, Q)` with `deg D <= 6`, for every
  nonzero `D` coprime to `F`, of any Y-parity (mixed included), with `Q` arbitrary.
  - `deg D <= 4`: the E1-kernel `V_d = {D in S_d : A D in J_{d+6}}` is zero.
  - `deg D = 5`: `V_5 = J_5 = <F_0..F_5>`; six exact rational certificates `lambda^T E2 = a_k^2`.
  - `deg D = 6`: `V_6 = J_6` (36-dim); 58 exact rational certificates in the J-relaxed weight-2 ideal
    (25 coordinate squares `a_k^2`, 33 quadratic forms), then a characteristic-zero Groebner basis (Singular)
    shows that the 35 coordinates of `V_6 / <F>` all vanish on the solution set, i.e. the only solutions are
    `D = c F` (not coprime).
* **Modular evidence only (NOT proved over Q):** the dimensions `dim V_d / J_d = 9+9, 27+27` at `d = 8, 9`;
  the `(J:A)` generator profile 6/3/3/456 in degrees 5/7/8/10; `Delta' = 4 t^3 L^6 M^3 + 27 f^2 g in J`;
  the single-prime fixed-denominator inconsistencies of wave 1 (all of degree 6 and hence superseded by the proof
  above).
* **Not excluded at all:** every degree `>= 7`.  In degree 7 the E1-admissible space is exactly
  `V_7 = J_7 (+) K_7` with `dim K_7 = 3` (Y-even, exact generators `nonj_d7.json`), `V_7^- = J_7^-`; the
  Hessian relaxation certifies none of the 123 coordinate squares mod p and its mod-p Groebner step timed out
  (`hpipe_d7.out`).  Nothing is known about degrees `>= 7` beyond E1 profiles.

## 1. Pinned environment

* Ubuntu, python 3.12, **python-flint 0.9.0** (`pip install python-flint==0.9.0`), **sympy 1.12**
  (only for parsing in the referee scripts), **Singular 4.3.2** (apt; one remote log line mentions 4.4.1 — any
  4.x works), **msolve 0.6.5** (apt), PyYAML.
* The lane scripts import the engine `hodge` from `/home/kejace/devel/sextics-work` (research branch copy:
  `research/engine/`) and read the instance record through `hodge.core.load_instance` from the upstream checkout
  (`ROOT` in `hodge/core.py`; the record is `data/rank-one-biquartic/instance.rec` of `hodge-computer/sextics`).
  A stranger should (i) check out `hodge-computer/sextics`, (ii) put `research/engine` on `PYTHONPATH` or edit the
  `sys.path.insert` line at the top of `w2lib.py`, (iii) point `ROOT` at the checkout.  The referee scripts parse
  the record themselves (`REC` constant at the top of `referee_replay.py`).
* Everything below runs on a laptop; the longest single step is `k7exact.py 7` (17 min) / `referee_replay.py e1`
  (11 min, dominated by a dense `8568 x 5544` rank at `d = 7`).

## 2. Referee replays (independent code; what each one establishes)

All commands are run inside this directory.  `referee_replay.py` uses its own polynomial arithmetic and its own
parser; from the packet it only takes the row-index -> monomial map and the basis-label order needed to read the
certificate JSON files.  Receipts: `referee_replay.out`, `referee_replay.json` (part d5),
`referee_replay_{d6,e1,k7}.json`, `referee_membership.json`, `referee_d6_charQ.{sing,out}`,
`referee_hessrelax6.g1.out`, `referee_fam_d5_p1000000009.{json,log.jsonl}`.

```
python3 referee_replay.py d5      # 12 s  : the 6 certificates of fam_d5_p1000000007_all.certs.json, checked against
                                  #         MY weight-2 E2 columns (21 a_k a_l + 540 a_k c + 126 Q): lambda^T E = a_k^2
python3 referee_replay.py d6      # 40 s  : the 58 certificates of hpipe_d6.blk1.certs.json(.gz): lambda kills all 2376
                                  #         Y-odd generators m F_i of J_12 and pairs with the TRUE E2 columns
                                  #         P_kl + P_lk (dense B_D(N), not the sparse Hessian representatives) to the
                                  #         quadratic form recorded in hpipe_d6.C.json
Singular -q referee_d6_charQ.sing # 0.1 s : char-0 Groebner basis of the 25 linear + 33 quadratic certified generators
                                  #         (regenerated from hpipe_d6.C.json): affine dim 1, 0 of 35 coordinates of
                                  #         V_6/<F> outside the radical  => V(I) = {D = cF}
python3 referee_replay.py e1      # 11 min: dim{D in S_d : A D in J_{d+6}} (all parities) at p = 1000000009 for d = 4..7:
                                  #         0, 6, 36, 129 with the J block at its regular-sequence rank in every degree
                                  #         => V_d = J_d over Q for d <= 6; dim_Q V_7 <= 129
python3 referee_replay.py k7      # ~2 min: the 3 generators of nonj_d7.json satisfy N(F) = A D exactly over Q and are
                                  #         independent modulo J_7 at p = 1000000009  => V_7 = J_7 (+) K_7, dim K_7 = 3
python3 referee_membership.py     # ~12 min: membership table at p = 1000000009 (non-membership rigorous over Q when
                                  #         the J block has its regular-sequence rank; membership mod p is evidence)
python3 fam_run.py 5 1000000009 --tag referee_fam_d5_p1000000009   # 8 s: lane code at an unused prime: rank 507,
                                  #         all 21 monomials a_k a_l in the row space
msolve -g 1 -t 1 -f ../../_aggregator/hess6/hessrelax6-rank-one-biquartic.ms -o referee_hessrelax6.g1.out  # 6 s
                                  #         aggregator's mixed-parity degree-6 relaxation at p = 10^9+7: 707 leading
                                  #         monomials, a pure power of each of the 35 variables (max exponent 3)
```

Unzip first if needed: `gunzip -k hpipe_d6.blk1.certs.json.gz` (the 30 MB original is excluded from the research
branch by size; the gzip copy is 12 MB).

## 3. Regenerating the packet's own receipts (lane code)

```
python3 fam_run.py 5 1000000007 --cert                      # 10 s : degree 5, exact certificates
python3 hess_pipeline.py 6 3 --sparse --tag hpipe_d6        # 7 min: degree 6 (steps A-D), writes hpipe_d6.*
python3 k7exact.py 7                                        # 17 min: exact K_7 generators -> nonj_d7.json
python3 hess_pipeline.py 7 3 --sparse --nonj nonj_d7.json --tag hpipe_d7   # 2 h; step B (Singular mod p) times out
python3 verify_e1kernel.py 1000000007                       # ~15 min: E1 profile d = 1..9, both parities
Singular -q quotJA_p65521.sing                              # 15 min: (J:A) profile at p = 65521
```

## 4. The logical chain (why the receipts prove the statement)

1. **Euler gauge.** For a certificate, `D (A D - N(F)) == 0 (mod F)`; `F` is irreducible and `D` coprime, so
   `N(F) = A D + F h`.  `N' = N - (h/6) Euler` leaves the cleared identity unchanged (the two `F`-terms cancel:
   `div(phi E) = (deg phi + 6) phi`, `(phi E)(D) = d phi D`) and has `N'(F) = A D`.  Hence `A D in J`
   (E1) and the identity becomes `E2: D div N' - N'(D) = F Q`.
2. **Koszul.** `N' = N_0(D) + K(c)` with `K(c)_i = sum_j c_ij F_j`, `c` antisymmetric of degree `d - 4`
   (first syzygies of a regular sequence are Koszul).  E2 is a weight-2 polynomial system in `(a, c, Q)`.
3. **E1-kernel over Q.** `dim_Q V_d <= dim_p V_d` whenever the `J_{d+6}` block has its regular-sequence rank mod p
   (then the syzygy part of the kernel is the same over both fields); `J_d subset V_d` always.  Receipts:
   `verify_e1kernel_1000000007.json`, `referee_replay_e1.json` (p = 1000000009).
4. **Certificates (Farkas duals).** A rational `lambda` on the degree-`2d` monomials with `lambda^T E2 = a_k^2`
   proves `a_k^2 in I`, so `a_k = 0` on every solution.  For `d = 6`, `lambda` is additionally required to kill
   every generator of `J_12` of the relevant parity: it then kills all Koszul and `Q` columns (which lie in
   `J_12`), so the same `lambda` is a certificate for the full E2 system.  Each `lambda` is verified exactly over
   Q (python-flint `fmpq`, Fractions).  Receipts: `fam_d5_*.certs.json`, `hpipe_d6.blk1.certs.json(.gz)`,
   `hpipe_d6.C.json`.
5. **Groebner over Q.** The ideal generated by the certified linear forms (`u` with `u^2 in I`) and certified
   quadratic forms has every coordinate of `V_6/<F>` in its radical (Singular, char 0): the solution set is the
   spurious line `D = cF`.  Receipts: `hpipe_d6.D.sing/.out`, `referee_d6_charQ.sing/.out`.
6. **Independent modular route (degree 6).** With `D = w(F)`, `w` a linear vector field (36 coefficients;
   `V_6 = J_6` makes this exhaustive), E2 mod `J_12` gives `A w^T Hess(F) w in J_12`, and `[A]` is injective on
   `(S/J)_6` (that is `V_6 = J_6` again), so `w^T Hess(F) w in J_6`: 419 quadrics in the 35 coordinates of the
   Euler slice `c_55 = 0`.  Their leading ideal mod `10^9+7` contains a pure power of every variable, so
   `(I_p)_71 = (S_p)_71`, a full-rank statement for a Macaulay matrix that persists over Q (rank can only grow),
   provided the mod-p quadrics are the reductions of the Q-quadrics; that holds because the standard monomials
   of `(S/J)_6` mod p form a Q-basis with p-integral coordinates when `dim J_6` mod p is the regular-sequence value
   (saturated-lattice argument).  This route was not needed for the proof but agrees with it.
7. **Non-membership transfer.** `g notin J` mod p with `rank_p J_k = rank_Q J_k` implies `g notin J` over Q;
   membership mod p transfers nothing.  Receipts: `referee_membership.json`, `quotJA_p65521.out` (Hilbert
   function of `S/J` printed and equal to the regular-sequence values), lane log `reverify-J-membership`.

## 5. Known gaps and caveats

* Degree 7 is open: `hpipe_d7` certified 0 of 123 coordinate squares mod p; the 12 MB Singular script
  `hpipe_d7.B.sing` (mod p, all 129 residual variables = 123 failing coordinates + 6 spurious centres, 924 quadratic
  forms) timed out after 3600 s.  A finer relaxation (keeping the Koszul unknowns, or working modulo `(F) + J^2`)
  or a remote Groebner run is needed; "degree 7 screened" must NOT be claimed.
* Degree 8 exact generators of `K_8` (modular dimension 9 + 9) were never computed (`remote_jobs.sh` not run).
* `Delta' in J` and the `(J:A)` profile are single/double-prime modular results.
* The wave-1 fixed-denominator "inconsistent" verdicts were obtained with engine v1/v2 at a single prime; they are
  now implied by the rigorous degree-6 exclusion and should not be cited independently.
* The exclusion says nothing about denominators of degree `>= 7`; the first unexcluded family is the 129-dimensional
  `V_7` (123 coordinates modulo the spurious `F S_1`).
