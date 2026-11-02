# Referee addendum to the wave-2 exclusion packet — `four-plus-two-residual` (q = 1, A = X0^3 Y0^3)

Written by the referee pass (Claude, Anthropic, acting for `kejace`), 2026-11-01/02. The agent's files in this directory are
unchanged; everything the referee added lives in `referee/` and in this file. Claim-by-claim verdicts: `referee/referee-report.md`.

## What this packet establishes (after refereeing)

Certificate (contract, q = 1): nonzero homogeneous D coprime to F, deg N_i = deg D + 1, deg Q = 2 deg D - 6, and
A D^2 - D N(F) + F (D div N - N(D)) = F^2 Q.  Taking this modulo F (F is irreducible: smooth) gives
E1: N(F) = A D + m F, i.e. [A][D] = 0 in R(F) = S/J(F); the rest is E2.

* **T1 [EXACT, VERIFIED]** R(F) = R(H) (x) R(B) (F = H + B in disjoint variables), and ker([A]: R(F)_d -> R(F)_{d+6}) =
  sum_k ker(a.: R(H)_{d-k} -> R(H)_{d-k+6}) (x) R(B)_k with a = X0^3 Y0^3 in the core.  Core kernel dimensions modulo J(H)_k
  ((Y0,Y1)-even / odd): k <= 5: 0/0; 6: 20/13; 7: 30/30; 8: 41/34; 9: 38/38; 10: 36/31; 11: 26/26; 12: 19/16.
  Consequences: **no E1-admissible D of degree <= 4 (any parity); every E1-admissible D of degree 5 is a polar F_v = v.grad F.**
* **T2 [EXACT, VERIFIED]** **No certificate with deg D = 5 exists over C**, so **deg D >= 6 for every certificate** of this instance.
  (For D = F_v: N = A v + (m/6) z + K, K Koszul; E2 becomes F_v A_v - A F_vv + F_v div K - K(F_v) = F Q; relaxing v(x)v and
  v(x)l to free tensors gives an exact linear system whose Sym^2-projection W is spanned by the eight mixed core x binary
  coordinates e_a e_b (a in {X0,X1,Y0,Y1}, b in {X2,Y2}); every v_i^2 vanishes on W, so W has no rank-one point.)
* **T3 [EXACT, VERIFIED, count corrected]** Among degree-6 Jacobian-type denominators D = sum_{a,c} w_ac z_c F_a (36 parameters)
  the G-eigenvectors (G = diagonal symmetry group, order 24) fall into **11** character classes (10 + 8 + 9 x 2 parameters; the
  draft's "13" double-counted two halves).  The sector-restricted Hessian relaxation of each class has no rank-one point except
  the Euler line D = 6 lambda F (not coprime to F).  **No G-eigenvector Jacobian-type denominator of degree 6 has a certificate.**
* **T4 [RANK, VERIFIED on samples]** Fixed rational D, sector-restricted full system: inconsistent mod p AND the exact rational
  nullity of M (n_m + Koszul nullity, valid for d <= 8) equals the modular nullity => inconsistent over Q, hence over C.
  Results per D in `remote/excl-symsweep-local-k6-8.jsonl` (see the referee report for the completed counts).
* **T5, T6: NOT RUN** (only the splitting prime of B was computed).  Nothing about the sum principle or generic degree 7-10
  members is established by this packet.
* **Structural D (referee addition)** B(X2,Y2) and H(X0,X1,Y0,Y1) (both in J(F), both coprime to F) are inconsistent over C
  [RANK]: rank_p M = 2298, rank_p [M|b] = 2299, exact nullity 304 = 226 (m) + 78 (Koszul), p = 10^9+9 (`referee/out/ref_t4.jsonl`).
  Note: the engine's dense solver reports the rank of the AUGMENTED matrix and `nullity = cols - rank_aug`; its
  "rank 2299 / nullity 303" for B and H (dossier A.1, aggregator log) therefore means rank M = 2298, nullity 304.

## Transfer arguments used (checked by the referee)
1. Hilbert functions: dim_p (S/J)_k >= dim_Q (S/J)_k >= CI_k, so equality with the complete-intersection series modulo one
   prime proves the series over Q (partials regular over Q and mod p).  Used for H and B.
2. Kernel dimensions: rank_p <= rank_Q for integer matrices.  "0 mod p" transfers to Q; for k >= 5 the map a. is surjective
   mod p, which forces rank_Q = dim R(H)_{k+6}, so dim ker_Q = dim R(H)_k - dim R(H)_{k+6} exactly (the fmpz computations agree).
3. Hessian relaxation: if (v, l, Q) is a certificate then (v(x)v, v(x)l, Q) solves the relaxed homogeneous linear system, so
   v(x)v lies in the projection W; the rank-one locus of W is the variety of the quadrics annihilating W.  "Ideal contains
   v_i^2 for all i" or "zero-dimensional (pure powers lead the Groebner basis)" => only v = 0.  Linear algebra over Q, Groebner
   bases over Q; conclusions hold over C.  Euler line D ~ F is excluded by relative primality.
4. Sector restriction (T3, T4): for a G-eigenvector D the certificate system is G-equivariant; averaging projects any solution
   to the forced isotypic component (characters chi(N_i) = chi(A) + chi(D) - chi(F) + wt(z_i), chi(Q) = chi(A) + 2 chi(D) - 2 chi(F)).
   The particular solution A w already lies in the sector; the Koszul map is equivariant.  Lossless.
5. [RANK]: inconsistent mod p means rank_p [M|b] = rank_p M + 1 <= rank_Q [M|b].  If an explicit rational solution space of
   dimension n with n = cols - rank_p M is exhibited (here: N = (m/6) z + K with k_ij of degree d - 4 <= 4, injective for
   d <= 8, and the exact rank of the (k, Q)-system), then rank_Q M = rank_p M and the system is inconsistent over Q, hence C.
   Farkas duals were not produced; the nullity argument is an equivalent rigorous route.
6. Not used here: msolve; the prescribed-support theorem (book ch. 1) — the results above are stated per denominator.

## Replaying the agent's packet (pinned: python-flint 0.9.0, Singular 4.3.2, sympy 1.12; msolve 0.6.5 is NOT used)
The scripts hard-code two paths: the engine `/home/kejace/devel/sextics-work` (synced to `research/engine/` on the research
branch) and the packet `/tmp/hodge/results/four-plus-two-residual/wave2-excl` (this directory); `hodge.core.load_instance`
reads `/home/kejace/devel/sextics/data/<id>/instance.rec` (the upstream checkout, `ROOT` in `engine/hodge/core.py`).
A stranger therefore does:

    git clone -b research https://github.com/kejace/sextics sx && cd sx/research
    mkdir -p /tmp/hodge/results/four-plus-two-residual && cp -r results/four-plus-two-residual/wave2-excl /tmp/hodge/results/four-plus-two-residual/
    # point the engine at your upstream checkout: edit ROOT in engine/hodge/core.py; put engine/ on PYTHONPATH
    export PYTHONPATH=$PWD/engine
    cd /tmp/hodge/results/four-plus-two-residual
    python3 wave2-excl/w2_kernels.py excl-kernels            # T1 (3 primes + exact; ~10 min)
    W2_RELAX_ONLY=1 python3 wave2-excl/w2_d5.py excl-d5        # T2 (exact; seconds; needs Singular)
    python3 wave2-excl/w2_d6jac2.py excl-d6jac2                # T3 (exact; all 15 (sub)families; sympy Groebner)
    python3 wave2-excl/w2_symsweep.py excl-symsweep 8 6 8     # T4 (8 workers, k = 6..8; ~1 h)
Outputs are appended to `remote/<job>.jsonl` (JSON lines; `receipts/` holds the copies the agent froze at 20:25).

## Independent referee replays (no engine, no w2lib; `referee/`)
    export HODGE_INSTANCE_REC=/path/to/data/four-plus-two-residual/instance.rec   # default: the kejace checkout
    python3 referee/ref_t1.py     # T1: Hilbert functions, core kernels (mod 10^9+9, exact at k = 5..8), 6-variable kernels d <= 6
    python3 referee/ref_t2.py     # T2: degree-5 relaxation, exact, own column construction, Singular rank-one test (+ control)
    python3 referee/ref_t3.py     # T3: character classes; Singular minAssGTZ / dim on the agent's saved annihilator quadrics
    python3 referee/ref_t4.py     # B, H, random members of the 10- and 8-parameter classes, degree-6 kernel reps, [RANK] at 10^9+9
    python3 referee/replay_d6jac2.py "MERGED Y-odd"   # re-runs the agent's exact degree-6 relaxation into referee/out/
Logs: `referee/out/*.jsonl` (one JSON line per experiment), Singular inputs `referee/out/*.sing`.
