# wave2-excl — `eight-plane-configuration` exclusion/rigor packet (2026-11-01)

Instance: q = 1, A = X0^3 Y0^3, F = the 236-term tau-even sextic of data/eight-plane-configuration/instance.rec.
Engine: ~/devel/sextics-work/hodge (core, fixed_d, jacobian); python-flint 0.9; msolve 0.6+ (grevlex, -g 1).
Primes: wave 1 used 2147483647 (local) and 2147483629 (remote); wave 2 uses 2147483587 everywhere, plus exact rationals where stated.

## Statements (what is proved vs screened)

S1 (rigorous). A notin J(F): mod 5, F = Fermat, J = (Z_i^5), and dim J_6 = 126 on both sides (saturation argument), so the class is nonzero.

S2 (rigorous). For d <= 9 and either Y-parity, {D in S_d : A*D in J(F)} = J(F)_d over Q.
   Proof: dim_Q ker <= nD - rank_p[A m | J_{d+6}] + dim_Q J_{d+6} = kerdim_p whenever rank_p J_{d+6} = dim_Q J_{d+6} (exact Koszul value,
   w2x_hilbert.py); kerdim_p = dim_Q J_d = 3|S_{d-5}| in every computed row; J_d is contained in the kernel. Computed independently with
   three primes and two code paths (wave-1 nullspace: e1kernel.jsonl; wave-2 normal forms: e1ladder.2147483587.jsonl).

S3 (rigorous, d <= 9). Every certificate denominator of degree d <= 9 is D = w(F), w_i in S_{d-5}, with gcd(D,F) = 1 iff w notin {h*Euler}.
   A certificate exists iff the Koszul-reduced system (R4 in the report) is consistent; a necessary condition is A*Hess_F(w,w) in J_{2d}.

S4 (rigorous). No certificate has deg D <= 6, for D of ANY Y-parity (pure sectors and mixed):
   d <= 4: S2 gives D in J_d = 0.  d = 5: Hess_F(w,w) = 0 is forced (J_4 = 0) and has only w = 0 (linearization rank 6 = dim Sym^2 per
   sector, 21 = dim Sym^2(C^6) mixed; hess.jsonl, hess_mixed.jsonl).
   d = 6: msolve on the quadric ideal (hess_d6_par*_p2147483587.ms and _p2147483647.ms, leading ideals .lm): Hilbert function [1,17,18,2,0] /
   [1,18,18,2,0]; mixed parity (hessmixed_d6_p2147483587.ms, 35 variables, 405 quadrics): quotient dimension 923. In all cases (S/I)_k = 0 for
   some k mod p, hence over Q (Macaulay rank can only drop mod p), hence V' = {0}; D = w(F) with w in V = T is an F-multiple.
   NOTE: the parity WLOG D -> D*tau(D) only shows that SOME certificate has even D; degree bounds therefore need the mixed sector too.

S5 (pending / partial). d = 7: 483 quadrics in 60 variables per parity (hess_d7_par*_p2147483587.ms); msolve job excl-epc-hess7.
   Subfamilies at d = 7 excluded (both parities; hess_sub.jsonl + hsub_d7_*.lm): grad (linearization), arr, blockX, blockY, coefX, coefB0 (msolve).
   d = 8 subfamilies excluded (hess_sub8_local.jsonl, wave-1 kernel K_10^+ with p = 2147483629): Hamiltonian fields C grad(G) for G in
   {K_i,H_i,diag_i,Q_iQ_j,Q_i^2} (linearization), coefB0 (msolve). Factor table (sum principle): exp c-factor-annihilator in the log.

S6 (screened). Individual Jacobian denominators (d = 5..9) tested with the Koszul-reduced E2 system: all inconsistent (two primes; d = 5
   X-partials inconsistent over exact rationals). Random and structured fields: see wave2-excl-remote-*.jsonl (jsample).

## Files
w2x_lib.py (helpers: exact parity-refined Hilbert function, J generators, normal-form reducer, logging)
w2x_hilbert.py (exact Hilbert function; consistency check of every wave-1 kernel row)
w2x_e1ladder.py (normal-form E1 ladder, third prime; receipts e1ladder.2147483587*.jsonl; bases annA_d*_par*.p*.pkl on the remote)
w2x_hess.py (Hessian relaxation on W_d/T_d; msolve inputs hess_d*_par*_p*.ms; receipts hess*.jsonl)
w2x_hess_sub.py (same on structured subfamilies; hsub_*.ms; hess_sub*.jsonl)
w2x_hess_mixed.py (mixed-parity Hessian relaxation; hessmixed_d*.ms; hess_mixed*.jsonl)
w2x_koszul.py (Koszul-reduced E2: modular and exact-rational; receipts in wave2-excl.log.jsonl / wave2-excl-remote-*.jsonl)
w2x_jsample.py (random + structured Jacobian fields, multiprocessing)
w2x_hess8_local.py (d = 8 Hessian relaxation on subfamilies via the degree-10 E1 kernel; hsub_d8_*.ms/.lm; hess_sub8_local.jsonl)
w2x_factors.py (annihilators and served spaces of low-degree factors)
remote_batch2a.sh / remote_batch2b.sh (remote batches)
run_msolve.sh (msolve wrapper writing JSON receipts)

## Replay
python3 w2x_hilbert.py; python3 w2x_e1ladder.py 5 6 7 8 9 10; python3 w2x_hess.py <d> <par>; msolve -g 1 -f hess_d6_par0_p2147483587.ms -o out.lm;
python3 w2x_koszul.py exact local partials; python3 w2x_koszul.py 2147483587 local J6cands; python3 w2x_jsample.py 8 8 <tag>.
