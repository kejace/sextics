# Wave 2 — `eight-plane-configuration` — lens EXCLUSION / RIGOR

STATUS: in progress — no witness. Rigorous so far: every certificate has deg D >= 7 (all D of degree <= 6 excluded over C); for deg D in {7,8,9} the denominator must be a Jacobian denominator D = w(F) satisfying the Hessian condition; d = 7 msolve run pending.

Log (append-as-you-go): /tmp/hodge/results/eight-plane-configuration/wave2-excl.log.jsonl
Packet: /tmp/hodge/results/eight-plane-configuration/wave2-excl/ (scripts w2x_*.py, run_msolve.sh, JSON receipts)
Remote jobs (lane eight-plane-configuration): excl-epc-ladder5to10, excl-epc-ladder11to12 (running), excl-epc-hess7 (queued, waiting for remote capacity).

## 0. Exact reductions used (all over Q; proofs in §5)

(R1) Parity: WLOG D has pure Y-parity (D -> D*tau(D)); E1 and E2 split by parity.
(R2) E1-kernel rigor: {D in S_d : A D in J} = J_d over Q iff the modular computation gives kerdim_p = dim_Q J_d AND rank_p J_{d+6} = dim_Q J_{d+6} (exact Koszul value). Both were checked for every wave-1 row (new script w2x_hilbert.py).
(R3) Jacobian-type reduction (d <= 9): every E1-admissible D of degree d is D = w(F) = sum_i w_i F_i, w_i in S_{d-5} with parity(w_{X_i}) = par(D), parity(w_{Y_i}) = par(D)+1. gcd(D,F)=1 <=> w not in T_d = {h*Euler : h in S_{d-6}}.
(R4) Koszul-reduced E2 (exact equivalence): for D = w(F), a certificate exists iff there are c_ij in S_{d-4} and Q with
     [D div N0 - N0(D)] + sum_{i<j} [D (F_j d_i c_ij - F_i d_j c_ij) - c_ij (F_j D_i - F_i D_j)] = F Q,  N0 = A w.
(R5) Hessian relaxation (necessary): A * Hess_F(w,w) in J_{2d}. For d <= 7 (mult. by A injective on R_{2d-6}): Hess_F(w,w) in J_{2d-6}. The condition is invariant under w -> w + h*Euler, so it is a condition on the class of w in W_d / T_d.

## 1. Re-verification of wave-1 rigorous claims (goal a)

| claim (wave 1) | wave-1 evidence | wave-2 check | verdict |
|---|---|---|---|
| A not in J(F) | rank 18 vs 19 mod 2^31-1 | mod-5 Fermat argument (exact, dim count J_6 = 126 both sides) | rigorous, agrees |
| Ann_R(A)_d = 0 for d <= 9, both parities | nullspace method, p = 2147483647 (local) and 2147483629 (remote, d = 9, 10) | (i) implied rank_p J_{d+6} equals exact Koszul dim for all 23 logged rows (w2x_hilbert.py) -> modular->rational argument airtight; (ii) independent normal-form method with p = 2147483587 for d = 1..4 locally (agree: 0), d = 5..12 running remotely | rigorous, agrees so far |
| dim Ann_R(A)_10 = 294 (even) / 261 (odd) = forced lower bound | p = 2147483629 | being recomputed (ladder job) | pending |
| six partials dF/dZ_i inconsistent | full system mod 2147483647, rank 1318 | Koszul-reduced system over exact rationals (flint fmpq): dF/dX0,dF/dX1,dF/dX2 inconsistent (rank 102 vs 103); Y-partials queued remotely | rigorous (3/6 fully exact), agrees |
| E_0(F), E_1(F), E_2(F), eps_X(F), X0 F_X0 inconsistent | full system mod 2147483647, rank 2359 | Koszul-reduced system mod 2147483587: all inconsistent (rank 370/371), plus Y0F_Y0, X-rotation, R_0(F), X0F_Y0, Y0F_X0, R_0+R_1 | agrees (two primes); superseded by §2 (all of J_6 excluded) |

## 2. Rigorous lower bound on deg D (goal b)

* d <= 4: no E1-admissible D at all (Ann = 0, J_d = 0). Rigorous.
* d = 5: D in J_5 = span(F_X0,F_X1,F_X2) (even) or span(F_Y0,F_Y1,F_Y2) (odd), D = w(F) with w constant. Hessian relaxation: Hess_F(w,w) in J_4 = 0, i.e. the second directional derivative of F vanishes identically; the linearization Sym^2(W) -> S_4 has rank 6 = dim Sym^2(W) (mod p, hence over Q), so w = 0. **All degree-5 denominators excluded** (not just the six basis elements).
* d = 6: D = w(F), w a linear vector field (18 parameters per parity), T_6 = span(Euler) (even sector only). Quadric ideal I = (coordinates of NF_J Hess_F(w',w')) on W' (17 / 18 variables): linearization rank 135 / 153 (kernel 18 in both sectors — the twisted-Euler tensors sum_j (Z_j d_i)(Z_l d_j)); msolve (p = 2147483587, grevlex, -g 1): Hilbert function of S/I = [1,17,18,2,0] (even) and [1,18,18,2,0] (odd), i.e. (S/I)_4 = 0, quotient dimension 38 / 39. Hence V' = {0} over F_p, so the degree-4 Macaulay matrix has full rank mod p, hence over Q, hence V'_C = {0}: **all coprime degree-6 denominators excluded** (both parities).
* d = 7: W' has dimension 60 per parity, 483 independent quadrics each (linearization rank 483 < 1830); semi-regular degree of regularity 4 -> msolve job excl-epc-hess7 (queued).
* d = 8, 9: Hessian relaxation leaves 468 quadrics in 156/159 variables (d = 8; dreg >= 8 by the Hilbert-series bound — not Groebner-feasible) and ~208 quadrics in ~350 variables (d = 9; positive-dimensional). Complete exclusion is out of reach by this route; see §3 for sampling/structured results.

Current rigorous statement: **every certificate has deg D >= 7; if deg D in {7,8,9} then D = w(F) is a Jacobian (tangency-locus) denominator with A*Hess_F(w,w) in J.**

## 3. Jacobian denominators d = 7, 8, 9 and the d = 10 kernel (goal c)
(pending)

## 4. Packet (goal d)
wave2-excl/w2x_hilbert.py (exact parity-refined Hilbert function + wave-1 consistency), w2x_lib.py, w2x_e1ladder.py (normal-form E1 ladder), w2x_hess.py (Hessian relaxation, msolve input), run_msolve.sh, w2x_koszul.py (Koszul-reduced E2, modular and exact), receipts: e1ladder.<p>.jsonl, hess.jsonl, hess_d*_par*_p*.ms/.lm, wave2-excl-remote-*.jsonl.

## 5. Proof sketches of the reductions
(R3) Ann_R(A)_d = 0 (R2) gives D in J_d; split w by parity. If w(F) = hF then (w - hE/6)(F) = 0 so w - hE/6 is a syzygy of the regular sequence (F_i) with coefficients of degree d-10 < 0, hence 0.
(R4) N(F) = AD + Fg forces N = Aw + gE/6 + K (Koszul exactness, c_ij in S_{d-4}); the g-terms cancel identically in D div N - N(D) - Dg; div(kappa_ij) = 0 gives the displayed formula.
(R5) Reduce (R4) modulo J: D, K(D), F Q, (grad_w w)(F) are in J, and w(w(F)) = (grad_w w)(F) + Hess_F(w,w). Hess_F(hE + w', hE + w') = 30 h^2 F + 10 h w'(F) + Hess_F(w',w') = Hess_F(w',w') mod J.
Rigor of msolve certificates: with rank_p J_{2d-6} = dim_Q J_{2d-6}, the mod-p quadric space is the reduction of the rational one; (S/I)_k = 0 mod p means the degree-k Macaulay matrix has full rank mod p, hence over Q, so V'_C = {0}.
