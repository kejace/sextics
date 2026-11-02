# Wave 2 theory group: DEFORMATION
Lanes: rank-one-biquartic, off-paired, one-reflection-skew, diagonal-permutation, dense-rational.
Work dir: /tmp/hodge/results/_cross/theory-deformation (log.jsonl). Engine v2.

## STATUS (updated after every result; newest first)
- [init 2026-11-02] Read all five dossiers + attack1 + wave2 reports + book ch.3,3a,4. Framework below.
  Box load ~540 at start (oversubscribed); heavy rank jobs deferred until load<140 (monitored).
- [RESULT A, analytic] rank-one-biquartic: first-order **E1** obstruction VANISHES trivially
  (D1=(2/3)L^2M), via Euler identity Sum_{i<3} X_i L_i = 2L. So the degree-6 support f CAN be
  corrected at first order in the E1 (Jacobian) sense; the genuine first-order obstruction is E2
  (the divergence/B-operator), computed below.

## Framework (first-order deformation of supports)
For F_t = F_0 + t C (C the symmetry-breaking even sextic) and a t=0 certificate (D0,N0,Q0) with
cleared identity E(F_0,D0,N0,Q0)=0, a fixed-degree certificate family F_t exists to first order iff
the linear system in (D1,N1,Q1) obtained from d/dt E(F_t,D_t,N_t,Q_t)=0 at t=0 is consistent:

  2 A D0 D1 - D1 N0(F0) - D0 N1(F0) - D0 N0(C) + C[D0 div N0 - N0(D0)]
    + F0[D1 div N0 + D0 div N1 - N1(D0) - N0(D1)] = F0^2 Q1 + 2 F0 C Q0.          (ORD1)

Two necessary sub-conditions:
  (E1') A D1 - Sum_i G_{i,0} C_i in J(F0)   [the derivative class [A C] must stay in [A].(R_{F0})]
  (E2') the divergence/B-defect must be absorbed mod F0^2.
Book 4.4 input (i): the IVHS class-derivative in direction C is multiplication by C in the Jacobian
ring, (R_{F0})_6 -> (R_{F0})_12, [A] -> [A C]; (E1') is exactly its solvability.

## Per-lane plan / state
- rank-one-biquartic: explicit t=0 cert (D0=f). E1' vanishes (RESULT A). Running E2' (ORD1) with D1
  free over S_6^+ : full fixed-degree-6 consistency w/ the C-source RHS. [queued, awaiting load]
- diagonal-permutation: continuation from sub-family O2=O3=0 (sep-even in Y) and O3-only (sep-even in X),
  each 2 coeffs away; need explicit cert at the sub-family point first. [setup]
- dense-rational: F_t=F0+eps F1 from Fermat; prior deform.py ran 4 of 12 subsets I, all first-order
  INCONSISTENT (rank_aug=rank+1). Finish remaining subsets + interpret the single obstruction functional.
- off-paired: t=0 is the holdout chain; holdout has NO explicit witness, so no cert to continue. State
  obstruction via the E1' class-derivative [A C] with C = X0^2X1X2Y0Y2/3. [setup]
- one-reflection-skew: neighbour one-reflection has an explicit degree-60 cert; first-order deform system
  ~1.4M x 4.1M (DiamondJM). Too large for a full run; compute E1'-class-derivative obstruction instead.
