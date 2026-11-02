# Attack plan: cyclic-triple-residual (q=1, A = X0^3 Y0^3)

F = Fermat sextic + (5/7)(X0^3X1Y0Y1 + X0^3X2^3 + X0^2Y0Y1Y2^2 + X0X1Y1^3Y2 + X0X2^3Y1Y2 + X1X2^3Y0Y1).
Z/3: g = (X2 -> zeta3 X2) fixes F.  Residue class of A Omega/F^2 has g-character zeta3 (A invariant, det g = zeta3).
The solved companion (A = X0 X2^2 Y0^3) lives in the g-INVARIANT sector = pullback from the Fano quotient
W = {Z^2 + (5/7) C Z + H = 0} in P(1,1,3,1,1,1) (double cover of P^4 branched on a sextic). Our class is in the
zeta3-eigenspace, which no uniruled intermediate quotient of the (Z/3)^2 x| Z/2 Galois closure serves (checked).
A certificate therefore needs a divisor T = X ∩ V(D) whose resolution carries (2,0)-forms in the zeta3-eigenspace
hitting this class; T must have non-normal or non-rational singularities (book §3.1).

Mechanisms to test (all exact, mod-p linear algebra via fixed_d):
 0. sanity: A not in J(F); character table; D=1 fails.           [done]
 1. Fermat calibration: on F0 = sum Z^6 the target X0^3Y0^3 is served by the cone hyperplanes X0 = lambda Y0 (lambda^6=-1):
    D = X0^6 + Y0^6.  Check which cone/sum denominators serve it on Fermat, to know what to deform.
 2. cones: D in k[Z_S], |S| = 3 (vertex P^2, base cubic/quartic curve of genus >= 1) and |S| = 4 (vertex line, base
    quartic/sextic surface with p_g > 0); definite Y-parity; coefficient scans.
 3. Fermat-type sums c_i Z_i^6 over subsets (incl. fqlx's excluded X2^6+Y0^6+Y1^6+Y2^6), pure-power part, F minus pieces.
 4. structure of F as a quadratic in Z = X2^3: U = X2^3 + (5/14) C, C, H, H' = H - t^2 C^2/4, discriminant; partials
    F_i and products; block radial derivatives; Hessian-related.
 5. E1-kernel sweep: dim ker([A]: R_d -> R_{d+6}) for d = 1..8 (R = Jacobian ring); candidates from any kernel.
 6. perturbation scans D0 + t*M around near-consistent D0 (rank-drop detection).
Each test is logged to log.jsonl.  Negative results with prime/rank/nullity are recorded.
