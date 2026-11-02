# cycle-q2 — CONSTRUCTION — part 3: covers creating a reflection (lens i)

Script `w2c_cover.py` (local, instant; `w2c_cover.json`).

## The §2.7 shape
A holdout-type cover needs, after a linear change of the Y-coordinates, F = f(rest) + y⁶ + ℓ1⁶ + ℓ2⁶ + y·(ℓ1 a1 + ℓ2 a2) with a1, a2 in the remaining variables; then T⁶ = ℓ1⁶ + ℓ2⁶, α = ℓ1/T, β = ℓ2/T, S = y² + T², P = yT, and y ↔ T is the created reflection, with the Fermat base α⁶ + β⁶ = 1 a *constant* coefficient curve.

## Non-coordinate linear forms are excluded structurally
The Y-sextic part of F is the Fermat form Y0⁶ + Y1⁶ + Y2⁶. If y, ℓ1, ℓ2 are linear forms with y⁶ + ℓ1⁶ + ℓ2⁶ = Y0⁶ + Y1⁶ + Y2⁶ then (y,ℓ1,ℓ2) is a linear automorphism of the Fermat sextic curve, and for degree ≥ 4 every such automorphism is monomial ((μ6)² ⋊ S3; Tzermias / Shioda). So the only candidates are (ζ-scaled) coordinate triples; the same holds for the X-sextic part, and a mixed choice (y an X-coordinate, ℓ's Y-coordinates) is covered by the enumeration below. A cover whose base depends on an active coordinate (e.g. T⁶ = Y0⁶ + Y2⁶ + ½X0X2³Y0Y2) is exactly the failure mode the book excludes in §2.7 ("a perturbation that makes the base depend on an active coordinate is a different failure that a new conic cannot repair").

## Enumeration (all 60 choices of y and {ℓ1,ℓ2} among the six coordinates)
Admissible: **0 of 60**. Minimal violators (one offending monomial each):
- y = Y0, ℓ = {X1,Y2} or {Y1,Y2}: X1³X2Y1Y2
- y = Y1, ℓ = {X2,Y0} or {Y0,Y2}: X0X2³Y0Y2
- y = Y2, ℓ = {X0,Y1} or {Y0,Y1}: X0³X1Y0Y1
Every other choice violates with ≥ 2 monomials. In each minimal case the obstruction is the single coupling monomial that closes the cycle — the "closing term" of why.md and of the dossier §1.7(3) — now checked for every coordinate pairing, not only the (Y0,Y2 | Y1) one.

## Consequence
Mechanism (i) is unavailable on this F, independently of the target: there is no finite cover of §2.7 type (coordinate or not) creating an exchange, and the §2.8 complementary quotient needs an order-6 diagonal element fixing X with Y-weights (5,1,5), which F does not admit (the closing monomial has Y-weight 10 ≡ 4; dossier §1.8). The cyclic structure has no endpoint, and every cover that symmetrises two couplings is broken by the third.
