---
name: hodge-sextics-math-insights
description: "Mathematical state of the open hodge.computer lanes as understood on 2026-11-01 (what is excluded, what mechanisms could work)"
metadata:
  node_type: memory
  type: project
  originSessionId: 2f213184-0ffd-40a3-8b2f-58f1332191ff
  modified: 2026-11-01T22:29:43.761Z
---

Facts established while attacking the 16 open lanes (2026-11-01):

- A certificate for (F,q,A) with fixed D is finite linear algebra in (N,Q). Engine: ~10-40 s per degree-6 denominator (compressed modular rref), ~400 s for the exact rational solve on dense cases. The (3,1)-support constraint (book ch.3 §3.1): the threefold X∩V(D) must be non-normal or have non-rational singularities; cones with vertex a line over surfaces of degree ≥4 (p_g>0) or vertex a plane over curves of genus ≥1 are the natural low-degree candidates; unions of hyperplane sections through a P^3 cannot support (3,1) classes.
- For generic (dense) sextics, E1 (A·D ∈ J(F)) forces D ∈ J(F) up to degree 9 (fqlx proved deg D ≥ 7 for dense-rational; our dossier agent extended the Jacobian-type classification to d=9, kernel dims 294/261 at d=10).
- Prior sweeps: DiamondJM swept ~1300 "natural degree-6 supports" over all q=1 open lanes (nothing). fqlx/benediamond/DiamondJM did deep geometric work on cycle-q2, one-reflection-skew, off-paired, dense-rational (see /tmp/hodge/discussions and dossiers).
- Character bookkeeping: residue character of monomial Z^e under diagonal weights w is Σ(e_j+1)w_j; quotient constructions serve only their averaging-ideal characters; the open "residual" lanes are exactly complementary characters.
- four-plus-two-residual: Thom–Sebastiani puts the target X0^3Y0^3 in the third-root (μ3) sector; the Fano companion W_H={H(X0,X1,Y0,Y1)+s^3+t^3} ⊂ P(1,1,1,1,2,2) has −K=O(2), but the binary block's triple-cover curve u^3=B(X2,Y2) (genus 4) has Prym of CM signature (2,1), so the elliptic contraction needs the specific B to be special (check J(T_B) splitting by point counts).
- diagonal-permutation: target is in the standard rep of S3; served on the sub-family where the O2,O3 coupling orbits vanish (paired-sign sector theorem); the hard part is continuation in the two coefficients 5/11, 5/13.

**Why:** avoids re-deriving the frontier analysis each session.
**How to apply:** consult before proposing a construction; update when a lane changes status. See [[hodge-sextics-project-setup]].
