STATUS: RUNNING — 2026-11-02 04:40 UTC — setting up the Fano variety of lines of W = X/mu3 in P(1^5,3); first enumerative counts
# Wave 4 — lane cyclic-triple-residual — route fano-lines

Work dir: /tmp/hodge/results/cyclic-triple-residual/wave4/ ; log: wave4/log.jsonl
Aside (user question "what is diamondjm"): `DiamondJM` is a GitHub contributor to the hodge.computer challenge (accepted solver of
`one-reflection-y2-sector`, author of the external hodgeHoldout sweep of ~1,300 degree-6 supports; see /tmp/hodge/notes/external-diamondjm.md).

## Setup (exact)
F = G(X0,X1,X2^3,Y0,Y1,Y2), G = Z^2 + L Z + M with (weighted P(1,1,1,1,1,3), Z = X2^3)
  L = (5/7)(X0^3 + X0 Y1 Y2 + X1 Y0 Y1),  M = X0^6+X1^6+Y0^6+Y1^6+Y2^6 + (5/7)(X0^3X1Y0Y1 + X0^2Y0Y1Y2^2 + X0X1Y1^3Y2).
W -> P^4 is the double cover branched along B' = {L^2 - 4M = 0}; the mu3-branch divisor B = W ∩ {Z=0} ≅ {M=0}.
A weighted line of W over a line l ⊂ P^4 is a binary cubic c with c^2 + L|_l c + M|_l ≡ 0, i.e. (2c+L|_l)^2 = (L^2-4M)|_l:
  lines of W <-> tritangent lines of the sextic 3-fold B' (two lines per l, swapped by the covering involution iota).
E_L = pi^{-1}(L) = {x2^3 = c(s,t)} is a plane cubic (j=0) in the plane <l, e_X2> ⊂ P^5; X ∩ <l,e_X2> = E_L ∪ E_{L'}.
Class bookkeeping in G(2,5): [𝓛] = a σ3 + b σ21; a = #lines through a general point, b = σ21-count; Plücker degree a+2b;
#lines meeting a general line and a general plane = a+b = degree of the X2-free sweep divisor of S_Σ.

## Results table
| step | what | prime | result | file |
|---|---|---|---|---|
