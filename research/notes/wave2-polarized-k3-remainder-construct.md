# Wave 2 — lane `polarized-k3-remainder`, lens CONSTRUCTION

STATUS: in progress — no witness (K3 sector provably out of reach of every quotient/cover mechanism on X; mixed sector construction being prepared).

Log: `/tmp/hodge/results/polarized-k3-remainder/wave2-construct.log.jsonl` (one JSON line per experiment).

## Results so far
1. [PROVEN, sympy, local] The linear automorphism group of F is larger than recorded: with U = X+Y, V = X−Y,
   h(U,V) = (U, ηV) (η primitive 6th root of unity) preserves F = ½(f(U)+f(V)) exactly, h^3 = the coordinate swap X ↔ Y,
   h^3 τ = ((X,Y) ↦ (−Y, X)). Hence Aut(X) ⊇ D_6 = ⟨h, τ⟩ of order 12 (and = D_6 given Aut(C_f) = 1 and uniqueness of the U/V splitting).
   The dossier's conditional claim "Aut = {1,τ}" is wrong (only the *diagonal-in-(X,Y)* group has order 2).
2. [PROVEN, sympy, local] X contains six planes Π_μ = {Y = μX}, (1+μ)^6 + (1−μ)^6 = 0, μ ∈ {±i, ±i(2±√3)}; equivalently {V = λU}, λ^6 = −1.
   They are the fixed planes of the involutions h^k τ, k odd.
