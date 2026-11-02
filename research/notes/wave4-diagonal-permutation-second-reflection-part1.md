# Wave 4 — diagonal-permutation — route second-reflection — Part 1: the exact argument

Notation. F is the instance sextic (data/diagonal-permutation/instance.rec), X = {F = 0} ⊂ P^5, s = (X0↔X1, Y0↔Y1),
r = X0 − X1, t = Y0 − Y1. The (01)-pencil is {r = μ t}; its member X_μ = X ∩ {r = μt} is a sextic threefold in P^4 with
coordinates (p, q, a, b, t): X0 = (p+μt)/2, X1 = (p−μt)/2, X2 = a, Y0 = (q+t)/2, Y1 = (q−t)/2, Y2 = b. F_μ is even in t, so s
acts on X_μ as the reflection t ↦ −t. Question (wave-3 §3.4, M-a): does X_μ (for generic μ, possibly after base change)
carry a SECOND reflection, i.e. a linear involution of P^4 with eigenvalues (1^4, −1) preserving F_μ, other than s?
All automorphisms of a smooth sextic threefold in P^4 are linear (Pic = Z·O(1) by Lefschetz, K = O(1)), and an involution
fixes a divisor iff it is a reflection of P^4; so the question is exactly: is Aut_lin(X_μ) ∋ a reflection ≠ s?

## Lemma A (Hessian-rank reduction; exact, any field k of characteristic ≠ 2,3,5)
Let G be a form of degree 6 in n variables over k̄ and g ∈ GL_n(k̄) with G∘g = G. Then
(i) rank Hess G(gv) = rank Hess G(v) for all v (Hess(G∘g)(v) = gᵀ Hess G(gv) g);
(ii) the polar hyperplane H_v = ker ∇G(v) satisfies H_{gv} = g(H_v) (∇G(gv) = g^{−T} ∇G(v)).
Hence g permutes the finite set R_1(G) = {v : rank Hess G(v) = 1} (when finite) and permutes their polar hyperplanes.

Application to F_μ (n = 5). If R_1(F_μ) = {e_a, e_b} (verified computationally, W4-1), then every g ∈ Aut_lin(F_μ)
permutes {e_a, e_b}, hence preserves L = span(e_a, e_b) and, since ∇F_μ(e_a) = 6·e_a^* (F has no monomial X2^5·(other)),
∇F_μ(e_b) = 6·e_b^*, also preserves M = H_{e_a} ∩ H_{e_b} = {a = b = 0}. As k̄^5 = L ⊕ M, g = (g_L, g_M) is block
diagonal. Restricting F_μ∘g = F_μ to L (where F_μ|_L = a^6 + b^6) gives g_L ∈ Aut(a^6 + b^6) = {monomial 2×2 matrices with
entries in μ_6} (the zero set of a^6 + b^6 is a regular hexagon in P^1 whose stabiliser in PGL_2 is dihedral of order 12).
Projectively there are 12 cases for g_L: a ↦ a, b ↦ ζ^j b and a ↦ b, b ↦ ζ^j a (j ∈ Z/6, ζ a primitive 6th root of unity).
For each case, g_M ∈ GL_3 is determined by the polynomial system F_μ∘(g_L, g_M) = F_μ (9 unknowns, ~140 equations of
degree ≤ 6). The solution scheme is reduced when H^0(T_{X_μ}) = 0 (its tangent space at a solution is {h ∈ gl_5 : ξ_h F_μ = 0}
= 0), so the degree of the zero-dimensional ideal equals the number of k̄-points.

Application to F (n = 6), Lemma A'. If R_1(F) = {e_0, …, e_5} (W4-2), every g ∈ Aut_lin(F) permutes the coordinate points,
i.e. is monomial: z_i ↦ λ_i z_{π(i)}. F∘g = F forces λ_i^6 = 1 (Fermat part) and c_m λ^m = c_{π(m)} for every monomial m of F;
since all c_m are rational and λ^m ∈ μ_6, the ratio c_{π(m)}/c_m must be ±1 (over F_p-bar: a 6th root of unity mod p — we
checked that no ratio of the distinct coefficient magnitudes {1, 5/7, 5/11, 5/13} is a 6th root of unity mod p).
The enumeration is a finite Z/6-linear problem (W4-6).

## Lemma B (good reduction and injectivity of specialisation for linear automorphisms; exact)
Let R be a complete DVR with fraction field K, uniformiser ϖ, residue field k, and let G ∈ R[z_1..z_n] be a form of degree
d ≥ 3, primitive, with {G = 0} smooth over k (hence over K). Then
(B1) every g ∈ PGL_n(K̄) with G∘g ∝ G lies in PGL_n(O_{K̄}) (i.e. extends to the smooth model), and
(B2) if {M ∈ gl_n(k̄) : ξ_M Ḡ ∈ k̄·Ḡ} = k̄·1 (equivalently H^0(X_k, T) = 0, via 0 → T_X → T_P|_X → O_X(d) → 0 and
H^0(T_P(−d)) = H^1(T_P(−d)) = 0), the reduction map Aut_lin(X_{K̄}) → Aut_lin(X_{k̄}) is an injective group homomorphism,
compatible with characteristic polynomials (so reflections map to reflections) and with any automorphism defined over R
(so s ↦ s).
Proof of (B1). Replace K by a finite extension L over which g is defined; O_L is a DVR. Scale g ∈ M_n(O_L) primitive and write
its Smith normal form g = u·diag(ϖ^{e_1}, …, ϖ^{e_n})·u' with u, u' ∈ GL_n(O_L), 0 = e_1 ≤ … ≤ e_n. Put G_1 = G∘u, G_2 = G∘u'^{−1};
both are primitive with smooth reduction. G∘g = cG becomes G_1(ϖ^{e_1}w_1, …, ϖ^{e_n}w_n) = c·G_2(w), i.e.
v(f^{(1)}_α) + e·α = m + v(f^{(2)}_α) for every exponent α, with m = v(c) ≥ 0. Smoothness of Ḡ_1 at the coordinate point e_1
forces a unit coefficient of G_1 on a monomial α_1 ∈ {d·ε_1, (d−1)ε_1 + ε_j}, whence m ≤ e·α_1 ≤ e_j ≤ e_n; smoothness of
Ḡ_2 at e_n forces a unit coefficient of G_2 on β_n ∈ {d·ε_n, (d−1)ε_n + ε_j}, whence (d−1)e_n ≤ e·β_n ≤ m. So (d−1)e_n ≤ e_n,
and d ≥ 3 gives e_n = 0: g ∈ GL_n(O_L).
Proof of (B2). If g ∈ GL_n(O_L) reduces to a scalar, write g = λ(1 + ϖ^r h) with h primitive, r ≥ 1, and (increasing r if h̄ is
scalar; the process stops because the group is finite and g is not scalar) h̄ non-scalar. Then G∘g = cG gives
(c−1)G ≡ ϖ^r ξ_h G (mod ϖ^{2r}), so ϖ^r | (c−1) and, writing c − 1 = ϖ^r c', c̄' Ḡ = ξ_{h̄} Ḡ with h̄ ∉ k̄·1: contradiction with
the hypothesis. Hence the kernel of reduction is trivial; compatibility statements are clear (p ≠ 2 for eigenvalues ±1).

## Consequences (all hypotheses verified computationally; see the results table)
(C1) R = Z_(p)-completion, p = 1073741827, G = F_{3/7} (n = 5): X_{3/7} smooth mod p (W4-3), H^0(T) = 0 mod p (W4-4),
Aut_lin(F̄_{3/7}) = {1, s} over F_p-bar (W4-1 + Lemma A + W4-5). By Lemma B: Aut_lin(X_{3/7}) over Q-bar = {1, s}.
(C2) R = Q̄[[μ − 3/7]], G = F_μ (coefficients in Q(μ), integral at μ = 3/7): special fibre X_{3/7} smooth over Q̄ (from W4-3,
a singular Q̄-point would reduce to a singular F̄_p-point), H^0(T_{X_{3/7}}) = 0 (W4-4 over Q). By Lemma B, Aut_lin of the
geometric generic fibre X_{η̄} of the (01)-pencil injects into Aut_lin(X_{3/7}) = {1, s} and contains s: Aut_lin(X_{η̄}) = {1, s}.
Every automorphism of X_{η̄} that could appear "after base change" is an automorphism of X_{η̄}; so there is NO hidden or
created second reflection on the (01)-pencil fibres — in fact no non-trivial fibre automorphism besides s at all.
(C3) R = Z_(p)-completion, G = F (n = 6): X smooth mod p (W4-3), H^0(T_X) = 0 (W4-4), R_1(F̄) = coordinate points (W4-2),
monomial automorphisms = S3 × ⟨τ⟩ (W4-6). By Lemma A' and Lemma B: Aut_lin(X) over Q-bar = S3 × ⟨τ⟩ (order 12), with no
reflection; this settles the wave-2 item J2 (hidden reflections of X: none) exactly, replacing the timed-out 27-unknown msolve
systems. A §2.7 partner coordinate (F = c y^6 + y G_5 + G_6) requires its pole v to have rank Hess F(v) = 1 (wave-2 §2.1), hence
v ∈ {e_i}; but every variable occurs with exponent 2 or 3 in some monomial of F, so no coordinate point is a §2.7 pole: no §2.7
partner coordinate (exact over F_p-bar; over Q-bar modulo the transfer of R_1, which Lemma B does not cover — rank-1 points of F
over Q-bar reduce to rank-≤1 points mod p, i.e. to coordinate points, so a non-coordinate Q-bar pole would have to be
p-adically close to a coordinate point; the exact Gröbner computation over Q is listed as the remaining check).

## What was NOT computed / caveats
- msolve (v0.6.5) correctness over F_p is assumed throughout, as elsewhere in the project.
- Nonlinear (birational) automorphisms of X_μ are excluded by Lefschetz only for SMOOTH fibres; singular members of the pencil are not
  the object of the §4.6 descent, whose cover is built from the generic fibre.
- Reflections on OTHER pencils/fibrations of X are not addressed by this route (the (01)-pencil is the one with a given first
  reflection); the full statement Aut_lin(X) = S3 × ⟨τ⟩ (C3) shows no second reflection can come from X itself on any pencil.
