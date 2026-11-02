# cycle-q2 — CONSTRUCTION — part 5: the one unexcluded mechanism, explicit equations, and where the book's recipe fails

## M*: S3-transport with base absorption (the mechanism that reaches the character)
Everything else in this lane kills χ_A before the geometry starts (parts 2–4). The one mechanism that reaches the
target's character is DiamondJM's: transport the second variation ∂_{s1}∂_{s2}(Ω/F_s) (2AΩ/F³ = ∂_{s1}∂_{s2}Ω/F_s|_0,
F_s = F + s1 m1 + s2 m2, m1 = X0X1X2Y0Y1Y2, m2 = X0X1X2²Y0Y1) through the paired-sign cubic and its Jacobian, keeping the
factor v² = t. Explicit auxiliary equations (`w2c_s3.py`, book R14 normalisation):

  cubic U = a x³ + b y³ + c z³ + d x²y + e y²z + f z²x,  (a,b,c) = (1+t³, 1+c1⁶t³, 1+c2⁶t³),  (d,e,f) = (t/2)(c1, c1c2, c2),  t = v² = (Y0/X0)²,
  s4 = c1c2 t³ (c1⁶c2 t³ + c1²c2⁶ t³ + c1c2² t³ + c1² + c1c2² + c2)/8,   s6 = 27a²b²c² + 4a²ce³ + 4ab²f³ + 6abcdef + 4bc²d³ − d²e²f²  (a polynomial in t³ of degree 18),
  S3(c1,c2):  η² = ξ³ + s4 ξ − s6/4  over P¹_t,   RES(c1,c2): the same over P¹_T with T = t³.

Verified hypotheses: (H1) s4, s6 are polynomials in t³ (t-exponents {3,6} and {0,3,…,18}) — the RES over the T-line exists, as cycle-three uses; (H2) over the t-line: deg s4 = 6 ≤ 12, deg s6 = 18 ⇒ n = 3, χ(S3) = 3, p_g = 2, q = 0; the model is minimal at t = 0 (ord Δ = 0) and t = ∞ (ord b = 0), discriminant of degree 36, j non-constant ⇒ κ(S3) = 1; (H3) the 2-forms are ω₀ = dt∧dξ/η and ω₁ = t·dt∧dξ/η, with characters ω, ω² under t ↦ ωt; the fibrewise element g_fib = (0,1,2,2,3,4) ∈ G_fib∖K (part 2) has χ_A(g_fib) = ω² and acts by t ↦ ω²t, so the target's (3,1)-part sits fibrewise in the sector of the **untwisted** form ω₀; (H4) [DiamondJM, cited] the fibrewise (0,2)-component of the transported source is nonzero (P ∉ J37) but lies in J37 + ∂_{c1}G·S19 + ∂_{c2}G·S19 — no obstruction at the fourfold level; (H5) [benediamond, cited] the vertical remainder along R(T) quadratic is nonzero (5-adic certificate), so no fibrewise primitive exists.

## Where the book's incidence recipe fails — precisely
The §3.3/3.4 recipe (polynomial curve family φ(τ) on the auxiliary surface, square incidence with invertible Jacobian
over F_p, arithmetic lifting, polynomial normal, termwise integration, trace) requires a **dominant** family of
polynomial (rational) curves on the auxiliary surface — that is the content of "every polynomial source has a
primitive". A surface dominated by a family of rational curves is uniruled (κ = −∞). S3 has κ = 1 (H2), the K3
fibres of the (λ,μ)-fibration have κ = 0, and X/P itself has κ = 0: **no finite-field incidence seed can exist for any
auxiliary surface that carries the ω²-sector**, because the failure is at dominance, not at seed-finding or lifting.
Parametrising curves on the RES by t with T = t³ does not help: termwise integration of t^{n+1/3} is fine, but the
trace back to the base needs a generically finite dominant map, which κ(S3) = 1 forbids. Hence
`HodgeComputer.Constructions.PolynomialCurves` and `HoldoutCover` cannot be the Lean engine for this lane.

## Why the support cannot be vertical either (new, small)
For a surface fibration π: X ⇢ B with generic fibre S of irregularity q(S) = 0 and injective Kodaira–Spencer map along
every curve (DiamondJM for S3 over (c1,c2)), the resolution of a vertical divisor π⁻¹(C) has p_g = 0 (Leray: R¹π_* = 0,
H⁰(C, R²)^{2,0} = flat (2,0)-sections = 0), so its Gysin image is of type (2,2) only; since α^{3,1} ≠ 0 (fqlx), α is not
supported on any vertical divisor of the (c1,c2)-fibration. Removing curves Γ_c from the fibres does not change the
fibrewise (0,2)-component either. So the primitive must be non-fibrewise: a rational current Ψ = Ψ₁dc1 + Ψ₂dc2 with
∇Ψ = σ, i.e. the class of σ in H²_dR(B°, V_{ω²}) must vanish, V_{ω²} the rank-12 (1,10,1) local system
(= H¹(P¹_T, R¹_RES ⊗ L_{ω²}), the Kummer-twisted cohomology of the rational elliptic surface — a reformulation that
makes the connection computable from the RES's twisted Picard–Fuchs data instead of DiamondJM's degree->130 frame).

## What would make M* constructive (not done here; the actual remaining work)
1. Compute the Gauss–Manin connection of V_{ω²} over C(c1,c2) in the twisted-de-Rham frame {t^j T^k dξ/η} of the
   RES (rank 12), and the class of σ in H²_dR(B°, V_{ω²}); if zero, the primitive Ψ is an explicit rational current and
   fqlx's pole-three clearing adapter turns it into a certificate with denominator = polar locus of Ψ (horizontal).
2. If nonzero on B° = C²∖discriminant, allow poles along a horizontal divisor Z (curves Γ_c ⊂ S3 with a Noether–Lefschetz
   jump in the ω₀-sector): P5 of the dossier, now with the precise requirement that Z be horizontal and non-normal.
Both are two-variable D-module computations, outside the polynomial-curve library; they are the only open door.
