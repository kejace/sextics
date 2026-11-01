# Dense rational · full coefficient support  (https://github.com/hodge-computer/sextics/discussions/59)


## COMMENT DC_kwDOVLOOJc4BIxgI by fqlx at 2026-10-20T16:38:14Z (updated 2026-10-20T16:38:14Z)

## Dense Rational search update — 20 October 2026

**No certificate has been found, and no Lean `HasCertificate` proof is claimed.** This AI-assisted report summarizes the completed search and exact obstruction calculations. The work includes separate implementations of the main arithmetic checks.

The input is the fixed 236-term Dense Rational polynomial, with $q=1$ and $A=X_0^3Y_0^3$, pinned to commit `d5881c13a7e785b6328cd823133f5245937b93e9`.

### 1. Every legal denominator has degree at least seven

The packet proves $\deg D\ge7$ for the full certificate contract over complex coefficients, allowing an arbitrary quotient $Q$. This is an unrestricted denominator-degree lower bound; it does not assume a particular denominator family or parity.

### 2. A nonlinear homogeneous polynomial Piola construction is excluded

Suppose $\Phi:\mathbb{C}^6\to\mathbb{C}^6$ is dominant, its components are homogeneous of positive degrees, and $P,b$ are polynomials satisfying

$$
F=P(\Phi),\qquad A=\det(d\Phi)\cdot b(\Phi).
$$

Then $\Phi$ must be invertible linear.

The arithmetic input uses the six **ambient partial derivatives taken before restriction** to $X_0=0$ or $Y_0=0$. There are no restricted-gradient syzygies with multiplier degree at most four, and the degree-five kernels consist precisely of the 15 constant Koszul vectors.

| Hyperplane | Multiplier degree | Matrix | Exact rank | Witness determinant mod 101 |
|---|---:|---:|---:|---:|
| $X_0=0$ | 4 | $715\times420$ | 420 | 32 |
| $X_0=0$ | 5 | $1001\times756$ | 741 | 48 |
| $Y_0=0$ | 4 | $715\times420$ | 420 | 94 |
| $Y_0=0$ | 5 | $1001\times756$ | 741 | 1 |

Exact kernels provide the characteristic-zero upper bounds; nonzero modular minors give matching lower bounds. An independent implementation rebuilt these matrices directly from the rational record.

The proof combines $J_\Phi\mid X_0^3Y_0^3$, successive divisibility along the two coordinate hyperplanes, adjugate identities, and the chain rule. Finiteness is unnecessary.

**Scope:** this excludes one polynomial Piola construction with the exact density identity. It does not exclude sums of densities, rational $b$, rational or inhomogeneous maps, or general certificates.

### 3. Degree seven reduces to a concrete quadratic system

Writing $G=9973F$, $A_0=9973A$, and $J=(\partial_iG)$, every Euler-normalized degree-seven solution of the first certificate equation has

$$
D=\nabla G\cdot w,\qquad N=A_0w+K\nabla G,
$$

with $w$ quadratic and $K$ skew-symmetric with cubic entries. The degree-eight calculation proves the injectivity needed to reduce a necessary condition to

$$
w^T\mathrm{Hess}(G)w\in J_8.
$$

After removing the six Euler gauge directions, this gives **951 homogeneous quadratic equations in 120 complex variables**.

Thirty-one numerical starts produced no candidate. Two bounded finite-field elimination attempts did not finish a Gröbner basis. These failures do not prove degree-seven nonexistence. A candidate $w$ would still require solving the full second equation for $K$ and checking denominator legality.

### 4. Quotient-pole searches have explicit limits

In coordinates $x_i=X_i$, $s=Y_0^2$, $a=Y_1/Y_0$, $b=Y_2/Y_0$, write $H=P_0+sP_1+s^2P_2+s^3P_3$.

The full weighted-polynomial numerator module with a simple $H_s$ pole is excluded even with arbitrary rational coefficients in $\mathbb{C}(a,b)$, all six vector components, and arbitrary $Q$. Extra fiber poles are outside this statement.

Further exact exclusions cover specified bounded coefficient spaces for $H_s^2$, the complementary discriminant factor, and the physical radial square $(sH_s)^2$. The packet records the precise bounds and enforces compatibility of coefficient values and derivatives.

### 5. The geometric construction still lacks a map to the target

A concrete auxiliary K3 surface from a genus-four trigonal slice falls within [Bolognesi–Laterveer’s algebraic Kuga–Satake results](https://arxiv.org/html/2310.11981v2). This does not supply a correspondence with verified action on the Dense target class. That connection remains missing.

### Reproduction and next step

The attached `Dense-Rational-search-packet.zip` contains the proofs, exact matrices, minor indices, scripts, search bounds, and replay logs. Start with `CONTINUATION.md`.

Core and full continuation replays passed; the Piola premises also passed a separate independent reconstruction. No Lean build was performed.

The immediate algebraic task is to find a non-Euler $w$ satisfying the degree-seven Hessian condition, then solve for $K$. No claim of unrestricted nonexistence is made.

[Dense-Rational-search-packet.zip](https://github.com/user-attachments/files/33728423/Dense-Rational-search-packet.zip)

