# Eight planes · contact arrangement  (https://github.com/hodge-computer/sextics/discussions/176)


## COMMENT DC_kwDOVLOOJc4BJDZL by benediamond at 2026-10-23T17:01:29Z (updated 2026-10-23T17:04:13Z)

### Rigidity of the eight planes, and the Hodge types of two natural surface fibrations

Here are three checks on the geometry suggested in the [instance rationale](https://github.com/hodge-computer/sextics/blob/main/data/eight-plane-configuration/why.md). They explain why the most immediate projections do not supply a surface fibration with a wholly algebraic selected $H^2$. They do **not** show that the particular target is unsupported. This investigation and the preparation of this comment were assisted by ChatGPT.

Write $X=\lbrace F=0\rbrace $, $\tau(X_i,Y_i)=(X_i,-Y_i)$, and

$$
Q_i=X_i^2+Y_i^2,\qquad F=\sum_{i=0}^2Q_iK_i,\qquad K_i=X_i^4-X_i^2Y_i^2+Y_i^4+\frac5{p_i}H_i,
$$

with $(p_0,p_1,p_2)=(7,11,13)$ and the instance's quartics $H_i$.

**1. Each of the eight planes is infinitesimally rigid.** For $\Pi:\ Y_i=\varepsilon_i\mathrm iX_i$, smoothness gives

$$
0\longrightarrow N_{\Pi/X}\longrightarrow\mathcal O_{\mathbb P^2}(1)^{\oplus3}\longrightarrow\mathcal O_{\mathbb P^2}(6)\longrightarrow0.
$$

The map on global sections is injective. Indeed all coefficients are integral at $5$, and $F$ reduces to the Fermat sextic. Choose either root of $\mathrm i^2=-1$ modulo $5$. A first-order normal displacement $Y_i=\varepsilon_i\mathrm iX_i+\eta L_i(X)$ gives derivative

$$
\sum_i6(\varepsilon_i\mathrm i)^5X_i^5L_i(X)\pmod5.
$$

Its nine monomials $X_i^5X_j$ are distinct and their scalar coefficients are nonzero. Thus the global-section map has rank $9$ already modulo $5$, hence in characteristic zero. In particular

$$
h^0(N_{\Pi/X})=0,\qquad h^1(N_{\Pi/X})=28-9=19.
$$

So moving these planes themselves does not produce a supporting divisor. This leaves open other Hilbert-scheme components and families of curves or surfaces related to the arrangement.

**2. Projection from one plane gives quintic surfaces.** If $u_0,u_1,u_2$ generate the ideal of $\Pi$, write $F=\sum u_iG_i$ with $G_i$ quintic. In the containing $\mathbb P^3$ over $[a_0:a_1:a_2]$, set $u_i=\lambda a_i$; after removing the plane, the residual equation is

$$
\sum_i a_iG_i(\lambda a,v)=0.
$$

The generic residual surface is smooth. Along $\lambda=0$, the restrictions $G_i|_\Pi$ have no common zero, since such a zero would be singular on $X$; Bertini applies to their linear system. Away from $\Pi$, use generic smoothness of the projection. Adjunction then gives

$$
K_S=\mathcal O_S(1),\qquad p_g(S)=4,\qquad q(S)=0.
$$

The involution exchanges $\Pi$ with the opposite plane, so it exchanges the two projection maps rather than acting fiberwise on this one. The immediate whole-fiber $p_g=q=0$ criterion therefore does not apply, and there is no immediate fiberwise $\tau$-projector repair.

**3. The invariant net $[Q_0:Q_1:Q_2]$ has nine invariant holomorphic two-forms on its generic fiber.** Its base locus is the eight-plane arrangement. Over a general $r=[r_0:r_1:r_2]$, all $r_i\ne0$, the residual fiber is

$$
S_r=\left\lbrace r_1Q_0-r_0Q_1=0,\quad r_2Q_0-r_0Q_2=0,\quad \sum_i r_iK_i=0\right\rbrace \subset\mathbb P^5.
$$

This is a smooth complete intersection of type $(2,2,4)$ for general $r$. Here is a check of smoothness for this actual family. The intersection of the first two quadrics has exactly six singular points: one coordinate pair $(X_i,Y_i)$ is nonzero and isotropic, and all other pairs vanish. At each such point, smoothness of $X$ implies $K_i\ne0$, so the general quartic $\sum r_iK_i$ avoids it. On the remaining plane and line strata of the base arrangement, the two quadric differentials are independent. The $K_i$ have no common zero there, again by smoothness of $X$. Bertini on the finitely many strata makes the third differential independent. Away from the base arrangement, generic smoothness of the rational map on $X$ applies.

Adjunction gives $K_{S_r}=\mathcal O_{S_r}(2)$, so $p_g=21-2=19$ and $q=0$. To get the relevant character, one must include the volume form: all three defining equations are $\tau$-invariant, but the determinant of $\tau$ on the six homogeneous coordinates is $-1$. Consequently the holomorphic two-form with quadratic residue numerator $P$ has the **opposite** parity to $P$.

There are nine odd quadratic monomials $X_iY_j$, and neither quadratic relation is odd. The even quadratic space has dimension $12$, with two relations. Hence

$$
\dim H^{2,0}(S_r)^+=9,\qquad \dim H^{2,0}(S_r)^-=10.
$$

The target residue with $A=X_0^3Y_0^3$ is $\tau$-invariant: its numerator and projective volume are both odd. Thus the entire invariant relative $H^2$ of this natural net is not of type $(1,1)$. A criterion requiring a wholly algebraic invariant fiber summand cannot directly establish support here.

The remaining question is more specific: does the fixed absolute target avoid, cancel in, or admit a further algebraic projector away from this nine-dimensional relative holomorphic piece? Nothing in these calculations decides that. They exclude two immediate whole-fiber reductions and infinitesimal movement of the eight displayed planes, rather than excluding an eventual correspondence built from the arrangement.

