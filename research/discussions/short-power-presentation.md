# Eight powers · linear-section remainder  (https://github.com/hodge-computer/sextics/discussions/174)


## COMMENT DC_kwDOVLOOJc4BJDb9 by benediamond at 2026-10-23T17:10:38Z (updated 2026-10-23T17:10:38Z)

### A product-of-curves model, with the actual quotient density

The eight-power presentation gives a compact birational model in which the remaining geometry is a coupled pair of plane-sextic pencils. This does not solve the instance, but it makes the effect of the global sign quotient and the transported source explicit.

Use the [fixed equation and numerator](https://github.com/hodge-computer/sextics/blob/faa9df2ec6406e4eacc135f677cf3711315d0e41/data/short-power-presentation/instance.rec)

$$
F=\sum_{i=0}^2X_i^6+\sum_{i=0}^2Y_i^6+c\big((L+M)^6+(L-M)^6\big),\qquad
c=\frac5{7^7},\qquad A=X_0^3Y_0^3,
$$

where $L=X_0+2X_1+3X_2$ and $M=Y_0+3Y_1+5Y_2$. On the chart $L=1$, write $m=M$ and, for $m\ne0$, put $y_i=Y_i/m$. Then

$$
X_0+2X_1+3X_2=1,\qquad y_0+3y_1+5y_2=1.
$$

Set $e=\sum X_i^6$ and $s=\sum y_i^6$. The hypersurface equation is exactly

$$
e+m^6s+2c(1+15m^2+15m^4+m^6)=0.
$$

Thus the geometric fiber over $(e,m)$ is an open part of $C_e\times D_s$, where the compact curves are

$$
C_e:\space \sum X_i^6=e(X_0+2X_1+3X_2)^6\subset\mathbf P^2,
$$

$$
D_s:\space \sum y_i^6=s(y_0+3y_1+5y_2)^6\subset\mathbf P^2,
\qquad s=\frac{-e-2c(1+15m^2+15m^4+m^6)}{m^6}.
$$

Both pencils have a smooth Fermat member at parameter zero, so their generic members are smooth plane sextics of genus $10$. The generic compact surface fiber therefore has $p_g=100$ and $q=20$; its middle cohomology retains $H^1(C_e)\otimes H^1(D_s)$.

The global sign involution sends $m$ to $-m$ and fixes the normalized $X_i,y_i$. Its quotient therefore has $n=m^2$ as a base coordinate and the **same** generic geometric surface fiber. The map from this base to the two pencil parameters is the cubic relation

$$
(s+2c)n^3+30cn^2+30cn+(e+2c)=0.
$$

In particular, the sign quotient does not project away the tensor-product cohomology. It removes the sign of the scale parameter, while the relation between the two varying curve parameters remains.

**The source can also be written without a hidden Jacobian.** Take $(x_0,x_1,y_0,y_1,m)$ as coordinates, with

$$
X=(x_0,x_1,(1-x_0-2x_1)/3),\qquad
Y=m\thinspace (y_0,y_1,(1-y_0-3y_1)/5).
$$

For the standard projective volume $\Omega=\sum_{j=0}^5(-1)^j Z_j\thinspace dZ_0\wedge\cdots\widehat{dZ_j}\cdots\wedge dZ_5$ in the order $(X_0,X_1,X_2,Y_0,Y_1,Y_2)$, restriction to $L=1$ gives

$$
\Omega=\frac{m^2}{15}\thinspace dx_0\wedge dx_1\wedge dy_0\wedge dy_1\wedge dm.
$$

Since $A=m^3x_0^3y_0^3$ and $dn=2m\thinspace dm$, the actual invariant form descends to

$$
\frac{A\Omega}{F^2}
=\frac{n^2x_0^3y_0^3}{30}\thinspace 
\frac{dx_0\wedge dx_1\wedge dy_0\wedge dy_1\wedge dn}{F_q^2},
$$

where

$$
F_q=e(x_0,x_1)+n^3s(y_0,y_1)+2c(1+15n+15n^2+n^3).
$$

The equality means that pulling the right-hand form back along $n=m^2$ gives the original form on this chart. The constants follow from the two determinants $1/3$ and $m^2/5$; the source carries $n^2$, not a constant coefficient on the parameter cover.

This gives a small explicit input for a correspondence or Gauss–Manin approach to the coupled tensor-product system. It does not show that the specified absolute residue has a nonzero component in every possible fiber summand, and the positive $p_g$ of this fibration does not exclude other fibrations or selected correspondences. What remains is an actual support or rational primitive for the displayed source, with its scale density and the cubic relation retained.

The original polynomial identity, parameter equations, involution action, and volume/source determinants were verified exactly. No solve, homogeneous certificate, or Lean proof is claimed. ChatGPT-assisted investigation for `benediamond`, with separate AI-agent checks; no independent human review is claimed.

