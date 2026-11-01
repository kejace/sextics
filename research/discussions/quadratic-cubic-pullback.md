# Quadratic pullback · residual density  (https://github.com/hodge-computer/sextics/discussions/175)


## COMMENT DC_kwDOVLOOJc4BJDaw by benediamond at 2026-10-23T17:06:33Z (updated 2026-10-23T17:06:33Z)

### The quadratic cover's invariant eigensheaves, and a conditional obstruction for pulled-back line cylinders

There is an exact calculation of the finite cover which narrows one natural use of the cubic's lines. It also measures how much of the invariant $H^{3,1}$ is missed by the known determinant density. This investigation and the preparation of this comment were assisted by ChatGPT. No solve is claimed.

Let $\phi:\mathbb P^5_Z\to\mathbb P^5_U$ be the map defined by the six quadrics of this instance, let $C:\sum U_i^3=0$, and write $X=\phi^{-1}(C)$ and $f=\phi|_X$. The quadrics are invariant under $\tau:(X_i,Y_i)\mapsto(X_i,-Y_i)$. Modulo $5$ they become the coordinate squares. Properness therefore excludes a characteristic-zero base point, and $\phi^*\mathcal O(1)=\mathcal O(2)$ excludes contracted curves. Thus $\phi$ is finite flat of degree $32$, as is its base change $f:X\to C$.

**The residual invariant space has dimension $207$.** The algebraic projector $p_{\mathrm{down}}=\frac1{32}f^*f_*$ has one-dimensional image on $H^{3,1}(X)$, represented by the determinant density $\delta$. As already checked in the [instance rationale](https://github.com/hodge-computer/sextics/blob/main/data/quadratic-cubic-pullback/why.md), the fixed $A=X_0^3Y_0^3$ is independent of $\delta$ in the sextic Jacobian ring: modulo $5$ the latter is $4X_0X_1X_2Y_0Y_1Y_2$, a distinct nonzero Fermat basis monomial.

The invariant residue space corresponds to odd $Y$-degree numerators, because $\det\tau=-1$. In degree $6$ there are

$$
3\binom72+\binom52^2+\binom72\cdot3=226
$$

such monomials. Exactly $18$ degree-six Jacobian relations $Z_i\partial_jF$ have this parity, with one index in each block. They are independent: the six quintic partials form a regular sequence and have no syzygy in degree $6$. Hence $h^{3,1}(X)^+=226-18=208$. Removing the cubic's line leaves a $207$-dimensional complement, and $[A]$ has a nonzero projection to it. This count alone gives no obstruction to divisor support.

**The useful finite-cover splitting is**

$$
(\phi_*\mathcal O_{\mathbb P^5})^+\cong\mathcal O\oplus6\mathcal O(-1)\oplus9\mathcal O(-2).\tag{1}
$$

Here is an algebraic derivation, valid for the actual quadrics. Put $B=\mathbb C[Z_0,\ldots,Z_5]/(Q_0,\ldots,Q_5)$. The $Q_i$ are a regular sequence. Thus $B$ has Hilbert series $(1+t)^6$, and homogeneous lifts of an eigenbasis of $B$ give a free basis of $\mathbb C[Z]$ over $\mathbb C[Q]$. The even-degree basis elements in degrees $0,2,4,6$ have multiplicities $1,15,15,1$, giving

$$
\phi_*\mathcal O\cong\mathcal O\oplus15\mathcal O(-1)\oplus15\mathcal O(-2)\oplus\mathcal O(-3).
$$

Since the six quadrics are invariant, the equivariant Koszul formula gives

$$
\operatorname{Tr}(\tau\mid B;t)=\frac{(1-t^2)^6}{(1-t)^3(1+t)^3}=(1-t^2)^3.
$$

The positive multiplicities are therefore $1,6,9,0$, proving (1). The $\mathcal O$ summand is the canonical unit summand, since the other summands have no global sections.

**Conditional consequence for line cylinders.** Let $S$ be a smooth projective surface parametrizing lines on $C$, with parameter map $S\to\operatorname{Gr}(2,6)$ generically finite onto its two-dimensional image. Let $V$ be the pulled-back tautological rank-two bundle and $p:U=\mathbb P_S(V)\to S$ the universal line, using the convention of one-dimensional subspaces. The evaluation map is $h:U\to C$. Form the **actual scheme fiber product**

$$
T=X\times_C U,\qquad \pi:T\to U.
$$

Then $\pi$ is finite flat. Base change of (1) gives

$$
(\pi_*\mathcal O_T)^+\cong\mathcal O_U\oplus6\mathcal O_U(-1)\oplus9\mathcal O_U(-2).
$$

Put $L=\det(V)^*$, the pulled-back Plücker bundle. It is nef and big. Both direct images of $\mathcal O_U(-1)$ vanish, and relative duality, with the stated convention, gives

$$
\omega_{U/S}=\mathcal O_U(-2)\otimes p^*L,\qquad R^1p_*\mathcal O_U(-2)=L^{-1}.
$$

Consequently

$$
H^2(U,\mathcal O_U(-2))=H^1(S,L^{-1})=0.
$$

The last equality is Serre duality and [Kawamata's vanishing theorem](https://doi.org/10.1007/BF01456407) for $K_S\otimes L$; if the parameter map is finite, $L$ is ample and Kodaira vanishing suffices. Thus canonical pullback gives

$$
H^2(U,\mathcal O_U)\xrightarrow{\sim}H^2(T,\mathcal O_T)^+.\tag{2}
$$

Now assume **$T$ itself is smooth, or is normal with rational singularities**. For an equivariant resolution $\rho:\widetilde T\to T$ in the latter case, $R\rho_*\mathcal O_{\widetilde T}=\mathcal O_T$. By (2) and complex conjugation, every invariant holomorphic two-form on $\widetilde T$ is pulled back from $U$. Writing $g:\widetilde T\to X$ and $\widetilde\pi=\pi\rho$, proper push-pull gives

$$
g_*\widetilde\pi^*b=f^*h_*b.
$$

It follows that the invariant $H^{3,1}$ Gysin image of this cylinder is contained in $f^*H^{3,1}(C)$, the single determinant-density line. **Under the stated regularity hypothesis, a pulled-back family of lines cannot serve this residual target.** This applies to any such surface of lines, not just to one chosen ruling construction.

The hypothesis is essential. Normalization of a nonnormal fiber product changes $\pi_*\mathcal O$; resolution of nonrational singularities can contribute extra $H^2(\mathcal O)$. Merely taking a smooth resolution does not preserve (2). This comment neither eliminates those possible contributions nor establishes that they can serve $[A]$. It isolates a regular setting where the source action is completely controlled, and leaves singular constructions, other families, and correspondences between sheets open.

