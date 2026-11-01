# Four plus two · suspension remainder  (https://github.com/hodge-computer/sextics/discussions/120)


## COMMENT DC_kwDOVLOOJc4BJDcS by benediamond at 2026-10-23T17:11:49Z (updated 2026-10-23T17:11:49Z)

### The binary contribution is a three-dimensional Prym, and the visible elliptic quotient misses the source

The $4+2$ splitting can be sharpened to an explicit character calculation. It identifies which part of the small binary factor a new correspondence would need to use.

Write the [fixed equation](https://github.com/hodge-computer/sextics/blob/faa9df2ec6406e4eacc135f677cf3711315d0e41/data/four-plus-two-residual/instance.rec) as $F=H(X_0,X_1,Y_0,Y_1)+B(x,y)$, where $(x,y)=(X_2,Y_2)$ and

$$
B(x,y)=x^6+\frac{15}{7}x^4y^2+\frac{10}{7}x^2y^4+y^6.
$$

The natural binary triple cover is $C_B:z^3=B(x,y)\subset\mathbf P(1,1,2)$. On $x=1$, with $t=y/x$ and $w=z/x^2$, it is

$$
w^3=b(t),\qquad b(t)=1+\frac{15}{7}t^2+\frac{10}{7}t^4+t^6.
$$

The polynomial is separable: modulo $5$ it is $1+t^6$, coprime to its derivative. The curve has genus four. A basis of its holomorphic differentials is

$$
\frac{dt}{w},\qquad \frac{dt}{w^2},\qquad \frac{t\thinspace dt}{w^2},\qquad \frac{t^2dt}{w^2}.
$$

At infinity $w$ has pole order two, so $t^i dt/w^j$ is regular precisely when $i\le2j-2$ for $j=1,2$; at the six finite branch points all four displayed forms are regular. This also gives the basis directly.

Let $\sigma(t,w)=(-t,w)$ and $\delta(t,w)=(t,\zeta_3w)$. The $\sigma$-invariant holomorphic line is $t\thinspace dt/w^2$. The anti-invariant space is

$$
H^{1,0}(C_B)^-=\left\langle\frac{dt}{w},\frac{dt}{w^2},\frac{t^2dt}{w^2}\right\rangle,
$$

with deck multiplicities $1$ for $\zeta_3^2$ and $2$ for $\zeta_3$. Therefore the anti-invariant abelian factor, equivalently the connected Prym of $C_B\to C_B/\sigma$, has dimension three and $\mathbf Q(\zeta_3)$ signature $(2,1)$, up to the ordering of embeddings.

The invariant quotient itself is the elliptic curve

$$
w^3=1+\frac{15}{7}v+\frac{10}{7}v^2+v^3,\qquad v=t^2.
$$

Its deck automorphism has order three and fixed points, so its $j$-invariant is $0$. But it is the wrong involution sector for the target. The prescribed numerator $A=X_0^3Y_0^3$ is unchanged by $Y_2\mapsto-Y_2$, while the projective volume changes sign. Hence the residue is anti-invariant under this reflection. Projection to the displayed elliptic quotient kills the binary component of this residue.

One can also track the exact deck character in the usual split-cover residue decomposition. In the ordinary sixth covers $U^6=H$ and $V^6=B$, the core $H^{2,1}$ residue has numerator $AU$, and the binary $H^{1,0}$ residue has numerator $V^3$. Their deck characters, including volume, are respectively $\zeta_6^2$ and $\zeta_6^4$. They descend under $u=U^2$ and $w=V^2$ to opposite order-three characters on

$$
T_H:u^3=H\subset\mathbf P(1,1,1,1,2),\qquad C_B.
$$

The binary holomorphic factor is the $\zeta_3^2$ line $dt/w$, which is indeed anti-invariant. The core factor is also anti-invariant under $(Y_0,Y_1)\mapsto(-Y_0,-Y_1)$: its volume is even and $A$ is odd. These two signs agree with the invariance of the original residue under the global $Y$ sign.

This explains concretely why the elementary ternary-split companion argument has a missing support contribution. The weighted companions are

$$
W_H:\space H+s^3+t^3=0\subset\mathbf P(1,1,1,1,2,2),\qquad
W_B:\space B+s^3+t^3=0\subset\mathbf P(1,1,2,2).
$$

Weighted adjunction gives a Fano fourfold for $W_H$, while the quasismooth model of $W_B$ resolves to a K3 surface. The $3+3$ construction has two Fano-supported factors and contracts a common elliptic factor; here the K3 side has no automatic corresponding support gain. Merely replacing its binary input by the visible elliptic quotient would remove the desired character.

The precise remaining question is whether this particular Prym/K3 component admits an additional algebraic correspondence that supplies the missing support. The signature calculation alone does not identify this Prym with a product of CM elliptic curves, or supply a definite-unitary finite-monodromy argument in a varying family. It also does **not** prove that this individual Prym is simple or exclude an exceptional splitting. No such splitting or serving correspondence is established here, and no nonexistence statement about certificates follows.

This uses the project's [split and K3 framework](https://github.com/hodge-computer/sextics/blob/faa9df2ec6406e4eacc135f677cf3711315d0e41/docs/03a-split-and-k3.md), with explicit differential and character calculations above. ChatGPT-assisted investigation for `benediamond`, including separate AI-agent checks; no independent human review, Lean proof, or solve is claimed.

