# Simultaneous permutations · standard sector  (https://github.com/hodge-computer/sextics/discussions/177)


## COMMENT DC_kwDOVLOOJc4BJDbi by benediamond at 2026-10-23T17:09:14Z (updated 2026-10-23T17:09:14Z)

### A source-preserving reduction to one transposition quotient

There is a small representation-theoretic reduction of this target which avoids averaging away its standard component. It does not yet give divisor support, but it gives a single quotient source on which to concentrate the geometric search.

Put $T_i=X_i^3Y_i^3$ and $B_i=T_0+T_1+T_2-3T_i$. Thus the prescribed numerator is

$$
A=T_0-T_1=\frac{B_1-B_0}{3}.
$$

Each $B_i$ is fixed by the transposition of the other two indices, and the three $B_i$ are conjugate under the simultaneous $S_3$ action. Consequently it suffices to support the one source

$$
B=B_2=T_0+T_1-2T_2
$$

on the quotient by $(01)$: pull the support back to $X$, apply the two appropriate permutations, and take the displayed linear combination. Simultaneous permutations of the $X$ and $Y$ blocks have determinant $1$, so these identities apply to the residue classes themselves, with no extra volume sign. In group-algebra terms,

$$
e_{\mathrm{std}}=\frac23\sum_{s\text{ a transposition}}\frac{1+s}{2}-2e_{S_3}.
$$

The full $S_3$ average kills $A$, whereas the proper-subgroup averages recover it. This reduction supplies the source identity; the missing mathematical step is to support the descended $B$.

**An explicit quotient chart, including the density.** Set

$$
X_0=\frac{p+r}{2},\quad X_1=\frac{p-r}{2},\quad X_2=a,
\qquad
Y_0=\frac{q+tr}{2},\quad Y_1=\frac{q-tr}{2},\quad Y_2=b,
\qquad u=r^2.
$$

The transposition negates $r$ and fixes $p,q,a,b,t,u$. Substitution in the [fixed equation](https://github.com/hodge-computer/sextics/blob/faa9df2ec6406e4eacc135f677cf3711315d0e41/data/diagonal-permutation/instance.rec) gives a polynomial $F_s(p,q,a,b,u,t)$, obtained simply by replacing each even power $r^{2j}$ by $u^j$. Its weighted degree is six, with $p,q,a,b$ of weight one, $u$ of weight two, and $t$ of weight zero. The descended numerator is

$$
B_s=\frac{(pq+tu)^3+3(pq+tu)(pt+q)^2u}{32}-2a^3b^3.
$$

On $a=1$, with the usual projective volume in the order $(X_0,X_1,X_2,Y_0,Y_1,Y_2)$, the affine volume is exactly

$$
dX_0\wedge dX_1\wedge dY_0\wedge dY_1\wedge db
=-\frac18\thinspace dp\wedge dq\wedge du\wedge dt\wedge db.
$$

Thus the quotient source is the residue of $-B_s\thinspace dp\wedge dq\wedge du\wedge dt\wedge db/(8F_s^2)$, rather than a numerator obtained by ignoring a variable-dependent Jacobian. The global $Y$ sign acts by $(q,b,t)\mapsto(-q,-b,-t)$; both $B_s$ and the five-volume are odd, so the residue remains invariant.

**The first surface fibration retains three holomorphic forms.** Fix $t$ and $\lambda=b/a$. The natural projective surface fiber is a sextic in $\mathbf P(1,1,1,2)$ with coordinates $(p,q,a,u)$. This is genuinely a smooth generic surface for the fixed coefficients: modulo $5$ all coupling terms and the mixed binomial coefficients disappear, leaving

$$
\frac{p^6+q^6+(1+t^6)u^3}{32}+(1+\lambda^6)a^6=0.
$$

At $t=\lambda=0$ its weighted cone is smooth away from the vertex and the surface misses the ambient singular point. Openness of smoothness therefore gives a smooth generic characteristic-zero fiber. Adjunction gives $K=\mathcal O(1)$ and $p_g=3$. The global $Y$ sign moves the base by $(t,\lambda)\mapsto(-t,-\lambda)$; it does not act as a further quotient on this generic geometric fiber.

So this particular surface fibration does not meet the whole-fiber Tate criterion. This is not an obstruction to the target: a smaller algebraically selected fiber component, another fibration, or a new correspondence could still support $B$. One also cannot get $B$ merely from the exceptional divisor over the fixed surface of the transposition: its $H^2(S)(-1)$ blowup summand is killed by pushforward to the original fourfold. An additional correspondence mixing it with the original primitive summand would need an explicit source-action calculation.

The polynomial substitution, numerator, projector identity, parity, and volume determinant were checked with exact rational arithmetic. This is a geometric reduction and a scoped limitation, with no certificate or Lean theorem claimed. ChatGPT-assisted investigation for `benediamond`, including separate AI-agent checks; no independent human review is claimed. The starting framework is the project's [geometric support chapter](https://github.com/hodge-computer/sextics/blob/faa9df2ec6406e4eacc135f677cf3711315d0e41/docs/03-geometric-support.md).

