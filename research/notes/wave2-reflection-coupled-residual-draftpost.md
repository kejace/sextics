# Draft discussion comment — `reflection-coupled-residual` (#154) — exclusion / rigor update (wave 2)

**Rigorous lower bound on the denominator degree, and what is only screened.** Same protected record
(`data/reflection-coupled-residual/instance.rec`, SHA256 `114c5b4f…38a2`), $A=X_0^3Y_0^3$, $q=1$, variables $(X_0,X_1,X_2,Y_0,Y_1,Y_2)$.
Write $J=(\partial_iF)$, $S/J$ the Jacobian ring, and grade everything by the Klein characters
$\chi(m)=(\deg_{Y_2}m \bmod 2,\ (\deg_{Y_0}m+\deg_{Y_1}m)\bmod 2)$ of the two sign symmetries $Y_2\mapsto -Y_2$, $(Y_0,Y_1)\mapsto -(Y_0,Y_1)$.

**Proven (over $\mathbb C$).**
1. *E1 by character.* A certificate denominator $D$ must satisfy $A D\in J$, and since $J$ is Klein-graded this holds
   for each character component of $D$ separately. The kernel of $[A]\colon (S/J)_d^{\chi}\to(S/J)_{d+6}^{\chi+(0,1)}$ is zero for
   $d\le 7$ in all four characters and, at $d=8$, one-dimensional in character $(0,0)$ and zero in the other three.
   (Computed mod $p=2147483587$ with the good-reduction check $\operatorname{rank}J_n^{\chi}=\dim S_n^{\chi}-$Fermat count in every
   degree used, so the modular kernel dimension bounds the rational one from above; agrees with the earlier computations mod
   $2^{31}-1$, $2^{31}-19$ and $32003$.) Hence: no denominator of degree $\le 4$; every denominator of degree $5,6,7$ lies in $J$;
   every denominator of degree 8 lies in $J_8\oplus\mathbb C\,K_8$, where $K_8$ is the rationally reconstructed class
   (37 terms, Klein-invariant) — and $A K_8\in J$ holds exactly over $\mathbb Q$ (rational numerators verified).
2. *Hessian relaxation of E2 on the Jacobian family.* For $D=g(F)=\sum g_i\partial_iF$ the numerators $N=A g$ solve E1, and
   modulo $J$ the second condition is independent of the choice of E1 solution (Koszul corrections lie in $J$):
   $\mathrm{ob}(g)=-[A\,g^{\mathsf T}\mathrm{Hess}(F)\,g]\in(S/J)_{2d}$ must vanish. It is invariant under $g\mapsto g+zG$ (these give
   $D\in(F)$), so it lives on $W=V/T$.
   * $d=5$: the 20 classes $[A\,\partial_i\partial_jF]$ ($\partial_{X_2}\partial_{Y_2}F=0$) are independent in $(S/J)_{10}$, so $g_ig_j=0$ for all pairs and $g=0$.
   * $d=6$: on $W$ (35 unknowns, all characters jointly, i.e. mixed-parity $D$ included) the 345 independent quadrics generate
     a zero-dimensional ideal mod $p$ (msolve, 709-element Gröbner basis; Singular: dim 0, multiplicity 1060, Hilbert function
     $1,35,285,477,225,36,1,0$). Since $\mathrm{HF}_{\mathbb Q}\le\mathrm{HF}_p$, the rational ideal contains all forms of degree 7, so
     $V(\mathrm{ob})=T$ over $\mathbb C$.
   * $d=7$, single character: the same computation on $W_\chi$ (dims 42/38/22/18, 274/281/132/117 independent quadrics) gives
     zero-dimensional ideals for $\chi=(0,1),(1,0),(1,1)$ [and $(0,0)$: PENDING].
   **Consequence: every certificate denominator has degree $\ge 7$; a degree-7 denominator must have components in at least two
   Klein characters [or lie in character $(0,0)$: PENDING].**
3. *The class $K_8$.* Fixed-$D$ test on the Klein-invariant sector (22016 equations, 4256 unknowns): inconsistent mod $p$, and the
   modular rank equals the rational rank because 476 explicit rational kernel vectors $(zG,0)$, $(F m e_i,\ D\partial_im-m\partial_iD)$
   match the modular nullity. So $K_8$ is not a certificate denominator (over $\mathbb C$). Second prime and the full 40301×7497
   system agree. [Singular locus of $X\cap V(K_8)$: PENDING.]

**Screened only (no claim of nonexistence).** [mixed-character degree 7 / degree-8 families / structured families d=7–12: PENDING]

**Reproduction.** Scripts and JSON receipts (primes, ranks, Hilbert functions) in the wave-2 packet
`wave2-excl/` (README with replay commands). Engine: exact polynomial arithmetic + python-flint; msolve 0.6, Singular.
No certificate, no Lean theorem is claimed; nothing here asserts that no certificate exists.
