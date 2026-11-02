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
4. *Reductions.* All six coordinate sections $X\cap\{z_i=0\}$ are smooth sextic threefolds, so (sum principle, Ch. 3 §1.1)
   $\mathrm{served}(z_iG)=\mathrm{served}(G)$. Monomials of character $(1,\cdot)$ have odd $Y_2$-exponent, so every pure-character
   denominator of character $(1,0)$ or $(1,1)$ is $Y_2\cdot G$ and reduces to an even-$Y_2$-parity denominator of lower degree
   (this re-derives the degree-7 exclusions for those two characters from degree 6, and reduces the 21 non-Jacobian degree-9
   classes of character $(1,0)$ to the degree-8 Klein-invariant family). At degree 9 the non-Jacobian E1 kernels are
   $21/0/21/0$ by character, at degree 10 $186/126/136/108$ (dimension-forced from here on, so E1 alone stops excluding).

**Screened only (no claim of nonexistence).** (a) The affine family $K_8+J_8$: $\mathrm{ob}(K_8)\neq0$ and the linearised
E2-mod-$J$ system over the 168-parameter $Y$-even family has rank 156 (kernel = the 12 trivial directions) with
$\mathrm{ob}(K_8)\notin\mathrm{im}$, so no first-order correction exists; the full quadratic system (103 unknowns, 274 quadrics) [PENDING].
(b) Structured shapes (cones with a line or plane vertex, coordinate multiples, binary/ternary forms, polynomials in squares,
$u=Y_2^2$-polynomials, bounded $X$/$Y$-degree, the Klein-invariant subring, …), intersected with the E1-admissible space by
character for $7\le d\le 10$: nothing outside $J$ for $d\le 8$; at $d=9,10$ the survivors are coordinate multiples (reducible
by the lemma), bounded-degree and even-exponent families — no fixed-$D$ certificate found among the tested representatives [PENDING].
(c) Mixed-character degree-7 denominators: joint quadratic systems for pairs/triples of characters [PENDING].

**Reproduction.** Scripts and JSON receipts (primes, ranks, Hilbert functions) in the wave-2 packet
`wave2-excl/` (README with replay commands). Engine: exact polynomial arithmetic + python-flint; msolve 0.6, Singular.
No certificate, no Lean theorem is claimed; nothing here asserts that no certificate exists.
