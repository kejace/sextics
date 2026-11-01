# Reflection · complementary character  (https://github.com/hodge-computer/sextics/discussions/154)


## COMMENT DC_kwDOVLOOJc4BIzT9 by fqlx at 2026-10-20T23:52:46Z (updated 2026-10-20T23:52:46Z)

Partial progress on the unchanged target: **a dominant genus-one construction and explicit partial primitive remove the infinity component of the actual source. A compact elliptic equation remains unsolved; this is not a certificate.**

I used $F$ verbatim from `data/reflection-coupled-residual/instance.rec`, $A=X_0^3Y_0^3$, $q=1$, and variable order $(X_0,X_1,X_2,Y_0,Y_1,Y_2)$. The record's SHA256 is `114c5b4f2c4c4ba36456a88fc19925d28e51399b4957f417ea5ca157904d38a2`. The companion quotient's Jacobian $2Y_2$ still has the opposite volume character and is not used as a solution.

Put $H(x,u)=F(x,z)$ with $u=z^2$. At the exact coupling $5/7$, a divided conic incidence imposing all thirteen coefficients of $H(\varphi(t),U(t))=0$ has a nonsingular algebraic characteristic-zero lift. A further slice makes $Y_0=\varphi_3$ linear and $U$ monic, with its four lower coefficients as independent parameters. Exact incidence/evaluation minors and the endpoint-preserving $\mathrm{GL}_2$ action establish dominance on the generic open of the **same smooth component**. The new evaluation determinant vanishes at the original seed; the nonzero derivative of $[t]Y_0$ and gauge determinant $-4[t]Y_0$ supply generic transversality. This is not a dominance assumption inferred from one modular curve.

The physical cover is $w^2=U(t)$, squarefree of genus one. Rebuilding the normal source on the new slice gives $p_A/(2w)$, retaining all parameter derivatives, determinant variations, source jet, and Hessian correction. Write

$$
U=t^4+u_3t^3+u_2t^2+u_1t+u_0,\qquad E_2=t^2+\tfrac12u_3t.
$$

Explicit rational currents remove the polynomial-exact and infinity terms. The infinity correction uses a normal-transgression identity $\rho=\sum_jD_jK_j$: moving infinity to the source-zero locus $Y_0=0$ kills the rational normal's remaining infinity contribution.

The exact unrestricted equation remaining **on this construction** is

$$
\sum_{\alpha\in\{u_0,u_1,u_2,u_3,t\}}
\left(2U D_\alpha C_\alpha-(D_\alpha U)C_\alpha\right)
=U(c_0+c_1E_2),\qquad C_\alpha\in E(t),
$$

where $E$ is the actual incidence coefficient field, finite separable over $\mathbf C(u_0,u_1,u_2,u_3)$. All parameter directions, arbitrary rational denominators, and moving time poles remain allowed. No denominator-degree bound or condition $Q=0$ is imposed. The actual compact coefficients are nonzero: $5c_0\equiv127\pmod{625}$ and $c_1/5\equiv16\pmod{625}$, with a proved precision budget. This excludes time-only completion of these particular currents, including after finite algebraic extension; it does not exclude other parameter currents or a target certificate. A sufficient compact-system reduction to one scalar differential equation is also supplied, but unsolved.

[Reproduction packet](https://github.com/user-attachments/files/33743911/reflection_residual_genus_one_continuation.zip), SHA256 `eaeeeee98815a5ded398e420f3afda8f81f752f7311ac295f358bfe51099adc7`. With SymPy 1.14.0, run `python3 verify_all.py` from the extracted continuation packet root. Local exact replay passed. The scalar correction rests on the reviewed geometric transgression argument; the symbolic partial-primitive checker verifies the subsequent subtraction conditional on that identity.

No residual certificate or Lean theorem was obtained. Protected sandbox/kernel verification was not run; no residual solve PR, CI submission, merge, or acceptance is claimed. This builds on the accepted holdout-bundle construction ([#161](https://github.com/hodge-computer/sextics/pull/161)), companion construction ([#150](https://github.com/hodge-computer/sextics/pull/150)), and repository PolynomialCurves machinery, with credit to their authors and contributors. Research and comment assisted by OpenAI ChatGPT GPT-6 Pro and Codex.



### REPLY DC_kwDOVLOOJc4BIz8q by fqlx at 2026-10-21T02:49:42Z

Continuation of the genus-one calculation: **the actual source now reduces to one holomorphic coefficient, but no certificate was obtained.**

Keep the same protected record, $A=X_0^3Y_0^3$, and pole parameter $q=1$. On $w^2=U(t)$, work over the actual finite incidence field $E/\mathbb C(u_0,u_1,u_2,u_3)$, retaining all four parameter derivatives. Deform the equation by $G_\varepsilon=G+\varepsilon A$ at fixed $U$ and gauges. For the full moving-normal determinant $\mathscr J$, a cofactor degree argument proves that $\mathscr J$ and $\dot{\mathscr J}$ are time-constant, also over dual numbers. Explicit cofactor currents leave

$$
g=-\frac12\dot{\mathscr J}.
$$

An explicit normal-change current connects the two Sylvester normals. For the $(G_2,G_4)$ normal, a boundary current, a closed odd moving-pole current, and a compact correction give the full identity

$$
\frac{p_A}{2w}
=\mathrm{div}_5\left(\frac{C^{(0)}+C^{(1)}}w\right)+\frac{R_{\ast}}w.
$$

All five current components, including the rational time correction and moving-pole derivatives, are explicit. The exact algebraic construction is accompanied by an independent replay of **64 substantive checks**. All nine coefficients of the cleared degree-at-most-eight identity vanish modulo $5^{143}$ (minimum absolute $5$-adic precision 143). On the established characteristic-zero branch, $v_5(R_{\ast})=0$ and $R_{\ast}\equiv529\pmod{625}$, so this selected combination is incomplete. This does not exclude other currents or certificates.

One concrete “rolling wheel” finishing condition is to find rational $h,k\in E$ with $g=T(h)+V(k)$. In depressed coordinates $U=x^4+px^2+qx+r$, $a=u_3$, $x=t+a/4$, set

$$
T=\partial_a,\qquad V=6q\partial_p+(8r-2p^2)\partial_q-pq\partial_r.
$$

Then $C=hT+k\bigl(V+(p+2x^2)\partial_x\bigr)$ gives $\mathrm{div}_5(C/w)=g/w$. Here $T$ fixes $p,q,r,x$; the displayed $q$ is a quartic coefficient, distinct from the protected pole parameter. No such decomposition was found; scaling leaves an Euler resonance. Other parameter currents and arbitrary additional moving poles remain allowed. [PR165](https://github.com/hodge-computer/sextics/pull/165) is useful source-transport machinery, but matching its polynomial-source hypothesis is still open here.

[Reproduction files and instructions](https://gist.github.com/fqlx/ae3df1cc3d0a3857b7b69f2ad80fc983) include the algebraic field/source definitions, full currents, proof, raw jets, verifier and receipts. No residual Lean/kernel verification is claimed. This continues the earlier work and draws on PRs #161, #150 and #165, with credit to their authors. Research and comment assisted by OpenAI ChatGPT GPT-6 Pro and Codex.



### REPLY DC_kwDOVLOOJc4BI-ZY by benediamond at 2026-10-22T18:13:05Z

### Radial audit of the source-root currents

Following `fqlx`’s [source-root construction](https://github.com/hodge-computer/sextics/discussions/154#discussioncomment-19087146), I obtain a simplification: its published $Y_0$-root Cartan forcing already has an explicit **compact** primitive. The actual instance remains open.

Use the same incidence field, all four physical parameter derivatives, and compact basis $(dt/w,(x^2+p/6)dt/w)$. For the actual radial field $S$, define

$$
\mathcal D_S f=S(f)+(\mathrm{div}S)f+\Gamma_Sf,
\qquad \Gamma_S=\mathrm{diag}(1,-1),
\qquad \mathrm{div}S=-10+T\nu.
$$

The covariant Euler identity is

$$
\mathcal D_S f=mf,\quad m\ne0
\quad\Longrightarrow\quad
\mathrm{Div}_{\nabla}\bigl(S^j f/m\bigr)=f.
$$

For the published current $\theta_Y=v^{-1}\iota_{\widetilde W}\mu$, the incidence identities give $d\theta_Y=0$, $\mathcal L_S\mu=0$, $[S,\widetilde W]=0$, and $S(v)=-2v$. Hence $\mathcal L_S\theta_Y=2\theta_Y$ and $\theta_Y=\frac12d(\iota_S\theta_Y)$. The normalized third-kind class is radially horizontal, so its compact forcing satisfies $\mathcal D_Sk_Y=2k_Y$; the compact primitive is simply $C_Y^j=S^jk_Y/2$.

I also constructed a closed Cartan current from the **two finite $X_0$ roots**, retaining their moving-root derivatives and normalized quadratic trace. It has the same weight and compact primitive. In contrast, the actual fixed-infinity compact source $\mathbf h$ satisfies $\mathcal D_S\mathbf h=0$. Thus constant combinations of these weight-two forcings cannot equal the nonzero source. Their previously observed arithmetic improvement does not by itself remove the weight-zero compact class.

A useful weight correction would require an even rational multiplier $f$ with $Sf=-2f$ and $\widetilde W(f)=0$. At the $Y_0$ root, writing $R=U(z)$, $r=R/B^4$, and $f=B^2\varphi$ on the normalized field $K=E^S$, this becomes

$$
D\varphi=\frac{Dr}{2r}\varphi,
$$

where $D$ is induced by $\widetilde W$. Equivalently, $\varphi/\sqrt r$ must be an odd first integral. The evident multiplier $\sqrt R$ loses the required parity.

Certified nonzero $5$-adic minors exclude constant combinations of $B^2$ and the four ordinary root coordinates for both currents; the $Y_0$ check uses one nearby Hensel point on the same component. This is only a five-dimensional exclusion, not an obstruction to arbitrary rational multipliers or certificates. The actual finishing equation $g=T(h)+V(k)$ remains unsolved.

The attachment contains the derivation, portable scripts, raw receipts, pinned target, and necessary public upstream dependencies. The weight identities are **structural deductions from the incidence identities**; zero congruences are implementation audits, not global polynomial certificates. No Lean proof is claimed.

[Derivation and reproducible checks](https://github.com/user-attachments/files/33835113/reflection-coupled-residual-radial-audit-2026-10-22.zip)

Research and writing assisted by ChatGPT (GPT-6 Astra Pro), for `benediamond`; credit to `fqlx` and the contributors underlying the incidence construction.



## COMMENT DC_kwDOVLOOJc4BJTkW by benediamond at 2026-10-27T12:25:35Z (updated 2026-10-28T15:41:55Z)

**Rational integration on the incidence curve, with independent compact coefficients and closed residue freedom.** For the fixed target `reflection-coupled-residual`, with source $A=X_0^3Y_0^3$ and certificate pole parameter $q_{\mathrm{cert}}=1$, this construction reduces the remaining source equation to rational integration on an algebraic curve. The zero-residue subproblem has a proved finite pole bound. The full construction retains three residue functions that can change the differential being integrated. **No rational primitive or certificate has been obtained.** The [research notes, supporting programs, data and checks (ZIP)](https://github.com/user-attachments/files/34087387/reflection-residual-research.zip) give the derivations and reproducible inputs.

The starting point is `fqlx`'s [incidence and actual-source construction](https://github.com/hodge-computer/sextics/discussions/154#discussioncomment-19084541) and its [source-root continuation](https://github.com/hodge-computer/sextics/discussions/154#discussioncomment-19087146). The four logarithmic compact currents have **independent rational coefficient functions**. Relations among their vector fields leave two effective compact directions; eliminating one coefficient gives

$$
du=\omega_g,\qquad g=(g_Y,g_0,g_1).
$$

Here $d$ differentiates along the incidence curve over the field of centered quartic parameters. This curve parametrizes compatible incidence coefficients; it is distinct from the genus-one physical cover $w^2=U(t)$. The functions $g_Y,g_0,g_1$ multiply the $Y_0$ residue current and the two channels associated with the quadratic $X_0$ block. Their closure equations, required weights and derivatives remain part of the system. They are independent unknowns subject to those equations, not three arbitrary constants. Formal local solvability holds at a verified regular germ; global rationality is the unresolved requirement. One rational solution that passes reconstruction would suffice for the certificate goal.

For the **baseline $g=0$**, let $z$ be the chosen incidence coordinate, $\delta$ its Jacobian minor, and $\mathcal T=\delta\partial_z$ the polynomial tangent derivation given by cofactors. Write the verified source as $\omega_0=(N/\delta^3)\,dz$. On the smooth projective normalization, every rational primitive satisfies

$$
f=\delta u\in L(65H_\infty),\qquad
\delta\mathcal T(f)-\mathcal T(\delta)f=N.
$$

The divisor $H_\infty$ records the largest pole order of the affine incidence coordinates at each point at infinity. The Riemann–Roch space $L(65H_\infty)$ consists of functions with those bounded poles and no others; it is finite-dimensional over the parameter field. Thus the displayed linear equation in that space is a complete rational-primitive test for the baseline on the selected component. The number 65 follows from the actual source and determinant degrees. It is not an imposed denominator-degree cutoff. The factor $\delta$ is an incidence Jacobian, not a proposed physical certificate denominator.

The pole-bound proof covers **every normalized place**, including higher ramification and branches above affine singularities. Its key step is to retain the source as $\omega_0=(H-L(V))\,dz$, with $V=B/\delta$ and $H=C/\delta$, and change the lift of the base derivation $L$ in a local uniformizer. The apparent high-order poles split into an exact derivative and a controlled remainder. This bounds a primitive if one exists; it does not assert that all residues or global obstructions vanish. See **One Jacobian factor suffices** in the notes.

Several structural results help organize the global calculation:

- The two block-sign symmetries give a quotient presentation with twelve coefficient ratios and eleven equations. Character averaging preserves solvability of the linear current system and retains all three residue functions. The baseline unknown $f=\delta u$ is invariant.
- A quartic degeneration certifies at least sixteen actual places at infinity across regular rank-three components, with cubic ramification and explicit source orders. At these places, parity excludes local logarithmic terms for the source-character system, and the baseline quotient pole bound for $f$ is 31. The general bounds retain the residue functions' actual poles. These places are neither a complete boundary list nor all assigned to one component.
- An explicit kernel classification excludes the described family of source-zero conic covers from the closure of every rank-three component. It covers all nonzero scales, all relevant sixth-root choices and all coprime quadratic factorizations. Other exceptional affine loci remain unclassified.

A regular component is selected by a rigorously isolated algebraic point. The resulting implicit description does not yet provide effective function-field arithmetic or a basis of the normalized section space. Likewise, the 91 certified points in a quotient fiber do not determine one component's degree. A complementary period calculation has an independently checked source evaluator and fourteen certified simple branch points, but its connecting paths and quadrature do not establish a nonzero period or rational exactness.

The baseline therefore faces a concrete **computational bottleneck**: construct effective arithmetic on one component, compute the complete normalized section space, and solve the displayed equation. Direct F4SAT saturation on a fiber at $(p,q,r,z)=(4,1,6,1)$ over $\mathbb F_{65537}$ was stopped without a finite quotient basis as successive sparse matrices grew to 675,752, 2,502,378 and 7,780,680 columns. A validated vectorized eight-thread replay reproduced that growth. This experiment did not reach integration; its failure is no obstruction to a certificate. **Direct saturation and the remaining computational gap** records the chart restrictions, controls and logs.

The **full residue problem also has a conceptual gap**. A proved finite primitive bound applies after specifying a closed residue tuple, using its actual denominators. No uniform pole bound for all unknown closed residue tuples has been established. Consequently, failure of the baseline would leave open whether residue currents cancel its differential class. Useful continuation requires either substantially more effective component arithmetic or an understanding of the closure equations and differential class that makes this cancellation tractable. Any candidate must still satisfy the weights, reconstruct the currents, descend, match the original source and admit legal denominator clearing.

This work builds on `fqlx`'s contributions and the collaborators credited in the linked packets. Research, computation and preparation assisted by Codex for `benediamond`.

