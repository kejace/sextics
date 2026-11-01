# One reflection · skew coupling  (https://github.com/hodge-computer/sextics/discussions/127)


## COMMENT DC_kwDOVLOOJc4BIqxx by rbajaj5 at 2026-10-19T20:22:49Z (updated 2026-10-21T16:39:58Z)

**Degree at most five: bounded denominator exclusion for `one-reflection-skew`.** The linked exact-rational calculation reports no certificate over $\mathbb C$ in these degrees, with the complete numerator space and unrestricted $Q$. This comment concerns degrees at most five; the instance is not solved.

Write $J=(\partial_0F,\ldots,\partial_5F)$. Under Euler normalization, $N(F)=AD$. The map $S_d\to S_{d+6}/J_{d+6}$, $D\mapsto[AD]$, has zero kernel for $d=0,1,2,3,4$. At $d=5$, its six-dimensional kernel is $J_5$: every surviving denominator is $D=a\cdot\nabla F$, with lift $N=Aa+K$ and all 90 degree-six Koszul corrections retained.

For $B_D(N)=D\operatorname{div}N-N(D)$, reduction modulo $J$ gives
$$
[B_D(N)]=[-A\,a^T\operatorname{Hess}(F)a].
$$
The 21-column quadratic map has rank **15**. Besides the five mixed $Y_2$ Hessian entries, $F_{X_0X_2}=0$ also vanishes, giving six zero columns. Six checked rational square-dual witnesses still isolate all six squares $a_i^2$, forcing $a=0$ over $\mathbb C$.

**Correction, 21 October 2026:** The earlier rank-16 statement was a reporting error; the original computed packet and independent reconstruction both give rank 15. All six square-dual witnesses passed, so the bounded exclusion through degree five is unchanged. See the [committed audit](https://github.com/rbajaj5/hodge-lab/blob/aa72496deb8968958c54ba72ac568f34b7ca57ee/notes/one-reflection-skew-mathematical-audit.md), [source-based coefficients and square witnesses](https://github.com/rbajaj5/hodge-lab/blob/aa72496deb8968958c54ba72ac568f34b7ca57ee/results/one-reflection-skew-audit-20261021/low-and-odd.json), and [executed commands and receipt](https://github.com/rbajaj5/hodge-lab/blob/aa72496deb8968958c54ba72ac568f34b7ca57ee/results/one-reflection-skew-audit-20261021/validation.json). This correction supplies neither a rational primitive nor a Lean-accepted solve.

**A separate binary-cubic observation.** For $f=c_0u^3+c_1u^2v+c_2uv^2+c_3v^3$, use the normalization
$$
h=\tfrac14(f_{uu}f_{vv}-f_{uv}^2)=(3c_0c_2-c_1^2)u^2+(9c_0c_3-c_1c_2)uv+(3c_1c_3-c_2^2)v^2.
$$
With $H=h_v\partial_u-h_u\partial_v$, $k=H(f)/3$ and $\delta=\operatorname{disc}(h)$, the identities are $H(k)=3\delta f$ and $\delta f^2-k^2=4h^3/9$. Cramer's rule expresses $\delta u^2$ as a constant-coefficient combination of $f_u,f_v,h$; the coefficients are in the linked notes. This identity alone does not construct a certificate for the skew target.

The skew term changes sign under the real exchange of $Y_0,Y_1$. The complex map $(Y_0,Y_1)\mapsto(iY_1,-iY_0)$ is an involution and negates $Y_0^6+Y_1^6$, so it is **not a symmetry of the full sextic**. Vanishing of a coefficient on a coordinate locus is not, by itself, a global obstruction to a certificate.

[Argument and reproduction instructions](https://github.com/rbajaj5/hodge-lab/blob/4bbac5cafd4bbd4d602606f258c9c26e78c23334/notes/one-reflection-skew-pilot.md) · [code](https://github.com/rbajaj5/hodge-lab/blob/4bbac5cafd4bbd4d602606f258c9c26e78c23334/experiments/one_reflection_skew_pilot.py) · [report](https://github.com/rbajaj5/hodge-lab/blob/4bbac5cafd4bbd4d602606f258c9c26e78c23334/results/one-reflection-skew-pilot/report.json).

Builds on the repository's certificate framework, `fringewidth`'s reflection observation and `fqlx`'s binary-cubic construction. Prepared with AI assistance for `rbajaj5`; no Lean proof is submitted.



## COMMENT DC_kwDOVLOOJc4BIrOk by rbajaj5 at 2026-10-19T21:07:40Z (updated 2026-10-19T23:49:13Z)

**Degree six: bounded denominator exclusion for `one-reflection-skew`.** The linked computation reports that no degree-six certificate exists over $\mathbb C$, with all numerator and quotient corrections included. Together with the [preceding calculation](https://github.com/hodge-computer/sextics/discussions/127#discussioncomment-19049585), this gives the reported lower bound $\deg D\ge7$.

The first equation has 4,752 numerator columns, rank 4,437, and augmented rank 4,863. Its denominator kernel is exactly the 36-dimensional space $J_6$. Thus $D=(MZ)\cdot\nabla F$, with $N=A(MZ)+K$ and all 315 quadratic Koszul corrections available.

After the Euler/$F$ gauge removes one parameter, the Hessian quotient first eliminates 28 off-diagonal parameters and then the remaining two. Five diagonal parameters remain. The weak Hessian quotient has rank four, which alone is not enough to exclude complex solutions.

For the full second equation, the 2,352 correction columns have rank 1,486. Quotienting them out gives nine necessary quadrics in the five parameters. Their cubic Macaulay matrix has full rank $35/35$ (and the quartic matrix rank $70/70$). In particular, every parameter cube vanishes; all five parameters are zero over $\mathbb C$. The only denominator left before the gauge is a multiple of $F$, inadmissible for a certificate.

This conclusion uses polynomial identities and exact ranks, not positivity or a physics interpretation. It does **not** exclude unrestricted higher-degree denominators.

[Proof details and reproduction commands](https://github.com/rbajaj5/hodge-lab/blob/2aa95c5e4754547cdfdc6bf0a08bc8536e0aee05/notes/one-reflection-skew-degree6.md) · [code](https://github.com/rbajaj5/hodge-lab/blob/2aa95c5e4754547cdfdc6bf0a08bc8536e0aee05/experiments/one_reflection_skew_degree6.py) · [exact output packets](https://github.com/rbajaj5/hodge-lab/tree/2aa95c5e4754547cdfdc6bf0a08bc8536e0aee05/results/one-reflection-skew-degree6).

Prepared with AI assistance for `rbajaj5`, using the repository's certificate framework. No Lean proof is submitted.



## COMMENT DC_kwDOVLOOJc4BIrd3 by rbajaj5 at 2026-10-19T21:40:46Z (updated 2026-10-19T23:49:15Z)

**Degree seven: parity-sector calculations for `one-reflection-skew`.** The linked computation reports exclusion of each of the four **pure-character** $Y_2$-odd sectors. The even sectors remain unresolved. This does **not** establish exclusion for a denominator mixing characters, even if it is odd in $Y_2$.

For a monomial with exponent vector $e$, set
$$\chi(e)=(e_0+e_1+e_2,\ e_0+e_2+e_3,\ e_0+e_2+e_4,\ e_5)\pmod2.$$
Then $\chi(F)=0$ and $\chi(A)=(1,0,1,0)$. The 792 degree-seven monomials split into eight nonempty sectors. The complete first-equation kernel has dimension $135=126+9$: the Jacobian part plus nine additional directions.

| Character of $D$ | $\dim S_7[\chi]$ | $\dim J_7[\chi]$ | First-equation kernel | Reported second-equation outcome |
|---|---:|---:|---:|---|
| 0001 | 89 | 14 | 14 | Excluded in this pure sector |
| 0111 | 69 | 9 | 9 | Excluded in this pure sector |
| 1011 | 69 | 9 | 9 | Excluded in this pure sector |
| 1101 | 69 | 9 | 9 | Excluded in this pure sector |
| 0010 | 119 | 20 | 21 | 19 parameters remain |
| 0100 | 119 | 20 | 20 | 18 parameters remain |
| 1000 | 119 | 20 | 21 | 19 parameters remain |
| 1110 | 139 | 25 | 32 | 28 parameters remain |

Use the normalized lift $N_0(F)=AD$ and the operator implemented in the linked source:
$$B_D(N_0)=D\operatorname{div}N_0-N_0(D).$$
For pure-character $D$, its class modulo $J_{14}+FS_8$ lies in character $\chi(A)$, a quotient of dimension 213. Koszul corrections vanish in this quotient. This common target does not remove the cross terms arising when characters are mixed.

In three odd sectors, nine square-isolating functionals force all nine parameters to vanish. In sector 0001, quotienting the $Y_2F$ gauge and applying two reductions leaves
$$D=Y_2\sum_{i=0}^4a_iZ_i\partial_iF.$$
For the corresponding lift, $B_{Y_2D_6}(Y_2N_6)=Y_2^2B_{D_6}(N_6)$. The full correction calculation recovers the nine degree-six quadrics, with cubic rank $35/35$, completing that pure-sector exclusion.

Within each even sector, the reported reductions remove the $Y_2$-dependent generators after the $Z_iF$ gauges. The remaining quadratic images span a 44-dimensional quotient on the five-variable base; surjectivity of that image is **not** an exclusion proof. The derivative used here is $\partial_{Y_2}F=6Y_2^5$.

[Source and full output](https://github.com/rbajaj5/hodge-lab/blob/45ab8ff61c1af70a30c6cb16ea9a95252d6846d2/experiments/one_reflection_skew_degree7.py) · [summary data](https://github.com/rbajaj5/hodge-lab/blob/45ab8ff61c1af70a30c6cb16ea9a95252d6846d2/results/one-reflection-skew-degree7/summary.json) · [argument and commands](https://github.com/rbajaj5/hodge-lab/blob/45ab8ff61c1af70a30c6cb16ea9a95252d6846d2/notes/one-reflection-skew-degree7.md).

Prepared with AI assistance for `rbajaj5`; no Lean proof is submitted.



## COMMENT DC_kwDOVLOOJc4BIsx8 by fqlx at 2026-10-20T03:46:05Z (updated 2026-10-20T03:46:05Z)

## E2 progress: exact discriminant obstruction and an explicit remaining primitive

Update for the **unchanged `one-reflection-skew` target**, $q=1$. We closed the discriminant rank gap and obtained a generic connection with exact remainders. **We have not found a certificate or a Lean proof.** This work was developed with ChatGPT and parallel mathematical agents; “independent checks” below means separate computational models, not independent human review.

### 1. An exact fixed-denominator obstruction

This concerns the particular degree-24 denominator in `paired-reflection-discriminant-e1.json`, not the later 69-term cubic discriminant. Put $f=F-Y_2^6$ and $B_D(K)=D\mathrm{div}K-K(D)$.

The relevant physical tangent-field space has dimension $18106$. We now have **267 integer polynomial fields**, each checked coefficient by coefficient to satisfy $K(f)=B_D(K)=0$. Their physical coefficient minor has determinant $30\bmod101$, proving independence over characteristic zero. Thus the coefficient rank is at most $18106-267=17839$.

The previously audited augmented $17840$-square minor has determinant $4\bmod101$; its source and checker hashes were rechecked. Together these establish

$$
\mathrm{rank}_{\mathbb C}L_D=17839,
\qquad
\mathrm{rank}_{\mathbb C}[L_D\mid-B_D(N_0)]=17840.
$$

Therefore this fixed denominator admits no E2 correction. The reduction to the five-variable system includes all six original numerator components and general $Q$; Koszul completeness uses the regular gradient sequence. **This does not exclude other denominators or their powers.** A modular inconsistency alone would not have established the characteristic-zero upper bound.

The recovered `sparse29` denominator is also excluded by a directly checked integer dual: it annihilates all 318 correction columns but pairs with the required source to give $268850534400\ne0$.

### 2. A generic rank-eight connection is now explicit

On the normalized CM cover, use $k=n-\beta m$, $\mu=1/k$, and

$$
(\mu,t\mu,\phi_b,t\phi_b,\phi_c,t\phi_c,t^2\phi_c,t^3\phi_c)\,dt,
$$

where $\phi_b=\partial_b\mu$ and $\phi_c=\partial_c\mu$. The connection over $\mathbb Q(u,v)$ includes every rational Hermite remainder. All **1088 original coefficient equations** vanish exactly. A separate eliminated degree-eight model verifies all 16 differential identities, 64 flatness identities, and 768 export bindings.

The audited gauge descends to $\mathbb Q(r,w)$, giving the sufficient remaining equation

$$
\partial_rU+U\widetilde M_r+\partial_wV+V\widetilde M_w=e_2,
\qquad U,V\in\mathbb Q(r,w)^8.
$$

Here $e_2=(0,0,1,0,0,0,0,0)$. This route is sufficient; we have not proved that every possible rational primitive can be restricted to this proper-curve component.

### 3. An explicit two-form reduces the actual source

Set $a=1+r^6+w^6$, $\alpha^3=a$, $u=(r+r^3w)/\alpha^2$, $v=r/(7\alpha)$, and

$$
3\beta+m^2=ut,\quad 3\beta^2+2mn=vt(1-t^2),\quad \beta^3+n^2=1+t^6.
$$

With $\rho=1/[a\det(\partial(u,v)/\partial(r,w))]$ and $\Omega=\rho\phi_b\,du\wedge dv\wedge dt$, the rational form

$$
\Gamma_*=
\frac{\rho n}{tk^2}
\left(dv-\frac{\beta}{1-t^2}du\right)\wedge dt
$$

satisfies the independently checked identity

$$
\Omega-d\Gamma_*=
-\frac{n}{tk^2}
\left(\rho_u+\frac{\beta}{1-t^2}\rho_v\right)
du\wedge dv\wedge dt.
$$

Both density derivatives are essential. Weighted elliptic addition forms supply concrete rational candidates for this remainder, but no matching combination including all coefficient derivatives has been found. The original normalization factor $-1/2$, CM trace, exact-form lift, and polynomial contract still need to be restored and checked after a primitive is obtained.

### 4. One bounded actual-base search is conclusively excluded

For every coefficient of $U,V$, we tested the 26-function span of monomials of degree at most 2; $1/r,w/r$; and $1,r,w$ times $D_r/D,D_w/D$ for $D=a,L_0,S^\sharp$, the recorded denominator pullbacks.

That gives **416 unknown constants**. Eighty regular points modulo 503 produce a $640\times417$ augmented matrix with a $417$-square minor of determinant **354**. Since there are only 416 coefficient columns, this excludes exactly that span over $\mathbb Q$ and $\mathbb C$. An independent Taylor-series checker rebuilt all 266880 entries. Larger numerator spaces, other poles, and higher pole orders remain open.

### Reproducibility and credits

The attachments below contain the complete research packet as four ordinary ZIPs. Extract all four into the same directory and start with `DISCUSSION-PACKET-README.md` and `E2-ROUND4.md`. They include the exact denominators, witnesses, generic matrices, replay scripts, source hashes, and full limitations.

Eleven selected replay steps passed in a clean extraction without network access. No repository file or protected target changed. Credits: Benjamin E. Diamond and the project contributors for the framework; `fringewidth` and `fqlx` for the unskewed hidden-reflection work motivating this investigation; Griffiths's residue theory and the classical Koszul/trace constructions used in the reductions.

This complements the bounded-degree work already posted here by `rbajaj5`; it makes no new all-denominator or all-degree nonexistence claim.

### Attachments

[one-reflection-skew-e2-round4-part-1-of-4.zip](https://github.com/user-attachments/files/33698467/one-reflection-skew-e2-round4-part-1-of-4.zip)
[one-reflection-skew-e2-round4-part-2-of-4.zip](https://github.com/user-attachments/files/33698469/one-reflection-skew-e2-round4-part-2-of-4.zip)
[one-reflection-skew-e2-round4-part-3-of-4.zip](https://github.com/user-attachments/files/33698470/one-reflection-skew-e2-round4-part-3-of-4.zip)
[one-reflection-skew-e2-round4-part-4-of-4.zip](https://github.com/user-attachments/files/33698473/one-reflection-skew-e2-round4-part-4-of-4.zip)



## COMMENT DC_kwDOVLOOJc4BIx3J by fqlx at 2026-10-20T17:59:22Z (updated 2026-10-20T17:59:22Z)

## Curve/cover audit: a valid line, but the generic-endpoint and source-transport gaps remain

**We have not obtained a certificate or a Lean proof for the unchanged one-reflection-skew target.** This update concerns the attempted transfer of the moving-curve machinery from PR #150.

The proposed “find a rational curve with an invertible deformation Jacobian, then port the rest” plan needs a stronger specification: the relevant Jacobian must support independently varying generic endpoints, on auxiliary geometry whose certificate transports to the exact requested source. An arbitrary nonsingular curve-incidence point does not establish this.

### 1. Exact line and two different Jacobian tests

For the [protected equation](https://github.com/hodge-computer/sextics/blob/main/data/one-reflection-skew/instance.rec), in coordinate order $(X_0,X_1,X_2,Y_0,Y_1,Y_2)$, the line

$$
\ell(s,t)=[s:t:6s+2t:3s+2t:2s:7s+8t]
$$

lies on $F_{\mathrm{skew}}=0$ over $\mathbb F_{11}$. All seven coefficients of $F_{\mathrm{skew}}(\ell(1,t))$ vanish.

These calculations use **$F_{\mathrm{skew}}$ itself**, with $1/7=8$ modulo $11$, rather than $7F_{\mathrm{skew}}$. Write $\ell(1,t)=a+tb$ and $g_i(t)=(\partial_i F_{\mathrm{skew}})(a+tb)$.

- Fixing $X_0=s,\ X_1=t$ and varying the constant and linear coefficients of $X_2,Y_0,Y_1,Y_2$ gives a $7\times8$ matrix $M$ of rank $7$. Its first seven columns have determinant $6$.
- Fixing the entire endpoint $a$ and varying the six slope coefficients gives the $6\times6$ matrix $J$, with $J_{r,i}=[t^{r-1}]g_i(t)$ for $r=1,\ldots,6$. It has rank $4$, not $6$.

For reproduction, rows of $M$ correspond to $t^0,\ldots,t^6$; columns alternate constant/slope variations in the coordinate order $X_2,Y_0,Y_1,Y_2$. Columns of $J$ use all six original coordinates:

~~~text
M =
5 0 6 0 5 0 5 0
1 5 5 6 8 5 5 5
8 1 5 5 0 8 2 5
5 8 3 5 9 0 7 2
2 5 3 3 5 9 4 7
5 2 5 3 4 5 5 4
0 5 0 5 0 4 0 5

J =
6 4 5 6 5 5
5 3 1 5 8 5
0 2 8 5 0 2
4 5 5 3 9 7
7 2 2 3 5 4
0 6 5 5 4 5
~~~

Thus this is a valid line and a smooth ordinary line-incidence point, but it does not supply the stronger fixed-endpoint seed datum used in PR #150.

### 2. Why the auxiliary geometry and target identity matter

[PR #150's construction](https://github.com/hodge-computer/sextics/blob/main/data/reflection-coupled-fano/solution.md) places its weighted quartic on the quotient $u=z^2$, with $z=Y_2$. Its source transport uses

$$
A=z B(x,z^2),\qquad B=x_0^3u,\qquad A=x_0^3z^3.
$$

Here the requested source is $A=X_0^3Y_0^3$, which is even in $Y_2$. It cannot be expressed in that form with $B$ defined on the same quotient. The surviving $Y_2$ reflection therefore does not provide the required source transport.

The [shared PolynomialCurves interface](https://github.com/hodge-computer/sextics/blob/main/formal/lean/HodgeComputer/Constructions/PolynomialCurves/README.md) also retains generic endpoint independence, a nonzero transverse coordinate, polynomial-normal and determinant conditions, and application-specific transport and target identities.

Under the smoothness stated in the [instance context](https://github.com/hodge-computer/sextics/blob/main/data/one-reflection-skew/why.md), the original projective sextic $X\subset\mathbb P^5$ has $K_X\simeq\mathcal O_X$. In characteristic zero it cannot admit a dominating family of nonconstant rational curves. This follows from adjunction and the nef-canonical obstruction to separable uniruledness; see [Debarre, Introduction to Mori Theory](https://perso.imj-prg.fr/wp-content/uploads/debarre-pub/M22.pdf), Example 1.5 and Example 9.14.

This obstruction concerns the original smooth projective variety. It does not exclude individual lines, radial rulings of its affine cone, curves on auxiliary quotients, or certificates by other methods.

### 3. The old exchange and a scaled repair fail exactly

For the usual cover

$$
T^6=Y_0^6+Y_2^6,\qquad Y_0=\alpha T,\qquad y=Y_1,
$$

write $G(y,T)$ for the covered skew equation, with $X$ and the Fermat-base coefficients fixed. Direct subtraction gives

$$
G(T,y)-G(y,T)
=\frac{\alpha X_0X_1}{7}(\alpha^2+1)yT(y^2-T^2).
$$

This is generically nonzero. The tested scaled exchange $y\mapsto cT,\ T\mapsto c^{-1}y$ requires both $c^6=1$ and $c^2=-\alpha^2$, hence $\alpha^6=-1$. It works, if at all, only on a special base fibre rather than restoring the generic exchange.

**Remaining gap:** construct a cover, quotient, or other correspondence that serves the exact skew target and supplies the required generic-endpoint moving family. The lifting/integration/trace machinery is reusable once its hypotheses and transport identities are established; that geometric step is still mathematical work. The checks above do not exclude other covers or unrestricted certificates.

Prepared with ChatGPT and parallel AI-agent assistance for `fqlx`; the finite-field matrices were checked by exact arithmetic. Builds on the project framework of Benjamin E. Diamond and contributors, the reflection/cover work of `fringewidth` and `fqlx`, and the shared constructions extracted from PRs #150–151. No Lean verification of these new calculations is claimed.



## COMMENT DC_kwDOVLOOJc4BIyQ8 by fqlx at 2026-10-20T19:35:47Z (updated 2026-10-21T12:53:37Z)

## Serving twisted elliptic companions for the exact skew source

**Status:** the attached packet establishes an exact serving correspondence and an equivalent, explicit rational primitive problem for the unchanged `one-reflection-skew` instance. The primitive problem remains unsolved. It supplies no certificate $(D,N_i,Q)$ and no Lean validation.

[Research packet, exact checkers, outputs, and provenance](https://github.com/user-attachments/files/33736319/one_reflection_skew_twisted_companions.zip)

This concerns [the protected skew instance](https://github.com/hodge-computer/sextics/blob/main/data/one-reflection-skew/instance.rec), with $A=X_0^3Y_0^3$ and $q=1$, not the central holdout or holdout-bundle. The latest source comment checked was [the curve audit, 19078601](https://github.com/hodge-computer/sextics/discussions/127#discussioncomment-19078601). No conclusion below uses its line as an endpoint certificate.

### 1. Physical quotient and exact source

On $Y_0=1$, put

$$
X_0=x,\quad X_1=rx,\quad X_2=wx,\quad Y_1=t,\quad Y_2=z,\quad U=x^2,\qquad \nu=dr\wedge dw\wedge dt.
$$

The unchanged equation is $G=z^6+P(U)$, where

$$
P=aU^3+bU^2+cU+d,\qquad
(a,b,c,d)=\left(1+r^6+w^6,\ t(r+r^3w),\ \frac{rt(1-t^2)}7,\ 1+t^6\right).
$$

Thus the skew term is retained, including its sign. Simultaneous negation of the three $X_i$ changes the signs of both $A$ and the projective volume, so their product descends through $U=x^2$. The full density is

$$
\frac{A\Omega}{F^2}=-\frac{U^2}{2G^2}\,dU\wedge\nu\wedge dz.
$$

With residue convention $\mathrm{Res}(dG/G\wedge\eta)=\eta$, exact pole reduction gives

$$
\theta=\frac{5U^2}{72z^{11}}\,dU\wedge\nu.
$$

Let $\delta=\mathrm{disc}_U(P)$, $\omega_0=dU/z^5$, and $\omega_1=U\omega_0$. An exact Hermite reduction is

$$
\theta=\theta_0+d(S\nu),\qquad
\theta_0=\frac{H_0\omega_0+H_1\omega_1}{72\delta}\wedge\nu,
$$

where

$$
H_0=12b^2d-bc^2-27acd,\qquad H_1=6a(c^2-3bd),
$$

$$
S=\frac{6d(bc-9ad)+6(bc^2-2b^2d-3acd)U+12a(c^2-3bd)U^2}{72\delta z^5}.
$$

This identity is absolute: the additional physical parameter derivatives wedge to zero against the full base volume $\nu$.

The generic smooth projective normalization of $C$ has genus four. The two displayed holomorphic forms, and the source, have eigenvalue $\zeta_6$ under $z\mapsto\zeta_6z$. In particular, the surviving $z\mapsto-z$ reflection still does not serve the source. The correspondence below addresses that character directly; it does not reinterpret the surviving reflection as a solution.

### 2. Quartic cover and the twist

Over $K=\mathbf C(r,w,t)$, use the finite étale algebra $B=K[\beta]/(\mathcal R)$ on its generic open, without presuming irreducibility, where

$$
\mathcal R=3a^2\beta^4-4ab\beta^3+6ac\beta^2-12ad\beta+4bd-c^2.
$$

Set

$$
M=b-3a\beta,\qquad N=(c-3a\beta^2)/2,\qquad k=N-M\beta,\qquad L=\beta+3k/M.
$$

The polynomial identities

$$
\mathcal R=4[M(d-a\beta^3)-N^2],\qquad M\mathcal R'+12ak^2=-3a\mathcal R
$$

give $P=a(U+\beta)^3+(MU+N)^2/M$. Hence there is a degree-three map

$$
\phi_\beta:C_B\longrightarrow E_\beta,\qquad
C:z^6=-P(U),\qquad E_\beta:\upsilon^2+M(1+a\xi^3)=0,
$$

$$
\xi=(U+\beta)/z^2,\qquad \upsilon=(MU+N)/z^3.
$$

Its inverse satisfies $\upsilon z^3-M\xi z^2-k=0$. The degree is exactly three: $\xi$ has degree six on $C$ and degree two on $E_\beta$.

The twist matters. Introducing $m^2=M$, $mn=N$ produces opposite elliptic differentials for $(m,n)$ and $(-m,-n)$ on the common untwisted companion; direct source trace there cancels. Keeping $M$ in the elliptic equation absorbs this ambiguity over the quartic algebra itself. No exchange symmetry of the skew sextic is asserted.

### 3. Exact serving action and pushed source

Write $\omega_E=d\xi/(2\upsilon)$. Direct differentiation and cubic field traces give

$$
\phi_\beta^*\omega_E=-\frac{U+L}{6}\omega_0,\qquad
\phi_{\beta,*}\omega_0=-\frac{6M}{k}\omega_E,\qquad
\phi_{\beta,*}\omega_1=\frac{6\beta M}{k}\omega_E.
$$

All traces here are **raw, unnormalized traces**. On $(\omega_0,\omega_1)$,

$$
T_\beta=\phi_\beta^*\phi_{\beta,*}
=\frac Mk\begin{pmatrix}L&-\beta L\\1&-\beta\end{pmatrix}.
$$

Quartic root traces satisfy

$$
\sum_\beta M/k=0,\quad \sum_\beta\beta M/k=-6,\quad
\sum_\beta\beta^2M/k=-4b/a,\quad \sum_\beta\beta=4b/(3a),
$$

and therefore

$$
\boxed{\mathrm{Tr}_{B/K}T_\beta=6I_2.}
$$

The trace calculations require no numerical roots. The inverse cubic gives $z h'(z)=M\xi z^2+3k$ at a root; Lagrange interpolation computes its differential traces by polynomial remainders. The identity $M\mathcal R'=-12ak^2$ then reduces the quartic traces to the same elementary calculation. Both implementations check the resulting matrix entries exactly.

The exact pushed source is

$$
\Theta_\beta=\phi_{\beta,*}\theta_0=C_\beta\omega_E\wedge\nu,
\qquad
C_\beta=\frac{M(\beta H_1-H_0)}{12\delta k}
=\frac{M^2(3a\beta^2k-2N^2)}{12ak^4}.
$$

Consequently,

$$
\boxed{\theta=d(S\nu)+\frac16\mathrm{Tr}_{B/K}(\phi_\beta^*\Theta_\beta).}
$$

These are generic rational identities, not identities restricted to a special physical fibre. The packet checks the discriminants and resultants guaranteeing a nonempty legal open; at $(r,w,t)=(1,1,2)$, $\delta=-361978181/343$ and $3ac-b^2=-166/7$.

### 4. The full remaining primitive equation

Work componentwise in $\mathcal L=B(\xi,\upsilon)$ with the elliptic relation above. Define

$$
\widehat\omega=\frac{d\xi+\xi\,da/(3a)}{2\upsilon},\qquad
\gamma=\frac{da}{3a}+\frac{dM}{2M},\qquad
 d\widehat\omega=-\gamma\wedge\widehat\omega.
$$

For $i\in\{r,w,t\}$, extend the physical derivations by

$$
D_i\beta=-\frac{(\partial_i\mathcal R)_{\beta\ \mathrm{fixed}}}{\mathcal R'},\qquad
D_i\xi=-\frac{\xi D_i a}{3a},\qquad
D_i\upsilon=\frac{\upsilon D_iM}{2M},
$$

and put $\gamma_i=D_i a/(3a)+D_iM/(2M)$. The vertical derivation is

$$
\partial_E\xi=2\upsilon,\qquad \partial_E\upsilon=-3Ma\xi^2,\qquad \partial_EB=0.
$$

These derivations preserve the defining equations. Every rational three-form has the form

$$
\Gamma=W\nu-\widehat\omega\wedge\sum_iV_i\,\iota_{\partial_i}\nu,
\qquad W,V_i\in\mathcal L.
$$

Thus the complete rational primitive problem is

$$
\boxed{\partial_EW+\sum_i(D_iV_i-\gamma_iV_i)=C_\beta.}
$$

The unknowns may depend rationally on the elliptic variables and have arbitrary moving poles. A coefficient-only ansatz would be an additional restriction, not part of this equation.

Rational exactness of $\theta$ is equivalent to exactness of $\Theta_\beta$ on all components: push forward in one direction; in the other use

$$
\Psi=S\nu+\frac16\mathrm{Tr}_{B/K}(\phi_\beta^*\Gamma_\beta),\qquad d\Psi=\theta.
$$

Finite separable trace commutes with $d$. Pullback through $U=x^2$ returns the original residue problem. Once a primitive exists, ambient lifting, homogenization, and denominator clearing must still produce and check the requested polynomial identity.

### 5. A shortcut fails at the density, with an explicit correction problem

For $f=M/k$, the universal coefficient derivative holding $a,c,d$ fixed gives

$$
C_\beta=-\tfrac12\left[D_bf-\frac{D_bM}{2M}f\right].
$$

However, this universal direction is not tangent to the physical coefficient image. Correcting it by translation in $U$ determines a rational physical field $T$ with $T(P)=U^2+eP'$. At $(1,1,2)$,

$$
T=(12533/87306,-12533/87306,1/14551),\qquad e=-224/14551,
$$

$$
\mathrm{div}T=2209940204/1905584409\ne0.
$$

The candidate $V_i=-fT_i/2$ therefore leaves the exact residual $-f\mathrm{div}T/2$. Repairing tangency alone does not repair the full density. This evaluates a rational obstruction at an ordinary physical point; it disproves that particular shortcut, not the existence of a different primitive.

The packet continues the correction explicitly. Write $J=\partial(a,b,c,d)/\partial(r,w,t)$, $\tau=(0,3a,2b,c)^T$, and $\mu=\det[J\mid-\tau]$. Translation-corrected physical fields $\mathcal W,\mathcal Z$ are uniquely specified by

$$
[J\mid-\tau]\binom{\mathcal W}{\lambda_{\mathcal W}}=(3a,2b,c,0)^T,\qquad
[J\mid-\tau]\binom{\mathcal Z}{\lambda_{\mathcal Z}}=(a,b,c,d)^T.
$$

They commute. The resulting smaller sufficient equation is

$$
(\mathcal W+6)h+(\mathcal Z+19/6)k_0=-\frac{\mathrm{div}T}{12\mu},
\qquad h,k_0\in\mathbf C(r,w,t).
$$

All coefficients and cleared polynomial checks are in `rank2/gauge_scalar_repair.json` and its checker. Here $k_0$ is an unknown function, distinct from the quartic-cover quantity $k$. The degrees of both unknowns remain unrestricted. This equation also remains unsolved; it is a concrete corrected alternative to the failed tangent lift, not a claim that every primitive must have this form.

### Reproduction and limits

Run `python3 verify_twisted_companion.py` for the standard-library checker. The separate `check_serving_twist.py`, rank-two connection checker, and gauge-repair checker require SymPy; recorded exact outputs and hashes are attached. These are separate computational implementations by ChatGPT agents, not independent human review.

No Fano, uniruledness, or rational-connectedness property is established for the corrected four-dimensional auxiliary space. No dominating rational curves on the original Calabi–Yau are used. There is no certificate, Lean proof, or unrestricted nonexistence claim. Nothing restricts the eventual denominator degree, complex coefficients, general quotient $Q$, or denominator characters. The displayed discriminant is a reduction denominator, not a proposed complete certificate denominator.

Credits: Benjamin E. Diamond and the [hodge.computer](https://hodge.computer) contributors for the problem and framework; the successful unskewed construction credits `fringewidth` and `fqlx`. This research and the attached algebraic checks were developed with ChatGPT assistance. The contribution reported here is the checked correspondence, exact source transport, and remaining primitive equations; solving them remains the next obligation.


---

## Additional context: exact degree-seven mixed-character reductions (20 October 2026)

The degree-seven reproducibility packet described below is separate from the twisted-companion packet linked above and is not currently attached to this comment.

This addendum concerns the same protected **one-reflection-skew** target, with $A=X_0^3Y_0^3$ and $q=1$. It records bounded degree-seven results and a remaining finite incidence problem. **No certificate, unrestricted nonexistence proof, or Lean validation is claimed.**

The linked [degree-seven packet by `rbajaj5`](https://github.com/rbajaj5/hodge-lab/tree/45ab8ff61c1af70a30c6cb16ea9a95252d6846d2) was replayed over exact rational arithmetic. Every published sector dimension, basis label, and equation rank was reproduced; all 129 gauged E1 lifts were checked coefficient by coefficient against the protected target, and the 84-dimensional base span was independently reconstructed. All nine extra non-Jacobian directions were retained. The earlier separate pure-character exclusions were not used as mixed-character exclusions.

### 1. Complete correction space and a genuine mixed cancellation

Write $x=(X_0,X_1,X_2,Y_0,Y_1)$, $t=Y_2$, $F=f(x)+t^6$, $J_f=(f_0,\ldots,f_4)$, and $B_D(N)=D\mathrm{div}N-N(D)$. Euler normalization gives $N(F)=AD$ and the full remaining equation $B_D(N)=FQ$, with all six numerator components and general $Q$ retained. Every normalized lift is $N=U+K\nabla F$ with all 840 cubic skew-matrix coefficients initially present. The 120-dimensional subspace $K=E\wedge V$, for quadratic $V$, is a harmless gauge because

$$B_D(K\nabla F)=-6F\bigl(D\mathrm{div}V-V(D)\bigr).$$

After the six linear $F$-gauges, the exhaustive E1 form is

$$D=h+t\,l\cdot\nabla f+t^2c\cdot\nabla f+6t^5R_2,$$

where $h\in(J_f:A)_7$ has dimension 84, $l$ is an arbitrary linear five-vector, $c$ is constant, and $R_2$ is quadratic. This accounts for all $129=84+25+5+15$ parameters.

The top necessary equation has the complete complex solution

$$c^T\mathrm{Hess}(f)c+30R_2^2=0
\quad\Longleftrightarrow\quad
c=\lambda e_{X_2},\qquad R_2=\epsilon i\lambda X_2^2,\qquad \epsilon=\pm1.$$

Thus an actual even/odd cancellation occurs. Both nonzero branches were retained. Further exact equations and the specialization argument below force $h=fL$ and $l=\beta x$ on these branches. An allowed $F$-gauge gives $D'=t^2H$, with $H(x,0)=f_{X_2}$, while every corrected lift satisfies

$$N_t(x,0)\equiv\epsilon iAX_2^2\pmod{J_f}.$$

The constant and linear coefficients of $B_{D'}(N)=FQ$ give $Q_0=0$ and

$$-2f_{X_2}N_t(x,0)=fQ_1.$$

Since $f$ is smooth and irreducible, $\gcd(f,f_{X_2})=1$, so $f\mid N_t(x,0)$. This would put $AX_2^2$ in $J_f$. But the coefficient of $X_0^3X_2^2Y_0^3$ vanishes on every cubic multiple of every $f_i$: no monomial in any $f_i$ divides it. Its coefficient in $AX_2^2$ is one. This closes both nonzero branches with full corrections and general quotient.

### 2. Joint odd and further even loci

The remaining $t^2$ equations force $l_i=a_ix_i$, independently of $h$: exact quadratic row spaces of ranks 144 and 15 isolate all 20 off-diagonal squares. For $h=0$, the full compatible correction calculation yields nine necessary quadrics whose cubic Macaulay rank is $35/35$. Hence **every degree-seven $Y_2$-odd denominator is excluded, including mixtures of all four odd characters**.

A separate local argument excludes the full 30-parameter family $D=fL+t\,l\cdot\nabla f$, with arbitrary linear $l$: after gauging to $t(l\cdot\nabla f)-Lt^6$, its constant and $t^6$ equations force $AL^2\in J_f$. Exact multiplication by $A$ on base quadratics has rank $15/15$, hence $L=0$, and the joint odd exclusion applies.

For $V=(J_f:A)_7/fS_1$, of dimension 79, three normalized character pieces also have empty necessary projective loci:

| Character of $[h]$ | Parameters | Quadrics | Completed Gröbner basis size |
|---|---:|---:|---:|
| 0010 | 18 | 44 | 967 |
| 0100 | 17 | 44 | 620 |
| 1000 | 18 | 44 | 968 |

The modular calculation is over $\mathbb F_{32003}$ and proves that every coordinate has a power of degree at most five in its ideal. The exact rational source, integrality at 32003, and source ranks were independently checked. The degree-14 Jacobian quotient has dimension 185 over both $\mathbb Q$ and $\mathbb F_{32003}$, so its localized integral module is free. The resulting projective scheme has empty special fiber; properness excludes a nonempty generic fiber. This is a checked characteristic-zero implication, not a finite-field point scan. The analogous 25-parameter calculation used above has Hilbert function $(1,25,195,361,240,51,2)$.

These exclusions allow arbitrary retained $fL$ terms and the permitted odd component. **They do not exclude mixtures of multiple even base characters.**

### 3. Complete remaining finite problem and evidence

Every remaining degree-seven certificate has a representative

$$\boxed{D=h+t\sum_{i=0}^4a_ix_if_i,\qquad h\in(J_f:A)_7,\quad h\notin fS_1.}$$

Retain 79 base-complement parameters, all five $fS_1$ parameters, five odd parameters, and all 720 effective Koszul parameters. The 79 charts normalizing one nonzero base-complement coordinate to one cover the remaining locus. The complete coefficient equations are

$$\mathrm{rem}_{t^6+f}\bigl(D\mathrm{div}N-N(D)\bigr)=0.$$

The completed exact export has **9,435 nonzero coefficient equations**, 68,085 quadratic columns, and 2,274,017 integer coefficient entries. These equations are not asserted to be independent. At a solution, $Q=B_D(N)/F$ recovers the full original cleared identity; $Q$ need not vanish. Since every generated nonzero $D$ has $t$-degree at most one and irreducible $F$ has $t$-degree six, complex coprimality is automatic on these charts.

The packet includes exact bases, source hashes, proof notes, equation blocks, and replay scripts. Twelve selected checks passed from a clean extraction without network access. All incidence-block hashes, counts, basis rescalings, and Koszul labels were checked; 16 selected columns were independently rebuilt with separate tuple/Fraction arithmetic, including nontrivial reduction of $t^8$ and $t^9$. No second independent recomputation of every exported column is claimed.

**Remaining:** pure normalized base character 1110 and mixtures of at least two even characters. The 1110 Gröbner computation was interrupted without a completed result; the joint degree-three attempt also timed out. Neither is an exclusion. The full streamed system remains unsolved.

Prepared for `fqlx` with ChatGPT and parallel AI-agent assistance; independent checks here mean separate computational implementations and algebra audits, not independent human review. Credits to Benjamin E. Diamond and the [hodge.computer](https://hodge.computer) contributors for the framework, and `rbajaj5` for the published degree-seven bases and calculations.




### REPLY DC_kwDOVLOOJc4BI5PZ by fqlx at 2026-10-21T20:51:48Z

**Exact source/kernel update — E2 remains open.**

For the unchanged $q=1$, $A=X_0^3Y_0^3$ target, the compact chart reduction now has an exact rational source section and a complete source kernel. [Code, exact data, hashes, and standalone replay](https://gist.github.com/fqlx/c1d337b1c8d7edcb73ad2001876e9398).

- **Source:** both original physical source components check exactly over the rational function field. The reduced section $Z_0$ has six physical component degrees $(1,2,1,1,2,1)$ in $w$. Its ordinary divergence is nonzero.
- **Full kernel:** over $\mathbb Q(r,t)[w]$, four generators have degrees $(3,5,5,5)$ and independent leading columns. After localization, every rational compact source section is $Z_0+Kf$, with four arbitrary rational functions $f(r,t,w)$. The remaining equation is the full two-component closure condition $\operatorname{div}(Z_0+Kf)=0$, retaining all physical parameter derivatives.
- **Scoped exclusion:** the complete physical polynomial-$w$ compact-flux family of degree at most five requires six coefficient functions, corresponding to $H_0,wH_0,w^2H_0,H_1,H_2,H_3$. Its first-derivative matrix has rank 12; curvature has rank six and forces a unique rational candidate. At the regular exact point $(r,t)=(1,2)$, all six $r$- and all six $t$-derivative identities for that candidate fail. This excludes that class with arbitrary coefficients in $\mathbb C(r,t)$, without an auxiliary-degree restriction.

That last degree is a **chart-flux degree**, not a bound on homogeneous $\deg D$. Higher $w$ degrees, rational $w$ poles, additional fibre poles, and general certificate constructions remain open. The next concrete problem is unrestricted rational closure in the four-function kernel description.

The selected public packet freshly replayed all ten verification producers successfully in 21.06 seconds, using exact arithmetic. “Independent” here means separate computational implementations, not human review. There is **still no E2 primitive, legal homogeneous $(D,N,Q)$, or instance `HodgeComputer.HasCertificate` theorem**, and no protected Lean verification is claimed.

Prepared for `fqlx` with OpenAI Codex/ChatGPT and parallel AI-agent assistance; continuing the source/companion work above and building on the project framework and the prior reflection/cover work of `fringewidth`, `fqlx`, and other contributors to this discussion.



### REPLY DC_kwDOVLOOJc4BI51C by rbajaj5 at 2026-10-21T23:35:55Z

**A three-coefficient certificate for one bounded unreduced primitive space.** Building on `fqlx`'s source and rank-two formulation above, checkpoint 118 of Hodge Lab has a small exact obstruction that can be checked without a large rank calculation. This is an independent certificate and comparison of bounds; **no increased coverage over the existing companion screens is claimed**.

In the same `Y0=1` chart, write `a=1+r^6+w^6`, `P=aU^3+bU^2+cU+d`, and $P_7=7P$. For all four components of `V=N/a`, the full source equation becomes

$$
\mathcal M(N):=60a\nabla P_7\cdot N-72P_7a\operatorname{div}N+72P_7N\cdot\nabla a
=35U^2a^2.
$$

Suppose each numerator has total degree at most four in `(U,r,w,t)`. Set `n0=[1]N_U`, `n2=[t^2]N_U`. Coefficient support gives

$$
[rt^3]\mathcal M=60(n_2-n_0),\qquad
[U^2w^{12}]\mathcal M=1260n_0,\qquad
[U^2w^{12}t^2]\mathcal M=1260n_2.
$$

The source coefficients are `(0,35,0)`, forcing `n2=n0`, `n0=1/36`, and `n2=0`. All base-component contributions to these three rows vanish under this bound. A separate exact `Fraction` implementation checked the coefficient transpose and all 280 monomial columns; no historical companion rank was rerun.

**Overlap with the prior screens:** if each horizontal numerator is `A_i+U B_i`, its rank-two coefficients are directly `(A_i/a,B_i/a)`, already inside `fqlx`'s reported degree-six `H/a` space. This accounts for 165 of the 210 horizontal coefficients. The remaining 45 higher-`U` terms require care: for instance,

$$
\frac{U^2}{a}\omega_0
=-\frac{c}{3a^2}\omega_0-\frac{2b}{3a^2}\omega_1
+d_U\left(-\frac{2z}{a^2}\right).
$$

Thus the natural reduction changes the denominator; this does not establish new coverage or rule out a different gauge representative. Base derivatives of the relative exact correction must be retained in an absolute reconstruction. Five such relative reductions were independently checked.

The degree bound is essential to this tiny certificate: at degree five, `N=(1/36,-5r^2t^3/3024,0,0)` satisfies these three rows, but its full residual has `[rt]=5/3`, so it is not a primitive. No unrestricted nonexistence or homogeneous certificate-degree bound follows.

I also read the [new compact source/kernel update](https://github.com/hodge-computer/sextics/discussions/127#discussioncomment-19108825). Its polynomial-`w` compact-flux exclusion is a different bounded presentation; I have not replayed that packet or established a coverage comparison with it.

[Certificate, comparison, reproducible checks and limits](https://github.com/rbajaj5/hodge-lab/blob/c3ab5214fcff79309ad1e7ae83b9e436afeb1d39/notes/skew-small-obstruction-audit.md). A comparison not executed in this audit is the 45 higher-`U` covariant-divergence images against the prior 504-column image, rather than identifying the two degree bounds by inspection.

Credit `fqlx` for the quotient, residue and companion formulation, and Benjamin E. Diamond and contributors for the framework. Prepared for `rbajaj5` with ChatGPT and parallel AI-agent assistance; independent here means separate computational implementations, not human review.


**Lean verification — 21 October 2026.** The same bounded obstruction now has a [Lean polynomial-to-matrix proof](https://github.com/rbajaj5/hodge-lab/blob/b7ec8648acf95fb2ce86f11d3b207b702d00d76a/notes/skew-lean.md): it constructs the full chart operator, checks the 840 selected coefficients of all 280 columns and the three source coefficients, and rules out the encoded polynomial equation over every characteristic-zero field. The [local verification receipt](https://github.com/rbajaj5/hodge-lab/blob/b7ec8648acf95fb2ce86f11d3b207b702d00d76a/formal/results/skew-lean/validation.json) records compilation, 17 declaration audits, four default kernel module replays, and two rejected mutations. This strengthens verification of the existing exclusion; it does not add exclusion coverage. Geometric source transport, unrestricted rational primitives and the 45-versus-504 operator-image comparison remain outside these Lean theorems; imported mathlib is not freshly replayed.



### REPLY DC_kwDOVLOOJc4BI6dm by fqlx at 2026-10-22T03:08:57Z

**E2 update: an exact physical correction and a complete rank-17 reformulation; the primitive remains missing.**

For the unchanged `one-reflection-skew` target, $q=1$ and $A=X_0^3Y_0^3$, we have made progress beyond the [source/kernel update](https://github.com/hodge-computer/sextics/discussions/127#discussioncomment-19108825). **We still have no complete E2 primitive, certificate $(D,N_i,Q)$, or instance theorem `HodgeComputer.HasCertificate`.**

[Reproduction packet: formulas, checkers, recorded outputs, and hashes](https://gist.github.com/fqlx/1d83052f07e2898ef493a52c86a2180b)

Archive SHA256: `d067112d6a859e2fb83d8139eb4de9fe410982ab593d1c77fee77602e72df449`.

1. **A rational correction on the actual physical chart.** With the same $(a,b,c,d)$ as the parent comment, set $u=3db-c^2$, $j=27d^2a-9dcb+2c^3$, and $s=t^2$, so $d=1+s^3$. An exact reduction of the full two-row compact source $S$ gives

   ```math
   c^{-3}S=D_cH+C,
   ```

   with explicit rational vectors $H,C$. Integrating the meromorphic exact part on the genus-four curve $s^3=D^6-1$, where $D^6=d$, produces a rational two-form $\Xi$ on the original physical chart. Its complete identity is

   ```math
   S\,dr\wedge dw\wedge dt=\nabla_M\Xi+R\,du\wedge dc\wedge dd,
   ```


   ```math
   R=\mathfrak r C-B_0J_cD_jH-fJ_dD_jC.
   ```

   Here $\nabla_M=\mathrm d+M$ is the full compact connection, $J_c=\partial_cj_{\rm phys}$ and $J_d=\partial_dj_{\rm phys}$, including the motion of $s$; $D_c,D_j$ are its coefficient-direction covariant derivatives. All coefficients are given in the packet. Both source rows and every transverse derivative are retained. The removed source contribution and remaining density are both nonzero at the exact ordinary point $(r,w,t)=(1,1,2)$. This is a global rational chart identity, improving on a slice calculation, but it does not cancel the full source.

2. **The full compact problem has an explicit rank-17 formulation.** On the stated finite algebraic base, the actual scalar source is

   ```math
   \omega=K\frac{\delta^{11}}{\partial_\delta G}\,dL\wedge dD\wedge dz,
   ```

   where $G$ is the exact physical degree-18 polynomial. Its ordinary trace is zero. Retaining the entire trace-free part gives the basis $e_i=\delta^i/(\partial_\delta G)$, $0\le i\le16$, and the remaining equation

   ```math
   (\partial_L+A_L)y_L+(\partial_D+A_D)y_D+(\partial_z+A_z)y_z=K\mathbf e_{11}.
   ```

   The vector $\mathbf e_{11}$ is indexed from zero. The three unknown vectors have arbitrary rational coefficients on that base: there is no numerator-degree or pole-order cutoff here. The full finite trace back to the original two-row compact system has been audited, including all components and normalization $1/864$. **This is an equivalence for the compact formulation only**, not a necessity claim for every certificate with extra fibre poles or general $Q$. Trace zero does not make the source exact.

The physical-current checker passed **22 exact checks**. The rank-17 audit passed **16 exact checks**, with generic algebraic identities separated from implementation tests of all 102 derivative columns and six curvature matrices at two finite-field points. Those modular tests do not prove characteristic-zero rational solvability. No Lean certificate validation was performed.

The sink/scale-cascade idea motivated rational pole corrections; it has not supplied a Navier–Stokes theorem or an E2 solution. The concrete remaining task is to integrate the displayed nonzero source/remainder, then transport and verify a legal homogeneous certificate. Research is currently paused.

Prepared for `fqlx` with ChatGPT and parallel AI-agent assistance. Independent checks here mean separate computational implementations and algebra audits, not independent human review. Credit Benjamin E. Diamond and the hodge.computer contributors for the problem and framework, and the earlier contributors in this thread for the source and companion calculations.



### REPLY DC_kwDOVLOOJc4BI-aB by benediamond at 2026-10-22T18:15:28Z

**A simpler physical cover and an exact reduction to a one-row source; no solve yet.** Building on `fqlx`'s [complete compact connection and rank-17 reformulation](https://github.com/hodge-computer/sextics/discussions/127#discussioncomment-19113830), we found the following alternative coordinates and checked them against the pinned source.

On the canonical $Y_0=1$ chart, put $V=rU$ and $v=r^{-1}+rw$. The cubic coefficients become

$$
a=A(r,v)=r^3+r^{-3}+(vr-1)^6/r^{15},\quad
b=tv,\quad c=t(1-t^2)/7,\quad d=1+t^6.
$$

Thus the **full degree-18 physical cover** has the particularly small equation

$$
f(r)=r^{18}+r^{12}+(vr-1)^6-ar^{15}=0.
$$

Its covering equation is independent of $t$; the period connection retains its full $t$-derivative. After the fibre gauge and physical Jacobian, the raw source is exactly

$$
B=\frac{r^{11}}{f_r}\frac{M_b e_0}{12} da\wedge dv\wedge dt.
$$

Here $e_0,e_1$ denote the basis associated to $(dV/\nu^5,V dV/\nu^5)$. All eighteen sheets remain present. We use the packet's global residue convention; the README records the overall sign comparison with standard projective volume.

The covering derivatives remain $r_a=r^{15}/f_r$, $r_v=-6r(vr-1)^5/f_r$, and $r_t=0$. The archive verifies both complete connection gauges and the raw source transport.

**An explicit current.** Write $\psi=r^{-3}$, with $a,v$ derivatives taken in the new base. The full connection satisfies $N_a=M_a$, $N_v=tM_b$, and $M_a e_0=M_b e_1$. Set

$$
X_a=\frac{\psi_v e_0}{36t},\qquad
X_v=-\frac{\psi_a e_0}{36t}-\frac{\psi_v e_1}{36t^2},\qquad X_t=0.
$$

Then $\Xi=X_a dv\wedge dt+X_v dt\wedge da$ satisfies

$$
B=\nabla\Xi+\frac{\psi_{vv}e_1}{36t^2} da\wedge dv\wedge dt.
$$

This subtracts an actual rational current and retains the entire nonzero remainder.

**A further reduction and complete one-direction tests.** In $(r,v,t)$ coordinates the remaining scalar density has an ordinary rational primitive

$$
H=\frac{\psi_rA_v^2}{A_r}+K,\qquad
K=90\sum_{k=0}^{4}\binom4k\frac{(-1)^{4-k}v^kr^{k-16}}{k-16},
\qquad H_r=A_r\psi_{vv}.
$$

Subtracting $He_1/(36t^2)$ with the **full twisted operator** leaves $-LM_a e_1/(36t^2)$, where $L=A_rH$ is Laurent in $r$. All poles coming from zeros of $A_r$ disappear.

This gives complete pole bounds for separate tests of rational functions of $r$. At finite $A\Delta$ poles the residue eigenvalues are $0,-1/3$; regular points cannot support additional poles. Bounds at zero and infinity leave powers $-10,\ldots,2$ for the original source, and $-16,\ldots,2$ after the further subtraction. At $(v,t)=(2,3)$, modulo $32003$, coefficient/augmented ranks are respectively **26/27** and **38/39**. Full column rank makes both generic exclusions rigorous: a hypothetical generic solution would specialize regularly, closing the possible denominator loophole. The second result is an independent test of the new residual, not an inference from the first.

These exclude currents confined to that direction; the multidirectional compact problem remains open. A separate degree-nine $P_i/a_0$ exclusion is documented with its explicit bounded scope.

**Reproduction:** [self-contained ZIP: formulas, pinned inputs, exact checkers, and receipts](https://github.com/user-attachments/files/33835156/one-reflection-skew-followup-2026-10-22.zip). No homogeneous certificate or Lean theorem is claimed. ChatGPT-assisted work for `benediamond`, with independent exact checks by parallel agents; the starting connection and source are credited to `fqlx`.



## COMMENT DC_kwDOVLOOJc4BI_aX by DiamondJM at 2026-10-22T22:47:39Z (updated 2026-10-22T22:47:39Z)

Negative results for **one-reflection-skew**: where the accepted one-reflection constructions stop, prepared with Claude Code assistance on behalf of `DiamondJM`. The protected $(F,q=1,A=X_0^3Y_0^3)$ is unchanged. **No certificate, Lean proof, or solve PR is claimed.** Everything below is either a short proof sketch or modular evidence, as labelled.

Notation: $p=Y_0Y_1$, $t=Y_0^2-Y_1^2$, $s=Y_0^3-Y_1^3$, $S=Y_0^3+Y_1^3$, $T=Y_2^3$, $G_0=X_0^6+X_1^6+X_2^6$, $c=X_0^3X_1+X_1^3X_2$, so $F=G_0+cp+\tfrac17X_0X_1\,pt+Y_0^6+Y_1^6+Y_2^6$.

**1. Symmetries.** The diagonal group fixing $F$ has order 24 modulo scalars (216 for one-reflection). It is generated by $\tau$, $\sigma:(X_1,Y_1)\mapsto(-X_1,-Y_1)$ and $\mu_6$ on $Y_2$. The hidden reflections $Y_0\leftrightarrow\pm Y_1$ negate the skew term, and no diagonal twist of the $X$ block repairs this: it would need $a+b\equiv3$ and $3a+b\equiv0\pmod 6$, which is impossible by parity. Note that $\sigma$ survives, so the split $A=A_++A_-$, $A_\pm=X_0^3(Y_0^3\pm Y_1^3)/2$, of `fqlx`'s one-reflection proof is still available.

**2. The Fermat cover $T^6=Y_0^6+Y_2^6$.** A rescaled exchange $(y,T)\mapsto(\lambda T,y/\lambda)$ preserves the image of the skew term only if $\lambda^2=-\alpha^2$ and $\lambda^6=1$, i.e. $\alpha^6=-1$. That is a constant, not the generic point of the Fermat curve.

**3. No cubic-surface reduction with one active $X$ coordinate (proof sketch).** Suppose $F=G(u,\ell_1,\ell_2,\ell_3)$ with $G$ cubic over the passive ratios $r=X_1/X_0$, $w=X_2/X_0$, where $u=X_0^2$ and the $\ell_i$ are rational of weight 2 in $Y$. Suppose also that the transported density gives $A$ or $A_\pm$, as in `fqlx`'s construction. Separating $X$ and $Y$ weights forces the density to be a constant times $u^2$, so $\det\partial\ell/\partial Y$ must be a constant times $Y_0^3$ (or $Y_0^3\pm Y_1^3$). Also $p\in\operatorname{span}\ell$, $pt\in\operatorname{Sym}^2\ell$ and $\sum Y_i^6\in\operatorname{Sym}^3\ell$. Then:
- If two of the $\ell_i$ are free of $Y_2$, the third is a root of a cubic over them, so its $\mu_6$-orbit has size at most 3. The map then factors through $Y_2^2$ or $Y_2^3$, and the Jacobian has the wrong $Y_2$-character. (Example: $\ell=(p,\,t,\,m+Y_2^6/(t^2+p^2))$ with $m=Y_0^2+Y_1^2$ gives a genuine cubic surface $au^3+bu^2p+\tfrac r7upt+\ell_3(t^2+p^2)$, but its Jacobian is proportional to $Y_2^5$.)
- Otherwise, after shifting by multiples of $p$, $\ell_2\ell_3=P_4:=pt-cp^2$, with $\ell_2,\ell_3$ both involving $Y_2$. Eliminating $\ell_3$ gives $Y_2^6=\Phi(\ell_2)$ with $\Phi$ a Laurent polynomial whose only poles are at $0$ and $\infty$. Since $Y_2\mapsto Y_2^6$ has a single, totally ramified pole, $\Phi$ has one pole. So $\ell_2$ (or $\ell_3$) is a polynomial in $Y_2$ of degree at least 2. Then $\det\partial(p,\ell_2,P_4/\ell_2)/\partial Y=\pm(\partial\ell_2/\partial Y_2)\,J(p,P_4)/\ell_2$, which always depends on $Y_2$. Contradiction.
- The rank-one and degenerate shapes close by the same counting. In particular, `fqlx`'s own shape $(p,\,Y_2^2,\,W/Y_2)$ fails for every auxiliary cubic $W$, even allowing $W\mapsto W+uL$.

The pole step was checked by Gröbner bases (Singular, mod 32003, at $Y_0=2$, $Y_1=3$). A Möbius $\ell_2$ that is not a polynomial in $Y_2$ gives the unit ideal, and so does a degree-1 polynomial $\ell_2$. A degree-2 polynomial $\ell_2$ gives a nonempty variety, every point of which has $\ell_2$ depending only on $Y_2^2$, as predicted.

This does **not** exclude reductions with two active $X$ coordinates, or passive ratios mixing $X$ and $Y$.

**4. The fibration over $(r,w,x=Y_0/Y_1)$.** Here the skew term is harmless. For both one-reflection and the skew sextic, the fibre is the plane sextic $g(X_0^2,Y_1^2)+Y_2^6=0$ with $g$ a binary cubic; the skew term only adds $\tfrac{r}{7}x(x^2-1)\,UV^2$ to $g$. The residue character of $A$ is $(+1,-1,\zeta)$ under ($X$-block sign, $(Y_0,Y_1)$ sign, $Y_2\mapsto\zeta Y_2$). This selects the 2-dimensional eigenspace $\{X_0Y_1^2,\,X_0^3\}$ of $H^{1,0}$ of the fibre, with no $(0,1)$ part. So the relevant variation is a signature-$(2,0)$ unitary piece of a $\mu_6$-cover of $\mathbb P^1$, not an elliptic curve, and the elliptic-companion machinery does not apply directly. The structure is identical for the solved one-reflection, so this fibration neither obstructs nor serves.

**5. Supports and deformation (modular, $p=268435399$).** None of the 55 natural degree-6 supports carries $A$: the 3-subsets and 2–4-subsets of sixth powers, and every coordinate sub-equation of $F$ containing a coupling. On the **solved** one-reflection, $s^2+T^2$, $S^2+T^2$ and $T^2-s^2/3$ carry none of $A$, $A_\pm$ either, while the accepted support is `fqlx`'s degree-30 $D_+$. So low-degree negatives on this family are weak evidence. I reproduced `fqlx`'s half-certificate mod a prime with $\sqrt3$; both identities hold exactly. A first-order deformation of it along the skew term is a linear system with about $1.4\times10^6$ unknowns against $4.1\times10^6$ degree-60 rows, because $D_+$ is semi-invariant only under $\tau$. I did not attempt it.

**Reproduction.** https://github.com/DiamondJM/hodgeHoldout, branch `claude/inspiring-turing-tx2pbu`: `skew/NOTES.md`, `skew/groebner/*.sing`, `skew/fibration.py`, `skew/fqlx_cert.py`, `skew/supp.py`, `sweep/fixed_d.py`.

Credit: `fqlx`'s one-reflection and holdout constructions; `fringewidth`'s identification of the hidden reflections; the repository authors' certificate formalism, prescribed-support theorem and character framework, with the certificate bridge attributed there to Benjamin E. Diamond; Griffiths's residue theory. Prepared with Claude Code assistance; no independent human review or Lean verification is claimed.

---
_Generated by [Claude Code](https://claude.ai/code)_



## COMMENT DC_kwDOVLOOJc4BJHdj by rbajaj5 at 2026-10-24T16:56:21Z (updated 2026-10-24T16:56:21Z)

**Full-current search update; no solve.** For the unchanged `one-reflection-skew` target, we tested all six scalar coefficients of three rank-two physical fluxes in

$$
B+B/A+B/A^2,\qquad
B=\operatorname{span}_{\mathbb Q}\{r^m v^n t^k:-10\le m\le2,\ 0\le n,k\le2\},
$$

where $A(r,v)=r^3+r^{-3}+(rv-1)^6/r^{15}$ is the physical chart coefficient, distinct from the protected source $X_0^3Y_0^3$. Both compact source rows, all physical $r,v,t$ derivatives, the packet's $1/12$ normalization and mixed-sector cancellations are retained.

**This declared finite family is excluded over $\mathbb Q$.** The highest-pole condition restricts only the $A^{-2}$ numerator to the bounded tangent module; all lower-sector coefficients remain free. This leaves 1,698 necessary parameters from the original 2,106. A 1,699-square augmented minor has determinant **14098 modulo 32003** at regular sample points. Separate implementations reconstructed the decisive minor; the linked audit explains the characteristic-zero implication and records each execution's scope.

[Proof details, reproducible checks, visualization and limits](https://github.com/rbajaj5/hodge-lab/blob/44c5789efc08a51b788844951690f3ec8dd6c162/notes/skew-A2-denominator-audit-20261024.md).

**No full primitive, homogeneous certificate, Lean verification of this exclusion, or accepted solve is claimed.** Other numerator supports, pole divisors and pole orders, and general rational or fibre-dependent currents remain open.

Credit **fqlx** for the [compact connection and source formulation](https://github.com/hodge-computer/sextics/discussions/127#discussioncomment-19113830), **benediamond** for the [physical cover and preceding reductions](https://github.com/hodge-computer/sextics/discussions/127#discussioncomment-19129985), and the hodge.computer contributors for the fixed problem and framework. Prepared for `rbajaj5` with ChatGPT Pro and Codex assistance; separate computational checks do not mean independent human review.

