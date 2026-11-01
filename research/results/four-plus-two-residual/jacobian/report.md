# Probe: Jacobian of the binary triple cover for `four-plus-two-residual`

Date: 2026-11-01. Instance `four-plus-two-residual` (q = 1, A = X0^3 Y0^3), F = H(X0,X1,Y0,Y1) + B(X2,Y2) with
B = X2^6 + (15/7) X2^4 Y2^2 + (10/7) X2^2 Y2^4 + Y2^6.
All computations: python3 (numpy 2.5.2, sympy 1.12, python-flint 0.9, mpmath) and Singular 4.x on this machine. Scripts, raw sums, and full output: `/tmp/hodge/results/four-plus-two-residual/jacobian/` (copied from `/tmp/hodge/scratch/jac42/`). Nothing in `/home/kejace/devel/sextics` was modified; nothing was posted.

Prior work: the discussion thread (#120, comment by `benediamond`, 2026-10-23, cached at `/tmp/hodge/discussions/four-plus-two-residual.md`) already contains the character bookkeeping of section 1 (same basis of differentials, same (2,1) signature, same observation that E' is the wrong sector) and explicitly leaves open whether the three-dimensional Prym is simple. Sections 2-3 below settle that.

## 0. Summary (one paragraph)

The bookkeeping in the task statement is correct (section 1). The binary factor relevant to A = X0^3 Y0^3 is the line spanned by the holomorphic form dt/w (weighted numerator u) on the genus-4 curve T_B: u^3 = B; it has deck character zeta^2, is anti-invariant under sigma: Y2 -> -Y2, and lies in the 3-dimensional Prym P of T_B -> E' = T_B/sigma (E' a j = 0 plane cubic whose H^{1,0} has the opposite character zeta). Point counts of T_B and E' over F_p and F_{p^2} for all good primes p <= 400 (36 split primes p = 1 mod 3, 38 inert primes), cross-checked against F_{p^3} counts at 12 primes and F_{p^4} counts at 4 primes, give the Frobenius characteristic polynomials of P. At 28 of the 36 split primes (first one p = 31) the characteristic polynomial of Frobenius on P is an irreducible sextic over Q with no ratio of roots a root of unity. **Hence P_{F-bar_p} is simple, so P over Q-bar is simple (specialisation of endomorphisms), so P contains no elliptic curve at all, in particular none with CM by Q(zeta_3).** Moreover the Frobenius fields Q(pi_p) differ between primes (different quadratic resolvents over Q(zeta_3)), so P has no CM by a sextic field, and End^0(P_{Q-bar}) = Q(zeta_3). Therefore the elliptic-contraction route of Chapter 3A, which needs an elliptic curve with CM by Q(zeta_3) inside P carrying the u-form line as its H^{1,0}, is not available for this B. This is an obstruction to that specific route only; it says nothing about the existence of certificates. The Fano companion W_H = {H + s^3 + t^3 = 0} in P(1,1,1,1,2,2) is quasi-smooth, has exactly three 1/2(1,1,1,1) singular points, -K = O(2), and Griffiths-Steenbrink Hodge numbers h^{4,0} = 0, h^{3,1} = 56, h^{2,2} = 307 (306 primitive), odd cohomology 0 (section 4).

Status labels used below: **[proved]** = mathematical fact with proof given or standard; **[computed]** = exact integer/finite-field computation performed here, reproducible from the scripts; **[heuristic]** = interpretation or expectation.

## 1. Bookkeeping: characters, differentials, Prym

### 1.1 The curves [proved, with computational checks]

In the chart X2 = 1, t = Y2/X2, w = u/X2^2:

- T_B: w^3 = f(t), f(t) = B(1,t) = 1 + (15/7) t^2 + (10/7) t^4 + t^6 = g~(t^2), g~(V) = 1 + (15/7) V + (10/7) V^2 + V^3 = p(1,V), where B(X2,Y2) = p(X2^2, Y2^2), p(U,V) = U^3 + (15/7) U^2 V + (10/7) U V^2 + V^3. [computed: identity checked by sympy.]
- Separability: disc(7 g~) = -11 * 2957; disc of 7 f = -2^6 7^2 11^2 2957^2. So B has six distinct roots (X is smooth; Singular confirms dim V(jacob F) = 0 with Milnor number 5^6 = 15625). Bad primes for everything below: 2, 3, 7, 11, 2957. [computed]
- Riemann-Hurwitz for T_B -> P^1 (degree 3, totally ramified over the six roots of B, unramified over infinity since deg B = 6 is divisible by 3 and the leading coefficient is 1): 2g - 2 = 3(-2) + 6*2 = 6, g(T_B) = 4. Independent check: Griffiths-Steenbrink for u^3 = B in P(1,1,2) gives h^{1,0} = dim R_2 = 4 (R = Jacobian ring of B - u^3, degree 6, sum of weights 4). [proved/computed]
- E' = T_B/sigma, sigma(t,w) = (-t,w): the plane cubic u^3 = p(U,V) (U = X2^2, V = Y2^2), i.e. w^3 = g~(V) with V = t^2. Genus 1, j = 0 (the order-3 automorphism u -> zeta u has fixed points, the three points u = 0). sigma has exactly 6 fixed points on T_B (Y2 = 0 or X2 = 0), consistent with 2*4 - 2 = 2(2*1 - 2) + 6. [proved]
- The sixth-root cover C_B: w^6 = B in P^2 (genus 10) is a double cover of T_B via u = w^2, branched at the six points w = 0. [proved]

### 1.2 Holomorphic differentials and characters [proved]

With rho(t,w) = (t, zeta w) and the pullback action rho^* on forms:

| form | regular? | rho^* eigenvalue | sigma^* |
|---|---|---|---|
| dt/w (weighted numerator u) | yes | zeta^{-1} = zeta^2 | -1 |
| dt/w^2 (numerator X2^2) | yes | zeta^{-2} = zeta | -1 |
| t dt/w^2 (numerator X2 Y2) | yes | zeta | +1 |
| t^2 dt/w^2 (numerator Y2^2) | yes | zeta | -1 |

(Regularity: at the six finite branch points ord(dt) = 2, ord(w) = 1; at the three points over infinity t has a simple pole, w a double pole, dt a double pole, so t^i dt/w^j is regular iff i <= 2j - 2.) Hence:

- H^{1,0}(T_B) = (zeta-part, dim 3) + (zeta^2-part, dim 1); signature of the Z[zeta_3]-action on J(T_B) is (3,1).
- H^{1,0}(E') = <t dt/w^2>, character zeta.
- H^{1,0}(P) = <dt/w, dt/w^2, t^2 dt/w^2>: characters (zeta^2, zeta, zeta), i.e. signature (2,1) with the single zeta^2-line being exactly the u-form. dim P = 3. [proved]
- Polarisation of P induced from J(T_B): for a ramified double cover of an elliptic curve with 6 branch points the Prym has polarisation type (1,1,2) (standard; not used below). [heuristic label only because I did not re-derive it]

### 1.3 Thom-Sebastiani / Shioda-Katsura placement of A = X0^3 Y0^3 [proved]

For F = H + B the primitive H^4 of X decomposes (after the Shioda-Katsura blow-up along the product of the two pure-block loci, plus joining terms) into the diagonal-mu_6-invariant part of H^3(Y_H) (x) H^1(C_B), with Y_H: w^6 = H in P^4 and C_B: w^6 = B in P^2. In Jacobian-ring terms:

- On Y_H (five variables, degree 6, sum of weights 5) the pole-order-m piece is R_{6m-5}; X0^3 Y0^3 w^a has degree 6 + a = 6m - 5 forces a = 1, m = 2: H^{2,1}(Y_H), character eta^{a+1} = eta^2 (eta a primitive 6th root, the volume form contributes one factor eta).
- On C_B (three variables, degree 6, sum of weights 3) the piece is R_{6m-3}; w^b has degree b = 6m - 3 forces b = 3, m = 1: H^{1,0}(C_B), character eta^4.
- eta^2 * eta^4 = 1: the pair is diagonal-invariant. Both characters have order 3, so the pair descends along w -> u = w^2 to (H^{2,1}(T_H)_{zeta}, H^{1,0}(T_B)_{zeta^2}) with T_H: u^3 = H in P(1,1,1,1,2). The curve-side form is w^3 Omega/(B - w^6) ~ w^3 dt/(6 w^5) = dt/(6 u): the u-form, character zeta^2 as in the table. Hodge type (2,1) + (1,0) = (3,1), consistent with q = 1 (pole order 2) on X.
- sigma-parity: A is even in Y2 and the volume form of P^5 is odd, so the residue of A Omega/F^2 is sigma-anti-invariant, matching the fact that dt/w is anti-invariant (lies in the Prym) and E' (the sigma-invariant quotient) does not see it. This is also why the reflection theorem of Chapter 3 does not serve this sector.
- A = X0^3 Y0^3 is its own normal form modulo the Jacobian ideal of F (Singular `reduce`), so the class is nonzero. [computed]

The companion T_H: u^3 = H in P(1,1,1,1,2) has K = O(6 - 6) = O (Calabi-Yau 3-fold) with Griffiths-Steenbrink numbers h^{3,0} = 1, h^{2,1} = 103 = 68 + 35, where R_{H,6} u^0 (dim 68) has deck character zeta and R_{H,4} u (dim 35) has deck character zeta^2 (convention, same as for T_B and as in Chapter 3A: u^k Omega carries zeta^{k+1}, the volume form contributing one zeta). Which piece carries the core class? Pulling back along Y_H -> T_H, u = U^2, the weighted volume form satisfies Omega_T |-> 2U Omega_Y, so the T_H-class A u^k Omega_T/(H - u^3)^2 pulls back to 2 A U^{2k+1} Omega_Y/(H - U^6)^2. The Y_H-class of the target is A U^1 (k = 0), so **the core factor is A u^0 in R_{H,6}, in H^{2,1}(T_H)_zeta** (deck character eta^2 = zeta, as it must be), not A u in R_{H,4} u. Likewise the curve factor V^3 Omega_C on C_B is the pullback of u^1 Omega_{T_B} (2k+1 = 3, k = 1), the u-form, deck character zeta^2. The sector of A is therefore

    H^{2,1}(T_H)_zeta  (x)  H^{1,0}(T_B)_{zeta^2}  =  (68-dimensional core piece)  (x)  (the u-line),

a (3,1)-class on X, consistent with q = 1. What is convention-independent and matters below: **the curve-side line is the u-form, whose deck character (zeta^2) is opposite to that of H^{1,0}(E') (zeta).** [proved]

### 1.4 What the 3A mechanism would need here [proved as a statement about the mechanism]

In Chapter 3A (3+3 split) the Fano companion W_f supports [H^2(T_f) (x) H^1(E)]^{mu_3} for the Fermat cubic E, and a translated diagonal in E x E contracts the two elliptic factors against each other, using that H^1(E)_zeta and H^1(E)_{zeta^2} pair nondegenerately. In the 4+2 situation the same construction would give support for [H^3(T_H) (x) H^1(E)]^{mu_3} on a resolution of W_H (section 4); the target's core factor H^{2,1}(T_H)_zeta appears there tensored with H^1(E)_{zeta^2} = H^{0,1}(E). One then needs to contract that H^1(E)_{zeta^2} against the curve side H^1(T_B)_{zeta^2} so that the image contains the u-form (type (1,0)). Algebraic correspondences preserve Hodge types, so this cannot be done by a correspondence from E itself (whose zeta^2-part is of type (0,1)); it needs the curve side to be, up to isogeny, an elliptic curve E'' with Q(zeta_3)-CM whose H^{1,0} is the u-line, after which the translated-diagonal pairing of 3A.3 between H^1(E)_{zeta^2} and H^1(E'')_{zeta}-type pieces applies. That requires a Z[zeta_3]-equivariant correspondence between E (or any j = 0 curve) and T_B whose image in H^1(T_B) contains the u-line, i.e. an elliptic curve E'' with CM by Q(zeta_3) inside J(T_B) with H^{1,0}(E'') = <dt/w>. Since dt/w is sigma-anti-invariant, E'' would have to sit inside P. Sections 2-3 show no such E'' exists.

Side remark [proved]: the visible elliptic quotient E' does serve the sector H^3(T_H)_{(opposite char)} (x) H^1(E'), i.e. curve-side forms t dt/w^2 (numerator X2 Y2). Those targets are odd in Y2, so they are already covered by the reflection theorem; the E'-contraction adds nothing beyond it. Everything in the Prym (curve numerators X2^2, Y2^2 and u) needs something new.

## 2. Point counts and Frobenius polynomials

### 2.1 Method [proved]

Let q = p^k with q = 1 mod 3, omega in F_q a primitive cube root of unity, chi(x) = x^{(q-1)/3} in {1, omega, omega^2}, iota: omega -> zeta (this fixes the prime pp = (p, zeta - omega) of Z[zeta] when p = 1 mod 3). Let H_+ be the eigenspace of H^1(T_B) on which rho acts (in the Lefschetz convention) by zeta, H_- its conjugate. Counting fixed points of rho^{-j} o Frob_q on the affine model w^3 = f(t) (for each t with f(t) != 0 either all three or none of the w satisfy w^{q-1} = omega^j, according to chi(f(t)) = omega^j; the three points at infinity behave like a point with f = 1) and averaging over j gives

    tr(Frob_q | H_+) = - S_q(f),   S_q(f) := sum_{t in P^1(F_q)} iota(chi(f(t)))  in Z[zeta],  with f(infinity) := 1,

and #T_B(F_q) = q + 1 + S_q + conj(S_q). The same formula with g~ in place of f gives E'. Subtracting gives the Prym: tr(Frob_q | H_+(P)) = -(S_q(f) - S_q(g~)) = -sum_V lambda(V) iota(chi(g~(V))) (lambda = quadratic character), a sextic character sum.

For p = 1 mod 3 the characteristic polynomial of Frob_p on H_+(P) is a cubic L_+(T) = 1 + b1 T + b2 T^2 + b3 T^3 with b_i in Z[zeta], and L_P(T) = L_+(T) conj(L_+)(T) in Z[T]. The functional equation (roots of conj(L_+) are p/alpha_i) gives b2 = conj(b1) b3 / p, so **b3 = p b2 / conj(b1)** whenever b1 != 0; thus F_p and F_{p^2} sums already determine L_P, and the F_{p^3} sum is an independent check (Newton's identity t3 = e1 t2 - e2 t1 + 3 e3). b1 was nonzero at all 36 split primes.

For p = 2 mod 3 Frobenius swaps H_+ and H_-, so L_P(T) = M(T^2) with M(S) = det(1 - Frob_{p^2} S | H_+(P)) in Z[S] of degree 3 and M(S) = 1 + b S + p b S^2 + p^3 S^3 (functional equation with leading coefficient p^3 forced by det Frob_p = p^3). **Lemma [proved]:** such an M always factors as (1 + pS)(1 + (b - p) S + p^2 S^2). (Three roots of absolute value p closed under conjugation with product -p^3 must include -p.) Consequently, at every inert prime, P mod p acquires a supersingular elliptic factor with L-factor 1 + pT^2, for *every* abelian 3-fold with Z[zeta_3]-action; this factor carries no information about elliptic factors over Q-bar. Over F_{p^2}, P ~ E_ss x B^2 with B an elliptic curve of trace p - b, ordinary iff p does not divide b - p. Inert primes therefore cannot exclude a geometric elliptic CM factor and were used only as consistency checks; the decisive data is at split primes. The F_{p^4} prediction #T_B(F_{p^4}) = p^4 + 1 - 2 (p - b)^2 was verified at p = 5, 17, 23, 29. [computed]

Implementation: `ffcount.py` (vectorised numpy arithmetic in F_p[x]/(m(x)), m irreducible found with python-flint; chi by binary exponentiation; all arithmetic exact integers mod p), `driver.py` (primes 5..400 minus bad primes; k = 1,2 for split primes, k = 3 additionally for split p <= 110, k = 2 for inert primes, k = 4 additionally for inert p <= 30), `analyze.py` (Z[zeta] arithmetic, reconstruction, factorisation over Q and Q(zeta_3) with sympy, exact root-of-unity test via Res_x(f(x), f(zx)) against cyclotomic polynomials Phi_m, m <= 200, plus an independent high-precision numerical angle test for m <= 2000, pp-adic Newton polygons, relative discriminant classes), `table.py`. Full per-prime output: `analysis_full.txt`, `analysis.json`, `tables.md`, `split_short.md`.

### 2.2 Internal consistency checks, all passed [computed]

- E' at split primes: pi_E' = -S_p(g~) has norm p at all 36 primes; -S_{p^2}(g~) = pi_E'^2 and (where computed) -S_{p^3}(g~) = pi_E'^3; #E'(F_p) = p + 1 - a_p with a_p = pi + conj(pi); a_p^2 - 4p = -3 m^2 at every prime (j = 0 CM shape).
- E' at inert primes: #E'(F_{p^2}) = (p+1)^2 and Frob_{p^2} = -p on H_+ (supersingular) at all 38 primes.
- b3 from the functional equation equals b3 from the F_{p^3} count at p = 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103, 109 (12/12).
- L_P has integer coefficients and all 6 roots of absolute value p^{-1/2} numerically at all 36 split primes; at inert primes the F_{p^2} trace of P is a rational integer at all 38 primes.
- pp-adic Newton polygon of L_+(P) is {0, 1, 2} (one unit root, two roots of valuation 1) at all 36 split primes; at pp-bar it is {0, 0, 1}; and v_pp(pi_E') = 1 at all 36. By the unit-root/Hodge comparison (unit roots correspond to H^1(O), i.e. H^{0,1}) this says h^{0,1}_+(P) = 1, h^{1,0}_+(P) = 2 and h^{1,0}_+(E') = 1, h^{0,1}_+(E') = 0: exactly the (2,1) signature of P and the zeta-character of H^{1,0}(E') from section 1, now confirmed independently at every split prime, and it pins H_+ = the zeta-eigenspace, so the u-form's conjugate is the unique unit-root direction in H_+(P). (The identification "unit root = H^{0,1}" is standard; the consistency of 36 primes with E' as calibration leaves no sign ambiguity.)

### 2.3 Results at split primes (p = 1 mod 3) [computed]

Condensed table (full data with b1, b2, b3 in Z[zeta], Newton polygons, factorisations over Q(zeta_3) in `tables.md`):

| p | #T_B(F_p) | #E'(F_p) | a_p(E') | c1, c2, c3 of L_P = 1 + c1 T + c2 T^2 + c3 T^3 + p c2 T^4 + p^2 c1 T^5 + p^3 T^6 | L_P over Q |
|---|---|---|---|---|---|
| 13 | 20 | 21 | -7 | -1, -1, 50 | reducible: (13*T**2 - 2*T + 1) * (169*T**4 + 13*T**3 - 12*T**2 + T + 1) |
| 19 | 30 | 27 | -7 | 3, -15, -142 | reducible: (19*T**2 + 8*T + 1) * (361*T**4 - 95*T**3 + 6*T**2 - 5*T + 1) |
| 31 | 24 | 25 | 7 | -1, -19, 194 | irreducible sextic |
| 37 | 30 | 27 | 11 | 3, 27, -158 | irreducible sextic |
| 43 | 54 | 57 | -13 | -3, 33, -322 | reducible: (43*T**2 + 8*T + 1) * (1849*T**4 - 473*T**3 + 78*T**2 - 11*T + 1) |
| 61 | 84 | 75 | -13 | 9, 135, 1150 | reducible: (61*T**2 + 13*T + 1) * (61*T**2 - 2*T + 1)^2 |
| 67 | 54 | 63 | 5 | -9, 129, -598 | irreducible sextic |
| 73 | 86 | 91 | -17 | -5, 71, -278 | irreducible sextic |
| 79 | 84 | 93 | -13 | -9, 45, -574 | reducible: (79*T**2 - 4*T + 1) * (6241*T**4 - 395*T**3 - 54*T**2 - 5*T + 1) |
| 97 | 96 | 93 | 5 | 3, -21, -542 | irreducible sextic |
| 103 | 86 | 91 | 13 | -5, 21, 1130 | irreducible sextic |
| 109 | 78 | 93 | 17 | -15, 147, -1082 | irreducible sextic |
| 127 | 156 | 148 | -20 | 8, -43, -2656 | irreducible sextic |
| 139 | 120 | 133 | 7 | -13, 17, 1010 | irreducible sextic |
| 151 | 158 | 175 | -23 | -17, 261, -4990 | irreducible sextic |
| 157 | 182 | 172 | -14 | 10, 35, 1660 | irreducible sextic |
| 163 | 182 | 181 | -17 | 1, -87, -154 | irreducible sextic |
| 181 | 176 | 163 | 19 | 13, -169, -5018 | reducible: (181*T**2 + 26*T + 1) * (32761*T**4 - 2353*T**3 - 12*T**2 - 13*T + 1) |
| 193 | 180 | 217 | -23 | -37, 723, -10958 | irreducible sextic |
| 199 | 168 | 172 | 28 | -4, 245, -1336 | irreducible sextic |
| 211 | 240 | 225 | -13 | 15, 297, 1418 | irreducible sextic |
| 223 | 234 | 201 | 23 | 33, 669, 11918 | irreducible sextic |
| 229 | 254 | 252 | -22 | 2, -85, 2252 | irreducible sextic |
| 241 | 218 | 228 | 14 | -10, -85, 5276 | irreducible sextic |
| 271 | 296 | 301 | -29 | -5, 261, -6886 | reducible: (271*T**2 - 29*T + 1) * (271*T**2 + 12*T + 1)^2 |
| 277 | 240 | 247 | 31 | -7, 159, -8978 | irreducible sextic |
| 283 | 260 | 316 | -32 | -56, 1665, -33040 | irreducible sextic |
| 307 | 342 | 327 | -19 | 15, -303, -9286 | irreducible sextic |
| 313 | 270 | 279 | 35 | -9, 183, -11678 | irreducible sextic |
| 331 | 384 | 364 | -32 | 20, 137, -4888 | irreducible sextic |
| 337 | 270 | 333 | 5 | -63, 1791, -35570 | irreducible sextic |
| 349 | 420 | 387 | -37 | 33, 903, 23182 | reducible: (349*T**2 + 37*T + 1) * (349*T**2 - 2*T + 1)^2 |
| 367 | 300 | 337 | 31 | -37, 725, -10582 | irreducible sextic |
| 373 | 422 | 412 | -38 | 10, 203, 14428 | irreducible sextic |
| 379 | 360 | 417 | -37 | -57, 2073, -47302 | irreducible sextic |
| 397 | 336 | 363 | 35 | -27, 639, -19802 | irreducible sextic |

Findings:

1. **P mod p is geometrically simple at 28 of the 36 split primes:** p = 31, 37, 67, 73, 97, 103, 109, 127, 139, 151, 157, 163, 193, 199, 211, 223, 229, 241, 277, 283, 307, 313, 331, 337, 367, 373, 379, 397. At each, L_P is an irreducible sextic over Q (so P_{F_p} is simple by Tate), and no ratio alpha_i/alpha_j (i != j) is a root of unity (exact resultant test against Phi_m, m <= 200, and the numerical test to m <= 2000 agree), so the characteristic polynomial of Frob_p^k is irreducible for every k, i.e. P is simple over F-bar_p.
2. At the 8 remaining split primes P mod p is not simple: at p = 13, 19, 43, 79, 181 it splits as (j = 0 elliptic curve, traces 2, -8, -8, 4, -26) x (simple abelian surface whose Frobenius eigenvalues have a ratio of order 3); at p = 61, 271, 349 it splits as E'_p x B_p^2 (the reduction of E' itself, or its twist at 271, times the square of an ordinary non-CM-by-Q(zeta_3) elliptic curve). In all 8 cases the elliptic eigenvalue in H_+ has pp-adic valuation 1, i.e. the sporadic mod-p elliptic factor carries the character zeta on H^0(Omega^1) and never the u-line's character zeta^2.
3. Frobenius fields differ between primes: the relative discriminant of the cubic L_+ over Q(zeta_3), reduced modulo squares, has odd exponent at the inert prime 1901 for p = 31, at 461 for p = 67, at 1319 for p = 73, at 23 for p = 97, etc. (28 distinct classes listed in `tables.md`). Hence the cubic extensions Q(zeta_3)(pi_p)/Q(zeta_3) are pairwise non-isomorphic for (at least) these primes.

### 2.4 Results at inert primes (p = 2 mod 3) [computed]

At all 38 inert primes 5 <= p <= 389 the structure of section 2.1 holds: M(S) = (1 + pS)(1 + c S + p^2 S^2) with c = b - p; c/p is not in {0, +-1, +-2} except at p = 5 (c = 2p: P mod 5 is entirely supersingular, Frob_{25} = -5), and the second factor is ordinary (p does not divide c) at all other 37 primes. As explained, these primes are consistent with any answer to the geometric question and carry no evidence either way; they only corroborate the Z[zeta_3]-structure and the counts.

## 3. Conclusions

### 3.1 Theorem-level statements [proved, given the exactness of the computed counts]

(a) P (the Prym of T_B -> E', equivalently the sigma-anti-invariant part of J(T_B)) is **absolutely simple**: its reduction at p = 31 (or any of the 27 other listed primes) is simple over F-bar_p, and an isogeny decomposition over Q-bar would specialise to one over F-bar_p (End(P_{Q-bar}) injects into End(P_{F-bar_p}) for good reduction). In particular P contains no elliptic curve, with or without CM.

(b) P is **not of CM type**, and End^0(P_{Q-bar}) = Q(zeta_3). Proof: for a simple abelian 3-fold the Albert classification allows End^0 to be Q, a totally real cubic field, or a CM field of degree 2 or 6; it contains Q(zeta_3), so it is Q(zeta_3) or a sextic CM field K. In the latter case, after a finite extension all endomorphisms are defined and the Frobenius at a prime above p lies in K, so Q(pi_p^f) = Q(pi_p) (no root-of-unity ratios) is a subfield of K of degree 6, hence equals K, for every p in the list of (a). But the fields Q(pi_31) and Q(pi_67) have different quadratic resolvents over Q(zeta_3) (section 2.3, item 3), contradiction.

(c) J(T_B) ~ E' x P with E' the j = 0 curve u^3 = p(U,V), and this is the complete isogeny decomposition over Q-bar. The only elliptic curve with CM by Q(zeta_3) in J(T_B) is E', whose H^{1,0} is the sigma-invariant line t dt/w^2 of character zeta; the u-form dt/w (character zeta^2) is not the H^{1,0} of any elliptic curve in J(T_B).

Caveats on rigour: the inputs are exact integer computations (finite-field enumeration, Z[zeta] arithmetic, sympy factorisation over Q, exact resultants) with the extensive redundancy of section 2.2; the mathematical steps are standard (Lefschetz trace formula for twisted Frobenius, Tate's theorem that an irreducible Frobenius characteristic polynomial implies simplicity over the finite field, specialisation of endomorphisms, Albert classification, Shimura-Taniyama for CM Frobenii). I did not run the computation through a second independent system (no Sage/Magma/PARI available), so "proved" above means "proved modulo the correctness of these scripts", which the twelve F_{p^3} cross-checks and the Newton-polygon pattern make very robust.

### 3.2 Answer to the route question [proved consequence + heuristic framing]

The elliptic-contraction mechanism of Chapter 3A, transplanted to the 4+2 split, needs an elliptic curve E'' with CM by Q(zeta_3) inside J(T_B) whose H^{1,0} is the u-line, so that the H^1(E) factor produced by the Fano companion W_H can be contracted against the curve side with the target's binary class in the image. By 3.1(a)-(c) no such E'' exists for this B: the u-line lies in the absolutely simple 3-fold P with End^0 = Q(zeta_3). **The elliptic-contraction route is not available for this instance.** This is an obstruction to that specific route only (any route that reduces the binary side to an elliptic curve with Q(zeta_3)-CM, including via any correspondence from any j = 0 curve); it is not a statement about the nonexistence of certificates, and it does not obstruct routes that work directly with the 3-dimensional P (e.g. a correspondence from a Picard-type curve or a surface whose H^2 maps onto H^1(T_H) (x) H^1(P)-type classes, or the "joining construction" mentioned in why.md).

[heuristic] P is a member of the 2-parameter family of abelian 3-folds with Z[zeta_3]-action of signature (2,1) (a Picard modular surface); the generic member has End^0 = Q(zeta_3), and 3.1 shows this particular B is generic in that sense. The eight split primes with a sporadic elliptic factor are the usual density-zero exceptional set and carry no geometric information; notably, even those sporadic factors have the wrong character.

[heuristic] A different-looking hope, that P might be isogenous to the Jacobian of a Picard curve y^3 = quartic (also signature (2,1)) and that such an identification could be exploited, is not excluded by this probe (P is principally polarisable up to isogeny only; the induced polarisation type is (1,1,2)). But even if true it would not reintroduce an elliptic factor.

## 4. The Fano companion W_H [computed/proved]

W_H = {H(X0,X1,Y0,Y1) + s^3 + t^3 = 0} in P(1,1,1,1,2,2), degree 6.

- Quasi-smooth: Singular, ring Q[X0,X1,Y0,Y1,s,t], `dim(std(jacob(W))) = 0` and `vdim = 2500 = 5^4 * 2 * 2` (the affine cone is smooth off the origin). Equivalent to smoothness of the sextic surface S_H = {H = 0} in P^3 (`dim V(jacob H) = 0`, Milnor number 625). [computed]
- Singular locus: the ambient singular stratum is the weight-2 line {X = Y = 0} = P(2,2) with isotropy mu_2 acting by (-1,-1,-1,-1) on the four weight-1 directions; W_H meets it in s^3 + t^3 = 0, i.e. exactly three points, each a cyclic quotient singularity 1/2(1,1,1,1) (Gorenstein since 4 * 1/2 is an integer, terminal since the age is 2 > 1). No other singularities (quasi-smooth and well-formed). [computed + proved]
- Canonical class: weighted adjunction, -K_W = O(1+1+1+1+2+2 - 6) = O(2); O(2) is Cartier (its local generator s or t is mu_2-invariant at the three singular points), ample, so W_H is a Gorenstein terminal Fano 4-fold with -K Cartier and the three singular points are the images of T_H x {c = 0} under the rational map T_H x E --> W_H, ([z:u],[a:b:c]) -> [cz : uac : ubc] (same formula as 3A.2; the three divisors T_H x {point} are contracted to the three singular points, and conversely the three divisors {H = 0} x {s = -zeta^k t} of W_H are contracted to S_H x {point} in T_H x E). The Miyaoka-Mori-type bound used in Chapter 3 gives covering rational curves of -K-degree <= 2n = 8, i.e. O(1)-degree <= 4. [proved]
- Hodge numbers (Griffiths-Steenbrink for quasi-smooth weighted hypersurfaces, h^{4-p,p}_prim = dim R_{6(p+1)-8}, R = C[X0,X1,Y0,Y1,s,t]/(jacob W) = R_H (x) C[s,t]/(s^2,t^2) with s, t of weight 2; Hilbert series of R_H is (1+z+z^2+z^3+z^4)^4 = 1,4,10,20,35,52,68,80,85,80,68,52,35,20,10,4,1): [computed]

      h^{4,0} = R_{-2} = 0,   h^{3,1} = R_4 = 56,   h^{2,2}_prim = R_10 = 306,   h^{1,3} = 56,   h^{0,4} = 0,
      h^{1,1} = 1, h^{2,2} = 307, all odd Betti numbers 0, b_4 = 419, e(W_H) = 423 (as a V-manifold / intersection cohomology; H^4 is pure).

  Decomposition of R_4 = R_{H,4} * 1 (+) R_{H,2} * {s,t} (+) R_{H,0} * st = 35 + 20 + 1: [heuristic reading] the 35 is H^{2,1}(T_H)_{zeta^2} (x) H^{1,0}(E)_zeta and the 1 is H^{3,0}(T_H)_zeta (x) H^{0,1}(E)_{zeta^2}, exactly the two invariant (3,1)-pieces of [H^3(T_H) (x) H^1(E)]^{mu_3} (E the Fermat cubic with H^{1,0}(E) of character zeta as in 3A.1); the target's core factor H^{2,1}(T_H)_zeta pairs with H^{0,1}(E)_{zeta^2} and sits in the invariant (2,2)-part of W_H, which the covering-curve support theorem covers equally; the extra 20 = 2 * p_g(S_H) (S_H sextic surface, p_g = 10) comes from the divisors {H = 0} x {s = -zeta^k t} contracted by the birational map to T_H x E, exactly parallel to the 3-fold case W_f of Chapter 3A where the same count gives R_5 = 18 + 20 + 3 = 41 = 21 + 2 * g(C_f). So W_H carries, beyond the T_H x E part, a copy of H^2(S_H)(-1)-type classes; none of this changes the sector analysis above.

## 5. Reproduction

    cd /tmp/hodge/results/four-plus-two-residual/jacobian
    python3 prelim2.py                   # identities, discriminants, Hilbert series, Hodge numbers
    Singular -q wh.sing; Singular -q wh2.sing   # quasi-smoothness, Milnor numbers, singular points, A mod J(F)
    python3 driver.py 5 400              # ~15 min on one core (parallelise by ranges as in the logs)
    python3 analyze.py > analysis_full.txt; python3 table.py > tables.md

Raw twisted sums are in `sums_*.jsonl` (per prime: S_q as (a,b) meaning a + b zeta, point counts, the chosen omega). Per-prime analysis in `analysis.json`.
