# Part 1 — Character transfer for a cyclic cover: Esnault–Viehweg, what a support must carry, and W-tower rigidity

Author: wave-2 trio-theory agent (restart). Receipts: /tmp/hodge/results/_cross/trio.log.jsonl, trio_*.py, trio-*.json.
Conventions: X ⊂ P^5 the smooth sextic 4-fold of the lane, ξ = [A Ω/F²] ∈ H^{3,1}(X) the target residue class (q = 1,
A = X0³Y0³ in all three lanes), G = ⟨g⟩ ≅ Z/n the cyclic symmetry used by the accepted companion, π: X → W = X/G,
B ⊂ W the branch divisor, χ the character of G on ξ (χ ≠ 1 is the whole point of the three lanes).

## 1.1 The three covers in one table

| lane | g | n | W (weighted hypersurface) | K_W | B | L (L^n = O(B)) | χ(ξ) | serving quotient (companion) |
|---|---|---|---|---|---|---|---|---|
| ctr | X2 ↦ ζ3 X2 | 3 | {H(x, Z) = 0} ⊂ P(1^5, 3), Z = X2³ | O(−2) | {Z = 0} ∈ \|O_W(3)\| | O_W(1) | ζ3 | W (Fano index 2) |
| rcr | Y2 ↦ −Y2 | 2 | {H(x, u) = 0} ⊂ P(1^5, 2), u = Y2² | O(−1) | {u = 0} ∈ \|O_W(2)\| | O_W(1) | −1 | W (Fano index 1) |
| qcr | s0: X0 ↔ Y0 | 2 | {F(w,v²,…) = 0} ⊂ P(1^5, 2), w = X0+Y0, t = (X0−Y0)² | O(−1) | {t = 0} ∈ \|O_W(2)\| | O_W(1) | −1 | X/s0 (Fano index 1) and X/⟨s0,σ0⟩ ⊂ P(1^4,2,2) (Fano index 2) |

[computed] The three branch threefolds B ≅ X ∩ {X2 = 0}, X ∩ {Y2 = 0}, X ∩ {X0 = Y0} are smooth: Singular, p = 32003,
`dim(std(jacob(f))) = 0` for each (receipts /tmp/hodge/results/_cross/trio_branch_{ctr,rcr,qcr}.sing); an m-primary Jacobian
ideal mod p is m-primary over Q, so B and hence W are smooth and the Esnault–Viehweg formulas below apply without any
normal-crossing or rounding corrections.

In every case π_*O_X = ⊕_{i=0}^{n−1} L^{−i} with g acting on the summand L^{−i} by ζ_n^{i} (the eigen-function is the
anti-invariant coordinate X2^i, Y2, resp. v = X0 − Y0), W is smooth exactly where B is smooth (local model x^n = z), and
χ(ξ) = (polynomial character of A)·det(g) is the nontrivial character because A is G-invariant while Ω picks up det g.

## 1.2 Esnault–Viehweg decomposition (the eigen-part in terms of W and B)

For a cyclic cover π: X → W of degree n with smooth (more generally normal-crossing, with the usual rounding) branch
divisor B and π_*O_X = ⊕ L^{(i)−1}, L^{(i)} = L^i ⊗ O(−⌊iB/n⌋) = L^i for i < n when B is reduced,
   π_*Ω^p_X = Ω^p_W ⊕ ⊕_{i=1}^{n−1} Ω^p_W(log B) ⊗ L^{−i},
the summand i being the ζ^i-eigensheaf [Esnault–Viehweg, Lectures on vanishing theorems, §3; also Esnault–Viehweg 1992
"Logarithmic de Rham complexes and vanishing theorems"]. Hence, writing H^{p,q}(X)_i for the ζ^i-eigenspace,

   H^{p,q}(X)_i = H^q(W, Ω^p_W(log B) ⊗ L^{−i})        (i = 1, …, n−1),     H^{p,q}(X)_0 = H^{p,q}(W).

Applied to the target sector (ζ^1 for ctr, the sign −1 for rcr/qcr):
* H^{4,0}(X)_1 = H^0(W, K_W + B − L) = H^0(W, O(0)) = C in all three lanes (K_W + B − L = (−2+3−1)H, (−1+2−1)H, (−1+2−1)H):
  the holomorphic 4-form Ω/F has the character χ (it is X2-twisted / Y2-odd / s0-odd).  For ctr the ζ3²-eigenspace has
  K_W + B − 2L = O(−1): no 4-form, consistent with the conjugate sector.
* H^{3,1}(X)_χ = H^1(W, Ω^3_W(log B) ⊗ L^{−1}),  H^{2,2}(X)_χ = H^2(W, Ω^2_W(log B) ⊗ L^{−1}),  H^{1,3}(X)_χ = H^3(W, Ω^1_W(log B)(−1)).
  The residue sequence 0 → Ω^p_W ⊗ L^{−1} → Ω^p_W(log B) ⊗ L^{−1} → Ω^{p−1}_B ⊗ L^{−1} → 0 relates these to the
  cohomology of the branch threefold B (a sextic threefold in P^4 in all three lanes: B ≅ X ∩ {X2 = 0}, X ∩ {Y2 = 0},
  X ∩ {X0 = Y0}).  Numerically the eigen-dimensions are those of the equivariant Jacobian ring (dossiers): the full
  χ-eigenspaces are (1; 107+τ-odd part; …) for ctr, (1, 126+126, …) for rcr (the r-odd part is the (−,−) ⊕ (−,+) Klein
  sectors), (1, 22+…, …) for qcr (the s0-odd part contains the ε-sector (0,22,51,22,0) and other Γ-characters).
  The lane targets lie in the subsectors cut out by the remaining symmetries (τ, ρ; r1, G_diag), which commute with G and
  act on the right-hand sides: the target sectors are the corresponding isotypic parts of H^1(W, Ω^3_W(log B)(−1)).

Why this is the right language for "character transfer": the companion construction serves H^4(W) = H^4(X)_0, which is
divisor-supported because W is a Fano (uniruled) variety (book ch. 3).  Nothing in that mechanism touches
H^*(W, Ω^•_W(log B) ⊗ L^{−1}), a twisted-coefficient cohomology on W that does not come from any variety dominated by
rational curves of W: a rational curve C ⊂ W meeting B transversally has π^{−1}(C) irreducible (a cyclic cover of P^1
branched at ≥ 2 points, genus (n−1)(deg B|_C − 2)/2 > 0 for deg B|_C ≥ 3), so the curve family that sweeps W lifts to a family of
positive-genus curves on X whose Jacobians carry the twisted classes (this is the construct lens's M5 / E5 picture in
EV language: the twisted sector is the H^1 of the family of hyperelliptic / j = 0 curves over the parameter space).

## 1.3 What a divisor support must carry (necessary conditions, W-language)

Let E ⊂ X be a divisor and Ẽ → E a resolution with Gysin map γ: H^2(Ẽ)(−1) → H^4(X).  WLOG E is G-invariant (replace
E by ∪_g gE, i.e. D by ∏ g·D: a legal denominator), so E = π^{−1}(E_W), and the resolution can be taken G-equivariant,
Ẽ → Ẽ_W := Ẽ/G a cyclic cover branched along B̃ (strict transform of B ∩ E_W plus exceptional components over which the
cover ramifies).  γ is G-equivariant, so ξ ∈ im γ forces
   (S1)  H^{2,0}(Ẽ)_χ ≠ 0, i.e.  H^0(Ẽ_W, Ω^2_{Ẽ_W}(log B̃) ⊗ L̃^{(1)−1}) ≠ 0   (EV on the resolved pair),
   (S2)  the Gysin image of H^{2,0}(Ẽ)_χ must be nonzero in H^{3,1}(X)_χ (not killed by the restriction of ξ to Ẽ: equivalently
         ∫_X ξ ∪ γ(η) = ∫_Ẽ ξ|_Ẽ ∪ η ≠ 0 for some η ∈ H^{0,2}(Ẽ)_χ̄ — the pairing with the (1,3)-part of the conjugate sector).
Three consequences, all exact:
(a) Smooth hypersurface sections do not work.  If E_W = W ∩ {d = 0} is smooth and meets B transversally (so Ẽ = E is
    the smooth divisor X ∩ {d(x, x_n^n) = 0}), then H^0(E_W, Ω^2_{E_W}(log B_E) ⊗ L^{−1}) sits between
    H^0(E_W, Ω^2_{E_W}(−1)) = 0 and H^0(B_E, Ω^1_{B_E}(−1)) = 0 (Bott/Lefschetz vanishing for complete intersections in
    weighted projective space, negative twist), hence vanishes: a smooth G-invariant divisor carries no twisted 2-forms.
    (In X-language: a smooth divisor of a sextic 4-fold has H^{2,0} = 0 by Lefschetz — the book's ch. 3 §3.1 statement; the
    EV form shows the same for the twisted part of singular-but-mild (rational singularities) E_W.)
(b) The twisted (2,0)-forms must be created by singularities of E_W (or non-normality of E) of a type whose resolution
    acquires twisted forms: cyclic-cover ramification along a singular curve of B ∩ E_W.  The two known mechanisms are
    exactly the two local models in which this happens:
      • rcr/qcr (n = 2): E_W swept by rational curves C_s tangent to B to even order at every intersection point, so
        π^{−1}(C_s) splits or has lower genus; the twisted forms live in H^0(S, Ω^1_S ⊗ f_*ω_{Ẽ/S}) for the slice S of the
        family (construct lens M5), with f_*ω of rank e − 1 for curves of weighted degree e (hyperelliptic lifts).
      • ctr (n = 3): a curve Γ of ordinary triple points of the surface B ∩ E_W, giving simple-elliptic Ẽ6 singularities of
        E along Γ whose exceptional j = 0 elliptic surface carries forms of character ζ3 (construct lens part 4).
    EV makes the common statement: in both cases the forms are sections of Ω^2_{Ẽ_W}(log B̃) ⊗ L̃^{(1)−1} coming from the
    exceptional divisors of Ẽ_W → E_W where B̃ fails to be the pull-back of B (the rounding ⌊B̃/n⌋ changes L̃^{(1)}).
(c) Quantitative budget (rigorous over Q via the infinitesimal variation of the sector along the symmetric family): a
    support that deforms with F along the G-symmetric family must carry the whole flat closure of ξ.  For ctr the construct
    lens proved Def_G·ξ = all 119 and Def_G²·ξ = all 26 ⇒ Σ_twisted h^{2,0}(Ẽ) ≥ 41 + 26 = 67 for very general members.
    The same test for rcr and qcr is run in Part 3 (new).

## 1.4 W-tower rigidity (theorem; closes "second covers", "fibre products", "changed branch sections" for all three lanes)

Setting.  π: X → W = X/G, G = ⟨g⟩ ≅ Z/n with n ∈ {2, 3}, ξ ∈ H^4(X, C) with g^*ξ = χ(g) ξ, χ ≠ 1 (hence χ faithful).
Let p: Z → W be ANY finite cover that factors through X (Z → X → W; e.g. Z = X ×_W X' for a second cyclic cover
X' → W branched along B' ∈ |O(B)| or any other member of any linear system, or Z = normalisation of a tower of cyclic
covers adjoined over W).  Let Ẑ → W be the Galois closure, M = Gal(Ẑ/W), N = Gal(Ẑ/X) ⊂ M (normal, M/N ≅ G).

Lemma 1 (how M acts on the pulled-back class).  Let ξ̂ = p_X^*ξ ∈ H^4(Ẑ).  For every h ∈ M, h^*ξ̂ = χ(h̄) ξ̂, where
h̄ ∈ M/N = G.  Proof: p_X ∘ h = h̄ ∘ p_X since h̄ is by definition the automorphism of X = Ẑ/N induced by h. ∎

Theorem 1 (rigidity).  For a subgroup K ⊂ M let e_K = |K|^{−1} Σ_{k∈K} k^* be the averaging projector, so that
H^4(Ẑ/K) = e_K H^4(Ẑ).  Then e_K ξ̂ ≠ 0 ⇔ K ⊂ N ⇔ Ẑ/K → X is a (branched) cover.  In particular the only intermediate
quotients of Ẑ whose cohomology sees ξ are the covers of X; W, X', X'' = Z/(diagonal), every "twisted" quotient and
every quotient of a fibre product over W, sees nothing.
Proof: e_K ξ̂ = (|K|^{−1} Σ_{k∈K} χ(k̄)) ξ̂, and Σ_{k∈K} χ(k̄) = |K ∩ N| Σ_{γ ∈ image(K)} χ(γ) vanishes unless χ is trivial on
image(K) ⊂ G, i.e. (χ faithful) image(K) = 1, i.e. K ⊂ N. ∎

Corollary 1 (fibre products / second branch sections).  Let X' → W be any cyclic cover branched along B' (any member
of |O(B)|, or of any other linear system, with any twist L'), Z = X ×_W X' (or its normalisation), with group
G × G'.  The ζ-eigenpart H^4(X)_χ pulls back to the (χ, 1)-isotypic part of H^4(Z); the only quotients of Z containing
it are Z and X.  The "diagonal" covers Z/{(g, g^{−1})} (branched along B + B', twist L ⊗ L') and Z/{(g,g)} carry the
(χ, χ)- and (χ, χ̄)-parts, which are orthogonal to the pull-back of H^4(X).  Hence no correspondence built from the
fibre product (compositions of push-forwards and pull-backs between intermediate quotients) moves a class of H^4(X')
or of any diagonal cover into H^4(X)_χ; the push-pull X' ← Z → X is zero on the χ-sector by character orthogonality.

Corollary 2 (abelian towers; the ctr construct lens's Lemma 1.4 in general form).  No tower of cyclic/abelian covers
over W (adjoining roots of any functions on W: the discriminant, the branch equation, the companion's auxiliary
polynomials, products thereof) can produce a uniruled, Fano, or otherwise supported intermediate variety that sees ξ:
every intermediate variety seeing ξ dominates X, and dominant images of uniruled varieties are uniruled while X is not
(p_g(X) = 1); more precisely a cover of X cannot be uniruled either, so the support problem is merely transported to a
cover of X with the same twisted sector (plus more).

What is NOT excluded by Theorem 1 (and is therefore the only place covers can help): covers Ẑ → P^4 (or → any base
Y ≠ W) whose Galois group does not contain the W-tower as a normal sub-tower — i.e. subgroups K ⊄ M.  These are exactly
the Galois closures of the NON-Galois maps X → P^4 of degree 6 (F even in a linear form, cubic in its square; or F
quadratic in Z = X2³) analysed in Part 2, where the class ξ̂ generates a genuinely non-abelian representation (the
3-dimensional reflection representation of B3 for rcr and qcr, the "class irrep" of the order-108 group for ctr) and
K ⊄ M can have fixed vectors.  Theorem 1 says the averaging-ideal test only has to be run on subgroups not contained in
M (and on K ⊂ N, which give covers of X): this is a large reduction of the enumeration and explains the shape of the
lenses' tables (the serving subgroups of ctr's 128 are either inside N·⟨twists⟩ or contain σ-conjugates).

## 1.5 Changing the branch section inside |O(B)| (why.md's "correspondence involving the branch locus")

Two natural readings, both negative:
(i) Deform B to B' ∈ |O_W(B)| and compare X with X' (both cyclic covers of W): Corollary 1 — the two twisted sectors
    H^4(X)_χ and H^4(X')_χ' are different local systems over the parameter space of B'; a correspondence between them
    through Z = X ×_W X' is zero on χ.  Over the FAMILY (B' moving in |O_W(B)|), the twisted sector forms a variation
    whose flat closure we compute in Part 3 (Def·ξ); a support for ξ at B must be an element of a flat family of
    supports only if it deforms, and then it carries the flat closure (quantitative budget §1.3(c)).
(ii) Use the branch threefold B itself as a correspondence: B ≅ X ∩ {x_n = 0} is a hyperplane section of X (fixed by G
    pointwise), and its H^3 (a weight-3 structure, level 3 for a sextic threefold: h^{3,0}(B) = 5) injects via Gysin
    into H^5 — the wrong degree; via the residue sequence its H^2 maps to H^1(W, Ω^3_W(log B)(−1)) only through
    H^1(B, Ω^2_B(−1)), which is the (2,1)-part of B twisted by O(−1) — nonzero in general, and this IS the twisted (3,1)
    sector's "boundary part"; but the classes of B are G-invariant (G acts trivially on B), so their Gysin images in H^4(X)
    lie in the invariant sector.  The branch locus contributes to the χ-sector only through log-forms with poles along B,
    never through algebraic cycles on B.

## 1.6 Summary of Part 1
* EV identifies the lane sectors as twisted cohomology H^q(W, Ω^p_W(log B) ⊗ O_W(−1)) of the Fano quotient W; the
  companion mechanisms live entirely in the untwisted summand.
* A support must be a G-invariant divisor E = π^{−1}(E_W) whose equivariant resolution has twisted (2,0)-forms; smooth
  (or rational-singularity) E_W never has them; they arise only from ramification of the cover along singular loci of
  B ∩ E_W — the hyperelliptic-slice (n = 2) and j = 0 elliptic-surface (n = 3) mechanisms are the two local models.
* Theorem 1 (rigidity): every cover of W containing X, in particular every fibre product X ×_W X' and every abelian
  tower over W, sees ξ only through covers of X.  Changing the branch section or adjoining roots over W cannot transfer
  support.  Only Galois closures over a different base (Part 2) and genuinely new divisor geometry (Part 3) remain.
