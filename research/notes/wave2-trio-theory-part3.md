# Part 3 — Decisive computations, receipts, honest verdicts

Receipts: /tmp/hodge/results/_cross/trio_sector_vhs.py (generalisation of the ctr construct lens's c5 test to lanes with
exchange symmetries), trio-qcr-sector-vhs.json, trio-rcr-sector-vhs.json (remote job trio-rcr-sectorvhs),
trio-qcr-galois.json, trio_qcr_ztsing.sing (remote job trio-qcr-ztsing).  Log: trio.log.jsonl.  Rigor convention: ranks
mod p = 1000003 and 998244353 (both agree); rank_p ≤ rank_Q and dim ker_p ≥ dim ker_Q are the directions used, and the
sector dimensions mod p equal the known complete-intersection values, so "rank_p = full" is rigorous over Q.

## 3.1 The sector-variation budget (what any support that deforms with F must carry)

Set-up (ctr construct lens part 6, now for all three lanes).  Let 𝓜 be the family of sextics with the lane's full visible
symmetry group Γ (ctr: 49 G-invariant monomials; rcr: Klein-invariant monomials; qcr: 26 Γ-orbit sums of G_diag-invariant
monomials).  The target's Γ-character sector S_A ⊂ H^4 is a sub-local system over 𝓜 with Hodge numbers
(0, h^{3,1}, h^{2,2}, h^{1,3}, 0).  Griffiths transversality: ∇_m ξ = m·ξ in R_12^χ, ∇_m∇_{m'} ξ = m m' ξ in R_18^χ.  Let 𝒲 be
the smallest sub-VHS of S_A defined over Q containing ξ.  Then 𝒲 ⊇ Def·ξ, Def²·ξ, their complex conjugates and
Def·(conjugates).  A divisor support for ξ that moves with F in 𝓜 (equivalently a support for very general members, by
the usual countable-union argument) supports 𝒲, hence its resolution carries at least dim 𝒲^{3,1} twisted (2,0)-forms.

| lane | n_Def | S_A = (h^{3,1}, h^{2,2}, h^{1,3}) | dim Def·ξ in (2,2) | dim Def²·ξ in (1,3) | K = {v ∈ S_A^{3,1}: Def·v = 0} | dim Def·S_A^{3,1} in (2,2) | 𝒲 | budget Σ h^{2,0}_tw(Ẽ) ≥ |
|---|---|---|---|---|---|---|---|---|
| ctr | 49 | (41, 119, 26) (not self-conjugate; S̄_A = (26,119,41)) | 41 | 26 = all | 0 | 119 = all | S_A ⊕ S̄_A | 41 + 26 = 67 (construct lens) |
| qcr | 26 | (22, 51, 22) (self-conjugate, ε real) | 22 | 22 = all | 0 | 51 = all | S_A | **22** (NEW) |
| rcr | see 3.2 | (126, 492, 126) (self-conjugate) | pending | pending | pending | pending | pending | pending |

qcr reasoning: Def²·ξ fills S_A^{1,3}, so 𝒲 ⊇ S_A^{1,3}; conjugating (𝒲 is defined over Q, the sector is real) 𝒲 ⊇ S_A^{3,1};
flatness then gives 𝒲 ⊇ Def·S_A^{3,1} = S_A^{2,2} (51 = all).  So 𝒲 = S_A: a support for the very general member of the
26-parameter Γ-family must carry the entire ε-sector, in particular 22 twisted (2,0)-forms.  (Def·ξ alone spans only
22 of the 51 (2,2)-classes — the first-order variation of ξ does not fill the sector, the second-order one does.)
Consequence for the lane's open degree-7/8 problems: a denominator D whose divisor E_W ⊂ W = X/s0 has only
rational-singularity-type degenerations (no twisted forms, §1.3(a)) can never work; a working D for F itself, if it is
not special to F, must produce ≥ 22 independent twisted 2-forms, e.g. ≥ 22 "slice" parameters in the hyperelliptic
mechanism — far above anything degree ≤ 8 offers (a degree-d hypersurface section of W carries ≤ O(d) such parameters).
This bound is for supports moving with F; special supports at the given F are not excluded.

## 3.2 rcr budget (remote job trio-rcr-sectorvhs)
(filled in below when the job returns; the remote is saturated — 356 cores allocated against the 200-core gate.)

## 3.3 qcr: klt check of the Galois closure (remote job trio-qcr-ztsing)
Z̃ = {e1(z²) = −a, e2(z²) = b, e3(z²) = −c} ⊂ P^7, CI (2,4,6), 8 variables; Singular computes the dimension and degree of
the singular locus (ideal + 3×3 minors of the Jacobian, p = 32003).  Expected: isolated or empty singular locus (then Z̃
has at worst isolated CI singularities; the quotients by B3-subgroups are klt wherever Z̃ is, and the BDPP step of Part 2
is rigorous there).  (result below)

## 3.4 The smallest decisive computations per lane for the most promising route (designs; not all run)
* ctr — route: isotrivial j = 0 elliptic surfaces over triple-point curves (EV: ramification of the μ3-cover along a
  curve where B ∩ E_W has ordinary triple points).  Budget 67 (part 6).  Decisive computation (construct lens item 14,
  unchanged by this note): the Fano variety 𝓕 of weighted lines on W and the twisted irregularity h^{1,0}_tw of the surface
  of lines tangent to B; Gröbner over F_p, ~hours on 16 cores; not run here.  A cheaper first test (new, from §1.3):
  for the 120 degree-9 plane families already enumerated (part 4), count the twisted 2-forms of the exceptional elliptic
  surfaces via the genus of the triple-point curve Γ (h^{2,0}_tw = h^{1,0}(Γ) = g(Γ)) — the sector-variation budget 67
  requires Σ g(Γ_i) ≥ 67, which no degree-9 family reaches (each Γ is a plane curve of degree ≤ 9, genus ≤ 28, and the
  families have ≤ 2 such curves): this re-derives the part-6 negative without the E2 tests.
* rcr — route: twisted hyperelliptic slices (M5; e = 2 conics lifting to genus-1 curves).  Budget: see 3.2 (if 126, the
  slice surface S must have h^0(Ω^1_S ⊗ λ) ≥ 126).  Decisive computation: the Hodge bundle λ = f_*ω on the surface S of
  conics with even contact to B through a fixed point / a pencil of hyperplane sections, and h^0(S, Ω^1_S ⊗ λ) by
  Riemann–Roch + vanishing; needs the explicit conic family of discussion #154 — not run here.  Status of the construct
  lens's gm_rank job (pulled 2026-11-02 04:31 UTC): the pencil-fibre sector dims (3, 90, 90, 3) reproduced at degrees
  7/13/19, final ∇-rank not yet written.
* qcr — route: the same hyperelliptic-slice mechanism on W = X/s0 (B = image of X ∩ {X0 = Y0}), or the fibration
  X/⟨σ0,σ1,σ2⟩ → P^1 × P^1 whose fibres are the double covers of the companion's cubic surfaces (χ(O) = 3, κ = 1, carrying
  the s0-odd 2-forms): supports of the form f^{−1}(C) for curves C in the base need h^0(C, (f_*ω_{E/C})^{odd} ⊗ ω_C) ≥ 22
  (budget), i.e. curves C of genus ≳ 11 or high degree in the (passive-ratio) base.  Decisive computation: the twisted
  Hodge bundle (rank 2, by χ(O) = 3 ⇒ p_g = 2 for q = 0) of the elliptic-surface fibration over P^1 × P^1 and its degree
  on curves — a single Riemann–Roch on the (6,6)-type discriminant; not run here.

## 3.5 Verdicts (per lane)
(filled in at the end)

### 3.5.1 cyclic-triple-residual (ctr) — verdict: NO ROUTE FOUND; theory closes all cover-type transfers
* EV: the target sector is the τ-invariant, G-sector part of H^1(W, Ω^3_W(log B)(−1)), W ⊂ P(1^5,3) the Fano quotient,
  B = {Z = 0} ≅ X ∩ {X2 = 0}.  The (4,0)-form shares the character; the sector S_A = (41, 119, 26) has no (4,0) part only
  because of the extra τ, g3 symmetries.
* Transfer by covers: Theorem 1 (rigidity) + the construct lens's exact enumeration of the order-108 closure: no uniruled
  quotient reaches ξ; fibre products X ×_W X' and all abelian towers over W are excluded without computation.
* Divisor supports must carry 67 twisted 2-forms for very general members of the 49-parameter G-family (construct lens
  part 6, reproduced in the table of §3.1); only the j = 0 elliptic-surface mechanism produces twisted forms, and it
  needs ≥ 7 triple-point curves of large genus ⇒ deg D ≳ 12 (dossier/lens estimate), outside every exact search run so far.
* Decisive next computation: Fano variety of lines of W and the twisted irregularity of the tangent-line surface (hours,
  16 cores); only worth running if someone is prepared to attack a degree ≳ 12 parametric family afterwards.

### 3.5.2 quadric-contact-residual (qcr) — verdict: NO ROUTE FOUND; the lane is structurally a second "reflection-coupled"
* NEW: the target is anti-invariant under all eight reflections of Γ⁻; for each of the four real ones F is a cubic in the
  square of the anti-invariant linear form, with Galois closure B3 (order 48) — exactly the rcr structure, four times.
  The averaging-ideal test over the closure (Part 2, formula (★)) gives no uniruled reaching quotient; the only K ≡ 0
  reaching quotients are X/H for residue-preserving H (κ = 0 by Reid–Tai, fqlx).  klt of Z̃: remote job pending.
* EV: target sector = ε-part of H^1(W, Ω^3_W(log B)(−1)) for W = X/s0 ⊂ P(1^5,2), B = image of X ∩ {X0 = Y0}; the companion's
  cubic-surface fibration lives in the untwisted summand of W' = X/⟨s0,σ0⟩.
* NEW budget: the Q-flat closure of ξ over the 26-parameter Γ-family is the whole ε-sector (22, 51, 22): a support moving
  with F must carry ≥ 22 twisted (2,0)-forms.  This rules out, for very general members, every denominator whose divisor
  on W has only mild singularities, and in particular every structured degree ≤ 8 family tested by the exclusion lens
  (consistent with their uniform negatives); at the specific F only special supports remain.
* Decisive next computation (not run): the twisted rank-2 Hodge bundle of the elliptic-surface fibration
  X/⟨σ0,σ1,σ2⟩ → P^1 × P^1 (fibres = double covers of the companion's cubic surfaces, χ(O) = 3, κ = 1) and the smallest
  base curve C with h^0(C, λ^{odd} ⊗ ω_C) ≥ 22; the corresponding f^{−1}(C) is the first candidate support family.
