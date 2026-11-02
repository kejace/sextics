# Wave 2 — `reflection-coupled-residual` — lens EXCLUSION / RIGOR

STATUS: no witness. RIGOROUS (over ℂ): every certificate denominator has degree ≥ 7 (d ≤ 4 by E1 in all four Klein characters; d = 5, 6 by E1 ⇒ D ∈ J_d plus the zero-dimensional Hessian relaxation on W = V/T, mixed parity included; E1 re-verified at a third prime); degree 7 is impossible for single-character denominators in characters (0,1), (1,0), (1,1) and for the mixed pair (1,0)+(1,1); the unique non-Jacobian degree-8 E1 class K8 is certified not a denominator; 7/7 tested degree-9 kernel members certified negative. OPEN: degree-7 character (0,0) block (8-hour GB run in progress) and the other mixed pairs/triples; degree-8 families (remote hess8); screens continue (hand-off section).

Started 2026-11-01 ~19:15 local. Work dir `/tmp/hodge/results/reflection-coupled-residual/wave2-excl/`,
log `/tmp/hodge/results/reflection-coupled-residual/wave2-excl.log.jsonl` (one JSON line per experiment).
Remote jobs named `excl-reflection-coupled-residual-*`.

## Operational note
Remote receipts are merged with `wave2-excl/merge_remote.py` after every `rjob pull`; `rjob submit` rsyncs the local lane dir
over the remote, which clobbered the remote `e1char.jsonl` once (degree-9 receipts for three characters lost; rerun queued with
W2LIB_SUFFIX so remote and local receipt files never share a name again).

## Plan
(a) re-verify wave-1/dossier rigorous claims with a different prime (p = 2147483587 / 2147483579; wave 1 used 2147483647, 2147483629; dossier used 32003), with explicit good-reduction (Hilbert function) checks;
(b) rigorous lower bound on deg D via E1 by Klein character + Hessian relaxation of E2 on the Jacobian family;
(c) non-Jacobian E1-kernel classes and structured families, degrees 7-12, character sectors, full E1+E2 (symmetric solver);
(d) packet; (e) draft post.

## Summary of the rigorous frontier (updated as results arrive)
* deg D ≤ 4: impossible (E1, all characters) — re-verified at a third prime.
* deg D = 5, 6: impossible (E1 ⇒ D ∈ J_d; Hessian relaxation zero-dimensional on W), mixed parity included.
* deg D = 7: impossible for single-character D in characters (0,1), (1,0), (1,1) [GB]; (1,0),(1,1) also by the Y2-parity
  reduction from degree 6; character (0,0) and mixed characters: computations running.
* deg D = 8: the E1-admissible space is J_8 ⊕ ℂK8 (all characters); K8 itself certified excluded; (1,0),(1,1) sectors reduce
  to degree-7 (0,·) sectors; (0,0),(0,1) sectors and the K8 family: GB computations running (remote hess8, local chain2).
* deg D ≥ 9: E1 kernels are dimension-forced (21/0/21/0 at d = 9, 186/126/136/108 at d = 10); only fixed-D tests of
  representatives/structured members (screens) and the reduction lemma apply.

## Results so far

### R1. Degree 5 (Hessian relaxation, local, p = 2147483587) — DONE
For D = Σ c_i F_i (c_i constants; this is the whole E1-admissible space at d = 5 once job e1char confirms
"E1 ⇒ D ∈ J_5" in all four characters), the E2-mod-J obstruction is ob(c) = −[A Σ c_i c_j F_ij] ∈ (S/J)_10.
The Hessian has exactly one vanishing off-diagonal entry, F_{X2 Y2} = 0; the other 20 classes [A F_ij] are
linearly independent in (S/J)_10 (rank 20 in each of the four target-character blocks: sizes 6/10/2/2;
good reduction verified: rank J_10^ψ = dim S_10^ψ − Fermat count for all ψ). Hence c_i c_j = 0 for all 20 pairs,
in particular c_i² = 0, so c = 0. **No denominator of degree 5 (any character, mixed parity included).**
Receipts: `wave2-excl/hess.jsonl` (exp phi-kernel/phi-verdict, d=5).

### R2. Degree 6 (Hessian relaxation, local, p = 2147483587) — DONE: deg-6 excluded
The linear criterion fails: ker Φ_6 has dimension 75 (per target character 17/26/13/19; trivial part only 1): the
Euler relations Σ_j z_j F_ij = 5F_i put every symmetrised tuple z⊗a + a⊗z (a ∈ S_1^6, 36 dims) into ker Φ. So the
genuine quadratic system is needed: on W = V/T (V = S_1^6 = 36, T = C·z, 35 unknowns, all four characters jointly,
so mixed-parity D included) the obstruction gives 963 nonzero quadrics of which 345 are independent (< 630 = dim S_2(W)),
so the degree-2 part alone does not decide. msolve (grevlex GB mod p, 709 elements, 1.8 s) + Singular: **dim 0, mult 1060,
minimal associated primes = {(w0,…,w34)}**, Hilbert function HF(0..7) = 1, 35, 285, 477, 225, 36, 1, **0**.
Hence I_p ⊇ S_7(W); since HF_Q(k) ≤ HF_p(k) (rank argument, good reduction in degree 12 verified), I_Q ⊇ S_7(W),
so V(ob) = T over ℂ: every g with ob(g) = 0 has g(F) ∈ (F). **No denominator of degree 6** (any character, mixed
parity included), given E1 ⇒ D ∈ J_6 in all characters (dossier mod 32003; re-checked by job e1char).
Files `wave2-excl/hess_d6_W00011011_p2147483587.{ms,json,lt,gb}`, `hess_d6_singular_minass.txt`.
Proof that B(T,·)=0 exactly: B(zG, g) = [A G Σ_j g_j Σ_i z_i F_ij] = 5[A G g(F)] = 0 since g(F) ∈ J.


### R3. Degree 7, pure Klein characters (Hessian relaxation + msolve/Singular, local, p = 2147483587) — 3 of 4 DONE
D ∈ J_7 (E1) of a single character χ ⇒ g ∈ V_χ, W_χ = V_χ/T_χ; the quadrics live in (S/J)_14^{(0,1)} (dim 426).
| χ | dim V | dim T | dim W | indep. quadrics | GB size | dim | mult | HF | verdict |
|---|---|---|---|---|---|---|---|---|---|
| (1,1) | 18 | 0 | 18 | 117 | 135 | 0 | 85 | 1,18,54,12,0 | V = {0} |
| (1,0) | 23 | 1 | 22 | 132 | 301 | 0 | 239 | 1,22,121,83,12,0 | V = {0} |
| (0,1) | 40 | 2 | 38 | 281 | 2987 | 0 | 1109 | 1,38,460,556,54,0 | V = {0} |
| (0,0) | 45 | 3 | 42 | 274 | — | — | — | — | msolve did not finish in 850 s / 3900 s (4 threads); 8-hour run started 21:02 |
Hence (modulo the (0,0) block) **no degree-7 denominator with a single Klein character**. Mixed-character degree-7
denominators: joint systems on ⊕W_χ — pairs (10,11): 40 unknowns/300 independent quadrics, (01,11): 56/470, (01,10): 60/479,
(00,11): 60/468, (00,10): 64/463, (00,01): 80/…; triples 78–102 unknowns; full 120 — exported and run by msolve with bounded
caps (`wave2-excl/chain2.sh`; receipts `gbcheck.jsonl`). **Pair (1,0)+(1,1): dim 0, mult 4399, HF 1,40,520,1543,1695,567,33,0
(msolve 166 s) — closed** (consistent with the Y2-parity reduction, R9).

### R4. The degree-8 E1-kernel class K8 (remote job excl-rcr-k8job) — DONE: K8 rigorously excluded
* K8 ∉ J(F): verified mod 2147483587 and 2147483579 (rank 168 vs 169); a rational polynomial outside J mod p is outside J over ℚ.
* A·K8 ∈ J(F) **exactly over ℚ**: the character-restricted E1 system (3276 × 3388, rank 2850) was solved mod two primes,
  CRT-lifted (62 bits, max denominator 18 bits) and the identity Σ N_i ∂_iF = A·K8 verified exactly (`wave2-excl/K8_E1_exact.json`,
  N-term counts 28/65/42/56/23/0). So the wave-1 class is a genuine rational E1-admissible non-Jacobian denominator.
* Fixed-D certificate test (rigorfix protocol): Klein-invariant-sector system 22016 × 4256 (of 7497 unrestricted columns,
  stabiliser order 4): **inconsistent mod p**, rank_p M = 3780, nullity 476; the 482 a-priori rational kernel vectors
  (zG, 0) and (F m e_i, D ∂_i m − m ∂_i D) have rank 476 = nullity, hence rank_ℚ M = rank_p M and the inconsistency holds
  over ℚ (and ℂ). Second prime agrees; the full unrestricted system (40301 × 7497) is also inconsistent mod p (rank 6703).
  **K8 itself is not a certificate denominator.** (427 s.) The affine family K8 + J_8 is treated in R6.
* Geometry of Z = X ∩ V(K8): the Jacobian ideal (F, K8, 2×2 minors) was not finished by Singular (30 min) nor msolve (50 min);
  whether Z is smooth remains open (not needed for the exclusion).

### R5. E1 by Klein character, independent re-verification (remote job excl-rcr-e1char, p = 2147483587; reps lifted with 2147483579)
Kernel of [A]: (S/J)_d^χ → (S/J)_{d+6}^{χ+(0,1)}; "ker" = dim{D ∈ S_d^χ : AD ∈ J}, "mod J" = its dimension modulo J_d^χ.
Good reduction (rank J_n^χ = dim S_n^χ − Fermat count) holds in every row, so each modular kernel dimension bounds the
rational one from above and the zeros are theorems.
| d | χ=(0,0): dim S, (S/J), ker, mod J | (0,1) | (1,0) | (1,1) |
|---|---|---|---|---|
| 1–4 | ker = 0 in all characters (dims 3/2/1/0, 10/6/3/2, 22/18/10/6, 48/38/22/18) | | | |
| 5 | 88, 85, 3, **0** | 78, 76, 2, **0** | 48, 47, 1, **0** | 38, 38, 0, **0** |
| 6 | 158, 144, 14, **0** | 138, 126, 12, **0** | 88, 82, 6, **0** | 78, 74, 4, **0** |
| 7 | 258, 213, 45, **0** | 238, 198, 40, **0** | 158, 135, 23, **0** | 138, 120, 18, **0** |
| 8 | 413, 301, 113, **1** (= K8, two primes, rational lift agrees with wave 1) | 378, 274, 104, **0** | 258, 194, 64, **0** | 238, 182, 56, **0** |
| 9 | 623, 381, 263, **21** | 360, 228, **0** | 263, 171, **21** | 242, 136, **0** |
| 10 | 460, 650, **186** | 426, 568, **126** | 318, 441, **136** | 302, 394, **108** |
(d = 9 rows re-run with suffix `.r2`; d ≥ 9 kernels are dimension-forced: (S/J)_9 → (S/J)_15 and (S/J)_10 → (S/J)_16 have
target dimensions 360/381/242/263 and 274/301/182/194 by the socle duality (S/J)_n ≅ (S/J)_{24−n}^∨.)
Agreement with wave 1 (Y-even sector, primes 2^31−1, 2^31−19) and the dossier (all characters mod 32003, d ≤ 6): complete.
New: the Y-odd characters at d = 7, 8 have no non-Jacobian E1-admissible class either, so at degree 8 the whole
E1-admissible space is J_8 ⊕ ℂ·K8 (mixed characters allowed).

### R7. Structured-family E1 screen (remote job excl-rcr-struct, p = 2147483587, d = 7..12) — d = 7 done
48 monomial-defined shapes (cones with a line or plane vertex, multiples of coordinates/squares, polynomials in squares,
X-only/Y-only, bounded X/Y-degree, the Klein-invariant subring, (X0,Y0)-ideal powers, …), each intersected with the
E1-admissible space by character. Degree 7: only `multiple-of-X0` (36), `multiple-of-X0^2` (6), `multiple-of-X1` (36),
`multiple-of-X1^2` (6) contain E1-admissible elements, all inside J_7 (consistent with R5); no cone, binary/ternary form,
polynomial in squares, u = Y2²-polynomial, … of degree 7 passes E1. Degree 8: only `cone-line: independent of X0,Y2` (1 element) and `multiple-of-X0` (126) meet the E1-admissible space, both
inside J_8: **no structured shape contains a member of the K8 family** (extends wave-1's 55-shape test). Degrees 9–12: running
(survivors outside J are saved as `struct_d<d>_<shape>_chi<ab>.json` for fixed-D tests).

### R6. The degree-8 family K8 + J_8 (remote job excl-rcr-k8job; wave-1 `e2lin.py` completed)
ob(K8) ≠ 0 in (S/J)_16^odd (dim 468). Linearisation over the Y-even family D = K8 + Σ g_i F_i (168 parameters):
rank L = 156 (kernel = the 12 trivial directions zG), rank [L | ob(K8)] = 157 ⇒ **no first-order correction inside the
E1 family removes the obstruction** (687 s, p = 2147483647). This is a screen (the obstruction is quadratic in g); the
full quadratic system on the Klein-invariant sector (1 + 102 unknowns, 274 quadrics) is exported by `hessK.py` and run by
msolve in the local chain (receipt pending). Fixed-D test of K8 itself: certified inconsistent (R4).

### R8. Reduction lemma (local, Singular p = 32003) and the degree-9 structured survivors
All six coordinate sections X ∩ {z_i = 0} are smooth (Jacobian ideal zero-dimensional mod p ⇒ smooth over ℚ). A smooth
threefold hypersurface has H^{2,0} = 0, so its (3,1)-served space is zero; by the sum principle (served(D1D2) =
served(D1) + served(D2), book ch. 3 §1.1 "finite sums: union of supports") **served(z_i·G) = served(G)**: any denominator
divisible by a coordinate serves exactly what its cofactor serves. Consequences: every degree-7 denominator divisible by a
coordinate is excluded (degree 6 is excluded outright); the degree-9 structured E1 survivors `multiple-of-X0/X1/X2` (one
non-Jacobian class each, = X_i·K8) reduce to K8 (excluded, R4); the 21 non-Jacobian `multiple-of-Y2` classes of character
(1,0) at degree 9 reduce to degree-8 cofactors (Y2·G with A·Y2·G ∈ J but A·G ∉ J in general — their cofactors are not
E1-admissible, so they cannot serve A: a certificate with denominator Y2·G would give one with denominator G by the lemma,
contradicting E1 for G unless G ∈ J_8 ⊕ ℂK8).
Remaining degree-9 structured survivors outside J (receipts `structscreen.jsonl`): 11 classes with even Y0,Y1,Y2
exponents and 3 + 11 classes inside (X0,X1,Y0,Y1)^4, all among the 21 + 21 kernel classes being tested by `excl-rcr-ker9`.
Degree 10 (kernels 186/126/136/108): the structured shapes meet the non-Jacobian E1 space in large families
(multiples of a coordinate: 292–296 classes each, reducible to degree 9 by the lemma; Y-degree ≤ 4: 184; X-degree ≤ 4: 132;
even exponents in Y0,Y1,Y2: 179 in (0,0)); no cone with a line or plane vertex, no binary/ternary form and no polynomial
in squares is E1-admissible outside J at degrees 8–10. These degree-10 families are screened only (no fixed-D tests).

### R9. Y2-parity reduction (structural)
Monomials of character (1,0) or (1,1) have odd Y2-exponent, hence a pure-character denominator D of such a character is
D = Y2·G with G of character (0,0) resp. (0,1) and degree d − 1. Since X ∩ {Y2 = 0} is smooth (R8), served(D) = served(G).
Iterating, **every pure-character denominator of odd Y2-parity reduces to one of even Y2-parity and lower degree**.
Consequences: (i) the degree-7 exclusions for (1,0), (1,1) are implied by the degree-6 exclusion (independent cross-check of
the GB results of R3); (ii) the degree-8 sectors (1,0), (1,1) are excluded once degree-7 characters (0,0), (0,1) are;
(iii) the 21 non-Jacobian degree-9 classes of character (1,0) are Y2·G with G ∈ S_8^{(0,0)}: a certificate with such a
denominator would give one with denominator G, which then must be E1-admissible, i.e. G ∈ J_8 ⊕ ℂK8 — so this whole
sector reduces to the degree-8 Klein-invariant family (R4, R6, hess8). The surviving "hard" sectors are (0,0) and (0,1)
(even Y2-parity) and mixed-character denominators.

### R10. Positive control (Fermat sextic, local)
On X0⁶+…+Y2⁶ with the same A, where D0 = X0⁶+X1⁶+X2⁶ ∈ J_6 is a known certificate denominator (DiamondJM's control), the
identical degree-6 pipeline gives only 52 independent quadrics on the 35-dimensional W and a solution cone of dimension 19
(Hilbert function 1, 35, 578, 6078, … never vanishing). So the relaxation detects known certificates and the zero-dimensional
results for the coupled F are not artefacts of the pipeline (`wave2-excl/control-fermat/`).

### R11. Degree-9 kernel representatives (remote job excl-rcr-ker9) and the degree-11 cones (excl-rcr-struct)
* The three rationally reconstructed representatives of the 21-dimensional non-Jacobian E1 kernel of character (0,0) at
  degree 9: restricted systems 32214 × 6804, all **inconsistent mod p and certified over ℚ** (rank 6006, nullity 798 = a-priori
  kernel rank), second prime agrees (822/299/249 s on the saturated box). The other 18 need more primes to reconstruct
  (`e1recon.py`, queued); four random integer combinations of the three: also certified inconsistent (**7/7 negative**).
* Structured screen at degree 11: **cones with a line vertex become E1-admissible for the first time** — polynomials in
  (X0,X1,X2,Y2) (vertex the line Y0Y1): a 46-dimensional E1-admissible space with 45 classes outside J_11; polynomials in
  (X0,X1,X2,Y0): 5. These are the book's natural (3,1)-support candidates (cone over a degree-11 surface with p_g > 0).
  The mod-J relaxation is useless here ((S/J)_22^{(0,1)} has dimension 6), so only fixed-D tests (screens) apply;
  restricted systems at degree 11 have ≈ 14 000 unknowns (local screen `chain3_cones.sh` queued behind the GB runs).
* Degree 12: cones with a line vertex exist for all 15 vertex lines (61 classes for X–X and X–Y0/Y1 lines, 185 for lines
  containing Y2, 255 for the (Y0,Y1) line; Klein-invariant: 14 / 105 / 126) — dimension-forced, screens only.

### Jobs submitted (remote)
* `excl-rcr-e1char`: E1 kernels by Klein character, d = 1..10, p = 2147483587 (+ 2147483579 for reconstruction).
* `excl-rcr-hess`: Φ-kernel at d = 6 (two primes); Macaulay-3 rank of the quadratic system per character at d = 7;
  msolve exports of the pair/full systems at d = 7.
* `excl-rcr-k8job`: DONE (R4, R6).
* `excl-rcr-struct`: DONE — structured-family E1 screen d = 7..12 (R7, R8, R11; 348 receipts in `structscreen.jsonl`).
* `excl-rcr-misc` (submission was still rsyncing at hand-off — verify with `python3 -m hodge.rjob status`; resubmit `bash wave2-excl/job_misc.sh` if absent): control on the same F with the served companion target X0³Y2³ (E1 kernels d = 1..9, receipts `control-companion/`), then 8-prime rational reconstruction of the degree-9 kernels (0,0) and (1,0) (`e1recon.r4.jsonl`, `e1reps_d9_chi*_Qk.json`).
* `excl-rcr-hess8`: degree-8 pure-character Jacobian families (χ = 11, 10, 01) and the K8 family (1 + 102 unknowns): msolve GB checks (running).
* `excl-rcr-ker9`: degree-9 E1-kernel representatives (21 in (0,0), plus (1,0)) fixed-D tests + random combinations; rerun of the lost degree-9 receipts.


## Methods (rigor protocol, in one place)
1. **E1 by character.** AD ∈ J ⇔ A·D_χ ∈ J for each Klein component (J is Klein-graded). Kernel of [A] on (S/J)_d^χ mod p with the
   good-reduction check rank J_n^χ = dim S_n^χ − Fermat count in degrees d and d+6; then dim_ℚ ker ≤ dim_p ker.
2. **Obstruction.** For D = g(F) with N = A·g, B_D(N) ≡ −A gᵀHess(F)g (mod J), independent of the E1 solution (Koszul corrections
   lie in J; verified numerically), and B(zG,·) = 5[A G g(F)] = 0. Certificate ⇒ ob(g) = 0 in (S/J)_{2d}.
3. **Modular ⇒ rational for the quadratic systems.** The modular quadrics are reductions of rational ones in the ideal (the
   annihilator lattice of J_{2d} is saturated when reduction is good), so HF_ℚ(k) ≤ HF_p(k) for the ideal on W = V/T;
   dim 0 mod p ⇒ HF_p(k0) = 0 ⇒ HF_ℚ(k0) = 0 ⇒ V(ob) = T over ℂ ⇒ every such D lies in (F).
4. **Fixed-D inconsistency.** rank_p M ≤ rank_ℚ M always; an exact rational kernel of dimension equal to the modular nullity
   (vectors (zG,0), (F m e_i, D∂_im − m∂_iD)) forces equality, so modular inconsistency is rational inconsistency.
5. **Reductions.** Sum principle (book ch. 3 §1.1) + smoothness of X ∩ {z_i = 0} (checked mod p) ⇒ served(z_iG) = served(G);
   odd-Y2-parity characters are Y2-multiples.
6. **Controls.** Fermat sextic (certificates exist): the pipeline returns a 19-dimensional solution cone at degree 6.

## State at hand-off (21:10 local) — where later results land
Local (brunos), still running, receipts appended automatically:
* `wave2-excl/gbcheck.jsonl`: the degree-7 character-(0,0) block (8-hour msolve run, 6 threads, started 21:02); then chain4
  (`chain4.sh`): exports of the (00,01,10) triple and the full 120-unknown system, the degree-8 K8-family export
  (`hessK_d8_chi00_p*.ms`, 1 + 102 unknowns) and its 30-minute GB run (`gbcheck_d8K.out`), two 20-minute pair runs.
* `wave2-excl/conescreen.jsonl` + `rigorfix.jsonl`: degree-11 Klein-invariant cone families (`chain3_cones.sh`).
Remote (hodge-sextics, heavily loaded): `excl-rcr-hess8` (degree-8 pure sectors (1,1), (1,0), (0,1) and the K8 family; receipts
`gbcheck.r3.jsonl`, `hess.r3.jsonl`), `excl-rcr-ker9` (degree-9 kernel tests; `rigorfix.r2.jsonl`), `excl-rcr-misc`
(companion control; degree-9 reconstructions). Pull with `python3 -m hodge.rjob pull --lane reflection-coupled-residual`
and merge with `wave2-excl/merge_remote.py` (union of JSONL lines; never lets a local copy clobber remote receipts).
Interpretation rules: a `gb-hilbert` receipt with `dim: 0` closes its block rigorously (see Methods); `msolve TIMEOUT`
means undecided; a `rigorfix` receipt with `certified: true` is a rigorous exclusion of that single denominator.
