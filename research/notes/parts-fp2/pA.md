## Appendix A. Computations (this session and the crashed 2026-11-01 session), reproducible; nothing in the repository was modified

Engine: `~/devel/sextics-work/hodge/{core.py,fixed_d.py}` (engine v2 — the F·N(D) sign bug was fixed at 17:25 local; every result below is tagged `engine: 2` or was computed independently). Logs: `/tmp/hodge/results/four-plus-two-residual/log.jsonl`, `e1kernel.jsonl`, `e1kernel-basis-d6-par0.json` (20 polys), `e1kernel-basis-d7-par0.json` (63 polys), `batch1.out`, `cone_vertex.out`; scripts `h.py` (helpers incl. `image_sector`), `batch1.py`, `batch2.py` (cubic cones — not run to completion), `cone_vertex.py`, `e1kernel.py`. New scripts: `/tmp/hodge/work/fp2/sym.py` (symmetry group, characters, sector dims, ages), `prym.py` (a_p table), `lpoly.py` (L-polynomials over F_{p^k}).

A.1 Fixed-denominator tests (p = 2147483647 for consistency; `e1` = A·D ∈ J(F) solvable):

| D | d | E1 | full system | rank | nullity | unknowns |
|---|---|---|---|---|---|---|
| 1 | 0 | fails (A ∉ J) | — | 18 | 0 | 18 |
| X0²+2Y0², X0²−2Y0², X0²+Y0², X0X1+Y0Y1 | 2 | fail | — | 168 | 0 | — |
| X0X1Y0Y1 | 4 | fails | — | 747 | 9 | — |
| K = (5/7)X0³X1Y0Y1+(10/7)X0X1Y0³Y1 | 6 | fails | — | 2214 | 162 | — |
| f(X0,X1) | 6 | fails | — | 2214 | 162 | — |
| g(Y0,Y1) | 6 | fails | — | 2214 | 162 | — |
| **B(X2,Y2)** | 6 | passes | **inconsistent** | 2299 | 303 | 2602 (rows 15088) |
| **H(X0,X1,Y0,Y1)** | 6 | passes | **inconsistent** | 2299 | 303 | 2602 (rows 16370) |
| Y0³+Y1³+ℓ²Y0 (ℓ = X0 − tX1, f(t,1)=0 mod p; Y-odd cubic cone, vertex plane through a point of S_H) | 3 | — | inconsistent (mod p) | 351 | 27 | 378 |

(The remaining `batch1.py` entries — f+g, F−K, pure powers, H−B, f+B, g+B, the six partials ∂_iF, five radial derivatives — and the `batch2.py` cubic cones and the other `cone_vertex.py` sextic cones were queued but the session crashed before they ran; the log contains only the rows above.)

A.2 Colon-space dimensions dim{D Y-even, deg d : A·D ∈ J(F)} (mod p): d=1:0 (3 monomials), 2:0 (12), 3:0 (28), 4:0 (66), 5:3 (=J_5; 126 monomials), 6:38 (J_6 = 18, quotient **20**; 236 monomials), 7:126 (J_7 = 63, quotient **63**; 396 monomials). Timings 0.1 s … 452 s.

A.3 Symmetry (exhaustive over (Z/12)^6, 288 affine solutions, 24 projective): generators and characters as in §1.3–1.4; sector dims (32, 80, 32) in degrees (6, 12, 18) of the Fermat ring with Y-odd parity and the χ_A character; ages as in §1.5.

A.4 Point counts and L-polynomials of C_B: w³ = 1 + (15/7)t² + (10/7)t⁴ + t⁶ (plus 3 or 1 points at infinity according to p ≡ 1 or 2 mod 3) and E: w³ = 1 + (15/7)v + (10/7)v² + v³ (same infinity rule; the earlier `prym.py` table under-counted E by 2 for p ≡ 1 mod 3 — the `lpoly.py` counts are the corrected ones):
* p=13: #C(F_{13^k}) = 20, 144, 2414, 28512; #E = 21, 147, 2268, 28371; L_E = 13T²+7T+1; L_P = (13T²−2T+1)(169T⁴+13T³−12T²+T+1).
* p=31: #C = 24, 936, 30624, 922896; #E = 25, 975, 30100, 925275; L_E = 31T²−7T+1; L_P = 29791T⁶−961T⁵−589T⁴+194T³−19T²−T+1 irreducible/Q, = N_{Q(√−3)/Q} of T³ + (−1/62 − 3√−3/62)T² + (−13/961 + 4√−3/961)T + 7/1922 + 5√−3/1922 (times 31³); |α_i| = √31; no α_i/α_j is a root of unity of order ≤ 2000.
* p=37: #C = 30, 1368, 49854, 1880016; #E = 27, 1323, 50544, 1874691; L_E = 37T²−11T+1; L_P = 50653T⁶+4107T⁵+999T⁴−158T³+27T²+3T+1 irreducible/Q, two conjugate cubics over Q(√−3).
* For p ≡ 2 mod 3 (11,17,23,29,41,47,53,59,71,83,89,101,107,113): #C = #E = p+1 (a_p = 0, supersingular, as expected for Z[ζ3]-CM-type Jacobians).
* gcd(b, b') = 1 mod 13, 31, 37 (good reduction). The functional-equation-completed L_C has all |roots| = √p (consistency check of the genus-4 count).

A.5 System-size table of §5 computed from exact monomial counts (`core.monomials`).

## Appendix B. All URLs appearing in the thread and in this dossier
* Discussion: https://github.com/hodge-computer/sextics/discussions/120
* Instance record permalink (thread): https://github.com/hodge-computer/sextics/blob/faa9df2ec6406e4eacc135f677cf3711315d0e41/data/four-plus-two-residual/instance.rec
* Split-and-K3 chapter permalink (thread): https://github.com/hodge-computer/sextics/blob/faa9df2ec6406e4eacc135f677cf3711315d0e41/docs/03a-split-and-k3.md
* Related accepted PRs cited by companions: #128 (holdout, `fqlx`), #134 (library, `benediamond`), #149 (crossed-ternary, `1337isnot1337`), #150 (reflection-coupled-fano, `fqlx`), #151 (holdout-ternary-deformation, `hodgepodgeconjecture`), #158 (higher-residue library), #159 (split-ternary-q2, `fqlx`) — all under https://github.com/hodge-computer/sextics/pull/<n>
* External references used by the book sections cited: Kollár–Larsen https://arxiv.org/pdf/math/0701466 ; Bolognesi–Laterveer https://arxiv.org/html/2310.11981v2 ; Shioda–Katsura (split blowup) https://www.jstage.jst.go.jp/article/tmj1949/31/1/31_1_97/_pdf ; Diamond https://arxiv.org/html/2603.29157v2
