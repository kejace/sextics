# Aggregator cross-checks of the degree-7 Hessian relaxations (receipts)

STATUS: 02:44 UTC — reflection-coupled-residual: the three closed Klein characters (1,1), (1,0), (0,1) re-derived with own code at p = 998244353 and found identical to the excl lens (quadric counts and colengths); (0,0) structurally identical (45/3/42 unknowns, 274 quadrics), Gröbner running locally. quadric-contact-residual: **all 28 degree-7 parity-class pairs re-derived, all m-primary** — the lens's items 3–4 (every pair of the eight (σ0, σ2, τ_Y) classes excluded ⇒ a degree-7 denominator needs character support on ≥ 3 classes) are confirmed by an independent implementation at a third prime.

## Method (common to both lanes)
E1 kernel mod J vanishes at d = 7 (both lanes, several primes), so a degree-7 denominator of a given character/class is `D = w(F)` with `w` a quadratic vector field whose components carry the matching characters. The certificate forces `A·Hess_F(w,w) ∈ J_14`; since the E1 kernel at d = 8 is one-dimensional (spanned by K8 resp. D8*, both in the trivial class), this is `Hess_F(w,w) ∈ J_8 ⊕ ℂ·K8` (resp. `⊕ ℂ·D8*`) — the Y-odd part of `Hess_F(w,w)` (cross terms of classes of odd Y-parity) must lie in `J_8^-` exactly. The Euler gauge `w ↦ w + hE` (`h` linear of the character) is removed by dropping the coefficient of `z_k²` in `w_k` for each `z_k` of that character. The quadrics are the coordinates of the normal forms in the engine's standard-monomial basis (`jacobian._Reducer`, both parities), after quotienting the even part by the K8/D8* direction. `msolve -g 1` (grevlex) gives the leading ideal; a pure power of every variable means the ideal is m-primary mod p, hence (Hilbert function over ℚ ≤ mod p, with the Jacobian ranks at their complete-intersection values) only the trivial solution over ℂ. Colengths (number of standard monomials) are order-independent and equal the lenses' "mult".

K8 was taken from the lane's rational expression (`/tmp/hodge/results/reflection-coupled-residual/k8_expr.txt`, the CRT-reconstructed class the lens verified exactly); D8* was recomputed from scratch as the exact quotient `(F_{X0}² − F_{Y0}²)/(X0² − Y0²)` (dossier identity, division verified by multiplication; 25 terms mod p).

## reflection-coupled-residual (`rcr7.py <chi>`; Klein characters `χ = (deg_Y2 mod 2, deg_Y0 + deg_Y1 mod 2)`)

| χ | dim V | dim T (Euler) | unknowns | independent quadrics | leading monomials | m-primary | colength | lens (p = 2147483587) |
|---|---|---|---|---|---|---|---|---|
| (1,1) | 18 | 0 | 18 | 117 | 135 | True | 85 | 117 quadrics, GB 135, mult 85, V = {0} — **agrees** |
| (1,0) | 23 | 1 | 22 | 132 | 272 | True | 239 | 132 quadrics, GB 301, mult 239, V = {0} — **agrees** |
| (0,1) | 40 | 2 | 38 | 281 | 2897 | True | 1109 | 281 quadrics, GB 2987, mult 1109, V = {0} — **agrees** |
| (0,0) | 45 | 3 | 42 | 274 | — | running (msolve timed out at 900 s on the saturated local box; relaunched with a 90-minute cap) | — | 274 quadrics, GB running |

## quadric-contact-residual (`qcr7.py <s> <t>`; classes `(b0,b2,y) = (e0+e3, e2+e5, e3+e4+e5) mod 2` of σ0, σ2, τ_Y; h1-characters mix freely)

| pair | type | unknowns | independent quadrics | leading monomials | m-primary | max pure power | seconds |
|---|---|---|---|---|---|---|---|
| {000,001} | edge | 32 | 217 | 1799 | True | 4 | 84.1 |
| {000,010} | edge | 32 | 220 | 1874 | True | 4 | 136.1 |
| {000,011} | face diag. | 32 | 220 | 1770 | True | 4 | 24.5 |
| {000,100} | edge | 32 | 235 | 1689 | True | 4 | 23.5 |
| {000,101} | face diag. | 32 | 235 | 1591 | True | 4 | 32.9 |
| {000,110} | face diag. | 28 | 215 | 812 | True | 4 | 9.5 |
| {000,111} | space diag. | 28 | 215 | 829 | True | 4 | 5.6 |
| {001,010} | face diag. | 32 | 220 | 1909 | True | 4 | 40.6 |
| {001,011} | edge | 32 | 220 | 1804 | True | 4 | 36.7 |
| {001,100} | face diag. | 32 | 235 | 1733 | True | 4 | 32.7 |
| {001,101} | edge | 32 | 235 | 1634 | True | 4 | 29.2 |
| {001,110} | space diag. | 28 | 215 | 859 | True | 4 | 4.3 |
| {001,111} | face diag. | 28 | 215 | 858 | True | 4 | 4.3 |
| {010,011} | edge | 32 | 208 | 2045 | True | 4 | 43.5 |
| {010,100} | face diag. | 32 | 229 | 1862 | True | 4 | 54.8 |
| {010,101} | space diag. | 32 | 229 | 1797 | True | 4 | 60.6 |
| {010,110} | edge | 28 | 209 | 889 | True | 4 | 3.0 |
| {010,111} | face diag. | 28 | 209 | 898 | True | 4 | 11.0 |
| {011,100} | space diag. | 32 | 229 | 1939 | True | 4 | 70.5 |
| {011,101} | face diag. | 32 | 229 | 1820 | True | 4 | 69.0 |
| {011,110} | face diag. | 28 | 209 | 896 | True | 4 | 2.9 |
| {011,111} | edge | 28 | 209 | 876 | True | 4 | 3.4 |
| {100,101} | edge | 32 | 233 | 1733 | True | 4 | 12.6 |
| {100,110} | edge | 28 | 221 | 757 | True | 4 | 2.3 |
| {100,111} | face diag. | 28 | 221 | 763 | True | 4 | 2.2 |
| {101,110} | face diag. | 28 | 221 | 788 | True | 4 | 2.2 |
| {101,111} | edge | 28 | 221 | 785 | True | 4 | 8.3 |
| {110,111} | edge | 24 | 174 | 434 | True | 3 | 2.0 |

(raw results in `/tmp/hodge/results/_aggregator/qcr7/qcr7_results.jsonl`, msolve inputs/outputs `qcr7_<s>_<t>.ms/.g1.out`; the lens quotes 227–258 "quadrics" for edges and 227–246 for diagonals — presumably raw coordinate counts — against 215–235 independent ones here; the conclusions coincide.)

## Verification-log entries
`rcr_d7_char11_relaxation`, `rcr_d7_char10_relaxation`, `rcr_d7_char01_relaxation`, `rcr_d7_char00_relaxation_structure`, `qcr_d7_pair_<s>_<t>` in `/tmp/hodge/notes/wave2-verification.log.jsonl`.

Replay: `cd /tmp/hodge/results/_aggregator/rcr7 && python3 rcr7.py 11 998244353 600`; `cd ../qcr7 && python3 qcr7.py 000 111 998244353 300` (each builds the `.ms`, runs msolve, prints a JSON line).
