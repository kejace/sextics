STATUS: RUNNING — 2026-11-02 05:09 UTC — DECIDED NEGATIVE (exact): Aut_lin(generic (01)-fibre) = {1, s}, Aut_lin(X) = S3×⟨τ⟩, no second reflection; sweep over 5 values of μ and 2 primes agrees; still running: rank-≤2 locus on P^5 (wave-2 J1, bonus) and Singular-over-Q cross-check of the rank-1 loci.

# Wave 4 — diagonal-permutation — route `second-reflection` (M-a of wave3 §3.4)

Question (J2/J2f). Is there a hidden or created second reflection on the (01)-pencil fibres
X_mu = X ∩ {r = mu t} (r = X0 − X1, t = Y0 − Y1), beyond the known reflection s (t ↦ −t)?
Positive ⇒ book §4.6 descent (Lean class L-B). Prior work: wave2 construct lens (jobs
`construct-dp-j1j2`, `construct-dp-fibre`), numerical screens J1-num/J2-num (nothing found).

## Harvest of the wave-2 remote jobs (04:30 UTC)
| job | item | status |
|---|---|---|
| construct-dp-j1j2 | J2 chart v0 (ℓ_v ∣ ∂_v³F, F(v)≠0; 57 eqs, 27 unknowns, msolve F_p, 4 threads) | TIMED OUT after 7000 s (no result) |
| construct-dp-j1j2 | J2 chart v1 | running since 22:20 EST (same system size; expected to time out) |
| construct-dp-fibre | J2f chart 0 (fibre hidden reflection, mu = 3/7; 27 unknowns) | TIMED OUT after 7000 s; charts 1–4 never ran (uncaught exception) |
| construct-dp-fibre | hessrank1 (rank-1 Hessian locus of F on P^5, 6 charts) | CRASHED immediately: NameError `RANK` used before assignment in w2c_hessrank.py (never ran) |
| construct-dp-fibre | J1f (rank ≤ 2 Hessian locus of F_mu, 5 charts, p = 1073741827) | DONE: charts p=1, q=1, t=1 EMPTY; charts a=1 (deg 5), b=1 (deg 19) zero-dimensional ⇒ locus finite, contained in the line L = {p = q = t = 0} = span(e_a, e_b) (modular, mu = 3/7) |

## Results table (running)
| id | what | result | rigor |
|---|---|---|---|
| W4-1 | rank-1 Hessian locus of F_mu (2x2 minors), mu = 3/7, p = 1073741827, 5 charts, u·z_j−1 tests | exactly {e_a, e_b} (both reduced, ideal degree 1); charts p,q,t empty | exact over F_p-bar (msolve) |
| W4-2 | rank-1 Hessian locus of F on P^5, same p, 6 charts | exactly the 6 coordinate points (multiplicity 1 at X_i, 3 at Y_i) | exact over F_p-bar |
| W4-3 | smoothness of X_{3/7} and of X mod p (Jacobian ideal, all charts) | smooth | exact over F_p-bar ⇒ smooth over Q-bar |
| W4-4 | h^0(T) = dim{M ∈ gl_n : ξ_M G ∈ k·G} − 1 for G = F_{3/7} (n=5) and F (n=6) | 0 over Q and over F_p | exact (linear algebra) |
| W4-5 | Aut_lin(F_{3/7}) over F_p-bar: block-diagonal reduction (Lemma A) + 12 msolve systems in 9 unknowns | order 2 = {1, s}; all 12 systems: diag j=0 gives g_M ∈ {1, diag(1,1,−1)} (both verified), the other 11 empty | exact over F_p-bar; ⇒ exact over Q-bar by Lemma B |
| W4-6 | monomial automorphisms of F (exact Z/6 enumeration over Q(ζ6)); Aut_lin(F) ⊂ monomial by W4-2 + Lemma A' | S3 × ⟨τ⟩, order 12, NO reflection, no §2.7 partner coordinate | exact (over F_p-bar via W4-2; over Q-bar by Lemma B) |
| W4-7 | generic fibre: Aut_lin(X_{η̄}) ↪ Aut_lin(X_{3/7}) (Lemma B over Q̄[[μ−3/7]]) | Aut_lin(X_{η̄}) = {1, s}: no hidden/created second reflection | exact |

## Robustness sweep (same pipeline W4-1/3/4/5 per (μ, p); all exact over F_p-bar, msolve)
| μ | p | rank-1 locus of F_μ | X_μ smooth mod p | h^0(T) (Q, F_p) | |Aut_lin(F_μ)| over F_p-bar | elements |
|---|---|---|---|---|---|---|
| 3/7 | 1073741827 | {e_a, e_b} | yes | 0, 0 | 2 | 1, s |
| 3/7 | 1000000009 | {e_a, e_b} | yes | 0, 0 | 2 | 1, s |
| 2/9 | 1073741827 | {e_a, e_b} | yes | 0, 0 | 2 | 1, s |
| −5/3 | 1073741827 | {e_a, e_b} | yes | 0, 0 | 2 | 1, s |
| −5/3 | 1000000009 | {e_a, e_b} | yes | 0, 0 | 2 | 1, s |
| 7/3 | 1000000009 | {e_a, e_b} | yes | 0, 0 | 2 | 1, s |
| 1 | 1073741827 | {e_a, e_b} | yes | 0, 0 | 2 | 1, s |
| 0 (control: τ-stable member) | 1073741827 | LARGER (extra rank-1 points in charts q, t) — Lemma A reduction not applicable | yes | 0, 0 | ≥ 4 (block-diagonal part: 1, s, τ = (q,t,b ↦ −q,−t,−b), sτ) — τ, sτ have eigenvalue pattern (1^3,−1^2), not reflections | lower bound only |
| 3/7 | 2147483647 | — | — | — | INVALID: msolve 0.6.5 rejects p = 2^31−1 for zero-dimensional systems (core dump); only its [-1] answers are valid | ignore |
Each valid row independently proves (via Lemma B) Aut_lin(X_μ) = {1, s} over Q-bar for that μ and hence (Lemma B over Q̄[[μ−μ0]]) for the geometric generic fibre.

## Decision (route second-reflection): NEGATIVE, exact
**Statement.** Let X_η̄ be the geometric generic fibre of the (01)-pencil {r = μt} on X. Then Aut_lin(X_η̄) = {1, s}
(s = the known reflection t ↦ −t). In particular there is no hidden second reflection, no reflection created by any base
change of the pencil, and no fibre automorphism of any other kind (no order-3 or order-6 elements, no second involution).
Also Aut_lin(X_{3/7}) = {1, s} for the specific smooth member μ = 3/7, and Aut_lin(X) = S3 × ⟨τ⟩ (order 12) on P^5 with no
reflection, settling wave-2 item J2 (hidden reflections / §2.7 partner coordinates of X: none).

**Which fibre automorphisms were excluded, and how.**
| class of automorphism of X_μ | excluded? | method | rigor |
|---|---|---|---|
| linear reflections of P^4 (eigenvalues (1^4,−1)) other than s | yes | Lemma A reduction (block-diagonal w.r.t. L ⊕ M) + 12 msolve systems (W4-5) | exact over F_p-bar (p = 1073741827), transferred to Q-bar and to the generic fibre by Lemma B (hypotheses W4-3, W4-4 verified) |
| any other linear automorphism (orders 3, 4, 6, …; involutions of type (1^3,−1^2)) | yes | same computation (it finds the FULL group, order 2) | same |
| automorphisms appearing only after base change of the pencil | yes | they are automorphisms of X_η̄, which inject into Aut_lin(X_{3/7}) (Lemma B over Q̄[[μ−3/7]]) | exact |
| non-linear automorphisms of a smooth fibre | yes | Lefschetz: Pic(X_μ) = Z·O(1), K = O(1) | standard |
| automorphisms of SINGULAR members of the pencil | not addressed | not needed for the §4.6 descent (built from the generic fibre) | — |
| birational self-maps of X_η̄ that are not automorphisms | not addressed (X_η̄ is of general type: Bir = Aut, so none beyond s) | standard | standard |

**Modular-only ingredients and their transfer.** Every Gröbner computation was done over F_p with msolve; each is used only
through Lemma B (the automorphism group over Q-bar injects into the one over F_p-bar, under the verified hypotheses "smooth mod p"
and "H^0(T) = 0 mod p"), so no "empty mod p ⇒ empty over Q" leap is made. Exact-over-Q (Singular) recomputation of the rank-1 loci
is running as an independent cross-check (not needed for the argument). msolve's correctness is assumed.

**Consequence for the construction programme (wave-3 §3.4 M-a).** The §4.6 descent through a second fibre reflection cannot be
realised on the (01)-pencil: W_t = X_μ/s stays the CY3 sextic in P(1,1,1,1,2) with no further quotient, and Z/⟨s̃,σ⟩ does not exist.
Nothing from this route can seed HoldoutCover/PolynomialCurves. The route is closed; P(nonempty) was estimated 0.02 in wave 3 and
is now 0 for linear fibre automorphisms (exact). Remaining open directions for the lane are M-b (served-space criterion for D8 /
(01)-stable degree-8/10 denominators) and non-quotient correspondences, both outside this route.
