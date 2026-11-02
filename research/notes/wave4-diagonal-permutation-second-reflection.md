STATUS: RUNNING — 2026-11-02 04:49 UTC — DECIDED NEGATIVE at mu = 3/7 (exact via Hessian-rank reduction + mod-p computation + specialization lemma): Aut_lin(X_mu) = {1, s}; robustness sweep over (mu, p) and write-up in progress

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
