# Wave 2 — off-paired — CONSTRUCTION lens — part 4: computations, verdict, and what would change it

## H. Computations of this session (all in /tmp/hodge/results/off-paired/wave2-construct/, log wave2-construct.log.jsonl)
| script | what | result |
|---|---|---|
| c1_identities.py (sympy, exact) | block form; (X2^3±Y2^3)^2 and (X2^6−Y2^6)^2 identities mod F; cubic factorisation of (H+Cu)^2−4u^6; disc of 2u^3−Cu−H; rational parametrisation of Σ; g-action on E_lambda | all TRUE |
| c2_sigma_n / c2b / c2c (Singular, Q(i)) | singular loci and local types of Σ_n ⊂ (P^1)^3, n=1,2,3,6 | Σ_1: 2A1 (rho=0,w=∞,s=±i) + 2A2 (rho=∞,w=0,s=±1); n ≥ 2: double curves along {rho=0,w=∞}, {rho=∞,w=0} (+ isolated A_k) |
| c3_sectors.py | character/Hodge count of H^4_prim sectors | h-inv Y-odd (56,144,56) all sigma2-anti; g-split 22/17/17 |
| c4_Skappa.sing (mod 32003) | quasi-smoothness of {C=k s^2, H=k s^3} ⊂ P(1,1,1,1,2), k=1,5; base curve | quasi-smooth ⇒ K=O(4), p_g=45, q=0; base curve = Fermat plane sextic (deg 6) ∪ irreducible deg-18 curve |
| c6_doubleplanes / c6b_types (Singular) | branch curves of the double-plane models of Σ_n, n=2..6 | n=2: mu=15 quadruple pt + A5; n=3: mu=19 quadruple pt; n=4: A5+A11 (K3); n=5 (octic): mu=15 + A13 + A5; n=6: A17 (K3) ⇒ p_g = 0,0,0,1,2,1 |
| Reid–Tai check (python) | pseudo-reflections in G; det weights | none; det ∈ {±1} |
| c5_served.py (remote job construct-off-paired-served) | served-space scan (mod 2^31−1, all 18 Y-odd sectors) of the structurally motivated denominators: X2^6−2Y2^6; X2^6+2X2^3Y2^3−3Y2^6; C^2−(X0^2+Y1^2)^4; C^2−(X0Y1+X1Y0)^4; C^2−X0^8; H^2−(X0^2+Y1^2)^6; (X2^3+Y2^3)(X2^3−Y2^3)H; [--heavy: 3-torsion product psi_3·tau^*psi_3 (deg 16), (C^2−mu^4)(X2^6−Y2^6)] | see §I (prediction: all trivial) |

## I. Served-space scan result (remote job construct-off-paired-served, engine v3, p = 2^31−1, all 18 Y-odd sectors)
| D | deg | char | E1 (A·D ∈ J) | nontrivial served classes (sum over sectors) | s |
|---|---|---|---|---|---|
| X2^6 − 2Y2^6 (pullback of a rational 3-section of Σ) | 6 | (0,0) | yes | 0 | 27 |
| X2^6 + 2X2^3Y2^3 − 3Y2^6 | 6 | mixed parity | — | rejected by the engine (mixed Y-parity); its τ-orbit product has degree 12 (not run) | — |
| C^2 − (X0^2+Y1^2)^4 (Λ_mu ∪ gΛ_mu, √kappa on Λ) | 8 | (0,0) | yes | 0 | 230 |
| C^2 − (X0Y1+X1Y0)^4 | 8 | (0,0) | yes | 0 | 227 |
| C^2 − X0^8 | 8 | (0,0) | yes | 0 | 223 |
| H^2 − (X0^2+Y1^2)^6 (cube root of kappa on Λ) | 12 | (0,0) | pending (job construct-off-paired-served30, sector (3,0) only) | pending | |
| (X2^3+Y2^3)(X2^3−Y2^3)·H (2-torsion of Y_± × the kappa=∞ fibre) | 12 | (0,0) | pending (served30) | pending | |
Interpretation: exactly as predicted in part 2 §C — the block multisection pullback (H1) and the vertical √kappa-twist with
linearly equivalent conjugate components (V1) serve only the Jacobian part.  Note that C^2 − mu^4 DOES pass the colon
condition E1, so it is a genuine candidate support that the geometry (not E1) kills.
Mordell–Weil screen (job construct-off-paired-mw + local primes): the Weierstrass fibrations Y_∓ : t^2 = ±2m^3 − Cm − H over
P^3 have NO polynomial section (m quadric, t cubic): msolve returns the empty variety mod 1073741827 (both signs),
2147483647 and 1000003 (eps = −1; for eps = +1 only the first prime: the local reruns timed out on the loaded box) — 30 unknowns, 84 cubic equations.  [SCREEN: affine system,
no properness; sections with denominators (m = M/d^2, t = T/d^3, deg d ≥ 1) not searched.]  Consequently the simplest
realisation of part 3 §E (reflection through a Mordell–Weil section) has no candidate; lead K.2(c) is weakened.
## J. Verdict for the five requested mechanisms
(i)  covers/exchange: FAILS by the no-hub theorem (part 3 §F); the descent-through-reflection lemma (part 3 §E) is the
     replacement criterion; its only known realisation (phi an automorphism of W) is reducible or has p_g ≥ 1.
(ii) Fano/surface-fibration companions: FAILS rigorously for every quotient X/K, K ≤ G, and for every elliptic/K3/abelian-twisted
     companion (Reid–Tai + bicanonical form, part 2 §D).  Open only: non-monomial p_g = q = 0 surface fibrations of X/G.
(iii) Galois closures: FAILS (det character of D6 only served by p_g = 1 rotation quotients; §4.7.2 for the cubic cover).
(iv) q-raising/deformation: no certificate on the sub-family to continue; first-order death of the holdout cover; nothing new.
(v)  N/A.
No finite-field incidence seed is produced because no auxiliary variety with the required support was found.

## K. What a construction must look like (sharpened by this session)
1. By part 2 §C the (3,1)-part of the sector involves H^{2,0} of the pencil S_kappa (p_g = 45) tensored with V^{1,0}; the
   genus-2 variation alone contributes nothing transcendental (Σ rational).  A serving divisor must therefore be of mixed
   type (M): fibrewise a correspondence between the genus-2 curve E_kappa and the surface S_kappa whose Jacobian/Albanese
   data couples J(E_kappa) to H^2_tr(S_kappa) — or a vertical divisor whose fibre components realise a Picard class of
   S̃_kappa with non-trivial monodromy (type V1 loophole), which requires NS(S̃_kappa^{gen}) ⊋ Z·h with monodromy.
2. Decisive next computations (not done here): (a) the Picard local system of the pencil {H = mu s^3} on the quartic double
   solid T = {s^2 = C}: compute ρ(S̃_kappa) at a few rational kappa by point counts / reduction mod p (van Luijk) — if ρ = 1
   generically, type V1 is dead and only (M) remains; (b) dossier E5 (restriction of the source to S_kappa and pairing with
   H^{2,0}(S̃_kappa)) decides V1 directly; (c) a search for birational involutions of W (part 3 §E) — e.g. Bertini/Geiser-type
   involutions of the genus-4 fibration W → P^3 or of the elliptic fibration Y_- → P^3 (translation by a section of Y_- that
   is NOT 2-torsion, lifted to X, is an automorphism of the fibre product that does not normalise sigma2 — the natural
   candidate for phi).  (c) is the most promising lead: sections of t^2 = 2m^3 − Cm − H over P^3 (Mordell–Weil of the
   Weierstrass fibration with (a4,a6) = (−C/2,−H/2)·const) give fibrewise translations; the 2-torsion ones are the mirrors.
3. If (c) succeeds, the book's recipe (polynomial curve family with invertible incidence Jacobian on the Fano quotient
   Z/<sigma,r>, polynomial normal, integration, trace, transport) applies verbatim and the shared Lean library suffices.
