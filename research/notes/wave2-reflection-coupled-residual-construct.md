# Wave 2 — `reflection-coupled-residual` — lens CONSTRUCTION

STATUS: in progress (started 2026-11-02 ~00:40 local; prior work read: dossier, attack1, sibling wave-2 notes excl/search/verify, book §2.7, §3.1–3.4, §4.6–4.7, companion solution; no experiments yet)

Target: q = 1, A = X0^3 Y0^3 on F = Σ Z^6 + (5/7) P (13 monomials), the same fourfold as the solved `reflection-coupled-fano`.
Klein group {1, r: Y2↦−Y2, ρ: (Y0,Y1)↦−(Y0,Y1), τ = rρ}. Residue of A: (r, ρ, τ) = (−, −, +); sector (h31,h22,h13) = (126, 492, 126).
Work dir: /tmp/hodge/results/reflection-coupled-residual/construct/ ; log: /tmp/hodge/results/reflection-coupled-residual/wave2-construct.log.jsonl

## Plan
(i) exchange-cover scan over all coordinate pairings (Fermat covers T^6 = Z_i^6 + Z_j^6, exchange T ↔ Z_k);
(ii) quotient companions: Klein quotients only (Aut(F) is Klein, dossier E5); Hodge numbers by sector;
(iii) Galois closure of the natural finite map X → P^4_x (degree 6, group B3 = Z/2 ≀ S3): averaging-ideal test on every intermediate quotient (reflection count ⇒ canonical class);
(iv) deformation/q-raising: one-reflection locus (§4.4), Fermat pencil;
(v) companion support modifications: branch divisor B_W, liftable rational curves (u|_C a square), lines on X (Fano scheme), six Fermat surfaces X ∩ {Z_i = Z_j = 0}.

## Experiments log
(none yet)
