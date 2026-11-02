# README-referee — lane `polarized-k3-remainder`, wave-2 exclusion packet (referee pass, 2026-11-02)

Referee: Claude (Anthropic) for `kejace`; no human review; nothing Lean-verified; no certificate claimed.
Inputs: draft post, excl report, this packet, construct and search reports, aggregator status + verification log
(19 lane entries), dossier, saved thread #115.  No adversarial verify-lens file exists for this lane.
Artefacts of this pass: `referee/` (scripts, receipts, logs) — all computed on the shared box under `ulimit -v 6 GB`,
niced, one process at a time; total compute < 10 CPU-minutes.

## 1. Independent replication (own code, tensor frame) — `referee/referee_tensor.py`, log `referee/referee.log.jsonl`
Frame: U = X+Y, V = X−Y; F = ½(f(U)+f(V)) with f = z0^6+z1^6+z2^6+(5/7)(z0^5z1+2z1^5z2+3z2^5z0+z0^2z1^3z2);
J(F) = J(f)(U)+J(f)(V), so R_F = R_f ⊗ R_f (bigraded).  A = 2^-6 ∏(U_i^2−V_i^2) = Σ_S (−1)^|S| m_{S^c}(U) m_S(V),
A_K3 = (η(U)−η(V))/64 with η = z0^2z1^2z2^2, A_mix = the six mixed terms.  E1 maps are built from the ternary ring only
(exact normal forms over Q, dims ≤ 19), so every rank below is exact over Q (fmpq_mat) and at four primes
(10^9+7, 2^31−1, 10^9+9, 999999893):
* C0 frame identity `F(instance) == ½(f(X+Y)+f(X−Y))` and `A == 2^-6 ∏(U_i^2−V_i^2)` exactly.
* C1 Hilb(R_f) = (1,3,6,10,15,18,19,18,15,10,6,3,1); η ≠ 0 in R_{f,6}.
* C2 A_K3 ≠ 0 and A_mix ≠ 0 modulo J(F) (bigraded pieces (6,0)+(0,6) vs (4,2)+(2,4); J_k = 0 for k ≤ 4).
* C5 K4 := ker(η: R_{f,4} → R_{f,10}) has dimension exactly 9 over Q (η onto R_10, rank 6).  C6 η: R_1 → R_7 injective.
* E1 for A: ker([A]: (R_F)_d → (R_F)_{d+6}) = 0 for d ≤ 9 EXACTLY over Q (domain dims 1,6,21,56,126,246,426,666,951,1246);
  d = 10: kernel 555 at all four primes (exact rank still running at hand-off; not load-bearing).
* E1 for A_K3: kernel 0 for d ≤ 7 exactly over Q; kernel exactly 90 at d = 8 and 306 at d = 9 over Q.
  With C5/C6 and the formal identities (η⊗1−1⊗η)(K4⊗K4) = 0 and (η⊗1−1⊗η)(ℓη⊗ℓ' + ℓ⊗ℓ'η) = 0 (deg ℓη^2 = 13 > 12)
  this proves the lens's R10–R11 description: ker_8 = K4⊗K4 (81) ⊕ {ℓ(U)ℓ'(V)(η(U)+η(V))} (9), exactly.

## 2. Packet-script replays (same code as the lens, re-run by the referee) — receipts in `referee/`
| replay | command (from the packet directory) | result | receipt |
|---|---|---|---|
| D = f(U) vs A_K3, d = 6 | `fixedD_koszul.py referee-fU-k3-d6 6 "[scale(U[i],Fr(1,42)) ...]" --target=k3` | char-0 block 966×167, rank_Q 136 / aug 137: inconsistent | `packet-replays.jsonl` |
| control D = G (illegal) vs A_K3 | `... referee-control-G-d6-k3 6 "[scale(z(i),Fr(1,6)) ...]" --target=k3` | 56/56 consistent | same |
| K3 relaxation d = 5 | `relax.py 5 65521 --target=k3 --exact-constraints --squares` | 6/6 parameters killed over Q | `relax-k3-d5-referee*` |
| K3 relaxation d = 6 | `relax.py 6 65521 --target=k3 --exact-constraints --squares` | 32/35 killed over Q; survivors c5,c10,c15 = coefficient of U_i in a_i; exact residual row space = span{(c5−c15)^2,(c5−c15)(c10−c15),(c10−c15)^2} (`residual-k3-d6-exact.out`) ⇒ c5=c10=c15 ⇒ D = 42c·f(U) mod gauge | `relax-k3-d6-referee*` |
| full relaxation d = 6 | `relax.py 6 65521 --exact-constraints --squares` | 35/35 killed over Q; block L-ranks 603/558/396/260/396/558 (mod p = Q) | `relax-full-d6-referee*` |
| D = f(U)U0 vs A_K3, d = 7 | `fixedD_koszul.py referee-fUU0-k3-d7 7 ...` | 1926×447, 377/378 exact | `packet-replays.jsonl` |
| D = f(U)U0^2 vs A_K3, d = 8 | `... referee-fUU0sq-k3-d8 8 ...` | 3366×888, 756/757 exact | same |
| D = f(U)U0^2 vs A, d = 8 | `... referee-fUU0sq-full-d8 8 ... --target=full` | blocks 756/757, 828/829, 516/517 exact (= thread) | same |
| served-space K3 projection | `servedK3.py referee-fU "<f(U)>" --exact`, same for `X0`, `Y0` | f(U): 4923×832, nullity 110, 17 served K3-monomial vectors, class dim 0 over Q; U0, V0: class dim 0 | `servedK3-referee.jsonl` |

## 3. Transfer and coverage checks (by hand)
* Normalization (thread REPORT.md) re-derived from the contract's cleared identity: N(F) = AD − FH, N' = N + Hz/6 gives
  N'(F) = AD and B_D(N') = B_D(N) + DH ∈ (F).  The Euler gauge D ↦ D + hF, N ↦ N + (hA/6)z maps certificates to
  certificates and legal denominators to legal ones (checked: B_{D+hF}(N') ≡ B_D(N) mod F), so gauge-fixing is sound.
* Koszul form: ∂_iF regular (X smooth) ⇒ M = Ta + S_q; B_D(S_q) ∈ J(G)^2 + (R_ijk) ⊆ L (div S = q_jP_k − q_kP_j and
  S(D) = Σ a_i R_ijk + J^2-terms) ⇒ the L-relaxation is a necessary condition.  Square elimination over Q is sound over C.
* Parameter coverage: a ∈ S_{d−5}^6 unrestricted (no parity ansatz); a ↦ D injective for d ≤ 9 (no Koszul syzygy of
  a-degree < 5); gauged dimension = dim J_d − dim (F)_d = 35, 120, 315 at d = 6, 7, 8 — matches the receipts.  Both
  Y-parities are therefore covered at d ≤ 7.  The degree-8 character sectors (`--char=b`) cover only μ6-homogeneous D.
* Degree 7 (thread p = 101, lens p = 65521; NOT replayed by the referee, ~30 min under the load cap): equal L-block ranks
  over Q and F_p ⇒ free cokernel over Z_(p) ⇒ the special fibre of the projective membership scheme is the mod-p system;
  nine empty charts + square reduction + homogeneous GB of affine dimension 0 ⇒ special fibre empty ⇒ generic fibre
  empty by properness.  The argument is correct as written; two independent implementations at two primes agree.
* Aggregator log (19 entries, 10^9+7 / 998244353) agrees on every overlapping item.  Its `hessrelax6_mixed_mprimary =
  False` means the plain mod-J Hessian relaxation does NOT close degree 6; the exclusion rests on the finer L-relaxation
  (thread + R1), which is what was replayed here.  No contradiction.

## 4. Gaps and wording changes (G = gap, W = wording)
* G1 Degree 7 for A_K3 is open (relaxation rank 935/7260, 19-target 1092, no squares).  Hence the sum-principle
  reduction only removes products of supports of degree ≤ 6, and it needs End_Hdg(T) = Q (thread) and η ∉ NS_prim⊗C
  (automatic if ρ = 1; only 1 ≤ ρ ≤ 3 is known).  Stated conditionally in the comment.
* G2/W Block cones g(U), g(V), g(U)h(V) "excluded in every degree" (R9, draft 2b) is a Hodge-theoretic sketch: the
  product-divisor part is fine (H^1(S_6) = 0), and the exceptional divisor over C×C contributes to the Künneth H^2⊗H^2
  part only the algebraic class [C]⊗[C] by the blow-up formula; the transport through the μ6-quotient and the book's
  projectors is not written.  Downgraded to "argument sketch, not a theorem".
* G3 Degree-8 results for A are mod-p and sector-restricted (char 0: cones g(U), g ∈ J(f)_8 with 6 residual quadrics;
  char 1: relaxation essentially empty; char 2: Σ U_i f_i(U) q_i(V) + g(V)); the full degree-8 problem is open.
* G4 All E2 tests of non-Jacobian K3-kernel elements are modular (2^31−1): screens only.
* G5 f(U)V0 (d = 7, A_K3) and the random J_7/J_8/J_9 members: exact in the lens's receipts, not replayed here.
* G6 The equivalence "HasCertificate(F,1,A) ⇔ GHC (coniveau 1) for Λ^2T ⊂ H^4(S_2×S_2)" (dossier §4.3) is a sketch
  with explicit hypotheses (End_Hdg(T) = Q, η ∉ NS_prim⊗C); the ⇐ direction also uses the book's existential supports
  for the mixed sector.  Posted as a reduction, not a result.
* G7 Construct Part B (mixed-sector recipe, F_5 seed, Hensel lift) and "Aut(X) = D_6 exactly" (needs Aut(C_f) = 1 from
  the probe) were not refereed; "Aut(X) ⊇ D_6" is trivial (f(ηV) = f(V)).
* W "served(f(U)) has zero K3 projection" is rigorous only in the `--exact` mode (lens's own correction); replayed exactly.

## 5. Verdict counts (claim table in `referee/referee-report.md`)
28 rows: VERIFIED 16, PLAUSIBLE 11, UNSUPPORTED 1 (block-cone statement as worded).  Strongest verified statement: the K3-sector
class A_K3 = P_0(A) admits no denominator of degree ≤ 6 (exact over Q), and E1 for the instance target is injective
through degree 9 exactly over Q; the thread's deg D ≥ 8 stands on two independent implementations.
