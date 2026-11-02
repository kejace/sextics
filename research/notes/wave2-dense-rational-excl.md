# Wave 2 — `dense-rational` — lens EXCLUSION / RIGOR

STATUS: in progress — no witness (this lens does not search for witnesses; nothing found by any test so far).
Started 2026-11-01 19:05. Log: `/tmp/hodge/results/dense-rational/wave2-excl.log.jsonl`.
Packet: `/tmp/hodge/results/dense-rational/wave2-excl/` (scripts `wave2_excl_*.py/.sh` in the lane dir, receipts in the packet dir).

## 0. Audit of wave-1 / dossier rigorous claims (so far)

| claim | wave-1 / dossier evidence | wave-2 re-check | verdict |
|---|---|---|---|
| `A = X0^3 Y0^3 not in J(F)_6` | p = 2147483647 (wave 1), mod 101 minor (packet) | p = 2147483587: rank J_6^- = 18 (= Koszul count), rank_aug 19 | agrees; rigorous (modular rank attains cap, augmented rank exceeds it) |
| `(J:A)_d = J_d` for d <= 4, both parities | p = 2147483647 (even), packet R1 (odd, mod 101) | p = 2147483587, both parities d = 1..4: K_full = 0 | agrees; rigorous (saturated-lattice argument) |
| degree 5 excluded | packet R2: "each single D = G_i fails E2" (six coordinate directions only, as summarized) | NEW exact computation: Sym^2(C^6) -> S_4, S -> sum S_ij F_ij has rank 21 (fmpq_mat, exact). With (J:A)_4 = 0 (both parities), E2 for D = v.grad F forces v^T Hess(F) v = 0 identically, hence v v^T = 0, v = 0: **every** direction v (mixed parity included) is excluded | strengthened; rigorous over C |
| degree-7 tau-pure exclusion via msolve (wave-1 §2.4/§4.7/§5) | logs `polar_d7_*_msolve.log` stop at F4 degree 3 (`exit 143`); placeholders never filled | **claim was unproven**; re-running on the remote with p = 1073741789, homogeneous input, `-g 1` leading-ideal certificate (jobs `excl-dense-rational-d7even/odd`) | pending |
| `(J:A)_d = J_d` for d = 5..9, kernels 294/261 at d = 10 | p = 2147483647 (even, direct), p = 101 and 10007 (dossier, dual) | p = 2147483587, both parities, direct nullspace (job `excl-dense-rational-small`) | pending |
| degree 6 excluded (packet R3, rank-one Hessian argument) | mod 101 minors | independent route: degree-6 Hessian relaxation as homogeneous quadratic systems in the even (17 unknowns, 218 quadrics), odd (18, 218) and mixed (35, 426) sectors; msolve `-g 1`, p = 1073741789, run locally (seconds): all three DRL leading ideals contain a pure power of every variable (even: 135 quadratic leading monomials; odd: 153; mixed: 495 of degree <= 3) | agrees; rigorous over C (projective reduction argument); **deg D >= 7 re-proved independently, mixed parity included** |

## 1. Rigorous results (wave 2, so far)

* R-a. A not in J_6 (two large primes + packet's mod-101 minor).
* R-b. (J:A)_d = J_d for d <= 4 in both parities (p = 2147483587; agrees with p = 2147483647).
* R-c. No degree-5 denominator at all: exact kernel computation above (new; closes the "only six directions" gap in the summary of packet R2).
* R-d. No degree-6 denominator in any parity sector (msolve leading ideals, three sectors): `deg D >= 7` for every certificate, re-proved by a route independent of the packet (relies on T1 at d = 6, both parities, and the projective-emptiness argument; receipts `wave2-excl/quadsys_d6_*_ltcheck.json`).

## 2. Running / planned

1. Degree-7 pure sectors (msolve, 32 threads each, remote). Certificate: leading ideal contains a pure power of every variable.
2. Small chain: degree-6 three sectors (msolve) + E1-kernel sweep d = 1..10, both parities, p = 2147483587, with kernel bases and E1 representatives at d = 10.
3. Then: d = 11, 12 kernel bases (`wave2_excl_chain_k1112.sh`); served-target spaces of ~30 named degree-5/6 denominators at a new prime (`wave2_excl_timage.py`); normalized E2 on degree-10 kernel representatives (`wave2_excl_e2reps.py`); kernel structure (`wave2_excl_kstruct.py`); all 12 first-order Fermat deformation seeds at a new prime (`wave2_excl_deform.py`).
