# wave2-excl — reproducible packet (lane polarized-k3-remainder, lens EXCLUSION/RIGOR, 2026-11-01/02)

Conventions. (U,V) frame: U = X+Y, V = X-Y (slots 0..2 = U, 3..5 = V); G = 7f(U) + 7f(V) = 14 F-hat;
T = prod_i (U_i^2 - V_i^2) = 64 A-hat (full target);  T_K3 = U0^2U1^2U2^2 - V0^2V1^2V2^2 = 64 x (K3-sector projection of A-hat).
In the XY frame: A_K3 = [(X0+Y0)^2(X1+Y1)^2(X2+Y2)^2 - (X0-Y0)^2(X1-Y1)^2(X2-Y2)^2]/64 (Y-odd, degree 6).
"Character" r = V-degree mod 6 (the mu6 grading V -> eta V, a symmetry of G, not of T).

Scripts (all read-only w.r.t. the repository; engine at ~/devel/sextics-work on PYTHONPATH):
* relax.py        — independent implementation of the L-relaxation (necessary quadratic system) for D in J(G):
                    D = sum P_i a_i, condition B_D(T a) in L = (G) + J(G)^2 + (R_ijk).  Per-character L-ranks mod p and
                    over Q (--exact-ranks), constraint rows (quotient coordinates), exact constraints (--exact-constraints,
                    d<=6), iterated square elimination (--squares).  Receipts relax-<target>-d<d>-p<p>.json (+ -exactQ.json, .npz).
* fixedD_koszul.py — exact rational fixed-D test in the Koszul normalization (unknowns: Koszul coefficients q_jk and
                    quotient Q).  For mu6-homogeneous D and target the system splits by character; ranks mod p and over Q.
                    An exact rank gap is a rigorous exclusion of that D for that target.  Receipts fixedD-koszul.jsonl.
* e1table.py      — kernel of multiplication by [A] and [A_K3] on the Jacobian ring (XY frame, engine), second prime.
* k3kernel.py     — non-Jacobian E1-kernel of [A_K3] at d=8 (XY frame) and fixed-D tests (mod p) of structured/random elements.
* k3kernel_uv.py  — same in the (U,V) frame with the engine's mu6 character-restricted solver (6x smaller systems);
                    exact basis of ker(eta: R_{f,4} -> R_{f,10}); structured pairings sum g_k(U) g_k(V).
* servedK3.py     — K3-block projection of the full served target space of a mu6-homogeneous D; --exact = rigorous
                    rational computation (fmpz nullspace).  Modular mode is evidence only (rank mod p <= rank over Q).
* relax_member.py — exact L-membership of B_D(Ta) for explicit a (certifies points of the relaxation over Q).
* survivors.py    — residual ideal of a relaxation after square elimination -> Singular dim/degree/minAssGTZ.
* tangent_probe.py — tangent hyperplane sections through {V=0} (mod-p probe of the t^2 + f(V) mechanism).
* relax_multi.py  — 19-target K3-block relaxation (same Koszul vector a for all classes B(U)-B(V)); zero-class check
                    with genuine J(f)_6 targets; per-character column grading.  Receipts multi-k3-d<d>-p<p>.json/.npz.
* zstruct.py      — tangent-space dimension of a relaxation zero set at points of the family f(U)·ℓ; degree-2 saturation
                    search for linear forms.  Receipt <prefix>-zstruct.json.
* k4_lll.py       — exact basis of K4 = ker(η: R_{f,4} -> R_{f,10}) (dim 9), membership verified exactly.  K4-basis-lll.json.
* struct_d7.py / struct_dn.py — exact Koszul fixed-D screens of random mu6-homogeneous members of J_d (all characters,
                    both targets).  Receipts fixedD-koszul-struct-d<d>.jsonl.
* relax.py --char=b — relaxation restricted to mu6-homogeneous D of character b (chain_char8.sh runs all sectors at d=8).
* member_family_d7.py — exact Q membership of the whole family f(U)·ℓ in the d=7 K3 relaxation (21 coefficient polys;
                    stopped under memory pressure, superseded by the sum-principle exclusion of that family).
* cands_d8d10_xy.jsonl + job_sweep_d8d10.sh — wave-1 structured degree-8/10 candidates (XY frame) for hodge.sweep (mod p).
* job_relax_full.sh / job_relax_k3.sh — remote drivers (degrees 6, 7; primes 65521 and 101).

Replay (examples):
  python3 relax.py 6 65521 --exact-constraints --squares --out=relax-full-d6-p65521
  python3 relax.py 6 65521 --target=k3 --exact-constraints --squares --out=relax-k3-d6-p65521
  python3 fixedD_koszul.py fU-k3-d6 6 "[scale(U[0],Fr(1,42)), scale(U[1],Fr(1,42)), scale(U[2],Fr(1,42)), {}, {}, {}]" --target=k3
  python3 fixedD_koszul.py control-G-d6-k3 6 "[scale(z(i),Fr(1,6)) for i in range(6)]" --target=k3     # must be consistent

Rigor conventions used in the receipts:
* "modular full rank => injective over Q" (rank mod p <= rank over Q) is used for all E1 injectivity statements.
* Exact statements say "over Q"/"exact"; everything else is "mod p" and is labelled as such.
* The L-relaxation is only a NECESSARY condition on D; its nonemptiness never implies a certificate.

## Replay recipe for the main new statement ("no denominator of degree <= 6 serves the K3-sector target")
1. E1 injectivity for A_K3, d <= 7 (mod 2^31-1, hence over Q):  `python3 e1table.py 2147483647 7 e1-AK3.json`  (reads targets A and A_K3;
   kernel_in_ring must be 0 for A_K3 in all cells d <= 7).  Receipt used: e1table-AK3-local-p2147483647.json.
2. Degree 5:  `python3 relax.py 5 65521 --target=k3 --exact-constraints --squares --out=relax-k3-d5`  -> all 6 parameters killed over Q.
3. Degree 6:  `python3 relax.py 6 65521 --target=k3 --exact-constraints --squares --out=relax-k3-d6`  -> 32/35 killed; survivors
   c5, c10, c15 (= a_i = c_i U_i, i < 3); residual quadrics (c0-c2)^2, (c0-c2)(c1-c2), (c1-c2)^2 (see survivors.py relax-k3-d6 --primes),
   so D ~ f(U) modulo the gauge D -> D + tG.
4. D = f(U):  `python3 fixedD_koszul.py fU-k3-d6 6 "[scale(U[0],Fr(1,42)), scale(U[1],Fr(1,42)), scale(U[2],Fr(1,42)), {}, {}, {}]" --target=k3`
   -> rank_Q 136, augmented 137 (inconsistent).  Control: `... control-G-d6-k3 6 "[scale(z(i),Fr(1,6)) for i in range(6)]" --target=k3` is consistent.
5. (Stronger, unconditional for these D)  `python3 servedK3.py fU "<f(U) in X-variables>" --exact` -> served K3-class dimension 0 over Q.
Hypotheses used downstream (not in steps 1-5): the thread's normalization lemma (REPORT.md), the book's prescribed-support theorem,
End_Hdg(T) = Q (thread) and η ∉ NS_prim⊗C for the dichotomy-based consequences.

## Glossary (conventions used in receipts)
* character r (or "char"): V-degree mod 6 of a (U,V)-monomial; the mu6 action V -> eta V is a symmetry of G = 7f(U)+7f(V)
  but not of the target T, so only mu6-homogeneous D (all monomials of one character) give block-split systems.
* Y-parity (engine, XY frame) = swap (U <-> V) symmetry; V-parity (UV frame) = character mod 2.  They are different splittings;
  totals agree (e.g. the degree-8 non-Jacobian K3 kernel: 51 + 39 (Y-parities) = 81 + 9 (V-parities, chars 4 and 1)).
* "gauged parameters": coefficients of a = (a_0..a_5) in D = sum P_i a_i after removing the Euler directions a -> a + h z/6
  (h of degree d-6), implemented as dropping the z5-divisible coefficients of a_5 (restricted to the sector when --char is used).
* "exact": computed over Q with fmpz/fmpq matrices (rigorous); "mod p": modular evidence unless the statement is of the form
  "full rank mod p => full rank over Q" (used for all E1 injectivity claims) or is accompanied by the flatness argument.
