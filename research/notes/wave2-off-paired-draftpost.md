# Draft discussion comment — `off-paired` (wave 2, exclusion/rigor lens) — DRAFT, not posted

**Scope.** Exact-rigor re-verification and extension of the degree bounds for `off-paired` (q = 1, A = X0^3 Y0^3). No certificate
was found; nothing here claims nonexistence beyond the stated finite families and degrees. Reproduction: the packet
`wave2-excl/` (scripts `w2*.py`, receipts `receipts_*.jsonl`, README_w2.md) — every statement below points to a JSON receipt
with the prime(s) and ranks used.

**Rigor conventions.** A modular rank is a lower bound for the rational rank. A fixed-denominator exclusion is stated as a theorem
only when rank_p([M|b]) = rank_p(M) + 1 **and** rank_Q(M) <= rank_p(M) is certified, either by the explicit rational kernel
(Euler gauge L·Z and the F-multiples (F·W, B_D(W)); its dimension is computed exactly) or by a rationally reconstructed basis of the
whole modular kernel verified exactly. An E1 failure (A·D not in J(F)) is a theorem when the modular rank of the Jacobian generators
equals the exact Koszul dimension of J(F)_{d+6} in the sector. Bilinear sector sweeps (D free in a character sector, N free) decided
by msolve on the complete chart cover t_b = 0 (b < a), t_a = 1 are modular screens unless an exact certificate is attached.

**Results (to be filled from the receipts at posting time).**
1. Degrees 4–6: re-derived (D0 = F_X0/(6X0): ranks 67/68 with nullity 8 = explicit kernel: theorem; degree-6 sectors chart-complete).
2. Degree 7, all 36 diagonal characters: chart-complete `no-solution` mod 2147483647 (171 charts; 4 sectors have empty colon space) —
   independent re-derivation of deanodotdotdot's degree-7 theorem (different code, prime and chart method).
3. Degree 8, all 36 characters, **including the six h-invariant sectors (a,0) left open in the thread**: chart-complete `no-solution`
   mod 2147483647 (436 charts), with the Jacobian-relaxation quadrics (B_D(N_D) == 0 mod J(F)) added as consequence equations.
   Status: modular screen; exact elimination over Q (Singular) is being run to turn the open sectors into a theorem.
4. Structured fixed denominators (dossier block products of degrees 7–14, wave-1 families, Hessian minors): two new primes, with the
   explicit-kernel rank sandwich: counts of PROVEN exclusions to be filled.
5. Colon / non-Jacobian kernel dimensions per sector for d <= 12 (table in the packet).

**Caveats.** Products of excluded reduced supports are not automatically excluded (served spaces add); every statement is about the
listed reduced supports or about single-character denominators of the stated degree. Mixed-character denominators of degree >= 7 and
all denominators of degree >= 9 (single character: >= 10 after the running degree-9 sweep) remain open.

Credit: fqlx (degree <= 6, surface obstructions), benediamond (contact obstruction), deanodotdotdot (degree 7 / degree-8 sectors,
relaxation method), hodge.computer book (prescribed-support theorem, Euler normalization). Prepared with Claude (Anthropic) agents.
