# Referee report — `cyclic-triple-residual` (#155), wave 2, exclusion/rigor lens

Referee pass of 2026-11-02 (Claude, Anthropic, for `kejace`). Inputs: the draft post, the exclusion report
`wave2-cyclic-triple-residual-excl.md`, the packet `wave2-excl/` (+ the remote mirror of its remote-only files), the
adversarial verify notes, the aggregator status, the dossier and the saved discussion thread. Independent code (no
engine imports): `research/results/cyclic-triple-residual/wave2-referee/ref_lib.py` and `ref_*.py`; receipts in
`wave2-referee/referee.log.jsonl`, `e1_table.json`, `hess_*`; copies of the remote-only receipts in `wave2-referee/mirror/`.
Fresh prime: 999999937 (not used by wave 1: 2147483647/2147483629; excl: 1000003/32003/2147483587; parallel excl
session: 1073741789/998244353; verify/aggregator: 1000000007, Singular 10007/32003). Exact over Q wherever feasible.

Verdict labels. VERIFIED: replayed by the referee (exact, or modular at the fresh prime with the transfer argument
checked); PLAUSIBLE: a receipt exists but it is modular-only or was not replayed here; UNSUPPORTED: no receipt in the
packet, or the draft's wording is stronger than the receipt.

## Claim table

| # | claim (draft section) | agent receipt | referee replay | verdict | wording for the post |
|---|---|---|---|---|---|
| 1 | `Ann([A])_d mod J_d = 0` for `d <= 7`, all parities and characters (§1) | kerverify (p = 1000003, 32003), wave 1 (2147483647); verify lens exact | `ref_e1.py`: exact over Q, all 18 sectors, d = 0..7; also mod 999999937; `dim J_n` = c.i. value for n <= 17 | VERIFIED | "exact over Q" |
| 2 | `Ann_8 = 1`, spanned by D8 (25 terms, sector (0,4)); `D8 notin J_8`, `A D8 in J_14`, explicit integer `P` (§1, §5) | e1_d8_exact.json (remote-only), kerverify | `ref_d8.py`: all three facts exact; own fmpz nullspace: `Ann_8 = 1`, packet D8 in its span; `ref_e1.py` exact d = 8 | VERIFIED | "exactly one-dimensional over Q" |
| 3 | `Ann_9 = 120` (116 forced + 4 sporadic), sector split (§1) | kerverify (3 primes); parallel session `a1_exact_8_12.out` (exact fmpz, 1124 s) | `ref_e1.py` mod 999999937: 120, identical sector split (upper bound); lower bound 116 from the verify lens's recount | VERIFIED | "120 (exact)" |
| 4 | `Ann_10 = 556`, `Ann_11 = 1020` (§1) | kerverify (3 primes); forced lower bounds 555, 1020 | `ref_e1.py` mod 999999937: 556 and 1020, sector splits identical to the agent's | PLAUSIBLE (exact value at d = 10); `Ann_11 = 1020` is exact (forced = modular) | "556 at four primes (555 <= Ann_10 <= 556), 1020 exact" |
| 5 | modular kernel dims are upper bounds for the Q-kernel when `dim J mod p` = c.i. value (§1) | stated in excl.md | argument re-derived: `rank_p J = rank_Q J` and `rank_p [J|A S] <= rank_Q [J|A S]`; c.i. equality checked at 999999937 in all degrees used | VERIFIED (logic) | keep |
| 6 | every certificate can be normalised to E1; no denominator of degree <= 4; degree 5-7 denominators are `D = G(F)`; degree 8: `c D8 + G(F)` (§1) | excl.md, attack1 | shift `N -> N + R Z/6` checked by hand (uses `D` coprime to the irreducible `F`); Koszul exactness for the syzygies | VERIFIED (logic) | keep |
| 7 | Hessian relaxation: E2 mod J reduces to `A G^T Hess(F) G in J_{2d}`; invariance under `G -> G + hE`; the slice is a complement of the Euler line (§2) | w2lib docstring | re-derived by hand (`N(D) = A G^T H G mod J`, `E^T H = 5 grad F`, `E^T H E = 30 F`, `h -> h Z_j` injective) | VERIFIED (logic) | keep |
| 8 | `d = 5`: gauge-fixed relaxation m-primary, no degree-5 denominator (§2) | relax_d5 (p = 1000003); verify: Singular char 0 | `ref_hess.py d5`: own quadrics exact over Q span all 21 quadratic forms; Singular char 0 `dim 0, vdim 7` | VERIFIED (unconditional) | "exact over Q" |
| 9 | `d = 6`, all parities/characters, 35 unknowns after gauge: m-primary, no degree-6 denominator (§2) | relax_d6_p1000003_g_gauge.json (35/35 pure powers); aggregator p = 1000000007; verify Singular GF(10007), GF(32003) | `ref_hess.py d6`: own quadrics (301 independent), p = 999999937, X0- and Y2-slice: Singular `dim 0, vdim 1068`, msolve 35/35 pure powers; `rank_p J_12 = 4437` = c.i. | VERIFIED | "rigorous (modular m-primary + transfer), replicated at three primes and two slices" |
| 10 | transfer "m-primary mod p => `V(I_Q) = {0}` on the slice" under `dim J_{2d} mod p` = c.i. value (§2) | excl.md | argument re-derived (primitive p-integral scaling; saturated lattice of `J_{2d}` reduces onto `J_{2d} (x) F_p` when the ranks agree) | VERIFIED (logic) | keep; state the hypothesis |
| 11 | `d = 7`, G-semi-invariant `D`: per-sector m-primary, "12 of 18 sectors" (§2 bullet) | hsector_d7_*.local.out: 13 sectors (0,0)-(0,5), (2,0)-(2,5), (4,0) | `ref_hess.py d7`: all 18 sectors at 999999937, X0-slice, Singular `dim 0` and msolve pure powers for every variable; `rank_p J_14 = 10122` = c.i. | VERIFIED and COMPLETED | "no G-semi-invariant denominator of degree 7 (all 18 characters)" |
| 12 | full `d = 7` relaxation (126 unknowns, 1441 quadrics) m-primary (§2 "[pending]") | relax7.out: quadrics built only; msolve GB unfinished after > 2 h on the loaded remote | not feasible locally (degree-4 Macaulay step ~ 10^7 columns) | UNSUPPORTED (pending) | "open" |
| 13 | full `d = 8` (337 unknowns) and the degree-8 sector chain (§2 "[pending]") | hsector_d8_00.out: built only | — | UNSUPPORTED (pending) | "open" |
| 14 | Fermat control (`t = 0`): same d = 6 pipeline is NOT m-primary (§2) | relax_d6_p1000003_fermat_gauge.json (12/35) | not replayed | PLAUSIBLE (negative control, modular is adequate) | keep as a control |
| 15 | degree-0 Koszul relaxation d = 5 m-primary; wspan diagnostics (§2) | relax0_d5_*; wspan_d6/d7.json | not replayed; d = 5 subsumed by row 8 | PLAUSIBLE (diagnostic) | drop or one clause |
| 16 | degree-9/10 "E1-admissible sector families satisfy the Koszul relaxation only on the illegal locus, sectors done: ..." (§3) | none at d >= 9 (only `bilin_d6_p1000003_s04`, d = 6, per-form tests = modular evidence) | — | UNSUPPORTED | delete |
| 17 | "every individual exact kernel element (120 at d = 9, 556 at d = 10) excluded ... sealed by exact null vectors" (§3) | remote kexact9: 120 exact elements (A K in J verified exactly by the agent); sector tests: 116 "inconsistent mod p, not sealed" (p = 1000003, 2147483587; nullity 69 unsealed), 4 exact (fmpz fallback: (0,0)#1, (0,3)#0, (2,5)#0, (4,1)#0); nothing at d = 10 | `ref_k9.py`: 8 of 8 elements checked are genuine, exactly over Q (the smallest-height element of 8 sectors, including all 4 of the agent's exactly sealed ones; the picks for sectors (4,0), (4,2), (4,4) were skipped, heights > 100 digits): each is sector-pure, `K notin J_9`, `A K in J_15` | PLAUSIBLE (screened) for 116/120; VERIFIED-by-agent-exact for 4; d = 10 UNSUPPORTED | "120 exact kernel elements screened mod p at two primes (4 exact)" |
| 18 | structured candidates up to degree 12 "all sealed inconsistent (table in the packet)" (§4) | no `families.out`; separate exact receipts exist elsewhere (verify lens; parallel session's Farkas duals; search log) | Farkas duals replayed exactly against all unrestricted columns: D8 (15015 cols), `E_X(F)`, `X0^6+X1^6+X2^6`, `X2^6+Y0^6+Y1^6+Y2^6` (`fqlx`), `X1^6+X2^6+Y1^6+Y2^6` | UNSUPPORTED as stated; the five Farkas duals VERIFIED | list only what has receipts |
| 19 | D8 excluded over C: sector system 2357 x 1188 sealed (33 exact null vectors, 1155/1156), exact fmpq rref (§5) | d8full2 log lines; verify lens exact 1155/1156 and 1041/1042 | `ref_d8.py`: own single-identity G-sector system (5506 x 1188 = 939 N + 249 Q), exact fmpz rank 1041, augmented 1042; Farkas dual (968 nonzeros) exact against all 15015 unrestricted columns, `y . rhs = 1` | VERIFIED (unconditional) | "excluded, exact Farkas dual" |
| 20 | `fqlx`'s exclusion of `X2^6+Y0^6+Y1^6+Y2^6` reproduced (§5) | dossier; verify lens exact 436/437 | Farkas replay exact (row 18) | VERIFIED | credit `fqlx` |
| 21 | dimension-forced lower bounds 0, 116, 555, 1020, 1325 (d = 8..12), sporadic 1, 4, 1, 0 (excl.md (a)) | forced_dims.json; verify lens own count | not re-run | PLAUSIBLE (exact combinatorics, independently recounted by the verify lens) | keep with attribution |
| 22 | structure of the sporadic degree-9 annihilators (`X0 D8, Y2 D8 in J`, `X2 D8` sporadic in (2,5), ...) (excl.md (a)) | "mod 1000003, local", no file | not replayed | PLAUSIBLE (modular) | mark "mod p" |
| 23 | NEW (referee attempt): pure-parity degree-7 families (Y-even / Y-odd `D`, 60 unknowns each after gauge) tested with the Hessian relaxation | parallel session built 720 quadrics (Y-even), msolve unfinished | `ref_hess.py d7yeven / d7yodd`: Y-even: 60 unknowns after gauge, 452 independent quadrics, p = 999999937; msolve exited rc=-11 (crash under a 12 GB virtual-memory cap) after 1380.0 s - INCONCLUSIVE; input saved as hess_d7_p999999937_X0_yeven.ms; relaunched with a larger cap at referee exit / Y-odd: msolve run not finished at referee time (no receipt) | Y-even: INCONCLUSIVE (msolve crashed); Y-odd: not run / unfinished | family statement only (a degree-7 `D` need not have pure parity) |

Counts (rows 1-22, one verdict per row): VERIFIED 12 (rows 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 19, 20); PLAUSIBLE 6
(rows 4, 14, 15, 17, 21, 22); UNSUPPORTED 4 (rows 12, 13, 16, 18). Row 23 is new work by the referee.

## Soundness of the transfer arguments used

* mod p -> characteristic 0 for E1 kernels: correct direction (modular kernel >= rational kernel) once
  `rank_p J_n = rank_Q J_n` in the two degrees used; `rank_Q J_n` is the complete-intersection value because `F` is
  smooth over Q (`F = Fermat mod 5`); the modular equalities are receipts (agent: p = 1000003, 32003; referee: 999999937).
* m-primary leading ideal -> no projective zero over `Q-bar`: the relaxation ideal is homogeneous and its mod-p generators
  are the reductions of the rational ones (saturated-lattice argument, needs the same rank equality in degree `2d`); a
  nonzero `Q-bar`-zero scaled to a primitive integral vector reduces to a nonzero zero mod a prime above p. Pure powers of
  every variable in the leading ideal of a Groebner basis <=> finite colength <=> `V = {0}` over `F_p-bar`. Singular's
  `dim(std(I)) = 0` was used as a second engine in every referee run. The gauge slice is a complement of the Euler line
  because `h -> h Z_j` is injective, and translation by the Euler line preserves `V(I)` (the quadrics change by elements of
  `J`), so `V(I) subset L <=> V(I) cap slice = {0}`.
* Farkas duals: a rational `y` with `y . col = 0` for every column of the unrestricted single cleared identity and
  `y . (-A D^2) = 1` proves inconsistency over every field of characteristic 0 (and the referee checked it against all
  columns with no parity or sector restriction). Sector reduction of fixed-D systems is lossless for G-semi-invariant `D`
  (project a solution onto the isotypic component of `A D^2`).
* Sealed modular ranks: `rank_p <= rank_Q`; `nullity_p` exactly verified rational null vectors force `rank_Q = rank_p`;
  with `rank_p(aug) = rank_p + 1` the system is inconsistent over Q. Correct. The agent's `sealed_exclusion` reports
  "SEALED" only when all null vectors are verified or when the exact fmpz fallback ran (4 degree-9 elements); the
  116 others are modular-only and must be called screened.
* Prescribed-support / power statements (verify lens rows 5, 13; `fqlx`'s "higher pole multiplicities with that same
  reduced support are also excluded") rest on the book's prescribed-support theorem at q = 1; not re-derived here.
* Sum principle at q = 1: not used by the excl lens; the search lens's "served-space" test did not run (no receipt).

## Reconciliation with the other lenses

* Verify lens (27 checks, written before the draft existed): no disagreement. Its two "weakened" items are confirmed
  here (degree-9 element exclusions are modular; D8 coset non-membership now exact). Its refutation of the attack1
  §4.2 plane-order prose does not affect the draft post (that text is not in it).
* Search lens vs. excl lens at d = 6, Y-even: the search lens's block-Macaulay bound (`HF_3 = 11, HF_4 = 1`, "at most
  one projective point modulo the F-gauge") is compatible with, and superseded by, the m-primary result (`HF_5 = 0`:
  three independent runs at p = 1000003, 1000000007, 999999937). The "one candidate point" is empty.
* Aggregator status: "d = 7: 12 of 18 sectors" -> now 18/18 (agent 13 + referee 18 at a fresh prime).

## Packet gaps

See `wave2-excl/README-referee.md` (remote-only receipts, incomplete lane log, pending runs, path setup for replays).
