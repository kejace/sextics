# README-referee — independent referee notes for packet `wave2-excl` (lane `cycle-q2`)

Written 2026-11-02 by the referee pass (Claude, Anthropic, for `kejace`); the agent's files were not modified.  Referee scripts,
receipts, regenerated certificates and logs live in the sibling directory `../referee/` (`/tmp/hodge/results/cycle-q2/referee/`,
mirrored to `research/results/cycle-q2/referee/` on the fork's `research` branch).

## Pinned environment used for the replays
python 3.12.3, python-flint 0.9.0 (import without `-I` for the agent's scripts; the referee verifier needs no flint),
Singular 4.3.2 (`sat` from `elim.lib` returns an ideal in this version; `sat(I, ideal)` saturates by an ideal), msolve 0.6.x
(not needed for any replay below; the 62-variable msolve run never finished).  Engine: `/home/kejace/devel/sextics-work`
(`hodge` v2, sign fix 2026-11-01); the engine is only needed to *regenerate* certificates, not to check them.

## What a stranger can replay, and how
1. Any stored Farkas certificate, with no engine and no sector machinery (full column set, exact integers):
   `cd wave2-excl && python3 -I ../referee/referee_farkas.py certs/B6.json certs/K10.json ...` (0.1-7 s each).
   The script checks `y^T M = 0` for every column of the FULL fixed-D system and `y^T (A D^2) != 0`.  Negative controls
   (`../referee/control_*.json`) fail as they should.  Results: `../referee/referee_farkas.results.jsonl`.
2. `Singular -q ann.sing` (minutes) — or `../referee/ann_ref.sing` at the fresh prime 999999937 (`../referee/ann_ref.out`).
3. The ten exact degree-6 sector relaxations: `for k in 001 002 004 005 022 024 040 050 302 310; do Singular -q sectors/exact_d6_k$k.sing; done`
   (each prints `DIM 0`; seconds).  The `.sing` files contain the exact quadrics over Q produced by `sector_exact.py`.
4. The trivial-sector cubic forms: `Singular -q polyfarkas/red_trivial6_exact_Q.sing` (5 s; `DIM_C_only 1`, `DIM_C_minus_Euler -1`).
5. The mixed-character unions (mod p): `python3 -I ../referee/referee_unions.py --size 2|3|4 --prime 1000000007` (3 s / 63 s / long);
   receipts `../referee/unions_size*_p*.json`.  This uses the correct saturation test (see gap G4).
6. Exact (over Q) union of the trivial sector with each nontrivial sector: `python3 ../referee/referee_pair_exact.py` (needs the engine;
   ~15 s per sector; receipts `../referee/pair_exact.results.jsonl`, Singular inputs `../referee/pair_exact_k*.sing`).
7. Exact repair of the nine singleton sectors: `python3 ../referee/referee_singletons.py` (needs the engine; ~1 min per sector),
   then step 1 on `../referee/certs_regen/exactQ_d6_k*.json`.

## Gaps found in the packet (none overturns the degree-5 bound; several change the degree-6 wording)
G1. 80 of the 447 certificate files store `y = null` (storage cap); they are receipts, not certificates.  Two of them were
    load-bearing for "every single-character degree-6 denominator" (`d6_..._k014_c0`, `d6_..._k303_c0`): regenerated in
    `../referee/certs_regen/` and verified; the rest are degree-7 elements and random points (screening only).
G2. The 60 exact polynomial-Farkas vectors of the trivial sector are not stored (only SHA-256 and the resulting cubic forms), so the
    identity `y(a)^T M(a) = 0` cannot be replayed from the packet; the zero-locus computation (step 4) can.  A sampled regeneration
    (3 vectors) is in `../referee/packet_shadow/polyfarkas/` with log `../referee/redlift_sample.log`.
G3. `sectors.py` tests the SYMMETRIC MOD-p LIFT of each kernel vector (line 169), which is in general not an E1-admissible denominator
    over Q.  Hence the 263 degree-7 and 48 degree-6 "element/random point" certificates certify inconsistency of integer polynomials that
    are (generically) not in the Q-kernel; they do not bear on the Q-families.  For the nine degree-6 singleton sectors the exact Q-kernel
    element has ~300-bit coefficients; the referee repair (step 7) decides them exactly.
G4. `subsets6.py` and `pair_exact.py` saturate by the PRODUCT of the nontrivial coordinates; that only proves "some nontrivial coordinate
    vanishes".  The needed statement ("all vanish", i.e. the union variety equals V0) requires saturation by the IDEAL of those
    coordinates (step 5/6).  With the correct test the mod-p results still hold at both primes for sizes 2, 3 (size 4: see receipts).
G5. The reduction argument T2 transfers EMPTINESS of a mod-p relaxation variety to Q-bar, so unions of nontrivial sectors are rigorously
    excluded.  For unions containing the trivial sector the mod-p variety is not empty (it is V0 mod p) and T2 does not exclude a
    denominator whose nontrivial components are p-adically small: those cases are rigorous only when computed exactly over Q (step 6
    does this for 27 pairs incl. all four order-2 characters; for the other pairs the kappa-kappa part of [H(D_a + D_c)] has its own
    character (2*kappa != 0) so the single-sector exact results already force D_c = 0; triples and quadruples remain mod-p screened).
G6. 18 ledger rows of `reverify.results.jsonl` are `ERROR ... mixed Y-parity` (the squared mixed-parity cubic 3-cones were never
    tested); the local packet has 101 named, 35 cone, 48 degree-6 and 263 degree-7 certificate files (the remote had 608 at hand-off).

## Referee artefact notes (added after the replays)
* `../referee/singletons_hessian.results.jsonl`: all nine singleton sectors have exact Q-kernel dimension 1 and a NONZERO exact Hessian class
  in `R_18` (`referee_singletons_hessian.py`, ~8 s per sector) — this alone excludes them over C.  `referee_singletons.py` additionally
  produces exact Farkas certificates for the true kernel elements (`../referee/certs_regen/exactQ_d6_k*.json`; entries of 10^4–10^5 bits,
  ~1–10 min per sector); `exactQ_d6_k013.json` (34 MB) exceeds the research-branch per-file cap and is mirrored as `.json.gz`.
* `../referee/pair_exact.results.jsonl`, `pair_exact_order2.results.jsonl`, `pair_exact_k302_pervar.sing`: exact over Q, union variety of
  {trivial + one nontrivial sector} equals V0 (ideal saturation) for 27 sectors including all four order-2 characters (0,0,3), (3,0,0),
  (3,0,3), (3,3,3); for a character kappa with 2*kappa != 0 the kappa-kappa part of [H(D_a + D_c)] has character chi_A + 2*kappa, distinct
  from chi_A and chi_A + kappa, so the single-sector exact result already forces D_c = 0 — the pair computation is only needed for order-2
  characters.
* `../referee/unions_size{2,3,4}_p*.json`: corrected (ideal-saturation) union test; sizes 2, 3 at both primes, size 4 at 10^9+7.
* `../referee/packet_shadow/polyfarkas/`: sampled regeneration (3 of 60) of the trivial-sector exact Farkas vectors; the regenerated cubic
  forms coincide bitwise with the stored cubics #0, #20, #40 of `polyfarkas/red_trivial6_exact.cubics_Q.txt`.
* `../referee/referee_farkas.results.jsonl`: 71 engine-free exact certificate checks (69 verified; 2 files without stored vector).
* Tools actually used: python 3.12.3, python-flint 0.9.0, Singular 4.3.2; primes 10^9+7, 998244353, 999999937 (msolve not needed).
