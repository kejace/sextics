# wave2-excl packet — cyclic-triple-residual (q = 1, A = X0^3 Y0^3) — exclusion / rigor lens

All computations use the exact instance from `data/cyclic-triple-residual/instance.rec` (F = Fermat sextic + (5/7)*(six
couplings), variables X0,X1,X2,Y0,Y1,Y2).  Engine: `/home/kejace/devel/sextics-work/hodge` (python-flint 0.9), msolve
0.6.x, Singular 4.3.2 on the remote box.  Every experiment appends a JSON line to `../wave2-excl.log.jsonl`
(keys: exp, D/family, d, prime, result, rank, seconds, job, note); receipts named `*.json` here.

## Scripts (this directory)
* `w2lib.py`        — shared library: sectors of the order-18 diagonal group, Jacobian reducers mod p, Hessian-relaxation
                       quadrics, sealed fixed-D exclusion (modular ranks + exact rational null vectors), exact E1-kernel
                       extraction (multi-prime CRT + rational reconstruction + exact E1 solve).
* `kerverify.py`    — (a) E1-kernel dimensions by sector at new primes + c.i. Hilbert-function check; exact facts about D8.
* `d8full.py`       — (a) D8 exclusion re-verified: unrestricted system mod p; sector system sealed; exact fmpq rref.
* `relax_run.py`    — (b) Hessian relaxation in degree d over the FULL E1-admissible family (all parities/characters):
                       quadrics mod p, msolve reduced GB, gauge-fixed m-primary test (rigorous), optional per-form tests.
* `relax0.py`       — (b,c) degree-0 Koszul relaxation (functionals killing F*S and all Koszul images), full or per sector.
* `wspan.py`        — diagnostic: rank of the Koszul span vs pi(J) per sector.
* `kernel_multi.py` — (c) exact E1-kernel elements at degree d in every sector + sealed fixed-D tests of each element.
* `families.py`     — (c) sealed exclusion of structured denominators (D8 products, polar products, pieces of F, ...).
* `sectors.sh`      — (c) chain of per-sector Koszul relaxations at degree d.

## How to replay
```
cd /tmp/hodge/results/cyclic-triple-residual/wave2-excl
python3 kerverify.py 1000003,32003 11          # E1 table + e1_d8_exact.json
SKIP_FORMS=1 python3 relax_run.py 6 1000003 8   # gauge-fixed m-primary test at d = 6 (seconds)
SKIP_FORMS=1 python3 relax_run.py 7 1000003 32  # d = 7 (minutes)
SKIP_FORMS=1 python3 relax_run.py 8 1000003 32  # d = 8 (needs e1_d8_exact.json)
python3 kernel_multi.py 9 6 16                  # exact degree-9 kernel elements + sealed tests
bash sectors.sh 9 1000003 8                     # per-sector Koszul relaxations at d = 9
python3 families.py 16 12                       # structured families
```
Statements, primes and ranks: see `/tmp/hodge/notes/wave2-cyclic-triple-residual-excl.md`.
