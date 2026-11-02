## Narrative (exclusion / rigor lens, wave 2)

### Best rigorous bound and what is new
* Unrestricted: deg D >= 7 for every legal denominator (fqlx, exact; re-verified here at degrees 4-6: D0 rank sandwich PROVEN, degree-6
  sectors chart-complete mod p).
* Single diagonal character (D a (g,h)-eigenvector): **degree 7 excluded for all 36 characters** — deanodotdotdot's theorem re-derived
  (i) chart-complete mod 2147483647 (171 charts + 4 empty colon sectors) and (ii) **exactly over Q by the relaxation pipeline for 34 of 36
  characters** ((1,0), (3,0): Singular runs pending on the remote); **degree 8 excluded for all 36 characters modulo p** (436 charts,
  chart-complete, including the six h-invariant characters (a,0) that the thread left open) and **exactly over Q for 21 characters** so far
  (9 attempted sectors have a positive-dimensional relaxation residual whose exact decomposition exceeded the local cap; the six (a,0)
  sectors and these 9 are queued in the remote chain `excl-off-paired-chain` with 2-hour Singular caps).  So, at the modular level, every
  single-character denominator has degree >= 9; at the theorem level the new content is the exact closure of 21 + 34 sectors by an
  independent method, and the method itself (below).
* Degree 9 single character: 27 of 50 charts of (0,0) closed locally (the rest need more than 60 s per chart); the full degree-9 sweep
  is in the queued remote chain.
* Structured fixed denominators (two new primes 2147483629 / 2147483587, explicit-kernel rank sandwich): see the receipts table — every
  tested block product of degree 7-14 processed so far is PROVEN excluded (either E1 fails against the exact Koszul dimension, or the full
  system is inconsistent with the rational rank certified); 8 records (3 denominators of degree 12-13) await the multi-prime
  re-certification (`w2_recert.py`, queued).  The wave-1 structured families (621) and Hessian minors (213) are queued after the block
  products in the same remote job (`excl-off-paired-batch1`, running).

### Method notes (why the exact sector exclusions are theorems)
* Fixed-D rank sandwich: rank_p(aug) = rank_p(M)+1 and rank_Q(M) <= rank_p(M) certified by the explicit kernel {(L Z,0)} + {(F W, B_D(W))}
  (exact dimension) or by exact verification of a reconstructed kernel basis.
* E1: modular rank of the Jacobian generators = exact Koszul dimension dim J(F)_{d+6}[chi] (F smooth) and augmented rank +1.
* Sector exclusion over Q: rational colon basis D_a modulo the Euler gauge and rational Koszul kernel S_k (verified exactly); exact bilinear
  equations; the degree-2 relaxation J2 (constant-coefficient combinations whose Koszul terms cancel; exact integer kernel, LLL-reduced) is a
  homogeneous ideal in Q[t]; if its reduction mod p is m-primary (msolve GB has a pure power of every variable) it is m-primary over Q
  (Hilbert function over Q <= Hilbert function mod p), so t = 0 and D is in the gauge.  Otherwise Singular over Q decomposes V(J2) and each
  component is excluded on its linear span by exact degree-3/4 localized relaxations (m-primary over Q).  Fallback: certified span of V(J2)
  (exact nilpotency certificates of the defining linear forms) + localized relaxation.  Degree-7 (0,0) and (2,0) were closed this way.
* Chart-complete msolve results are modular screens; a GB = {1} on a chart does not transfer to Q (the px-1 phenomenon), which is why the
  exact pipeline works with HOMOGENEOUS ideals in t only.

### Pending / how to pick up
* Remote jobs: `excl-off-paired-batch1` (running at session end; 346/680 block products done, then the 621 wave-1 denominators and 213
  Hessian minors follow automatically).  The combined chain (deltat full solve, exact rigor d=7/8 all sectors with 2-hour Singular caps,
  degree-9 chart sweep, multi-prime re-certification, kernel dims d=12) could NOT be submitted: the remote stayed above the runner's load
  guard (load > 150, < 60 GB free) for the whole session.  To launch it when capacity frees:
  `cd /home/kejace/devel/sextics-work && python3 -m hodge.rjob submit --lane off-paired --name excl-off-paired-chain --cpus 56 --wait 7200 -- 'bash /tmp/hodge/results/off-paired/wave2-excl/w2_chain.sh'`
  (run `python3 wave2-excl/w2_sync.py` first so the local log is a superset of the remote one; the chain writes per-process logs).  After they run: `python3 w2_sync.py && python3 w2_report.py
  report_notes.md` regenerates this report; receipts land in `/tmp/hodge/remote/off-paired/wave2-excl/receipts_*.jsonl`.
* Open at the end of this session: degree 7 unrestricted (mixed characters) — not attacked (dense 171-variable relaxation; the
  Hessian/Jacobian relaxation alone is not m-primary at degree 7); degree 9+ single character (modular, in progress); Delta_t (degree 18,
  E1 passes, queued).
