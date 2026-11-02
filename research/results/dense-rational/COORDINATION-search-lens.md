# Coordination note (search lens, dense-rational) — written 2026-11-02 00:50 by the second search agent (scripts w2lib.py, w2_*.py)

Two search-lens agents were spawned on this lane. To avoid duplicate CPU on the overloaded box, the second agent
(w2_*.py, jobs /tmp/hodge/jobs/search-dense-rational-{dims,relax*,vp*,struct*}.log) covers ONLY:
  * Hessian relaxation at d = 7 even (already running, cheap), d = 8 even/odd/mixed, d = 9 even/odd  (w2_relax.py)
  * full normalized-E2 varpro at d = 7 MIXED parity, d = 8 even/odd, d = 9 even/odd               (w2_varpro.py)
  * exact mod-p E1-kernel dims d = 4..12 both parities via hodge.jacobian.mult_kernel (1 core)    (w2_dims.py)
  * structured linear families vs E1 at d = 10, 11, 12 (cones, bidegree forms, block products)    (w2_struct_e1.py)
The first agent (w2s_*.py: nvp7e/nvp7o, nvp10ker/nvp10full, w2s_hessrelax, w2s_ms7) keeps d = 7 pure sectors, d = 10 kernel
families, numerical E1 spaces and the msolve runs.  Both append batch lines to wave2-search.log.jsonl (distinct "batch"
labels: second agent uses "hessian-relaxation", "full-E2-varpro", "e1-kernel-dims", "structured-family-E1").
Report: the second agent writes its section into /tmp/hodge/notes/wave2-dense-rational-search.md between the markers
<!-- AGENT-B BEGIN --> / <!-- AGENT-B END --> (merge-preserving rewrite); please keep that block when rewriting.
