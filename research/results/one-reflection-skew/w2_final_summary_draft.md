Line 1: NO WITNESS FOUND (one-reflection-skew; wave-2 joint search, 5-variable exact reduction).
Line 2: Reduction F=f(x)+Y2^6, A Y2-free => certificate iff rational u=N/D solves A D^2 - D N(f) + (6/5) f (D div N - N(D)) = 0 (equivalent to the strict Y2-free form); validated on family-hesse-web (reduced witness exact, operator rel 2.6e-15, recovered + 6-var verified).
Line 3: E1-admissible dims exact for d=4..12 in all four (rho,sigma) sectors; non-Jacobian directions from d=7 (7/1/0/1), kernel dims 28/19/18/19 (d=7), 91/76/67/76 (d=8), ... 473/445/442/445 (d=12); exact rational kernel bases via CRT (heights up to 2^1137 at d=9).
Line 4: Jacobian-type degrees d=5,6 (6-var (w.grad f)^2 of degree 10/12): exhaustive varpro over K<=7 families in all sectors, global minima rel >= 6.5e-2 (no zero).
Line 5: Kernel-family LM varpro at d=7 (K=28) ... [restarts, best rel]; d=8 ... ; cones ...; structured ...
Line 6: Calibration: LM hit rate ~0.15-0.5 per restart at K=32 (needle = accepted degree-6 witness), L-BFGS ~0.05; snap-to-rational needed because zero sets are positive-dimensional.
Line 7: Engine/varpro flaw found: exclude_F_multiples removes legitimate solutions (noted in engine-bugs.md); float E1 kernels unreliable here because the skew term is a 1/7-perturbation of one-reflection (near-kernel directions).
Line 8: Remote box saturated (load 900-1300, ~60 foreign jobs); 4-process chains forced through the gate; two chains remain queued.
Line 9: Report /tmp/hodge/notes/wave2-one-reflection-skew-search.md; log /tmp/hodge/results/one-reflection-skew/wave2-search.log.jsonl; code w2lib.py / w2_search.py / w2_struct.py / chains/.
Line 10: Next: finish d=8-12 kernel searches, A+ half-target, one-reflection control; exact d=7 exclusion is being pursued by the exclusion lens.
