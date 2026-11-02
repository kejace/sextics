Note for the sibling capacity agent (shared remote dir), 2026-11-02 ~02:00 UTC:
- I measured kernel sparse-list arithmetic (your Reflect.lean approach, my Core0/1/2 variants): ~280 us and ~25 KB per merge
  step in the kernel => memory-bound around 3e5 merge steps per declaration; not a game changer vs ring.
- Alternative that scales: Kronecker big-integer evaluation, src/Reflect/Kron.lean (soundness: homogeneity + l1 bound +
  mixed-radix injectivity), generators gen/kron_gen.py (synthetic), gen/kron_submit.py (full proposal), gen/bind_test.py
  (HodgeWitness binding by rfl + matchW).  Remote results: K_m3_d6_k40 14.8 s user (template 205 s), K_m3_d8_k150 with
  e2 raw 435k products: 58 s user / 5.5 GB.  IMPORTANT: mark big data defs `noncomputable` (the code generator blows up on
  10k-leaf trees: that was the maxRecDepth failure of K_m3_d10_k*).
- Template ceiling seen so far: e2m3k150 (185k raw products) fails at 2M heartbeats after 605 s / 12.5 GB.
- Report (local to my machine, synced to the research branch): /tmp/hodge/notes/lean-capacity.md
