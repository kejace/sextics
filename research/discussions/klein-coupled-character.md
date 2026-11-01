# Coupled Klein · nontrivial character  (https://github.com/hodge-computer/sextics/discussions/116)


## COMMENT DC_kwDOVLOOJc4BJMgE by fqlx at 2026-10-26T03:27:03Z (updated 2026-10-26T03:27:03Z)

Partial update for `klein-coupled-character`: the original $q=1$, $A=X_0^3Y_0^3$ target is unchanged, and **step 23 remains unsolved**.

Two construction families now have explicit obstructions:

1. On the source-adapted degree-21 support at $\epsilon=1$, the three fixed $t/r/s$ Hamiltonians, **individually balanced at the 13 ordinary residue orbits**, give distinct compact controls. Cofactor elimination forces the only possible rational base-amplitude candidate. After clearing its denominator, the necessary identity has residual $(39,22,26,84,83)\ne0$ modulo 101 at the specified good point. This excludes that three-control family even with arbitrary rational amplitude poles, retaining the amplitude-dependent fiber corrections. Pointwise first-jet surjectivity did not give a rational solution.

2. A **single normalized variable-trigonal trace**, with arbitrary rational $w(a,b)$ and one common seed, cannot produce the required boundary source $K_\eta$. On a proper resolution of the degree-seven density cover, every polar curve is rational except one deck-stable curve. The seven target norm-pole curves are Klein curves cyclically permuted by the deck action; properness, residues and equivariance give a contradiction. The weighted source is regular at infinity; that boundary was retained.

[Derivations, scripts, recorded outputs, reproduction commands and SHA-256 manifest](https://gist.github.com/fqlx/a66456607d1bb86dee4294fa2128f692/06d3d8e2bd59cb6496b1a31e9b9399fdebf82d0a).

These conclusions exclude only the stated families. Collective residue cancellation **before individual balancing**, independent seeds/weights, additional compact controls and fiber-dependent amplitudes remain open. No full original $D/N/Q$ certificate or Lean proof has been obtained.

AI-assisted research and update prepared with Codex/GPT-6 for `fqlx`; prior human references are credited in the packet. No additional mathematical replay was run solely for this post.

