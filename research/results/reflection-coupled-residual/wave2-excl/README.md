# Wave-2 exclusion packet — `reflection-coupled-residual` (q = 1, A = X0^3 Y0^3)

Everything here is reproducible with the engine at `/home/kejace/devel/sextics-work/hodge` (python-flint 0.9,
numpy, scipy), msolve 0.6+ and Singular 4.x. Variables X0,X1,X2,Y0,Y1,Y2. Klein character of a monomial:
chi(m) = (deg_{Y2} m mod 2, (deg_{Y0} m + deg_{Y1} m) mod 2). F has chi = (0,0), A has chi = (0,1).

## Statements (what is proven; receipts in the JSONL files; primes recorded in every receipt)

S1. For d <= 4 no denominator of degree d satisfies E1 (A D in J(F)) in any character  [receipts: e1char.jsonl]
S2. For d = 5, 6, 7 the E1-admissible space equals J_d in every character (kernel of [A] on (S/J)_d^chi is 0)
    [e1char.jsonl; modular kernel dimension bounds the rational one from above because the Hilbert function of
    S/J mod p equals the Fermat count in the relevant degrees (field `good_*` in the receipts)].
S3. Degree 5: for D = sum c_i F_i the E2-mod-J obstruction ob = -[A sum c_i c_j F_ij] has the 20 nonzero classes
    [A F_ij] (F_{X2 Y2} = 0) linearly independent in (S/J)_10 => c = 0.  [hess.jsonl, exp phi-kernel/phi-verdict d=5]
S4. Degree 6: on W = (S_1^6)/(C z) (35 unknowns, all characters jointly) the quadratic system ob(g) = 0 is
    zero-dimensional mod p (dim 0, mult 1060, HF(0..7) = 1,35,285,477,225,36,1,0) => I_Q contains S_7(W) =>
    V(ob) = T over C => D = g(F) in (F).  [gbcheck / hess_d6_* files, hess_d6_singular_minass.txt]
S5. Degree 7, pure characters: see gbcheck.jsonl (hess_d7_W<chi>_* files).
(further statements appended as experiments finish)

## Why the modular computations are proofs
* kernel / rank statements: for an integer matrix M, rank_p(M) <= rank_Q(M); when the Jacobian ideal has good
  reduction in the degree (rank J_n^chi mod p = dim S_n^chi - Fermat count), the mod-p quotient (S/J)_n^chi is the
  reduction of the rational one and the modular kernel dimension of [A] bounds the rational one from above.
* quadratic systems: HF_T(k) <= HF_Q(k) <= HF_p(k); HF_p(k0) = 0 for the system on W = V/T implies that the rational
  ideal contains all forms of degree k0 on W, i.e. the only solutions over C are the trivial ones (g(F) in (F)).
* the obstruction: for D = g(F) and any E1 solution N, B_D(N) is independent of N modulo J (Koszul corrections lie in J)
  and equals -[A sum g_i g_j F_ij]; B(zG, .) = 5[A G g(F)] = 0 exactly.
* fixed-D inconsistency (rigorfix.jsonl): modular inconsistency plus an exact rational kernel of dimension equal to
  the modular nullity (a-priori vectors (zG,0) and (F m e_i, D d_i m - m d_i D)) certifies rank_Q M = rank_p M.

## Replay
    python3 e1char.py 1 10 2147483587 2147483579
    python3 hess.py phi 5 2147483587 ; python3 hess.py export 6 2147483587 00 01 10 11 ; python3 gbcheck.py hess_d6_W00011011_p2147483587.ms 2
    python3 hess.py export 7 2147483587 00 ; python3 gbcheck.py hess_d7_W00_p2147483587.ms 2
    python3 rigorfix.py K8 "$(python3 -c 'import json;print(json.load(open("../kernel_d8_X03Y03.json"))[0])')" --full
