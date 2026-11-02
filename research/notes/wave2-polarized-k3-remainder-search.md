# Wave 2 — `polarized-k3-remainder` joint search (remote, 2026-11-01)

**STATUS: NO WITNESS FOUND (search running; this file is rewritten after every batch).**

Lane: q = 1, A = X0X1X2Y0Y1Y2, F = ½(f(U)+f(V)), U = X+Y, V = X−Y. Diagonal symmetry group = {1, τ} only, so the
"character sectors" are the two Y-parities of D. Wave-1 facts used: [A] has maximal rank on (R/J)_d in every degree
(injective d ≤ 9 ⇒ every E1-admissible D of degree ≤ 9 is Jacobian-type D = Σ a_i F_i; surjective d ≥ 10 with
kernel 555/1020/1325 at d = 10/11/12); degrees ≤ 7 excluded (thread); sum principle ⇒ only D with an irreducible
factor of degree ≥ 8 matter.

Method (new in wave 2, `w2lib.py`): reduced Koszul E2 system of book §5.2.  For Jacobian-type D = Σ a_i F_i take
N0_i = A a_i (exact E1 lift; the slack dF(N) = AD + FG is removable, G drops out of E2).  E2 ⟺
B_D(N0) + T_D(c) ∈ F·S_{2d−6}, c antisymmetric of degree d−4.  Q-block eliminated by the X0-remainder map (F monic
in X0).  Variable projection + Levenberg–Marquardt (Kaufman Jacobian) on the unit sphere of a, F-multiples
projected out.  Hessian relaxation = the same system modulo J(F): A·Hess_F(a,a) ∈ J(F)_{2d} (necessary), with
deflation against the known spurious family D = g·h, g = f(U)−f(V).

Log: `/tmp/hodge/results/polarized-k3-remainder/wave2-search.log.jsonl` (merged from remote `w2-*.jsonl`).

## Jobs
(see below; updated per batch)
