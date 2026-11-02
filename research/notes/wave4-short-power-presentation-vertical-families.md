STATUS: RUNNING — 2026-11-02 04:55 UTC — E1 vertical d=7..13 = J-part only (exact); D_res fails E1 (exact); E2-mod-J relaxation + fixed-D solves of u*w running

# wave 4 — short-power-presentation — route: vertical families (E1 on C[e(X), e(Y), L, M] and D_res)

Owner: wave-4 agent (route vertical-families). Prior work: wave-2 construct lens (`w2c_vert.py`, never ran: remote gate), wave-3 ranking §3.15.

## Method
`J(F) = I' + (F0, F3)` with `I' = (X1^5 - 2X0^5, X2^5 - 3X0^5, Y1^5 - 3Y0^5, Y2^5 - 5Y0^5)` (`F1 - 2F0 = 6(X1^5 - 2X0^5)` etc., exact).
`I'` is a Groebner basis with pure-power leading terms, so `T = S/I'` is free over `C[X0, Y0]` with basis `X1^a X2^b Y1^c Y2^e`, `a,b,c,e <= 4`.
E1 (`A D in J`) <=> `A D mod I'` lies in the span of `{F0 m, F3 m}` in `T_{d+6}` (parity sector of the target).
Transfer: rank_p(J-block) <= rank_Q(J-block) <= generic (complete-intersection) value; equality mod p ("Hilbert check") gives
dim ker_p >= dim ker_Q, so a zero kernel mod p is a zero kernel over Q (exact).
Observation: `X0 F0 + X1 F1 + X2 F2 = 6 e(X) + L h_L`, so `e(X) ≡ -L h_L/6`, `e(Y) ≡ -M h_M/6` mod J: the whole vertical ring
collapses onto `C[L, M]` modulo J.

## Results table
(filled as computations finish)
