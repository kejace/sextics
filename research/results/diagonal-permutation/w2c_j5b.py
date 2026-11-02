"""J5b: E1-kernel (both parities) on R_6 of: the sign-type quadric Q_orb = u_X wedge u_Y (sweep of the S3-orbit planes of
transposition-fixed points), and the cone over S_01 with a random (non-equivariant) vertex line."""
import sys, json, time, random
sys.path.insert(0, '/home/kejace/devel/sextics-work')
exec(open('w2c_randcone.py').read().split("log = open")[0])  # reuse kernel_full and setup
Fs = sp.sympify(core.to_str(F))
log = open('wave2-construct.log.jsonl','a')
random.seed(7)
Qorb = sp.expand((X0-X1)*(Y1-Y2) - (X1-X2)*(Y0-Y1))
cands = {'Q_orb = (X0-X1)(Y1-Y2)-(X1-X2)(Y0-Y1) [sign rep, sweep of S3-orbit planes]': Qorb}
# cone over S_01 with random vertex line: projection to P^3_0 = {r=t=0} along a random 2-plane W: x -> x - (components along W)
# choose W spanned by w1,w2 random; the cone equation: F0(pi(x)) where pi kills W and fixes the complement {r=t=0}... simplest:
# coordinates (p,q,a,b,r,t); P^3_0 = {r=t=0}; a vertex line L = span(w1,w2) with w_i = (alpha_i, r=1/0, t=0/1) generic;
# pi(x) = x - r(x)*w1 - t(x)*w2 where w1,w2 have r,t-components (1,0),(0,1): then pi(x) has r=t=0.
p_,q_,a_,b_,r_,t_ = sp.symbols('p q a b r t')
sub_fw = {X0:(p_+r_)/2, X1:(p_-r_)/2, X2:a_, Y0:(q_+t_)/2, Y1:(q_-t_)/2, Y2:b_}
Fpq = sp.expand(Fs.subs(sub_fw, simultaneous=True))
F0 = Fpq.subs({r_:0, t_:0})
for trial in range(2):
    w1 = [random.randint(-3,3) for _ in range(4)]; w2 = [random.randint(-3,3) for _ in range(4)]
    pi = {p_: p_ - r_*w1[0] - t_*w2[0], q_: q_ - r_*w1[1] - t_*w2[1], a_: a_ - r_*w1[2] - t_*w2[2], b_: b_ - r_*w1[3] - t_*w2[3]}
    Dpq = sp.expand(F0.subs(pi, simultaneous=True))
    back = {p_: X0+X1, r_: X0-X1, a_: X2, q_: Y0+Y1, t_: Y0-Y1, b_: Y2}
    cands['cone over S_01, vertex line w1=%s w2=%s'%(w1,w2)] = sp.expand(Dpq.subs(back, simultaneous=True))
pp,aa,tt,rr,qq,bb = X0+X1, X2, Y0-Y1, X0-X1, Y0+Y1, Y2
cands['cone over Fix(tau s)_+ curve, vertex E_-'] = sp.expand(Fs.subs({X0:pp/2, X1:pp/2, X2:aa, Y0:tt/2, Y1:-tt/2, Y2:0}, simultaneous=True))
cands['cone over Fix(tau s)_- curve, vertex E_+'] = sp.expand(Fs.subs({X0:rr/2, X1:-rr/2, X2:0, Y0:qq/2, Y1:qq/2, Y2:bb}, simultaneous=True))
for name, D in cands.items():
    Dd = core.parse(str(D))
    t0=time.time(); k, rk = kernel_full(F, Dd, core.degree(Dd), p); dt=time.time()-t0
    r = {"exp":"J5b-divisorE1-kernel","D":name,"d":core.degree(Dd),"prime":p,"result":"ker dim %d (both parities)"%k,"seconds":round(dt,1),"job":"local","note":"0 => V(D) serves no primitive (3,1) class"}
    print(r, flush=True); log.write(json.dumps(r)+"\n"); log.flush()
