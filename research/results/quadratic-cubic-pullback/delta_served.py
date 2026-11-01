"""Is delta (the determinant density) served by the pulled-back join denominators?  Fixed-D solve with target A := delta."""
import sys, time, json
sys.path.insert(0,'/tmp/hodge/results/quadratic-cubic-pullback')
import qcp
from qcp import core, fixed_d, P
delta = P(open(qcp.WD + '/delta.txt').read())
names = sys.argv[1:]
Q = qcp.Q
Ds = {'join012': core.add(core.add(core.power(Q[0],3),core.power(Q[1],3)),core.power(Q[2],3)),
      'join345': core.add(core.add(core.power(Q[3],3),core.power(Q[4],3)),core.power(Q[5],3)),
      'join034': core.add(core.add(core.power(Q[0],3),core.power(Q[3],3)),core.power(Q[4],3)),
      'Fpart6': {e:c for e,c in qcp.F.items() if e[0]+e[1]+e[2]==6},
      'fermatX': P('X0**6 + X1**6 + X2**6')}
for nm in names:
    D = Ds[nm]
    t0 = time.time()
    ok, info = fixed_d.e1_solvable(qcp.F, delta, 1, D)
    print(nm, 'E1(delta):', ok, {k:v for k,v in info.items() if k in ('rank','cols')}, flush=True)
    out = fixed_d.solve_fixed_D_fast(qcp.F, delta, 1, D, verbose=True, want_exact=False)
    rec = {"D": nm, "d": 6, "family": "target=delta", "e1": ok, "consistent": out['consistent'], "rank": out['rank'], "nullity": out['nullity'], "seconds": round(time.time()-t0,1), "note": f"target A:=delta with D={nm}", "engine": 2}
    print(json.dumps(rec), flush=True)
    with open(qcp.LOG,'a') as f: f.write(json.dumps(rec)+'\n')
