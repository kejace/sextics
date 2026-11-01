import sys, pickle, json
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
from five import *
for d in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    for par in (0, 1):
        info, basis = e1_kernel5(d, par=par)
        print(json.dumps(info), flush=True)
        with open(f'e1ker5_d{d}_par{par}.pkl', 'wb') as fh:
            pickle.dump({'p': P1, 'basis': basis, 'info': info}, fh)
        log(f'<5var E1 kernel sweep degree {d} parity {par}>', d, '5var:e1-kernel-sweep', info['kerD'] > 0, None,
            info['kerD'], info['ker_mod_J'], info['seconds'],
            f'5var: dim{{D rho-par {par} deg {d}: AD in J_f}}={info["kerD"]}, dim J_f,{d}^par={info["dimJ"]}, kernel mod J={info["ker_mod_J"]}, matrix {info["rows"]}x{info["cols"]}')
