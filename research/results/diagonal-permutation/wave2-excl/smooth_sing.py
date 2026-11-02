"""Smoothness of F modulo the engine primes: vdim of J(F) over Z/p (Singular). 15625 = 5^6 <=> zero-dimensional
Jacobian ideal <=> the six partials are a regular sequence mod p <=> Hilbert function of S/J(F) mod p is Koszul."""
import subprocess, sys, time
from common import *
F, A, q, rec = core.load_instance(IID)
primes = [int(x) for x in sys.argv[1:]] or [2147483647, 2147483629, 1000003]
for p in primes:
    script = f"""
ring r = {p}, (X0,X1,X2,Y0,Y1,Y2), dp;
poly F = {singular_poly(F)};
ideal J = jacob(F);
option(redSB);
ideal G = std(J);
int v = vdim(G);
print("VDIM " + string(v));
intvec h = hilb(G, 1);
print("HILB1 " + string(h));
quit;
"""
    open(f'smooth_{p}.sing', 'w').write(script)
    t0 = time.time()
    out = subprocess.run(['Singular', '-q', f'smooth_{p}.sing'], capture_output=True, text=True, timeout=3600).stdout
    vd = [l for l in out.splitlines() if l.startswith('VDIM')]
    vdim = int(vd[0].split()[1]) if vd else None
    log({'exp': 'smooth-modp-vdim', 'p': p, 'vdim': vdim, 'koszul_expected': 15625, 'smooth_mod_p': vdim == 15625,
         'seconds': round(time.time() - t0, 1), 'tool': 'Singular 4.3.2 std/vdim', 'raw': out[-400:]})
